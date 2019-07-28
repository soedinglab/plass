#include <stdlib.h>
#include <mmseqs/lib/omptl/omptl_algorithm>
#include "Timer.h"
#include "NucleotideMatrix.h"
#include "ReducedMatrix.h"
#include "Indexer.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "simd.h"
#include "Debug.h"
#include "AminoAcidLookupTables.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "LocalParameters.h"

#ifdef OPENMP
#include <omp.h>
#endif

#define HI_NIBBLE(b) (((b) >> 4) & 0x0F)
#define LO_NIBBLE(b) ((b) & 0x0F)

struct Kmer {
    struct TwoLetters {
        unsigned int first : 4;
        unsigned int last : 4;
    };

    uint64_t kmer;
    unsigned int id;
    TwoLetters firstAndLastLetter;
    short pos;
    bool isReverse;
    Kmer() {}

    static bool compareRepSequenceAndIdAndPos(const Kmer &first, const Kmer &second) {
        if (first.kmer < second.kmer)
            return true;
        if (second.kmer < first.kmer)
            return false;
        if (first.id < second.id)
            return true;
        if (second.id < first.id)
            return false;
        if (first.pos < second.pos)
            return true;
        if (second.pos < first.pos)
            return false;
        return false;
    }
};

std::pair<size_t, char> getSubstituion(const char lastLetter, const size_t currKmerIndex, const bool isReverse,
                                       const char * reverseLookup, Kmer *pKmer, const size_t currPos, const size_t maxSize);


void printKmer(size_t idx, int size);

// Compute reverse complement of k-mer in 2-bit-per-nucleotide encoding (A: 00, C: 01, G: 10, T: 11)
uint64_t revComplement(const uint64_t kmer, const int k) {
    // broadcast 64bit to 128 bit
    __m128i x = _mm_cvtsi64_si128(kmer);

    // create lookup (set 16 bytes in 128 bit)
    // a lookup entry at the index of two nucleotids (4 bit) describes the reverse
    // complement of these two nucleotid in the higher 4 bits (lookup1) or in the
    // lower 4 bits (lookup2)
    __m128i lookup1 = _mm_set_epi8(0x50,0x10,0xD0,0x90,0x40,0x00,0xC0,0x80,0x70,
                                   0x30,0xF0,0xB0,0x60,0x20,0xE0,0xA0);
    __m128i lookup2 = _mm_set_epi8(0x05,0x01,0x0D,0x09,0x04,0x00,0x0C,0x08,0x07,
                                   0x03,0x0F,0x0B,0x06,0x02,0x0E,0x0A);

    // _mm_set1_epi8: create 128 bit with all bytes set to given value
    // here: 0x0F (00001111) and 0xF0 (11110000)
    // _mm_and_si128: bitwise AND
    __m128i kmer1 = _mm_and_si128(x, _mm_set1_epi8(0x0F)); // get lower 4 bits
    __m128i kmer2 = _mm_and_si128(x, _mm_set1_epi8(0xF0)); // get higher 4 bits

    // shift right by 2 nucleotids
    kmer2 >>= 4;

    // use _mm_shuffle_epi8 to look up reverse complement
    kmer1 = _mm_shuffle_epi8(lookup1, kmer1);
    kmer2 = _mm_shuffle_epi8(lookup2, kmer2);

    // _mm_or_si128: bitwise OR
    x = _mm_or_si128(kmer1, kmer2);

    // set upper 8 bytes to 0 and revert order of lower 8 bytes
    x = _mm_shuffle_epi8(x,
                         _mm_set_epi8(0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0,1,2,3,4,5,6,7));

    // shift out the unused nucleotide positions (1 <= k <=32 )
    // broadcast 128 bit to 64 bit
    return (((uint64_t)_mm_cvtsi128_si64(x)) >> (uint64_t)(64-2*k));

}

int correctreads(int argc, const char **argv, const Command& command)  {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.kmerSize = 5;
    par.parseParameters(argc, argv, command, true, 0, 0);

    Debug(Debug::INFO) << "Sequence database: " << par.db1 << "\n";
    DBReader<unsigned int> seqDb (par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    seqDb.open(DBReader<unsigned int>::NOSORT);
    NucleotideMatrix subMat(par.scoringMatrixFile.nucleotides, 1.0, 0.0);
    Debug(Debug::INFO) << "Output database: " << par.db2 << "\n";
    DBWriter dbw(par.db2.c_str(), par.db2Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_NUCLEOTIDES);
    dbw.open();
    const unsigned int BUFFER_SIZE = 1024;

    //ACTG 0123 TGAC
    char reverseLookup[4] = {2, 3, 0, 1};
    Timer timer;
    timer.reset();
    size_t offset = 0;
    Kmer * allKmers = new Kmer[seqDb.getAminoAcidDBSize()*2];
    // Create a 1D Tensor on length 20 for input data.
    Debug(Debug::INFO) << "Extract kmers\n";
#pragma omp parallel
    {
        Kmer * threadKmerBuffer = new Kmer[BUFFER_SIZE];
        size_t bufferPos = 0;
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        Sequence seq(par.maxSeqLen, seqDb.getDbtype(), &subMat,  par.kmerSize, false, false);

#pragma omp for schedule(static)
        for (size_t id = 0; id < seqDb.getSize(); id++) {
            char *seqData = seqDb.getData(id, thread_idx);
            unsigned int seqLen = seqDb.getSeqLens(id)-2;

            unsigned int dbKey = seqDb.getDbKey(id);
            unsigned int pos = 0;
            while (pos < (seqLen - (par.kmerSize -1))) {
                const char *kmer = seqData+pos;
                uint64_t kmerIdx = 0;
                for(size_t kmerPos = 0; kmerPos < par.kmerSize; kmerPos++){
                    kmerIdx = kmerIdx << 2;
                    kmerIdx = kmerIdx |  (kmer[kmerPos]>>1)&3;
                }
                threadKmerBuffer[bufferPos].kmer = kmerIdx;
                threadKmerBuffer[bufferPos].pos = static_cast<short>(pos);
                char firstLetter = static_cast<char>((kmer[0]>>1)&3);
                char lastLetter = static_cast<char>((kmer[par.kmerSize-1]>>1)&3);
                threadKmerBuffer[bufferPos].firstAndLastLetter.first = firstLetter;
                threadKmerBuffer[bufferPos].firstAndLastLetter.last  = lastLetter;
                threadKmerBuffer[bufferPos].isReverse=false;
                bufferPos++;

                size_t revKmerIdx = revComplement(kmerIdx, par.kmerSize);
//                printKmer(kmerIdx, par.kmerSize);
//                printKmer(revKmerIdx, par.kmerSize);
                threadKmerBuffer[bufferPos].kmer = revKmerIdx;
                threadKmerBuffer[bufferPos].firstAndLastLetter.first=reverseLookup[lastLetter];
                threadKmerBuffer[bufferPos].firstAndLastLetter.last =reverseLookup[firstLetter];
                threadKmerBuffer[bufferPos].isReverse=true;

                if (bufferPos+1 >= BUFFER_SIZE) {
                    size_t writeOffset = __sync_fetch_and_add(&offset, bufferPos);
                    memcpy(allKmers + writeOffset, threadKmerBuffer, sizeof(Kmer) * bufferPos);
                    bufferPos = 0;
                }
                bufferPos++;
                pos++;
            }
        }
        if(bufferPos > 0){
            size_t writeOffset = __sync_fetch_and_add(&offset, bufferPos);
            memcpy(allKmers + writeOffset, threadKmerBuffer, sizeof(Kmer) * bufferPos);
        }

        delete [] threadKmerBuffer;
    }
    Debug(Debug::INFO) << "Done\n";
    Debug(Debug::INFO) << "Time for extracting kmers: " << timer.lap() << "\n";

    Debug(Debug::INFO) << "Sort kmer ... ";
    timer.reset();
    omptl::sort(allKmers, allKmers + offset, Kmer::compareRepSequenceAndIdAndPos);

    Debug(Debug::INFO) << "Done\n";
    Debug(Debug::INFO) << "Time for sort: " << timer.lap() << "\n";

    if(allKmers[0].kmer != allKmers[1].kmer){
        std::pair<size_t, char> ret = getSubstituion(allKmers[0].firstAndLastLetter.last, allKmers[0].kmer,
                                                     allKmers[0].isReverse, reverseLookup, allKmers, 0, offset);
        allKmers[0].kmer  = ret.first;
    }
    for(size_t pos = 1; pos < offset; pos++){
        Kmer & currKmer = allKmers[pos];

        // correct if it is a singleton
        if(allKmers[pos-1].kmer != currKmer.kmer &&
           allKmers[pos+1].kmer != currKmer.kmer){
            std::pair<size_t, char> ret =getSubstituion(allKmers[pos].firstAndLastLetter.last, allKmers[pos].kmer,
                                                        allKmers[pos].isReverse, reverseLookup, allKmers, pos, offset);
            allKmers[pos].kmer = ret.first;
        }
    }

    dbw.close();
    seqDb.close();

    return EXIT_SUCCESS;
}

void printKmer(size_t idx, int size) {
    char output[32];
    char nuclCode[4] = {'A','C','T','G'};
    int temp = idx;
    for (int i=size-1; i>=0; i--)
    {
        output[i] = nuclCode[ idx&3 ];
        idx = idx>>2;
    }
    output[size]='\0';
    std::cout << output << std::endl;
}

std::pair<size_t, char> getSubstituion(const char lastLetter, const size_t currKmerIndex, const bool isReverse,
                                       const char * reverseLookup, Kmer *kmerArray, const size_t currPos, size_t maxSize) {
    //    ATTGA 0
    //    ATTGT 3
    //    ATTGT 3
    //    ATTTA <- A + 3
    //    ATTTC <- C + 2, C - 1
    //    ATTTG <- C + 1, C - 2
    //    ATTTG
    //    ATTTT <- T - 3
    //    ATTTT




    size_t startIndex = 0;
    size_t endIndex = 0;
    switch (lastLetter){
        case 0:
            startIndex = currKmerIndex + 1;
            endIndex   = currKmerIndex + 3;
            break;
        case 1:
            startIndex = currKmerIndex - 1;
            endIndex = currKmerIndex + 2;
            break;
        case 2:
            startIndex = currKmerIndex - 2;
            endIndex   = currKmerIndex + 1;
            break;
        case 3:
            startIndex = currKmerIndex - 3;
            endIndex   = currKmerIndex - 1;
            break;
        default:
            std::cout << "this should not happen" << std::endl;
            break;
    }
    // backward search
    size_t pos = (currPos != 0) ? currPos - 1 : 0;
    bool foundKmer = false;
    char corrLastLetter=-1;
    size_t prevKmer = SIZE_MAX;
    size_t corrKmer= SIZE_MAX;
    while(pos > 0 && kmerArray[pos].kmer >= startIndex ){
        if(prevKmer == kmerArray[pos].kmer){
            corrKmer = kmerArray[pos].kmer;
            corrLastLetter = kmerArray[pos].firstAndLastLetter.last;
            foundKmer = true;
            break;
        }
        prevKmer = kmerArray[pos].kmer;
        pos--;
    }

    prevKmer = SIZE_MAX;
    pos = currPos + 1;
    while(foundKmer == false && pos < maxSize && kmerArray[pos].kmer <= endIndex ){
        if(prevKmer == kmerArray[pos].kmer){
            corrKmer = kmerArray[pos].kmer;
            corrLastLetter = kmerArray[pos].firstAndLastLetter.last;
            foundKmer = true;
            break;
        }
        prevKmer = kmerArray[pos].kmer;
        pos++;
    }

    if(foundKmer==true){

        char replaceChar = (isReverse) ? reverseLookup[corrLastLetter] : corrLastLetter;
        return std::make_pair(corrKmer, replaceChar);
    }

    return std::make_pair(currKmerIndex, lastLetter);
}

#undef HI_NIBBLE
#undef LO_NIBBLE