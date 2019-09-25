#include "NucleotideMatrix.h"
#include "Sequence.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "LocalParameters.h"
#include "simd.h"

#ifdef OPENMP
#include <omp.h>
#endif

#include <omptl/omptl_algorithm>
#include <cstdlib>

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

    Debug::Progress progress(seqDb.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        Kmer* threadKmerBuffer = new Kmer[BUFFER_SIZE];
        size_t bufferPos = 0;
        Sequence seq(par.maxSeqLen, seqDb.getDbtype(), &subMat,  par.kmerSize, false, false);

#pragma omp for schedule(static)
        for (size_t id = 0; id < seqDb.getSize(); id++) {
            progress.updateProgress();

            char *seqData = seqDb.getData(id, thread_idx);
            unsigned int seqLen = seqDb.getSeqLen(id);

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

                size_t revKmerIdx = Util::revComplement(kmerIdx, par.kmerSize);
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

        delete[] threadKmerBuffer;
    }
    Debug(Debug::INFO) << "Time for extracting kmers: " << timer.lap() << "\n";
    timer.reset();

    Debug(Debug::INFO) << "Sort kmer ... ";
    omptl::sort(allKmers, allKmers + offset, Kmer::compareRepSequenceAndIdAndPos);
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
