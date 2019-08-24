/*
 * Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
 * detect circular fragments
 */

#include "DBReader.h"
#include "DBWriter.h"
#include "Indexer.h"
#include "LocalParameters.h"
#include "NucleotideMatrix.h"

#include <algorithm>

#define HIT_RATE_THRESHOLD 0.24
// threshold to distinguish cyclic/terminal redundant genomes from random hits on linear genomes
// chosen based on analysis on ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ (last modified 7/11/19)


int checkcycle(int argc, const char **argv, const Command& command) {

    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> *seqDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str(),  par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    seqDbr->open(DBReader<unsigned int>::NOSORT);
    //TODO: DBReader<unsigned int> *seqDbr ?

    DBWriter cycleResultWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_NUCLEOTIDES);
    cycleResultWriter.open();

    DBWriter linearResultWriter(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_NUCLEOTIDES);
    linearResultWriter.open();

    const size_t kmerSize = par.kmerSize;

    //TODO: openmp, threads
    //TODO: splits
    unsigned int thread_idx = 0;

    int seqType  =  seqDbr->getDbtype();
    BaseMatrix *subMat;

    if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.nucleotides, 1.0, 0.0);
    }else {
        Debug(Debug::ERROR) << "Module checkcyle only supports nucleotide input database" << "\n";
        EXIT(EXIT_FAILURE);
    }

    struct kmerSeqPos {
        size_t kmer;
        unsigned int pos;

        static bool compareByKmer(const kmerSeqPos &first, const kmerSeqPos &second) {
            if (first.kmer < second.kmer)
                return true;
            if (second.kmer < first.kmer)
                return false;
            if (first.pos < second.pos)
                return true;
            if (second.pos < first.pos)
                return false;
            return false;
        }
    };

    Indexer indexer(subMat->alphabetSize - 1, kmerSize);
    Sequence seq(par.maxSeqLen, seqType, subMat, kmerSize, false, false);
    kmerSeqPos *frontKmers = new(std::nothrow) kmerSeqPos[par.maxSeqLen/2+1];
    Util::checkAllocation(frontKmers, "Can not allocate memory");
    kmerSeqPos *backKmers = new(std::nothrow) kmerSeqPos[par.maxSeqLen/2+1];
    Util::checkAllocation(backKmers, "Can not allocate memory");

    unsigned int *diagHits = new unsigned int[par.maxSeqLen/2+1];
    for (size_t id = 0; id < seqDbr->getSize(); id++) {

        char *nuclSeq = seqDbr->getData(id, thread_idx);
        unsigned int seqLen = seqDbr->getSeqLens(id) - 2;

        if (seqLen >= par.maxSeqLen) {
            Debug(Debug::WARNING) << "Sequence too long: " << seqDbr->getDbKey(id) << ". It will be skipped. "
                                     "Max length allowed is " << par.maxSeqLen << "\n";
            continue;
        }
        seq.mapSequence(id, seqDbr->getDbKey(id), nuclSeq);

        //TODO: try spaced kmers?
        //TODO: limit the number of kmers in the first half of the sequence? only first 15%?

        /* extract front and back kmers */
        unsigned int frontKmersCount = 0;
        while (seq.hasNextKmer() && frontKmersCount < seqLen/2+1) {

            int *kmer = (int *) seq.nextKmer();

            uint64_t kmerIdx = indexer.int2index(kmer, 0, kmerSize);
            (frontKmers + frontKmersCount)->kmer = kmerIdx;
            (frontKmers + frontKmersCount)->pos = seq.getCurrentPosition();

            frontKmersCount++;
        }

        unsigned int backKmersCount = 0;
        while (seq.hasNextKmer()) {

            int *kmer = (int *) seq.nextKmer();

            uint64_t kmerIdx = indexer.int2index(kmer, 0, kmerSize);
            (backKmers + backKmersCount)->kmer = kmerIdx;
            (backKmers + backKmersCount)->pos = seq.getCurrentPosition();

            backKmersCount++;
        }

        std::sort(frontKmers, frontKmers + frontKmersCount, kmerSeqPos::compareByKmer);
        std::sort(backKmers, backKmers + backKmersCount, kmerSeqPos::compareByKmer);

        /* calculate front-back-kmermatches */
        unsigned int kmermatches = 0;
        std::fill(diagHits, diagHits + seqLen/2 +1 , 0);

        unsigned int idx = 0;
        unsigned int jdx = 0;
        while (idx < frontKmersCount && jdx < backKmersCount) {

            if (frontKmers[idx].kmer < backKmers[jdx].kmer)
                idx++;
            else if (frontKmers[idx].kmer > backKmers[jdx].kmer)
                jdx++;
            else {
                size_t kmerIdx= frontKmers[idx].kmer;
                unsigned int pos = frontKmers[idx].pos;
                while (jdx < backKmersCount && kmerIdx == backKmers[jdx].kmer) {

                    int diag = backKmers[jdx].pos - pos;
                    if (diag >= static_cast<int>(seqLen/2)) {
                        //diag = abs(diag);
                        diagHits[diag-seqLen/2]++;
                        kmermatches++;
                    }
                    jdx++;
                }
                while(idx < frontKmersCount && kmerIdx == frontKmers[idx].kmer) {
                    idx++;
                }
            }

        }
        //std:: cout << "number of kmermatches " << kmermatches << std::endl;
        /*for (size_t i=0; i < seqLen/2; i++) {
            if (diagHits[i] != 0)
                std:: cout << id << "\t" << seqLen << "\t"  << i+seqLen/2 << "\t" << diagHits[i] << std::endl;
        }*/

        /* calculate maximal hit rate on diagonal bands */
        int splitDiagonal = -1;
        float maxDiagbandHitRate = 0.0;

        if (kmermatches > 0) {
            for (unsigned int d = 0; d < seqLen / 2; d++) {
                if (diagHits[d] != 0) {
                    unsigned int diag = d + seqLen / 2;
                    unsigned int diaglen = seqLen - diag;
                    unsigned int gapwindow = diaglen * 0.01;
                    unsigned int lower = std::max(0, static_cast<int>(d - gapwindow));
                    unsigned int upper = std::min(d + gapwindow, seqLen / 2);
                    unsigned int diagbandHits = 0;

                    for (size_t i = lower; i <= upper; i++) {
                        diagbandHits += diagHits[i];
                    }

                    float diagbandHitRate = static_cast<float>(diagbandHits) / (diaglen - kmerSize + 1);
                    if (diagbandHitRate > maxDiagbandHitRate) {
                        maxDiagbandHitRate = diagbandHitRate;
                        splitDiagonal = diag;
                    }
                }

            }
        }

        if (maxDiagbandHitRate >= HIT_RATE_THRESHOLD) {
            //std::cout << "cyclic" << std::endl;

            unsigned int len = seqDbr->getSeqLens(id) - 1; //skip null byte
            std::string seq;
            if (par.chopCycle) {
                seq = std::string(nuclSeq, splitDiagonal);
                nuclSeq = (char *) seq.c_str();
                len = splitDiagonal;
            }

            cycleResultWriter.writeData(nuclSeq, len, seqDbr->getDbKey(id), 0);//thread_idx);
        }
        else {
            //std::cout << "linear" << std::endl;
            unsigned int len = seqDbr->getSeqLens(id) - 1; //skip null byte
            linearResultWriter.writeData(nuclSeq, len, seqDbr->getDbKey(id), 0);
        }
    }
    delete diagHits;
    delete frontKmers;
    delete backKmers;


    cycleResultWriter.close(true);
    linearResultWriter.close(true);


    //TODO: split main in functions
    //TODO: remap?
    //TODO: chop cycle
}
