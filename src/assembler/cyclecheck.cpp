/*
 * Written by Annika Jochheim <annika.jochheim@mpibpc.mpg.de>
 * detect circular fragments
 *
 * constraint: fragments contain redundant parts at most 3 times
 */

#include "DBReader.h"
#include "DBWriter.h"
#include "Indexer.h"
#include "LocalParameters.h"
#include "NucleotideMatrix.h"

#include <algorithm>
#ifdef OPENMP
#include <omp.h>
#endif


#define HIT_RATE_THRESHOLD 0.24
// threshold to distinguish cyclic/terminal redundant genomes from random hits on linear genomes
// chosen based on analysis on ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ (last modified 7/11/19)
// verified for kmerSize = 22

void setCycleCheckDefaults(LocalParameters *p) {
    p->kmerSize = 22;
    p->chopCycle = false;
}

int cyclecheck(int argc, const char **argv, const Command& command) {

    LocalParameters &par = LocalParameters::getLocalInstance();
    setCycleCheckDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> *seqDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str(),  par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    seqDbr->open(DBReader<unsigned int>::NOSORT);

    DBWriter cycleResultWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_NUCLEOTIDES);
    cycleResultWriter.open();

    const size_t kmerSize = par.kmerSize;
    int seqType  =  seqDbr->getDbtype();
    BaseMatrix *subMat;

    if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.nucleotides, 1.0, 0.0);
    }else {
        Debug(Debug::ERROR) << "Module cyclecheck only supports nucleotide input database" << "\n";
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

    Debug::Progress progress(seqDbr->getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        /* 1. split sequence in 3 parts and extract kmers from each part (frontKmers, middleKmers and backKmers)
         * 2. find kmermatches between 1 third and 2 third, 1 third and 3 third, 2 third and 3 third
         * 3. sum up kmermatches for diagonalbands
         * 4. find longest diagonal (= smallest diag index) which fullfill threeshold
         */

        Indexer indexer(subMat->alphabetSize - 1, kmerSize);
        Sequence seq(par.maxSeqLen, seqType, subMat, kmerSize, false, false);
        kmerSeqPos *frontKmers = new(std::nothrow) kmerSeqPos[par.maxSeqLen / 3 + 1];
        Util::checkAllocation(frontKmers, "Can not allocate memory");
        kmerSeqPos *middleKmers = new(std::nothrow) kmerSeqPos[par.maxSeqLen / 3 + 1];
        Util::checkAllocation(middleKmers, "Can not allocate memory");
        kmerSeqPos *backKmers = new(std::nothrow) kmerSeqPos[par.maxSeqLen / 3 + 1];
        Util::checkAllocation(backKmers, "Can not allocate memory");
        unsigned int *diagHits = new unsigned int[(par.maxSeqLen / 3) * 2 + 1];

#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < seqDbr->getSize(); id++) {
            progress.updateProgress();

            char *nuclSeq = seqDbr->getData(id, thread_idx);
            unsigned int seqLen = seqDbr->getSeqLen(id);

            if (seqLen >= par.maxSeqLen) {
                Debug(Debug::WARNING) << "Sequence " << seqDbr->getDbKey(id) << " too long. It will be skipped. "
                                                                                          "Max length allowed is "
                                      << par.maxSeqLen << "\n";
                continue;
            }
            seq.mapSequence(id, seqDbr->getDbKey(id), nuclSeq, seqLen);

            //TODO: try spaced kmers?
            //TODO: limit the number of kmers in the first half of the sequence? only first 15%?

            /* extract kmers */
            unsigned int frontKmersCount = 0, middleKmersCount = 0, backKmersCount = 0;
            unsigned int thirdSeqLen = seqLen / 3;
            while (seq.hasNextKmer()) {

                unsigned int pos =  seq.getCurrentPosition();
                const unsigned char *kmer = seq.nextKmer();
                uint64_t kmerIdx = indexer.int2index(kmer, 0, kmerSize);

                if (pos  < thirdSeqLen + 1) {
                    (frontKmers + frontKmersCount)->kmer = kmerIdx;
                    (frontKmers + frontKmersCount)->pos = seq.getCurrentPosition();
                    frontKmersCount++;
                }
                else if (pos < 2 * thirdSeqLen + 1) {
                    (middleKmers + middleKmersCount)->kmer = kmerIdx;
                    (middleKmers + middleKmersCount)->pos = seq.getCurrentPosition();
                    middleKmersCount++;
                }
                else // pos > 2* thirdSeqLen
                {
                    (backKmers + backKmersCount)->kmer = kmerIdx;
                    (backKmers + backKmersCount)->pos = seq.getCurrentPosition();
                    backKmersCount++;
                }
            }

            std::sort(frontKmers, frontKmers + frontKmersCount, kmerSeqPos::compareByKmer);
            std::sort(middleKmers, middleKmers + middleKmersCount, kmerSeqPos::compareByKmer);
            std::sort(backKmers, backKmers + backKmersCount, kmerSeqPos::compareByKmer);

            /* calculate front-back-kmermatches and front-middle-kmermatches */
            unsigned int kmermatches = 0;
            std::fill(diagHits, diagHits + 2*thirdSeqLen + 1, 0);

            unsigned int idx = 0;
            unsigned int jdx = 0;
            unsigned int kdx = 0;

            while(idx < frontKmersCount && (jdx < backKmersCount || kdx < middleKmersCount) ) {

                size_t kmerIdx = frontKmers[idx].kmer;
                unsigned int pos = frontKmers[idx].pos;

                while (jdx < backKmersCount && backKmers[jdx].kmer < kmerIdx) {
                    jdx++;
                }
                while (kdx < middleKmersCount && middleKmers[kdx].kmer < kmerIdx) {
                    kdx++;
                }

                while (jdx < backKmersCount && kmerIdx == backKmers[jdx].kmer) {
                    int diag = backKmers[jdx].pos - pos;
                    if (diag >= static_cast<int>(seqLen / 3)) {
                        diagHits[diag - seqLen / 3]++;
                        kmermatches++;
                    }
                    jdx++;
                }

                while (kdx < middleKmersCount && kmerIdx == middleKmers[kdx].kmer) {
                    int diag = middleKmers[kdx].pos - pos;
                    if (diag >= static_cast<int>(seqLen / 3)) {
                        diagHits[diag - seqLen / 3]++;
                        kmermatches++;
                    }
                    kdx++;
                }

                idx++;
                while (idx < frontKmersCount && kmerIdx == frontKmers[idx].kmer) {
                    idx++;
                }
            }

            /* calculate middle-back-kmermatches */
            jdx = 0, kdx = 0;
            while (kdx < middleKmersCount && jdx < backKmersCount) {

                if (middleKmers[kdx].kmer < backKmers[jdx].kmer)
                    kdx++;
                else if (middleKmers[kdx].kmer > backKmers[jdx].kmer)
                    jdx++;
                else {
                    size_t kmerIdx = middleKmers[kdx].kmer;
                    unsigned int pos = middleKmers[kdx].pos;
                    while (jdx < backKmersCount && kmerIdx == backKmers[jdx].kmer) {

                        int diag = backKmers[jdx].pos - pos;
                        if (diag >= static_cast<int>(seqLen / 3)) {
                            //diag = abs(diag);
                            diagHits[diag - seqLen / 3]++;
                            kmermatches++;
                        }
                        jdx++;
                    }
                    while (kdx < middleKmersCount && kmerIdx == middleKmers[kdx].kmer) {
                        kdx++;
                    }
                }

            }

            /* calculate hit rate on diagonal bands */
            unsigned int splitDiagonal = 0;

            if (kmermatches > 0) {
                for (unsigned int d = 0; d < 2 * thirdSeqLen; d++) {
                    if (diagHits[d] != 0) {
                        unsigned int diag = d + thirdSeqLen;
                        unsigned int diaglen = seqLen - diag;
                        unsigned int gapwindow = diaglen * 0.01;
                        unsigned int lower = std::max(0, static_cast<int>(d - gapwindow));
                        unsigned int upper = std::min(d + gapwindow, 2 * thirdSeqLen);
                        unsigned int diagbandHits = 0;

                        for (size_t i = lower; i <= upper; i++) {
                            if (diagHits[i] <= diagHits[d])
                                diagbandHits += diagHits[i];
                        }

                        float diagbandHitRate = static_cast<float>(diagbandHits) / (diaglen - kmerSize + 1);
                        if (diagbandHitRate > HIT_RATE_THRESHOLD) {
                            splitDiagonal = diag;
                            break;
                        }
                    }

                }
            }

            if (splitDiagonal != 0) {

                unsigned int len = seqDbr->getEntryLen(id)-1;
                std::string seq;
                if (par.chopCycle) {
                    seq = std::string(nuclSeq, splitDiagonal);
                    seq.push_back('\n');
                    nuclSeq = (char *) seq.c_str();
                    len = seq.size();
                }
                cycleResultWriter.writeData(nuclSeq, len, seqDbr->getDbKey(id), thread_idx);

            }
        }
        delete[] diagHits;
        delete[] frontKmers;
        delete[] middleKmers;
        delete[] backKmers;
    }

    cycleResultWriter.close(true);

    seqDbr->close();
    delete seqDbr;

    return EXIT_SUCCESS;

}
