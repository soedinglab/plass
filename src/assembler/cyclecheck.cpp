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

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        Indexer indexer(subMat->alphabetSize - 1, kmerSize);
        Sequence seq(par.maxSeqLen, seqType, subMat, kmerSize, false, false);
        kmerSeqPos *frontKmers = new(std::nothrow) kmerSeqPos[par.maxSeqLen / 2 + 1];
        Util::checkAllocation(frontKmers, "Can not allocate memory");
        kmerSeqPos *backKmers = new(std::nothrow) kmerSeqPos[par.maxSeqLen / 2 + 1];
        Util::checkAllocation(backKmers, "Can not allocate memory");
        unsigned int *diagHits = new unsigned int[par.maxSeqLen / 2 + 1];

#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < seqDbr->getSize(); id++) {

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

            /* extract front and back kmers */
            unsigned int frontKmersCount = 0;
            while (seq.hasNextKmer() && frontKmersCount < seqLen / 2 + 1) {

                const unsigned char *kmer = seq.nextKmer();

                uint64_t kmerIdx = indexer.int2index(kmer, 0, kmerSize);
                (frontKmers + frontKmersCount)->kmer = kmerIdx;
                (frontKmers + frontKmersCount)->pos = seq.getCurrentPosition();

                frontKmersCount++;
            }

            unsigned int backKmersCount = 0;
            while (seq.hasNextKmer()) {

                const unsigned char *kmer = seq.nextKmer();

                uint64_t kmerIdx = indexer.int2index(kmer, 0, kmerSize);
                (backKmers + backKmersCount)->kmer = kmerIdx;
                (backKmers + backKmersCount)->pos = seq.getCurrentPosition();

                backKmersCount++;
            }

            std::sort(frontKmers, frontKmers + frontKmersCount, kmerSeqPos::compareByKmer);
            std::sort(backKmers, backKmers + backKmersCount, kmerSeqPos::compareByKmer);

            /* calculate front-back-kmermatches */
            unsigned int kmermatches = 0;
            std::fill(diagHits, diagHits + seqLen / 2 + 1, 0);

            unsigned int idx = 0;
            unsigned int jdx = 0;
            while (idx < frontKmersCount && jdx < backKmersCount) {

                if (frontKmers[idx].kmer < backKmers[jdx].kmer)
                    idx++;
                else if (frontKmers[idx].kmer > backKmers[jdx].kmer)
                    jdx++;
                else {
                    size_t kmerIdx = frontKmers[idx].kmer;
                    unsigned int pos = frontKmers[idx].pos;
                    while (jdx < backKmersCount && kmerIdx == backKmers[jdx].kmer) {

                        int diag = backKmers[jdx].pos - pos;
                        if (diag >= static_cast<int>(seqLen / 2)) {
                            //diag = abs(diag);
                            diagHits[diag - seqLen / 2]++;
                            kmermatches++;
                        }
                        jdx++;
                    }
                    while (idx < frontKmersCount && kmerIdx == frontKmers[idx].kmer) {
                        idx++;
                    }
                }

            }

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
        delete[] backKmers;
    }

    cycleResultWriter.close(true);

    seqDbr->close();
    delete seqDbr;

    return EXIT_SUCCESS;

}
