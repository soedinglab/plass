//
// Created by Annika Seidel on 2019-08-09
//

//#include "checkcycle.h"

#include "DBReader.h"
#include "DBWriter.h"
#include "Indexer.h"
#include "LocalParameters.h"
#include "NucleotideMatrix.h"

#include <algorithm>

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


    int seqType  =  seqDbr->getDbtype();
    BaseMatrix *subMat;

    if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.nucleotides, 1.0, 0.0);
    }else {
        Debug(Debug::ERROR) << "Module checkcyle only supports nucleotide input database" << "\n";
        EXIT(EXIT_FAILURE);
    }

    std::cout << "subMat->alphabetSize - 1: " << subMat->alphabetSize - 1 << std::endl;
    Indexer indexer(subMat->alphabetSize - 1, kmerSize);

    //TODO: better solution than maxseqlen? at least warning
    Sequence seq(par.maxSeqLen, seqType, subMat, kmerSize, false, false);
    kmerSeqPos *frontKmers = new kmerSeqPos[par.maxSeqLen+1];
    kmerSeqPos *backKmers = new kmerSeqPos[par.maxSeqLen+1];

    for (size_t id = 0; id < seqDbr->getSize(); id++) {

        char *nuclSeq = seqDbr->getData(id, thread_idx);
        unsigned int seqLen = seqDbr->getSeqLens(id) - 2;
        seq.mapSequence(id, seqDbr->getDbKey(id), nuclSeq);

        //TODO: try spaced kmers?
        //TODO: limit the number of kmers in the first half of the sequence? only first 15%?
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

        unsigned int kmermatches = 0;
        //TODO: better solution than diag array?
        unsigned int *diagKmerCounter = new unsigned int[seqLen/2+1];
        std::fill(diagKmerCounter, diagKmerCounter + seqLen/2 +1 , 0);

        unsigned int idx = 0;
        unsigned int jdx = 0;
        while (idx < frontKmersCount && jdx < backKmersCount) {

            if (frontKmers[idx].kmer < backKmers[jdx].kmer)
                idx++;
            else if (frontKmers[idx].kmer > backKmers[jdx].kmer)
                jdx++;
            else {
                size_t kmerIdx= frontKmers[idx].kmer;
                unsigned int rangeStart = jdx;
                while (jdx < backKmersCount && kmerIdx == backKmers[jdx].kmer) {
                    jdx++;
                }
                unsigned int rangeEnd = jdx;

                unsigned int pos = frontKmers[idx].pos;
                for (size_t rangeIdx=rangeStart; rangeIdx < rangeEnd; rangeIdx++){
                    int diag = backKmers[rangeIdx].pos - pos;
                    if (diag >= static_cast<int>(seqLen/2)) {
                        //diag = abs(diag);
                        diagKmerCounter[diag-seqLen/2] ++;
                        kmermatches++;
                    }
                }
                while(idx < frontKmersCount && kmerIdx == frontKmers[idx].kmer) {
                    idx++;
                }
            }

        }

        //std:: cout << "number of kmermatches " << kmermatches << std::endl;
        for (size_t i=0; i < seqLen/2; i++) {
            if (diagKmerCounter[i] != 0)
                std:: cout << id << "\t" << seqLen << "\t"  << i+seqLen/2 << "\t" << diagKmerCounter[i] << std::endl;
        }
        //TODO: sort by diag to find diag with most matches

    }


    //TODO: remap?
}
