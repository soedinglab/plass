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

    Indexer indexer(subMat->alphabetSize - 1, kmerSize);

    //TODO: better solution than maxseqlen? at least warning
    Sequence seq(par.maxSeqLen, seqType, subMat, kmerSize, false, false);
    kmerSeqPos * kmers = new kmerSeqPos[par.maxSeqLen+1];
    for (size_t id = 0; id < seqDbr->getSize(); id++) {

        char *nuclSeq = seqDbr->getData(id, thread_idx);
        unsigned int seqLen = seqDbr->getSeqLens(id) - 2;
        std::cout<< "seqLen: " << seqLen << std::endl;

        seq.mapSequence(id, seqDbr->getDbKey(id), nuclSeq);

        unsigned int seqKmerCount = 0;
        while (seq.hasNextKmer()) { // && seq.getCurrentPosition() < seqLen/2+1) {

            int *kmer = (int *) seq.nextKmer();

            uint64_t kmerIdx = indexer.int2index(kmer, 0, kmerSize);
            (kmers + seqKmerCount)->kmer = kmerIdx;
            (kmers + seqKmerCount)->pos = seq.getCurrentPosition();
            seqKmerCount++;
        }

        std:: cout << "number of extracted kmers " << seqKmerCount << std::endl;
        std::sort(kmers, kmers + seqKmerCount, kmerSeqPos::compareByKmer);

        unsigned int kmermatches = 0;

        size_t kmerIdx = kmers[0].kmer;
        unsigned int pos = kmers[0].pos;
        int diag;
        //TODO: better solution than diag array?
        unsigned int *diagKmerCounter = new unsigned int[seqLen/2+1];
        std::fill(diagKmerCounter, diagKmerCounter + seqLen/2 +1 , 0);
        for (size_t jdx = 1; jdx < seqKmerCount; jdx++) {

            //TODO: use hashSeqPair?

            size_t currKmerIdx = kmers[jdx].kmer;

            if (kmerIdx != currKmerIdx) {
                kmerIdx = currKmerIdx;
                pos = kmers[jdx].pos;
            }
            else {
                diag = kmers[jdx].pos - pos;


                if (diag > static_cast<int>(seqLen/2)) {
                    //diag = abs(diag);
                    diagKmerCounter[diag-seqLen/2] ++;
                    kmermatches++;
                }
            }


        }

        std:: cout << "number of kmermatches " << kmermatches << std::endl;

        for (size_t i=0; i < seqLen/2; i++) {
            std:: cout << i << ": " << diagKmerCounter[i] << std::endl;
        }

        //TODO: sort by diag to find diag with most matches

    }

    //TODO: use two kmer arrays to limit the number of kmers in the first half of the sequence
    //TODO: remap?
}
