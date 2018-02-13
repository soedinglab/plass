#include "ReducedMatrix.h"
#include "Indexer.h"
#include "Parameters.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "DBWriter.h"
#include "Debug.h"
#include "DBReader.h"
#include "keras_model.h"
#include "predict_coding_acc9260_56x96.model.h"
#ifdef OPENMP
#include <omp.h>
#endif

int filternonecoding(int argc, const char **argv, const Command& command)  {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);

    Debug(Debug::INFO) << "Sequence database: " << par.db1 << "\n";
    DBReader<unsigned int> seqDb (par.db1.c_str(), par.db1Index.c_str());
    seqDb.open(DBReader<unsigned int>::NOSORT);

    Debug(Debug::INFO) << "Output  file: " << par.db2 << "\n";

    DBWriter dbw(par.db2.c_str(), par.db2Index.c_str(), static_cast<unsigned int>(par.threads));
    dbw.open();

    // Initialize model.
    KerasModel model;
    std::string submat((const char *)predict_coding_acc9260_56x96_model, predict_coding_acc9260_56x96_model_len);

    model.LoadModel(submat);
    SubstitutionMatrix subMat("blosum62.out", 2.0, 0.0);
    ReducedMatrix redMat(subMat.probMatrix, subMat.subMatrixPseudoCounts, 7, subMat.getBitFactor());

    // Create a 1D Tensor on length 20 for input data.
#pragma omp parallel
    {
        Tensor in(56);
        float counter[255];
        std::fill(counter, counter + 255, 1.0);
        Sequence seq(par.maxSeqLen, Sequence::AMINO_ACIDS, &subMat,  par.kmerSize, false, false);
        Sequence rseq(par.maxSeqLen, Sequence::AMINO_ACIDS, &redMat, 2, false, false);
        Indexer indexer(redMat.alphabetSize, 2);
        float diAACnt[redMat.alphabetSize * redMat.alphabetSize];
        std::fill(diAACnt, diAACnt + redMat.alphabetSize * redMat.alphabetSize, 1.0);
#pragma omp for schedule(static)
        for (size_t id = 0; id < seqDb.getSize(); id++) {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            std::vector<float> data;
            char * seqData = seqDb.getData(id);
            unsigned int dbKey = seqDb.getDbKey(id);
            seq.mapSequence(id, dbKey, seqData);

            //printf("%5d ", seq.L);
            float totalAACnt = 0;
            for (size_t pos = 0; pos < seq.L; pos++) {
                if (seq.int_sequence[pos] < subMat.alphabetSize - 1) {
                    counter[seq.int_sequence[pos]] += 1.0;
                    totalAACnt += 1.0;
                }
            }
            for (size_t aa = 0; aa < subMat.alphabetSize - 1; aa++) {
                data.push_back(counter[aa] / (totalAACnt + subMat.alphabetSize - 1));
                counter[aa] = 1.0;
                //printf("%.4f ", data.back());
            }

            rseq.mapSequence(id, dbKey, seqData);
            float totalDiAACnt = 0;
            while (rseq.hasNextKmer()) {
                const int *kmer = rseq.nextKmer();
                // ignore x
                if (kmer[0] == redMat.alphabetSize - 1 || kmer[1] == redMat.alphabetSize - 1) {
                    continue;
                }
                size_t index = indexer.int2index(kmer);
                diAACnt[index] += 1.0;
                totalDiAACnt += 1.0;
            }
            size_t kmer[2];
            for (size_t raa = 0; raa < (redMat.alphabetSize * redMat.alphabetSize); raa++) {
                indexer.index2int(kmer, raa, 2);
                if (kmer[0] == redMat.alphabetSize - 1 || kmer[1] == redMat.alphabetSize - 1) {
                    continue;
                }
                data.push_back(diAACnt[raa] / (totalDiAACnt + (redMat.alphabetSize - 1) * (redMat.alphabetSize - 1)));
                //printf("%.4f ", data.back());
                diAACnt[raa] = 1.0;
            }
            //printf("\n");
            //printf("%d\n", data.size());
            in.data_ = data;
            // Run prediction.
            Tensor out;
            model.Apply(&in, &out);
            if (out.data_[0] > 0.2) {
                // -1 dont write \0 byte
                dbw.writeData(seqData, seqDb.getSeqLens(id) - 1, dbKey, thread_idx);
            }else{
                dbw.writeData("\n",  1, dbKey, thread_idx);
            }
//        out.Print();
        }
    }
//    std::cout << "Filtered: " << static_cast<float>(cnt)/ static_cast<float>(seqDb.getSize()) << std::endl;
    seqDb.close();
    dbw.close(DBReader<unsigned int>::DBTYPE_AA);
    return 0;
}

