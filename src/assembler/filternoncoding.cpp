#include "ReducedMatrix.h"
#include "Indexer.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"

#include "Debug.h"
#include "AminoAcidLookupTables.h"
#include "DBReader.h"
#include "DBWriter.h"

#include "LocalParameters.h"

#include "kerasify/keras_model.h"
//#include "predict_coding_acc9540_57x32x64.model.h"
//#include "predict_coding_acc9642_57x32x64.model.h"
//#include "predict_coding_acc9598_57x32x64.model.h"
#include "predict_coding_acc9743_57x32x64.model.h"

//#include "predict_coding_acc9260_56x96.model.h"
//#include "predict_coding_acc9623_57x32x64.model.h"

#ifdef OPENMP
#include <omp.h>
#endif

int filternoncoding(int argc, const char **argv, const Command& command)  {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    Debug(Debug::INFO) << "Sequence database: " << par.db1 << "\n";
    DBReader<unsigned int> seqDb (par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    seqDb.open(DBReader<unsigned int>::NOSORT);

    Debug(Debug::INFO) << "Output file: " << par.db2 << "\n";
    DBWriter dbw(par.db2.c_str(), par.db2Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, seqDb.getDbtype());
    dbw.open();

    // Initialize model.
    KerasModel model;
//    model.LoadModel(std::string((const char *)predict_coding_acc9623_57x32x64_model, predict_coding_acc9623_57x32x64_model_len));
    model.LoadModel(std::string((const char *)predict_coding_acc9743_57x32x64_model, predict_coding_acc9743_57x32x64_model_len));

    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0, 0.0);
    ReducedMatrix redMat7(subMat.probMatrix, subMat.subMatrixPseudoCounts, subMat.aa2num, subMat.num2aa, subMat.alphabetSize, 7, subMat.getBitFactor());
//    ReducedMatrix redMat3(subMat.probMatrix, subMat.subMatrixPseudoCounts, subMat.aa2int, subMat.int2aa, subMat.alphabetSize, 3, subMat.getBitFactor());

    // Create a 1D Tensor on length 20 for input data.
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        Tensor in(56);
        float counter[255];
        std::fill(counter, counter + 255, 1.0);
        Sequence seq(par.maxSeqLen, seqDb.getDbtype(), &subMat,  par.kmerSize, false, false);
        Sequence rseq2mer(par.maxSeqLen, seqDb.getDbtype(), &redMat7, 2, false, false);
//        Sequence rseq5mer(par.maxSeqLen, Sequence::AMINO_ACIDS, &redMat3, 5, false, false);
        Indexer indexerDi(redMat7.alphabetSize, 2);
//        Indexer indexerPenta(redMat3.alphabetSize, 5);
        float *diAACnt = new float[redMat7.alphabetSize * redMat7.alphabetSize];
        std::fill(diAACnt, diAACnt + redMat7.alphabetSize * redMat7.alphabetSize, 1.0);
//        int pentaIdxRange = MathUtil::ipow<size_t>(redMat3.alphabetSize, 5);
//        float *pentaAACnt = new float[pentaIdxRange];
//        std::fill(pentaAACnt, pentaAACnt + pentaIdxRange, 1.0);
        Charges charge;
        Doolittle doolittle;

#pragma omp for schedule(static)
        for (size_t id = 0; id < seqDb.getSize(); id++) {
            std::vector<float> data;
            char *seqData = seqDb.getData(id, thread_idx);
            unsigned int dbKey = seqDb.getDbKey(id);
            unsigned int seqLen = seqDb.getSeqLen(id);
            seq.mapSequence(id, dbKey, seqData, seqLen);
//            printf("%5d ", seq.L);
//            float chargeFlt = Util::averageValueOnAminoAcids(charge.values, seqData);
//            float doolittleFlt = Util::averageValueOnAminoAcids(doolittle.values, seqData);
//            printf("%3f ", chargeFlt);
//            printf("%3f ", doolittleFlt);
            float totalAACnt = 0;
            data.push_back(static_cast<float>(seq.L));
            for (int pos = 0; pos < seq.L; pos++) {
                if (seq.numSequence[pos] < subMat.alphabetSize - 1) {
                    counter[seq.numSequence[pos]] += 1.0;
                    totalAACnt += 1.0;
                }
            }
            for (int aa = 0; aa < subMat.alphabetSize - 1; aa++) {
                data.push_back(counter[aa] / (totalAACnt + subMat.alphabetSize - 1));
                counter[aa] = 1.0;
//                printf("%.4f ", data.back());
            }
            // di matrix
            {
                rseq2mer.mapSequence(id, dbKey, seqData, seqLen);
                float totalDiAACnt = 0;
                while (rseq2mer.hasNextKmer()) {
                    const unsigned char *kmer = rseq2mer.nextKmer();
                    // ignore x
                    if (static_cast<int>(kmer[0]) == redMat7.alphabetSize - 1 ||
                        static_cast<int>(kmer[1]) == redMat7.alphabetSize - 1) {
                        continue;
                    }
                    size_t index = indexerDi.int2index(kmer);
                    diAACnt[index] += 1.0;
                    totalDiAACnt += 1.0;
                }
                size_t kmer[2];
                float diRealRange = static_cast<float>((redMat7.alphabetSize - 1) * (redMat7.alphabetSize - 1));
                for (int raa = 0; raa < (redMat7.alphabetSize * redMat7.alphabetSize); raa++) {
                    indexerDi.index2int(kmer, raa, 2);
                    if (static_cast<int>(kmer[0]) == redMat7.alphabetSize - 1 ||
                        static_cast<int>(kmer[1]) == redMat7.alphabetSize - 1) {
                        continue;
                    }
                    data.push_back(diAACnt[raa] / (totalDiAACnt + diRealRange));
//                    printf("%.4f ", data.back());
                    diAACnt[raa] = 1.0;
                }
            }
            // 5mer
//            {
//                rseq5mer.mapSequence(id, dbKey, seqData);
//                float totalPentaAACnt = 0;
//                while (rseq5mer.hasNextKmer()) {
//                    const int *kmer = rseq5mer.nextKmer();
//                    // ignore x
//                    if (kmer[0] == redMat3.alphabetSize - 1 ||
//                        kmer[1] == redMat3.alphabetSize - 1 ||
//                        kmer[2] == redMat3.alphabetSize - 1 ||
//                        kmer[3] == redMat3.alphabetSize - 1 ||
//                        kmer[4] == redMat3.alphabetSize - 1) {
//                        continue;
//                    }
//                    size_t index = indexerPenta.int2index(kmer);
//                    pentaAACnt[index] += 1.0;
//                    totalPentaAACnt += 1.0;
//                }
//                size_t kmer5[5];
//                int pentaRealRange = static_cast<float>(MathUtil::ipow<size_t>(redMat3.alphabetSize - 1, 5));
//                for (int raa = 0; raa < pentaIdxRange; raa++) {
//                    indexerPenta.index2int(kmer5, raa, 5);
//                    if (static_cast<int>(kmer5[0]) == redMat3.alphabetSize - 1 ||
//                        static_cast<int>(kmer5[1]) == redMat3.alphabetSize - 1 ||
//                        static_cast<int>(kmer5[2]) == redMat3.alphabetSize - 1 ||
//                        static_cast<int>(kmer5[3]) == redMat3.alphabetSize - 1 ||
//                        static_cast<int>(kmer5[4]) == redMat3.alphabetSize - 1) {
//                        continue;
//                    }
//                    data.push_back(pentaAACnt[raa] / (totalPentaAACnt + pentaRealRange));
//                    printf("%.4f ", data.back());
//                    pentaAACnt[raa] = 1.0;
//                }
//            }
//            printf("\n");
            //printf("%d\n", data.size());
            in.data_ = data;
            // Run prediction.
            Tensor out;
            model.Apply(&in, &out);
            if (out.data_[0] > par.proteinFilterThreshold) {
                // -1 dont write \0 byte
                dbw.writeData(seqData, seqDb.getEntryLen(id)-1, dbKey, thread_idx);
            } else {
                dbw.writeData("\n",  1, dbKey, thread_idx);
            }
//        out.Print();
        }

        delete[] diAACnt;
//        delete[] pentaAACnt;
    }
//    std::cout << "Filtered: " << static_cast<float>(cnt)/ static_cast<float>(seqDb.getSize()) << std::endl;
    dbw.close(true);
    seqDb.close();

    return EXIT_SUCCESS;
}

