#include "SubstitutionMatrix.h"
#include "MultipleAlignment.h"

#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"

#ifdef OPENMP
#include <omp.h>
#endif

int findPosOfM(const char *seq) {
    int stopPos = -1;
    int pos = 0;
    while(seq[pos] != '\0'){
        if (seq[pos] == 'M') {
            stopPos = pos;
            break;
        }
        pos++;
    }
    return stopPos;
}

int findassemblystart(int argn, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argn, argv, command, true, 0, 0);

    DBReader<unsigned int> qDbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    qDbr.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> *tDbr = &qDbr;

    DBReader<unsigned int> resultReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter resultWriter(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
    resultWriter.open();

    // + 1 for query
    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0f, 0.0f);
    Debug(Debug::INFO) << "Start computing start in sequences.\n";
    int *addStopAtPosition = new int[qDbr.getSize()];
    std::fill(addStopAtPosition, addStopAtPosition + qDbr.getSize(), -1);
    Debug::Progress progress(resultReader.getSize());
    const float threshold = 0.2;

#pragma omp parallel
    {
    unsigned int thread_idx = 0;
#ifdef OPENMP
    thread_idx = (unsigned int) omp_get_thread_num();
#endif

#pragma omp for schedule(dynamic, 100)
    for (size_t id = 0; id < resultReader.getSize(); id++) {
        progress.updateProgress();
        // Get the sequence from the queryDB
        unsigned int queryKey = resultReader.getDbKey(id);
        const size_t qId = tDbr->getId(queryKey);
        char *querySeqData = tDbr->getData(qId, thread_idx);
        int queryPosOfM = findPosOfM(querySeqData);
        if (queryPosOfM == -1){
            continue;
        }
        bool hasStopM = false;
        if (queryPosOfM > 0){
            hasStopM = querySeqData[queryPosOfM-1] == '*';
        }

        struct PositionOfM {
            unsigned int id; int mPos; bool hasM; bool hasStopM;
            PositionOfM(unsigned int id, int mPos, bool hasM, bool hasStopM)
                    : id(id), mPos(mPos), hasM(hasM), hasStopM(hasStopM) {}
        };
        std::vector<PositionOfM> stopPositions;
        stopPositions.emplace_back(qId,queryPosOfM, true, hasStopM);

        char *results = resultReader.getData(id, thread_idx);
        while (*results != '\0') {
            char dbKey[255 + 1];
            Util::parseKey(results, dbKey);
            const unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
            const char *entry[255];
            const size_t columns = Util::getWordsOfLine(results, entry, 255);
            Matcher::result_t res;
            if (columns >= Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
                res = Matcher::parseAlignmentRecord(results);
            } else {
                Debug(Debug::ERROR) << "ERROR: Backtrace is missing for at result: " << id  << "\n";
                EXIT(EXIT_FAILURE);
            }

            const size_t edgeId = tDbr->getId(key);
            if (edgeId == qId){
                results = Util::skipLine(results);
                continue;
            }
            char *dbSeqData = tDbr->getData(edgeId, thread_idx);
            int posOfM = -1;
            bool hasM = false;
            bool hasStopM = false;
            if (res.qStartPos >= queryPosOfM  && queryPosOfM <= res.qEndPos){
                int queryMoffset = queryPosOfM - res.qStartPos;
                int dbMPos = res.dbStartPos + queryMoffset;
                posOfM = dbMPos;
                hasM = dbMPos >= 0 && (dbSeqData[dbMPos] == 'M');
                if (hasM == false){
                    goto endOfLoop;
                }
                if (dbMPos > 0 && hasM){
                    hasStopM = dbSeqData[dbMPos - 1] == '*';
                }
            } else {
                goto endOfLoop;
            }
            stopPositions.emplace_back(edgeId,posOfM, hasM, hasStopM);
            endOfLoop:
            results = Util::skipLine(results);
        }
        int stopMCount = 0;
        int mCount = 0;
        for (size_t seqIdx = 0; seqIdx < stopPositions.size(); seqIdx++){
            stopMCount += stopPositions[seqIdx].hasStopM;
            mCount += stopPositions[seqIdx].hasM;
        }
        if (stopPositions.size() > 1){
            const float frequency = static_cast<float>(stopMCount) / static_cast<float>(stopPositions.size());
            if (frequency >= threshold){
                for (size_t seqIdx = 0; seqIdx < stopPositions.size(); seqIdx++){
                    int target;

                    int curVal = stopPositions[seqIdx].mPos;
                    __atomic_load(&addStopAtPosition[stopPositions[seqIdx].id], &target ,__ATOMIC_RELAXED);
                    do {
                        if (target >= curVal) break;
                    } while (!__atomic_compare_exchange(&addStopAtPosition[stopPositions[seqIdx].id],  &target,  &curVal , false,  __ATOMIC_RELAXED, __ATOMIC_RELAXED));
                }
            }
        }
    }

    }
#pragma omp parallel
    {
        std::string str;
        str.reserve(100000);

        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

#pragma omp for schedule(dynamic, 100)
        for(size_t id = 0; id < qDbr.getSize(); id++){
            unsigned int queryKey = qDbr.getDbKey(id);
            char *querySeqData = tDbr->getData(id, thread_idx);
            int mPos = addStopAtPosition[id];
            if (mPos == -1){
                resultWriter.writeData(querySeqData, strlen(querySeqData), queryKey, thread_idx);
            } else {
                str.append("*");
                str.append(querySeqData + mPos);
                resultWriter.writeData(str.c_str(), str.length(), queryKey, thread_idx);
                str.clear();
            }
        }
    }

    // cleanup
    resultWriter.close();
    resultReader.close();
    qDbr.close();
    delete [] addStopAtPosition;

    Debug(Debug::INFO) << "\nDone.\n";

    return EXIT_SUCCESS;
}



