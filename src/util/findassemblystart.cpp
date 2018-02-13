/*
 * findassemblystart
 * written by Martin Steinegger <martin.steinegger@mpibpc.mpg.de>.
 */

#include <cstdio>
#include <map>
#include <fstream>
#include <unistd.h>
#include <math.h>

#ifdef OPENMP
#include <omp.h>
#endif

#include "SubstitutionMatrix.h"
#include "MultipleAlignment.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Parameters.h"
#include "Util.h"

int findPosOfM(char *seq) {
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
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argn, argv, command, 3, true, true);

    DBReader<unsigned int> qDbr(par.db1.c_str(), par.db1Index.c_str());
    qDbr.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> *tDbr = &qDbr;

    DBReader<unsigned int> resultReader(par.db2.c_str(), par.db2Index.c_str());
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    const float threshold = 0.2;
    DBWriter resultWriter(par.db3.c_str(), par.db3Index.c_str(), par.threads);
    resultWriter.open();

    // + 1 for query
    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0f, 0.0f);
    Debug(Debug::INFO) << "Start computing start in sequences.\n";
    int * addStopAtPosition = new int[qDbr.getSize()];
    std::fill(addStopAtPosition, addStopAtPosition + qDbr.getSize(), -1);

#pragma omp parallel for schedule(dynamic, 100)
    for (size_t id = 0; id < resultReader.getSize(); id++) {
        Debug::printProgress(id);
        // Get the sequence from the queryDB
        unsigned int queryKey = resultReader.getDbKey(id);
        const size_t qId = tDbr->getId(queryKey);
        char *querySeqData = tDbr->getData(qId);
        int queryPosOfM = findPosOfM(querySeqData);
        if(queryPosOfM==-1){
            continue;
        }
        bool hasStopM = false;
        if(queryPosOfM > 0){
            hasStopM = querySeqData[queryPosOfM-1] == '*';
        }

        struct PositionOfM {
            unsigned int id; int mPos; bool hasM; bool hasStopM;
            PositionOfM(unsigned int id, int mPos, bool hasM, bool hasStopM)
                    : id(id), mPos(mPos), hasM(hasM), hasStopM(hasStopM) {}
        };
        std::vector<PositionOfM> stopPositions;
        stopPositions.push_back(PositionOfM(qId,queryPosOfM, true, hasStopM));

        char *results = resultReader.getData(id);
        while (*results != '\0') {
            char dbKey[255 + 1];
            Util::parseKey(results, dbKey);
            const unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
            char *entry[255];
            const size_t columns = Util::getWordsOfLine(results, entry, 255);
            Matcher::result_t res;
            if (columns >= Matcher::ALN_RES_WITH_OUT_BT_COL_CNT) {
                res = Matcher::parseAlignmentRecord(results);
            }else{
                Debug(Debug::ERROR) << "ERROR: Backtrace is missing for at result: " << id  << "\n";
                EXIT(EXIT_FAILURE);
            }

            const size_t edgeId = tDbr->getId(key);
            if(edgeId == qId){
                results = Util::skipLine(results);
                continue;
            }
            char *dbSeqData = tDbr->getData(edgeId);
            int posOfM = -1;
            bool hasM = false;
            bool hasStopM = false;
            if(res.qStartPos >= queryPosOfM  && queryPosOfM <= res.qEndPos){
                int queryMoffset = queryPosOfM - res.qStartPos;
                int dbMPos = res.dbStartPos + queryMoffset;
                posOfM = dbMPos;
                hasM = (dbSeqData[dbMPos] == 'M');
                if(hasM == false){
                    goto endOfLoop;
                }
                if(dbMPos > 0 && hasM){
                    hasStopM = dbSeqData[dbMPos - 1] == '*';
                }
            }else{
                goto endOfLoop;
            }
            stopPositions.push_back(PositionOfM(edgeId,posOfM, hasM, hasStopM));
            endOfLoop:
            results = Util::skipLine(results);
        }
        int stopMCount = 0;
        int mCount = 0;
        for(size_t seqIdx = 0; seqIdx < stopPositions.size(); seqIdx++){
            stopMCount += stopPositions[seqIdx].hasStopM;
            mCount += stopPositions[seqIdx].hasM;
        }
        if(stopPositions.size() > 1 ){
            float frequence = static_cast<float>(stopMCount) / static_cast<float>(stopPositions.size());
            if(frequence >= threshold){
                for(size_t seqIdx = 0; seqIdx < stopPositions.size(); seqIdx++){
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

#pragma omp parallel
    {
        std::string str;
        str.reserve(100000);
#pragma  omp for schedule(dynamic, 100)
        for(size_t id = 0; id < qDbr.getSize(); id++){
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif
            unsigned int queryKey = qDbr.getDbKey(id);
            char *querySeqData = tDbr->getData(id);
            int mPos = addStopAtPosition[id];
            if(mPos == -1){
                resultWriter.writeData(querySeqData, strlen(querySeqData), queryKey, thread_idx);
            }else{
                str.append("*");
                str.append(querySeqData + mPos);
                resultWriter.writeData(str.c_str(), str.length(), queryKey, thread_idx);
                str.clear();
            }
        }
    }
    // cleanup
    delete [] addStopAtPosition;
    resultWriter.close(DBReader<unsigned int>::DBTYPE_AA);
    resultReader.close();
    qDbr.close();

    Debug(Debug::INFO) << "\nDone.\n";

    return EXIT_SUCCESS;
}



