#include "LocalParameters.h"
#include "DistanceCalculator.h"
#include "Matcher.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "MathUtil.h"

#include <limits>
#include <cstdint>
#include <queue>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

class CompareResultBySeqId {
public:
    bool operator() (const Matcher::result_t & r1,const Matcher::result_t & r2) {
//        if(r1.seqId < r2.seqId )
//            return true;
//        if(r2.seqId < r1.seqId )
//            return false;
        int score1 = BIT_CLEAR(r1.score, 31);
        int score2 = BIT_CLEAR(r2.score, 31);
        if(score1 < score2 )
            return true;
        if(score2 < score1 )
            return false;
        if(r1.dbKey > r2.dbKey )
            return true;
        if(r2.dbKey > r1.dbKey )
            return false;
        /*  int seqLen1 = r1.qEndPos - r1.qStartPos;
          int seqLen2 = r2.qEndPos - r2.qStartPos;
          if(seqLen1 < seqLen2)
              return true;
          if(seqLen2 < seqLen1 )
              return false;*/
        return false;
    }
};

typedef std::priority_queue<Matcher::result_t, std::vector<Matcher::result_t> , CompareResultBySeqId> QueueBySeqId;
Matcher::result_t selectFragmentToExtend(QueueBySeqId &alignments,
                                             unsigned int queryKey) {
    // results are ordered by score
    while (alignments.empty() == false){
        Matcher::result_t res = alignments.top();
        alignments.pop();
        size_t dbKey = res.dbKey;
        const bool isReverse = res.score<0;
        const bool notRightStartAndLeftStart = !(res.dbStartPos == 0 &&  res.qStartPos == 0 );
        const bool rightStart = res.dbStartPos == 0 && (res.dbEndPos != res.dbLen-1);
        const bool leftStart = res.qStartPos == 0   && (res.qEndPos != res.qLen-1);
        const bool isNotIdentity = (dbKey != queryKey);

        if ((rightStart || leftStart) && notRightStartAndLeftStart && isNotIdentity){
            return res;
        }
    }
    return Matcher::result_t(UINT_MAX,0,0,0,0,0,0,0,0,0,0,0,0,"");
}


int doassembly(LocalParameters &par) {
    DBReader<unsigned int> *sequenceDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    sequenceDbr->open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> * alnReader = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    alnReader->open(DBReader<unsigned int>::NOSORT);

    DBWriter resultWriter(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, sequenceDbr->getDbtype());
    resultWriter.open();

    bool reverseResult = false;
    int seqType = sequenceDbr->getDbtype();
    BaseMatrix *subMat;
    if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.c_str(), 1.0, 0.0);
        reverseResult = true;
    } else {
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0, 0.0);
    }

    SubstitutionMatrix::FastMatrix fastMatrix = SubstitutionMatrix::createAsciiSubMat(*subMat);

    unsigned char * wasExtended = new unsigned char[sequenceDbr->getSize()];
    std::fill(wasExtended, wasExtended+sequenceDbr->getSize(), 0);
    Debug::Progress progress(sequenceDbr->getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

        std::vector<Matcher::result_t> alignments;
        alignments.reserve(300);
        char *queryRevSeq = NULL;
        int queryRevSeqLen = par.maxSeqLen;
        if (reverseResult == true) {
            queryRevSeq = new char[queryRevSeqLen];
        }
#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < sequenceDbr->getSize(); id++) {
            progress.updateProgress();
            unsigned int queryId = sequenceDbr->getDbKey(id);
            char *querySeq = sequenceDbr->getData(id, thread_idx);
            char *querySeqToUse;
            unsigned int querySeqLen = sequenceDbr->getSeqLens(id) - 2;
            unsigned int leftQueryOffset = 0;
            unsigned int rightQueryOffset = 0;
            std::string query(querySeq, querySeqLen); // no /n/0

            if (reverseResult == true) {
                NucleotideMatrix *nuclMatrix = (NucleotideMatrix *) subMat;
                for (int pos = querySeqLen - 1; pos > -1; pos--) {
                    int res = subMat->aa2int[static_cast<int>(querySeq[pos])];
                    queryRevSeq[(querySeqLen - 1) - pos] = subMat->int2aa[nuclMatrix->reverseResidue(res)];
                }
            }
            std::string queryRev(queryRevSeq,querySeqLen);

            char *alnData = alnReader->getDataByDBKey(queryId, thread_idx);
            alignments.clear();
            Matcher::readAlignmentResults(alignments, alnData);
            QueueBySeqId alnQueue;
            bool queryCouldBeExtended = false;
            bool isReverse;
            while(alignments.size() > 1){
                bool queryCouldBeExtendedLeft = false;
                bool queryCouldBeExtendedRight = false;
                for (size_t alnIdx = 0; alnIdx < alignments.size(); alnIdx++) {
                    float scorePerCol = static_cast<float>(BIT_CLEAR(alignments[alnIdx].score,31)) / static_cast<float>(alignments[alnIdx].alnLength);
                    float alnLen = static_cast<float>(alignments[alnIdx].alnLength);
                    float ids = static_cast<float>(alignments[alnIdx].seqId) * alnLen;
                    alignments[alnIdx].seqId = ids / (alnLen + 0.5);
                    alignments[alnIdx].score = static_cast<int>(scorePerCol*100);
                    alignments[alnIdx].score = BIT_SET(alignments[alnIdx].score, 31);
                    if(alignments[alnIdx].qStartPos > alignments[alnIdx].qEndPos  && reverseResult){
                        // alignment is on reverse of query sequence
                        alignments[alnIdx].qStartPos = alignments[alnIdx].qLen-alignments[alnIdx].qStartPos-1;
                        alignments[alnIdx].qEndPos = alignments[alnIdx].qLen-alignments[alnIdx].qEndPos-1;
                        alignments[alnIdx].score = BIT_CLEAR(alignments[alnIdx].score, 31);
                    }

                    alnQueue.push(alignments[alnIdx]);
                    if (alignments.size() > 1)
                        __sync_or_and_fetch(&wasExtended[sequenceDbr->getId(alignments[alnIdx].dbKey)],
                                            static_cast<unsigned char>(0x40));
                }
                std::vector<Matcher::result_t> tmpAlignments;
                Matcher::result_t besttHitToExtend;
                while ((besttHitToExtend = selectFragmentToExtend(alnQueue, queryId)).dbKey != UINT_MAX) {
                    querySeqToUse = (char *) query.c_str();
                    isReverse = false;
                    querySeqLen = query.size();
                    if (reverseResult && BIT_CHECK(besttHitToExtend.score,31)==false){
                        querySeqToUse = (char *) queryRev.c_str();
                        size_t temp = leftQueryOffset;
                        leftQueryOffset = rightQueryOffset;
                        rightQueryOffset = temp;
                        isReverse = true;
                    }
//                querySeq.mapSequence(id, queryKey, query.c_str());
                    unsigned int targetId = sequenceDbr->getId(besttHitToExtend.dbKey);
                    if (targetId == UINT_MAX) {
                        Debug(Debug::ERROR) << "Could not find targetId  " << besttHitToExtend.dbKey
                                            << " in database " << sequenceDbr->getDataFileName() << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                    char *targetSeq = sequenceDbr->getData(targetId, thread_idx);
                    unsigned int targetSeqLen = sequenceDbr->getSeqLens(targetId) - 2;
                    // check if alignment still make sense (can extend the query)
                    if (besttHitToExtend.dbStartPos == 0) {
                        if ((targetSeqLen - (besttHitToExtend.dbEndPos + 1)) <= rightQueryOffset) {
                            continue;
                        }
                    } else if (besttHitToExtend.qStartPos == 0) {
                        if (besttHitToExtend.dbStartPos <= leftQueryOffset) {
                            continue;
                        }
                    }
                    __sync_or_and_fetch(&wasExtended[targetId], static_cast<unsigned char>(0x10));
                    int qStartPos, qEndPos, dbStartPos, dbEndPos, score;
                    int diagonal = (leftQueryOffset + besttHitToExtend.qStartPos) - besttHitToExtend.dbStartPos;

                    int dist = std::max(abs(diagonal), 0);
                    if (diagonal >= 0) {
//                    targetSeq.mapSequence(targetId, besttHitToExtend.dbKey, dbSeq);
                        size_t diagonalLen = std::min(targetSeqLen, querySeqLen - abs(diagonal));
                        DistanceCalculator::LocalAlignment alignment = DistanceCalculator::computeSubstitutionStartEndDistance(
                                querySeqToUse + abs(diagonal),
                                targetSeq, diagonalLen, fastMatrix.matrix);
                        qStartPos = alignment.startPos + dist;
                        qEndPos = alignment.endPos + dist;
                        dbStartPos = alignment.startPos;
                        dbEndPos = alignment.endPos;
                        score = alignment.score;
                    } else {
                        size_t diagonalLen = std::min(targetSeqLen - abs(diagonal), querySeqLen);
                        DistanceCalculator::LocalAlignment alignment = DistanceCalculator::computeSubstitutionStartEndDistance(
                                querySeqToUse,
                                targetSeq + abs(diagonal),
                                diagonalLen, fastMatrix.matrix);
                        qStartPos = alignment.startPos;
                        qEndPos = alignment.endPos;
                        dbStartPos = alignment.startPos + dist;
                        dbEndPos = alignment.endPos + dist;
                        score = alignment.score;
                    }

                    if (dbStartPos == 0 && qEndPos == (querySeqLen - 1) ) {
                        if((!isReverse && queryCouldBeExtendedRight == true) || (isReverse && queryCouldBeExtendedLeft == true)) {
                            float alnLen = qEndPos - qStartPos;
                            float scorePerCol = static_cast<float>(score) / (alnLen+0.5);
                            besttHitToExtend.score = static_cast<int>(scorePerCol*100);
                            if(!isReverse){
                                besttHitToExtend.score = BIT_SET(besttHitToExtend.score,31);
                            }
                            tmpAlignments.push_back(besttHitToExtend);
                            continue;
                        }
                        size_t dbFragLen = (targetSeqLen - dbEndPos) - 1; // -1 get not aligned element
                        std::string fragment = std::string(targetSeq + dbEndPos + 1, dbFragLen);
                        if (fragment.size() + query.size() >= par.maxSeqLen) {
                            Debug(Debug::WARNING) << "Sequence too long in query id: " << queryId << ". "
                                    "Max length allowed would is " << par.maxSeqLen << "\n";
                            break;
                        }
                        //update that dbKey was used in assembly
                        __sync_or_and_fetch(&wasExtended[targetId], static_cast<unsigned char>(0x80));
                        if (reverseResult) {
                            NucleotideMatrix *nuclMatrix = (NucleotideMatrix *) subMat;
                            const char *c_fragment = fragment.c_str();
                            char *c_fragmentRev = new char[dbFragLen];
                            for (int pos = dbFragLen - 1; pos > -1; pos--) {
                                int res = subMat->aa2int[static_cast<int>(c_fragment[pos])];
                                c_fragmentRev[(dbFragLen - 1) - pos] = subMat->int2aa[nuclMatrix->reverseResidue(res)];
                            }
                            std::string fragmentRev = std::string(c_fragmentRev,dbFragLen);
                            if (!isReverse) {
                                queryCouldBeExtendedRight = true;
                                query = query + fragment;
                                queryRev = fragmentRev + queryRev;
                            } else {
                                queryCouldBeExtendedLeft = true;
                                queryRev = queryRev + fragment;
                                query = fragmentRev + query;
                            }
                        }
                        else {
                            query += fragment;
                            queryCouldBeExtendedRight = true;
                        }
                        rightQueryOffset += dbFragLen;

                    } else if (qStartPos == 0 && dbEndPos == (targetSeqLen - 1)) {
                        if ((!isReverse && queryCouldBeExtendedLeft == true)|| (isReverse && queryCouldBeExtendedRight == true)) {
                            float alnLen = qEndPos - qStartPos;
                            float scorePerCol = static_cast<float>(score) / (alnLen+0.5);
                            besttHitToExtend.score = static_cast<int>(scorePerCol*100);
                            if(!isReverse){
                                besttHitToExtend.score = BIT_SET(besttHitToExtend.score,31);
                            }
                            tmpAlignments.push_back(besttHitToExtend);
                            continue;
                        }
                        size_t dbFragLen = dbStartPos;
                        std::string fragment = std::string(targetSeq, dbFragLen); // +1 get not aligned element
                        if (fragment.size() + query.size() >= par.maxSeqLen) {
                            Debug(Debug::WARNING) << "Sequence too long in query id: " << queryId << ". "
                                    "Max length allowed would is " << par.maxSeqLen << "\n";
                            break;
                        }
                        // update that dbKey was used in assembly
                        __sync_or_and_fetch(&wasExtended[targetId], static_cast<unsigned char>(0x80));
                        if (reverseResult) {
                            NucleotideMatrix *nuclMatrix = (NucleotideMatrix *) subMat;
                            const char *c_fragment = fragment.c_str();
                            char *c_fragmentRev = new char[dbFragLen];
                            for (int pos = dbFragLen - 1; pos > -1; pos--) {
                                int res = subMat->aa2int[static_cast<int>(c_fragment[pos])];
                                c_fragmentRev[(dbFragLen - 1) - pos] = subMat->int2aa[nuclMatrix->reverseResidue(res)];
                            }
                            std::string fragmentRev = std::string(c_fragmentRev,dbFragLen);
                            if (!isReverse) {
                                queryCouldBeExtendedLeft = true;
                                query = fragment + query;
                                queryRev = queryRev + fragmentRev;
                            } else {
                                queryCouldBeExtendedRight = true;
                                queryRev = fragment + queryRev;
                                query = query + fragmentRev;
                            }
                        }
                        else {
                            queryCouldBeExtendedLeft = true;
                            query = fragment + query;
                        }
                        leftQueryOffset += dbStartPos;
                    }
                    if (isReverse){
                        size_t temp = leftQueryOffset;
                        leftQueryOffset = rightQueryOffset;
                        rightQueryOffset = temp;
                    }

                }
                if (queryCouldBeExtendedRight || queryCouldBeExtendedLeft){
                    queryCouldBeExtended = true;
                }
                alignments.clear();
                break;
                querySeqLen = query.size();
                querySeq = (char *) query.c_str();
                for(size_t alnIdx = 0; alnIdx < tmpAlignments.size(); alnIdx++){
                    int idCnt = 0;
                    int qStartPos = tmpAlignments[alnIdx].qStartPos;
                    int qEndPos = tmpAlignments[alnIdx].qEndPos;
                    int dbStartPos = tmpAlignments[alnIdx].dbStartPos;
                    int diagonal = (leftQueryOffset + besttHitToExtend.qStartPos) - besttHitToExtend.dbStartPos;
                    int dist = std::max(abs(diagonal), 0);
                    if (diagonal >= 0) {
                        qStartPos+=dist;
                        qEndPos+=dist;
                    }else{
                        dbStartPos+=dist;
                    }
                    unsigned int targetId = sequenceDbr->getId(tmpAlignments[alnIdx].dbKey);
                    char *targetSeq = sequenceDbr->getData(targetId, thread_idx);
                    for(int i = qStartPos; i <= qEndPos; i++){
                        int targetRes = static_cast<int>(targetSeq[dbStartPos+(i-qStartPos)]);
                        int queryRes = static_cast<int>(querySeq[i]);
                        idCnt += (queryRes == targetRes) ? 1 : 0;
                    }
                    float seqId =  static_cast<float>(idCnt) / (static_cast<float>(qEndPos) - static_cast<float>(qStartPos) + 0.5);
                    tmpAlignments[alnIdx].seqId = seqId;
                    if(seqId >= par.seqIdThr){
                        alignments.push_back(tmpAlignments[alnIdx]);
                    }
                }
            }
            if (queryCouldBeExtended == true) {
                query.push_back('\n');
                __sync_or_and_fetch(&wasExtended[id], static_cast<unsigned char>(0x20));
                resultWriter.writeData(query.c_str(), query.size(), queryId, thread_idx);
            }
        }
    } // end parallel

// add sequences that are not yet assembled
#pragma omp parallel for schedule(dynamic, 10000)
    for (size_t id = 0; id < sequenceDbr->getSize(); id++) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        //   bool couldExtend =  (wasExtended[id] & 0x10);
        bool isNotContig =  !(wasExtended[id] & 0x20);
//        bool wasNotUsed =  !(wasExtended[id] & 0x40);
//        bool wasNotExtended =  !(wasExtended[id] & 0x80);
        //    bool wasUsed    =  (wasExtended[id] & 0x40);
        //if(isNotContig && wasNotExtended ){
        if (isNotContig){
            char *querySeqData = sequenceDbr->getData(id, thread_idx);
            unsigned int queryLen = sequenceDbr->getSeqLens(id) - 1; //skip null byte
            resultWriter.writeData(querySeqData, queryLen, sequenceDbr->getDbKey(id), thread_idx);
        }
    }

    // cleanup

    resultWriter.close(true);
    alnReader->close();
    delete [] wasExtended;
    delete alnReader;
    delete [] fastMatrix.matrix;
    delete [] fastMatrix.matrixData;
    sequenceDbr->close();
    delete sequenceDbr;
    Debug(Debug::INFO) << "\nDone.\n";

    return EXIT_SUCCESS;
}

int assembleresult(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, 3);

    MMseqsMPI::init(argc, argv);

    // never allow deletions
    par.allowDeletion = false;
    Debug(Debug::INFO) << "Compute assembly.\n";
    return doassembly(par);
}

