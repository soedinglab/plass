
#include "NucleotideMatrix.h"
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
#include <string>
#include <vector>
#include <sstream>

#ifdef OPENMP
#include <omp.h>
#endif

class CompareResultBySeqId {
public:
    /*bool operator() (const Matcher::result_t & r1,const Matcher::result_t & r2) {
        if(r1.seqId < r2.seqId )
            return true;
        if(r2.seqId < r1.seqId )
            return false;
        if(r1.alnLength < r2.alnLength )
            return true;
        if(r2.alnLength < r1.alnLength )
            return false;
        if(r1.dbKey > r2.dbKey )
            return true;
        if(r2.dbKey > r1.dbKey )
            return false;

        return false;
    }*/
    bool operator() (const Matcher::result_t & r1,const Matcher::result_t & r2) {
        unsigned int mm_count1 = (1 - r1.seqId) * r1.alnLength + 0.5;
        unsigned int mm_count2 = (1 - r2.seqId) * r2.alnLength + 0.5;

        unsigned int alpha1 = mm_count1 + 1;
        unsigned int alpha2 = mm_count2 + 1;
        unsigned int beta1 = r1.alnLength - mm_count1 + 1;
        unsigned int beta2 = r2.alnLength - mm_count2 + 1;

        //double c=(std::tgamma(beta1+beta2)*std::tgamma(alpha1+beta2))/(std::tgamma(alpha1+beta1+beta2)*std::tgamma(beta1));
        double log_c = (std::lgamma(beta1+beta2)+std::lgamma(alpha1+beta1))-(std::lgamma(alpha1+beta1+beta2)+std::lgamma(beta1));

        //double r = 1.0; // r_0 =1
        double log_r = 0.0;
        double p = 0.0;
        for (size_t idx = 0; idx < alpha2; idx++) {

            p += exp(log_r + log_c);
            //r *= ((alpha1+idx)*(beta2+idx))/((idx+1)*(idx+alpha1+beta1+beta2));
            log_r = log(alpha1+idx)+log(beta2+idx)-(log(idx+1) + log(idx+alpha1+beta1+beta2)) + log_r;
        }
        //p *= c;

        if (p < 0.45)
            return true;
        if (p  > 0.55)
            return false;
        if (r1.dbLen - r1.alnLength < r2.dbLen - r2.alnLength)
            return true;
        if (r1.dbLen - r1.alnLength > r2.dbLen - r2.alnLength)
            return false;

        return true;
    }
};

typedef std::priority_queue<Matcher::result_t, std::vector<Matcher::result_t> , CompareResultBySeqId> QueueBySeqId;
Matcher::result_t selectBestFragmentToExtend(QueueBySeqId &alignments,
                                             unsigned int queryKey) {
    // results are ordered by seqid
    while (alignments.empty() == false){
        Matcher::result_t res = alignments.top();
        alignments.pop();
        size_t dbKey = res.dbKey;
        const bool notRightStartAndLeftStart = !(res.dbStartPos == 0 && res.qStartPos == 0);
        const bool rightStart = res.dbStartPos == 0 && (res.dbEndPos != static_cast<int>(res.dbLen)-1);
        const bool leftStart = res.qStartPos == 0   && (res.qEndPos != static_cast<int>(res.qLen)-1);
        const bool isNotIdentity = (dbKey != queryKey);
        if ((rightStart || leftStart) && notRightStartAndLeftStart && isNotIdentity){
            return res;
        }
    }
    return Matcher::result_t(UINT_MAX,0,0,0,0,0,0,0,0,0,0,0,0,"");
}


int dohybridassembleresult(LocalParameters &par) {
    DBReader<unsigned int> *nuclSequenceDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str(),  par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    nuclSequenceDbr->open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> *aaSequenceDbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(),  par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    aaSequenceDbr->open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> * nuclAlnReader = new DBReader<unsigned int>(par.db3.c_str(), par.db3Index.c_str(),  par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    nuclAlnReader->open(DBReader<unsigned int>::NOSORT);

    DBWriter nuclResultWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_NUCLEOTIDES);
    nuclResultWriter.open();

    DBWriter aaResultWriter(par.db5.c_str(), par.db5Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
    aaResultWriter.open();

    NucleotideMatrix subMat(par.scoringMatrixFile.nucleotides, 1.0f, 0.0f);
    SubstitutionMatrix::FastMatrix fastMatrix = SubstitutionMatrix::createAsciiSubMat(subMat);

    unsigned char * wasExtended = new unsigned char[nuclSequenceDbr->getSize()];
    std::fill(wasExtended, wasExtended+nuclSequenceDbr->getSize(), 0);
    Debug::Progress progress(nuclSequenceDbr->getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        std::vector<Matcher::result_t> nuclAlignments;
        nuclAlignments.reserve(300);

        #pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < nuclSequenceDbr->getSize(); id++) {
            progress.updateProgress();
            unsigned int queryKey = nuclSequenceDbr->getDbKey(id);

            char *nuclQuerySeq = nuclSequenceDbr->getData(id, thread_idx);
            unsigned int nuclQuerySeqLen = nuclSequenceDbr->getSeqLen(id);

            unsigned int aaQueryId = aaSequenceDbr->getId(queryKey);
            char *aaQuerySeq = aaSequenceDbr->getData(aaQueryId, thread_idx);
            unsigned int aaQuerySeqLen = aaSequenceDbr->getSeqLen(aaQueryId);

            unsigned int nuclLeftQueryOffset = 0;
            unsigned int nuclRightQueryOffset = 0;
            std::string nuclQuery(nuclQuerySeq, nuclQuerySeqLen); // no /n/0
            std::string aaQuery(aaQuerySeq, aaQuerySeqLen); // no /n/0

            bool excludeLeftExtension = (aaQuery[0] == '*');
            bool excludeRightExtension = (aaQuery[aaQuerySeqLen-1] == '*');

            char *nuclAlnData = nuclAlnReader->getDataByDBKey(queryKey, thread_idx);

            nuclAlignments.clear();
            Matcher::readAlignmentResults(nuclAlignments, nuclAlnData, true);

            QueueBySeqId alnQueue;
            bool queryCouldBeExtended = false;
            while(nuclAlignments.size() > 1){
                bool queryCouldBeExtendedLeft = false;
                bool queryCouldBeExtendedRight = false;
                for (size_t alnIdx = 0; alnIdx < nuclAlignments.size(); alnIdx++) {
                    if(nuclAlignments[alnIdx].seqId < par.seqIdThr)
                        continue;
                    alnQueue.push(nuclAlignments[alnIdx]);
                    if (nuclAlignments.size() > 1) {
                        size_t id = nuclSequenceDbr->getId(nuclAlignments[alnIdx].dbKey);
                        __sync_or_and_fetch(&wasExtended[id],
                                            static_cast<unsigned char>(0x40));
                    }
                }
                std::vector<Matcher::result_t> tmpNuclAlignments;

                Matcher::result_t nuclBesttHitToExtend;

                while ((nuclBesttHitToExtend = selectBestFragmentToExtend(alnQueue, queryKey)).dbKey != UINT_MAX) {
                    nuclQuerySeqLen = nuclQuery.size();
                    nuclQuerySeq = (char *) nuclQuery.c_str();

//                nuclQuerySeq.mapSequence(id, queryKey, nuclQuery.c_str());
                    unsigned int nuclTargetId = nuclSequenceDbr->getId(nuclBesttHitToExtend.dbKey);
                    if (nuclTargetId == UINT_MAX) {
                        Debug(Debug::ERROR) << "Could not find nuclTargetId  " << nuclBesttHitToExtend.dbKey
                                            << " in database " << nuclSequenceDbr->getDataFileName() << "\n";
                        EXIT(EXIT_FAILURE);
                    }

                    char *nuclTargetSeq = nuclSequenceDbr->getData(nuclTargetId, thread_idx);
                    unsigned int nuclTargetSeqLen = nuclSequenceDbr->getSeqLen(nuclTargetId);
                    
                    unsigned int aaTargetId = aaSequenceDbr->getId(nuclBesttHitToExtend.dbKey);
                    char *aaTargetSeq = aaSequenceDbr->getData(aaTargetId, thread_idx);
                    unsigned int aaTargetSeqLen = aaSequenceDbr->getSeqLen(aaTargetId) ;

                    // check if alignment still make sense (can extend the nuclQuery)
                    // avoid extension over start/stoppcodons
                    if (nuclBesttHitToExtend.dbStartPos == 0) {
                        if (((nuclTargetSeqLen - (nuclBesttHitToExtend.dbEndPos + 1)) <= nuclRightQueryOffset) || excludeRightExtension ||
                              aaTargetSeq[0] == '*') {
                            continue;
                        }
                    } else if (nuclBesttHitToExtend.qStartPos == 0) {
                        if ((nuclBesttHitToExtend.dbStartPos <= static_cast<int>(nuclLeftQueryOffset)) || excludeLeftExtension ||
                           aaTargetSeq[aaTargetSeqLen-1] == '*') {
                            continue;
                        }
                    }
                    __sync_or_and_fetch(&wasExtended[nuclTargetId], static_cast<unsigned char>(0x10));
                    int qStartPos, qEndPos, nuclDbStartPos, nuclDbEndPos;
                    int diagonal = (nuclLeftQueryOffset + nuclBesttHitToExtend.qStartPos) - nuclBesttHitToExtend.dbStartPos;
                    int dist = std::max(abs(diagonal), 0);
                    DistanceCalculator::LocalAlignment alignment = DistanceCalculator::ungappedAlignmentByDiagonal(
                            nuclQuerySeq, nuclQuerySeqLen,
                            nuclTargetSeq, nuclTargetSeqLen,
                            diagonal, fastMatrix.matrix, par.rescoreMode);
                    if (diagonal >= 0) {

//                    nuclTargetSeq.mapSequence(nuclTargetId, nuclBesttHitToExtend.dbKey, dbSeq);

                        qStartPos = alignment.startPos + dist;
                        qEndPos = alignment.endPos + dist;
                        nuclDbStartPos = alignment.startPos;
                        nuclDbEndPos = alignment.endPos;
                    } else {

                        qStartPos = alignment.startPos;
                        qEndPos = alignment.endPos;
                        nuclDbStartPos = alignment.startPos + dist;
                        nuclDbEndPos = alignment.endPos + dist;
                    }

                    if (nuclDbStartPos == 0 && qEndPos == (static_cast<int>(nuclQuerySeqLen) - 1) ) {
                        if(queryCouldBeExtendedRight == true) {
                            tmpNuclAlignments.push_back(nuclBesttHitToExtend);
                            continue;
                        }
                        size_t nuclDbFragLen = (nuclTargetSeqLen - nuclDbEndPos) - 1; // -1 get not aligned element
                        size_t aaDbFragLen = (nuclTargetSeqLen/3 - nuclDbEndPos/3) - 1; // -1 get not aligned element

                        std::string fragment = std::string(nuclTargetSeq + nuclDbEndPos + 1, nuclDbFragLen);
                        std::string aaFragment = std::string(aaTargetSeq + nuclDbEndPos/3 + 1, aaDbFragLen);

                        if (fragment.size() + nuclQuery.size() >= par.maxSeqLen) {
                            Debug(Debug::WARNING) << "Sequence too long in nuclQuery id: " << queryKey << ". "
                                    "Max length allowed would is " << par.maxSeqLen << "\n";
                            break;
                        }
                        //update that dbKey was used in assembly
                        __sync_or_and_fetch(&wasExtended[nuclTargetId], static_cast<unsigned char>(0x80));
                        queryCouldBeExtendedRight = true;
                        nuclQuery += fragment;
                        aaQuery += aaFragment;

                        nuclRightQueryOffset += nuclDbFragLen;

                    } else if (qStartPos == 0 && nuclDbEndPos == (static_cast<int>(nuclTargetSeqLen) - 1)) {
                        if (queryCouldBeExtendedLeft == true) {
                            tmpNuclAlignments.push_back(nuclBesttHitToExtend);
                            continue;
                        }
                        int hasStart = (aaTargetSeq[0] == '*')? 1:0;
                        std::string fragment = std::string(nuclTargetSeq, nuclDbStartPos); // +1 get not aligned element
                        std::string aaFragment = std::string(aaTargetSeq, nuclDbStartPos/3 + hasStart); // +1 get not aligned element

                        if (fragment.size() + nuclQuery.size() >= par.maxSeqLen) {
                            Debug(Debug::WARNING) << "Sequence too long in nuclQuery id: " << queryKey << ". "
                                    "Max length allowed would is " << par.maxSeqLen << "\n";
                            break;
                        }
                        // update that dbKey was used in assembly
                        __sync_or_and_fetch(&wasExtended[nuclTargetId], static_cast<unsigned char>(0x80));
                        queryCouldBeExtendedLeft = true;
                        nuclQuery = fragment + nuclQuery;
                        aaQuery = aaFragment + aaQuery;
                        nuclLeftQueryOffset += nuclDbStartPos;
                    }

                }
                if (queryCouldBeExtendedRight || queryCouldBeExtendedLeft){
                    queryCouldBeExtended = true;
                }
                nuclAlignments.clear();
                nuclQuerySeq = (char *) nuclQuery.c_str();
                break;
                for(size_t alnIdx = 0; alnIdx < tmpNuclAlignments.size(); alnIdx++){
                    int idCnt = 0;
                    int qStartPos = tmpNuclAlignments[alnIdx].qStartPos;
                    int qEndPos = tmpNuclAlignments[alnIdx].qEndPos;
                    int dbStartPos = tmpNuclAlignments[alnIdx].dbStartPos;

                    int diagonal = (nuclLeftQueryOffset + nuclBesttHitToExtend.qStartPos) - nuclBesttHitToExtend.dbStartPos;
                    int dist = std::max(abs(diagonal), 0);
                    if (diagonal >= 0) {
                        qStartPos+=dist;
                        qEndPos+=dist;
                    }else{
                        dbStartPos+=dist;
                    }
                    unsigned int targetId = nuclSequenceDbr->getId(tmpNuclAlignments[alnIdx].dbKey);
                    char *nuclTargetSeq = nuclSequenceDbr->getData(targetId, thread_idx);
                    for(int i = qStartPos; i < qEndPos; i++){
                        idCnt += (nuclQuerySeq[i] == nuclTargetSeq[dbStartPos+(i-qStartPos)]) ? 1 : 0;
                    }
                    float seqId =  static_cast<float>(idCnt) / (static_cast<float>(qEndPos) - static_cast<float>(qStartPos));
                    tmpNuclAlignments[alnIdx].seqId = seqId;
                    if(seqId >= par.seqIdThr){
                        nuclAlignments.push_back(tmpNuclAlignments[alnIdx]);
                    }
                }
            }
            if (queryCouldBeExtended == true) {
                nuclQuery.push_back('\n');
                aaQuery.push_back('\n');
                __sync_or_and_fetch(&wasExtended[id], static_cast<unsigned char>(0x20));
                nuclResultWriter.writeData(nuclQuery.c_str(), nuclQuery.size(), queryKey, thread_idx);
                aaResultWriter.writeData(aaQuery.c_str(), aaQuery.size(), queryKey, thread_idx);
            }
        }
    } // end parallel

// add sequences that are not yet assembled
#pragma omp parallel for schedule(dynamic, 10000)
    for (size_t id = 0; id < nuclSequenceDbr->getSize(); id++) {
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
            char *querySeqData = nuclSequenceDbr->getData(id, thread_idx);
            unsigned int queryLen = nuclSequenceDbr->getEntryLen(id) - 1; //skip null byte
            nuclResultWriter.writeData(querySeqData, queryLen, nuclSequenceDbr->getDbKey(id), thread_idx);
            char *queryAASeqData = aaSequenceDbr->getData(id, thread_idx);
            unsigned int queryAALen = aaSequenceDbr->getEntryLen(id) - 1; //skip null byte
            aaResultWriter.writeData(queryAASeqData, queryAALen, aaSequenceDbr->getDbKey(id), thread_idx);
        }
    }

    // cleanup
    aaResultWriter.close(aaSequenceDbr->getDbtype());
    nuclResultWriter.close(nuclSequenceDbr->getDbtype());
    nuclAlnReader->close();
    delete [] wasExtended;
    delete nuclAlnReader;

    delete [] fastMatrix.matrix;
    delete [] fastMatrix.matrixData;
    aaSequenceDbr->close();
    delete aaSequenceDbr;
    nuclSequenceDbr->close();
    delete nuclSequenceDbr;
    Debug(Debug::INFO) << "\nDone.\n";

    return EXIT_SUCCESS;
}

int hybridassembleresults(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    MMseqsMPI::init(argc, argv);

    // never allow deletions
    par.allowDeletion = false;
    Debug(Debug::INFO) << "Compute assembly.\n";
    return dohybridassembleresult(par);
}

