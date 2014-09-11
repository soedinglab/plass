//
//  QueryScoreLocal.h
//  kClust2_xcode
//
//  Created by Martin Steinegger on 03.06.13.
//  Copyright (c) 2013 Martin Steinegger. All rights reserved.
//

#ifndef QUERYSCORELOCAL_H
#define QUERYSCORELOCAL_H
#include "QueryScore.h"
#include <vector>
#include <map>



class QueryScoreLocal : public QueryScore {
    
public:
    QueryScoreLocal(int dbSize, unsigned short * seqLens, int k, short kmerThr, double kmerMatchProb, float zscoreThr);
   
    ~QueryScoreLocal();
    
    std::pair<hit_t *, size_t> getResult (int querySeqLen, unsigned int identityId);

    void reset();
    
    // NOT needed for Local scoring
    void setPrefilteringThresholds();
    
    // Sort local results by sequence id, i and j
    static bool compareLocalResult(LocalResult first, LocalResult second);
private:
    unsigned short gapOpenPenalty;
    unsigned short gapExtendPenalty;
    unsigned short minKmerMatch; // the minimal score from kmer scores



};
#endif /* defined(QUERYSCORESEMILOCAL_H) */