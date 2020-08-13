#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include <Parameters.h>
#include <MultiParam.h>
#include <algorithm>
#include <cfloat>

class LocalParameters : public Parameters {
public:
    static void initInstance() {
        new LocalParameters;
    }

    static LocalParameters& getLocalInstance() {
        if (instance == NULL) {
            initInstance();
        }
        return static_cast<LocalParameters&>(LocalParameters::getInstance());
    }

    std::vector<MMseqsParameter *> assemblerworkflow;
    std::vector<MMseqsParameter *> nuclassemblerworkflow;
    std::vector<MMseqsParameter *> hybridassemblerworkflow;

    std::vector<MMseqsParameter *> assembleDBworkflow;
    std::vector<MMseqsParameter *> nuclassembleDBworkflow;
    std::vector<MMseqsParameter *> hybridassembleDBworkflow;

    std::vector<MMseqsParameter *> assembleresults;
    std::vector<MMseqsParameter *> cyclecheck;
    std::vector<MMseqsParameter *> createhdb;
    std::vector<MMseqsParameter *> extractorfssubset;
    std::vector<MMseqsParameter *> filternoncoding;
    std::vector<MMseqsParameter *> hybridassembleresults;
    std::vector<MMseqsParameter *> reduceredundancy;


    int filterProteins;
    int deleteFilesInc;
    int minContigLen;
    float clustSeqIdThr;
    float clustCovThr;
    float proteinFilterThreshold;
    bool cycleCheck;
    bool chopCycle;

    MultiParam<int> multiNumIterations;
    MultiParam<int> multiKmerSize;
    MultiParam<int> multiAlnLenThr;
    MultiParam<float> multiSeqIdThr;

    PARAMETER(PARAM_FILTER_PROTEINS)
    PARAMETER(PARAM_PROTEIN_FILTER_THRESHOLD)
    PARAMETER(PARAM_DELETE_TMP_INC)
    PARAMETER(PARAM_MIN_CONTIG_LEN)
    PARAMETER(PARAM_CLUST_MIN_SEQ_ID_THR)
    PARAMETER(PARAM_CLUST_C)
    PARAMETER(PARAM_CYCLE_CHECK)
    PARAMETER(PARAM_CHOP_CYCLE)
    PARAMETER(PARAM_MULTI_NUM_ITERATIONS)
    PARAMETER(PARAM_MULTI_K)
    PARAMETER(PARAM_MULTI_MIN_SEQ_ID)
    PARAMETER(PARAM_MULTI_MIN_ALN_LEN)


private:
    LocalParameters() :
            Parameters(),
            multiNumIterations(INT_MAX,INT_MAX),
            multiKmerSize(INT_MAX,INT_MAX),
            multiAlnLenThr(INT_MAX,INT_MAX),
            multiSeqIdThr(FLT_MAX,FLT_MAX),
            PARAM_FILTER_PROTEINS(PARAM_FILTER_PROTEINS_ID,"--filter-proteins", "Filter Proteins", "filter proteins by a neural network [0,1]",typeid(int), (void *) &filterProteins, "^[0-1]{1}$"),
            PARAM_PROTEIN_FILTER_THRESHOLD(PARAM_PROTEIN_FILTER_THRESHOLD_ID,"--protein-filter-threshold", "Protein Filter Threshold", "filter proteins lower than threshold [0.0,1.0]",typeid(float), (void *) &proteinFilterThreshold, "^0(\\.[0-9]+)?|1(\\.0+)?$"),
            PARAM_DELETE_TMP_INC(PARAM_DELETE_TMP_INC_ID,"--delete-tmp-inc", "Delete temporary files incremental", "Delete temporary files incremental [0,1]",typeid(int), (void *) &deleteFilesInc, "^[0-1]{1}$", MMseqsParameter::COMMAND_COMMON | MMseqsParameter::COMMAND_EXPERT),
            PARAM_MIN_CONTIG_LEN(PARAM_MIN_CONTIG_LEN_ID, "--min-contig-len", "Minimum contig length", "Minimum length of assembled contig to output", typeid(int), (void *) &minContigLen, "^[1-9]{1}[0-9]*$"),
            PARAM_CLUST_MIN_SEQ_ID_THR(PARAM_CLUST_MIN_SEQ_ID_THR_ID,"--clust-min-seq-id", "Clustering seq. id. threshold","Seq. id. threshold passed to linclust algorithm to reduce redundancy in assembly (range 0.0-1.0)",typeid(float), (void *) &clustSeqIdThr, "^0(\\.[0-9]+)?|1(\\.0+)?$", MMseqsParameter::COMMAND_CLUST),
            PARAM_CLUST_C(PARAM_CLUST_C_ID,"--clust-min-cov", "Clustering coverage threshold","Coverage threshold passed to linclust algorithm to reduce redundancy in assembly (range 0.0-1.0)",typeid(float), (void *) &clustCovThr, "^0(\\.[0-9]+)?|1(\\.0+)?$", MMseqsParameter::COMMAND_CLUST),
            PARAM_CYCLE_CHECK(PARAM_CYCLE_CHECK_ID,"--cycle-check", "Check for circular sequences", "Check for circular sequences (avoid infinite extension of circular or long repeated regions) ",typeid(bool), (void *) &cycleCheck, "", MMseqsParameter::COMMAND_MISC | MMseqsParameter::COMMAND_EXPERT),
            PARAM_CHOP_CYCLE(PARAM_CHOP_CYCLE_ID,"--chop-cycle", "Chop Cycle", "Remove superfluous part of circular fragments (see --cycle-check)",typeid(bool), (void *) &chopCycle, "", MMseqsParameter::COMMAND_MISC | MMseqsParameter::COMMAND_EXPERT),
            PARAM_MULTI_NUM_ITERATIONS(PARAM_MULTI_NUM_ITERATIONS_ID, "--num-iterations", "Number of assembly iterations","Number of assembly iterations performed on nucleotide level,protein level (range 1-inf)",typeid(MultiParam<int>),(void *) &multiNumIterations, ""),
            PARAM_MULTI_K(PARAM_MULTI_K_ID, "-k", "k-mer length", "k-mer length (0: automatically set to optimum)", typeid(MultiParam<int>), (void *) &multiKmerSize, "", MMseqsParameter::COMMAND_CLUSTLINEAR | MMseqsParameter::COMMAND_EXPERT),
            PARAM_MULTI_MIN_SEQ_ID(PARAM_MULTI_MIN_SEQ_ID_ID, "--min-seq-id", "Seq. id. threshold", "Overlap sequence identity threshold [0.0, 1.0]", typeid(MultiParam<float>), (void *) &multiSeqIdThr, "", MMseqsParameter::COMMAND_ALIGN),
            PARAM_MULTI_MIN_ALN_LEN(PARAM_MULTI_MIN_ALN_LEN_ID, "--min-aln-len", "Min alignment length", "Minimum alignment length (range 0-INT_MAX)", typeid(MultiParam<int>), (void *) &multiAlnLenThr, "", MMseqsParameter::COMMAND_ALIGN)


    {
        // assembleresult
        assembleresults.push_back(&PARAM_MIN_SEQ_ID);
        assembleresults.push_back(&PARAM_MAX_SEQ_LEN);
        assembleresults.push_back(&PARAM_THREADS);
        assembleresults.push_back(&PARAM_V);
        assembleresults.push_back(&PARAM_RESCORE_MODE); //temporary added until assemble and nuclassemble use same rescoremode

        extractorfssubset.push_back(&PARAM_TRANSLATION_TABLE);
        extractorfssubset.push_back(&PARAM_USE_ALL_TABLE_STARTS);
        extractorfssubset.push_back(&PARAM_THREADS);
        extractorfssubset.push_back(&PARAM_V);

        filternoncoding.push_back(&PARAM_PROTEIN_FILTER_THRESHOLD);
        filternoncoding.push_back(&PARAM_THREADS);
        filternoncoding.push_back(&PARAM_V);

        //cyclecheck
        cyclecheck.push_back(&PARAM_MAX_SEQ_LEN);
        cyclecheck.push_back(&PARAM_CHOP_CYCLE);
        cyclecheck.push_back(&PARAM_THREADS);
        cyclecheck.push_back(&PARAM_V);

        //createhdb
        createhdb.push_back(&PARAM_COMPRESSED);
        createhdb.push_back(&PARAM_V);

        //reduceredundancy (subset of clustering parameters which have to be adjusted)
        reduceredundancy.push_back(&PARAM_ALPH_SIZE);
        reduceredundancy.push_back(&PARAM_CLUSTER_MODE);
        reduceredundancy.push_back(&PARAM_K);
        reduceredundancy.push_back(&PARAM_KMER_PER_SEQ);
        reduceredundancy.push_back(&PARAM_KMER_PER_SEQ_SCALE);
        reduceredundancy.push_back(&PARAM_IGNORE_MULTI_KMER);
        reduceredundancy.push_back(&PARAM_MIN_SEQ_ID);
        reduceredundancy.push_back(&PARAM_COV_MODE);
        reduceredundancy.push_back(&PARAM_C);
        reduceredundancy.push_back(&PARAM_MAX_SEQ_LEN);
        reduceredundancy.push_back(&PARAM_WRAPPED_SCORING);
        reduceredundancy.push_back(&PARAM_GAP_OPEN);
        reduceredundancy.push_back(&PARAM_GAP_EXTEND);
        reduceredundancy.push_back(&PARAM_ZDROP);
        reduceredundancy.push_back(&PARAM_THREADS);


        // assembledb workflow
        assembleDBworkflow = combineList(rescorediagonal, kmermatcher);
        assembleDBworkflow = combineList(assembleDBworkflow, extractorfs);
        assembleDBworkflow = combineList(assembleDBworkflow, assembleresults);
        assembleDBworkflow = combineList(assembleDBworkflow, filternoncoding);

        assembleDBworkflow.push_back(&PARAM_FILTER_PROTEINS);
        assembleDBworkflow.push_back(&PARAM_NUM_ITERATIONS);
        assembleDBworkflow.push_back(&PARAM_DELETE_TMP_INC);
        assembleDBworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
        assembleDBworkflow.push_back(&PARAM_RUNNER);

        // easyassembleworkflow
        assemblerworkflow = combineList(assembleDBworkflow, createdb);
        
        // nucl assembledb workflow
        nuclassembleDBworkflow = combineList(rescorediagonal, kmermatcher);
        nuclassembleDBworkflow = combineList(nuclassembleDBworkflow, assembleresults);
        nuclassembleDBworkflow = combineList(nuclassembleDBworkflow, cyclecheck);

        nuclassembleDBworkflow.push_back(&PARAM_CYCLE_CHECK);
        nuclassembleDBworkflow.push_back(&PARAM_MIN_CONTIG_LEN);
        nuclassembleDBworkflow.push_back(&PARAM_NUM_ITERATIONS);
        nuclassembleDBworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
        nuclassembleDBworkflow.push_back(&PARAM_DELETE_TMP_INC);
        nuclassembleDBworkflow.push_back(&PARAM_RUNNER);

        // easynuclassembleworkflow
        nuclassemblerworkflow = combineList(nuclassembleDBworkflow, createdb);

        // hybridassembleresults
        hybridassembleresults.push_back(&PARAM_MIN_SEQ_ID);
        hybridassembleresults.push_back(&PARAM_MAX_SEQ_LEN);
        hybridassembleresults.push_back(&PARAM_RESCORE_MODE);
        hybridassembleresults.push_back(&PARAM_THREADS);
        hybridassembleresults.push_back(&PARAM_V);

        // hybridassembledbworkflow
        hybridassembleDBworkflow = combineList(extractorfs, hybridassembleresults);
        hybridassembleDBworkflow = combineList(hybridassembleDBworkflow, nuclassembleDBworkflow);
        hybridassembleDBworkflow = combineList(hybridassembleDBworkflow, reduceredundancy);
        hybridassembleDBworkflow.push_back(&PARAM_CLUST_MIN_SEQ_ID_THR);
        hybridassembleDBworkflow.push_back(&PARAM_CLUST_C);

        // hybridassembledbworkflow special parameter: replace with MultiParam to make aa and nucl values independent
        hybridassembleDBworkflow = removeParameter(hybridassembleDBworkflow, PARAM_K);
        hybridassembleDBworkflow.push_back(&PARAM_MULTI_K);
        hybridassembleDBworkflow = removeParameter(hybridassembleDBworkflow, PARAM_MIN_SEQ_ID);
        hybridassembleDBworkflow.push_back(&PARAM_MULTI_MIN_SEQ_ID);
        hybridassembleDBworkflow = removeParameter(hybridassembleDBworkflow, PARAM_MIN_ALN_LEN);
        hybridassembleDBworkflow.push_back(&PARAM_MULTI_MIN_ALN_LEN);
        hybridassembleDBworkflow = removeParameter(hybridassembleDBworkflow, PARAM_NUM_ITERATIONS);
        hybridassembleDBworkflow.push_back(&PARAM_MULTI_NUM_ITERATIONS);

        // easynuclassembleworkflow
        hybridassemblerworkflow = combineList(hybridassembleDBworkflow, createdb);

        filterProteins = 1;
        deleteFilesInc = 1;
        proteinFilterThreshold = 0.2;
        clustSeqIdThr = 0.97;
        clustCovThr = 0.99;
        minContigLen = 1000;
        chopCycle = true;
        cycleCheck = true;

        multiNumIterations = MultiParam<int>(12,20);
        multiKmerSize = MultiParam<int>(14,22);
        multiAlnLenThr = MultiParam<int>(0,0);
        multiSeqIdThr = MultiParam<float>(0.97,0.97);

    }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
