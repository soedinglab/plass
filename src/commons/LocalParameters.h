#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include <Parameters.h>
#include <algorithm>

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

    PARAMETER(PARAM_FILTER_PROTEINS)
    PARAMETER(PARAM_PROTEIN_FILTER_THRESHOLD)
    PARAMETER(PARAM_DELETE_TMP_INC)
    PARAMETER(PARAM_NUM_PROT_ITERATIONS)
    PARAMETER(PARAM_NUM_NUCL_ITERATIONS)
    PARAMETER(PARAM_MIN_CONTIG_LEN)
    PARAMETER(PARAM_CLUST_THR)
    PARAMETER(PARAM_CYCLE_CHECK)
    PARAMETER(PARAM_CHOP_CYCLE)
    int filterProteins;
    int deleteFilesInc;
    int numProtIterations;
    int numNuclIterations;
    int minContigLen;
    float clustThr;
    float proteinFilterThreshold;
    bool cycleCheck;
    bool chopCycle;

private:
    LocalParameters() :
            Parameters(),
            PARAM_FILTER_PROTEINS(PARAM_FILTER_PROTEINS_ID,"--filter-proteins", "Filter Proteins", "filter proteins by a neural network [0,1]",typeid(int), (void *) &filterProteins, "^[0-1]{1}$"),
            PARAM_PROTEIN_FILTER_THRESHOLD(PARAM_PROTEIN_FILTER_THRESHOLD_ID,"--protein-filter-threshold", "Protein Filter Threshold", "filter proteins lower than threshold [0.0,1.0]",typeid(float), (void *) &proteinFilterThreshold, "^0(\\.[0-9]+)?|1(\\.0+)?$"),
            PARAM_DELETE_TMP_INC(PARAM_DELETE_TMP_INC_ID,"--delete-tmp-inc", "Delete temporary files incremental", "delete temporary files incremental [0,1]",typeid(int), (void *) &deleteFilesInc, "^[0-1]{1}$"),
            PARAM_NUM_PROT_ITERATIONS(PARAM_NUM_PROT_ITERATIONS_ID, "--num-prot-iterations", "Number of assembly prot iteration","Number of assembly iterations performed on amino acid level [1, inf]",typeid(int),(void *) &numProtIterations, "^[1-9]{1}[0-9]*$"),
            PARAM_NUM_NUCL_ITERATIONS(PARAM_NUM_NUCL_ITERATIONS_ID, "--num-nucl-iterations", "Number of assembly nucl iteration","Number of assembly iterations performed on nucleotide level [1, inf]",typeid(int),(void *) &numNuclIterations, "^[1-9]{1}[0-9]*$"),
            PARAM_MIN_CONTIG_LEN(PARAM_MIN_CONTIG_LEN_ID, "--min-contig-len", "Minimum contig length", "Minimum length of assembled contig to output", typeid(int), (void *) &minContigLen, "^[1-9]{1}[0-9]*$"),
            PARAM_CLUST_THR(PARAM_CLUST_THR_ID,"--clust-thr", "Clustering threshold","Threshold to reduce redundancy in assembly by applying the linclust algorithm (clustering threshold) (range 0.0-1.0)",typeid(float), (void *) &clustThr, "^0(\\.[0-9]+)?|1(\\.0+)?$"),
            PARAM_CYCLE_CHECK(PARAM_CYCLE_CHECK_ID,"--cycle-check", "Check for circular sequences", "Check for circular sequences (avoid infinite extension of circular or long repeated regions) ",typeid(bool), (void *) &cycleCheck, "", MMseqsParameter::COMMAND_MISC | MMseqsParameter::COMMAND_EXPERT),
            PARAM_CHOP_CYCLE(PARAM_CHOP_CYCLE_ID,"--chop-cycle", "Chop Cycle", "Remove superfluous part of circular fragments",typeid(bool), (void *) &chopCycle, "", MMseqsParameter::COMMAND_MISC | MMseqsParameter::COMMAND_EXPERT)



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
        reduceredundancy.push_back(&PARAM_CLUSTER_MODE);

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
        nuclassembleDBworkflow.push_back(&PARAM_CYCLE_CHECK);
        nuclassembleDBworkflow.push_back(&PARAM_MIN_CONTIG_LEN);
        nuclassembleDBworkflow.push_back(&PARAM_CLUST_THR);
        nuclassembleDBworkflow.push_back(&PARAM_NUM_ITERATIONS);
        nuclassembleDBworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
        nuclassembleDBworkflow.push_back(&PARAM_DELETE_TMP_INC);
        nuclassembleDBworkflow.push_back(&PARAM_RUNNER);
        nuclassembleDBworkflow = combineList(nuclassembleDBworkflow, kmermatcher);
        nuclassembleDBworkflow = combineList(nuclassembleDBworkflow, rescorediagonal);
        nuclassembleDBworkflow = combineList(nuclassembleDBworkflow, assembleresults);
        nuclassembleDBworkflow = combineList(nuclassembleDBworkflow, cyclecheck);

        // easynuclassembleworkflow
        nuclassemblerworkflow = combineList(nuclassembleDBworkflow, createdb);

        // hybridassembleresults
        hybridassembleresults.push_back(&PARAM_MIN_SEQ_ID);
        hybridassembleresults.push_back(&PARAM_MAX_SEQ_LEN);
        hybridassembleresults.push_back(&PARAM_THREADS);
        hybridassembleresults.push_back(&PARAM_V);
        hybridassembleresults.push_back(&PARAM_RESCORE_MODE);


        // hybridassemblerworkflow
        hybridassembleDBworkflow.push_back(&PARAM_CYCLE_CHECK);
        hybridassembleDBworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
        hybridassembleDBworkflow.push_back(&PARAM_RUNNER);
        hybridassembleDBworkflow = combineList(hybridassembleDBworkflow, extractorfs);
        hybridassembleDBworkflow = combineList(hybridassembleDBworkflow, hybridassembleresults);
        hybridassembleDBworkflow = combineList(hybridassembleDBworkflow, nuclassembleDBworkflow);
        hybridassembleDBworkflow = combineList(hybridassembleDBworkflow, reduceredundancy);

        std::remove(hybridassembleDBworkflow.begin(), hybridassembleDBworkflow.end(), &PARAM_NUM_ITERATIONS);
        hybridassembleDBworkflow.push_back(&PARAM_NUM_PROT_ITERATIONS);
        hybridassembleDBworkflow.push_back(&PARAM_NUM_NUCL_ITERATIONS);

        // easynuclassembleworkflow
        hybridassemblerworkflow = combineList(hybridassembleDBworkflow, createdb);

        filterProteins = 1;
        deleteFilesInc = 1;
        proteinFilterThreshold = 0.2;
        clustThr = 0.97;
        minContigLen = 1000;
        chopCycle = false;
        cycleCheck = true;

    }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
