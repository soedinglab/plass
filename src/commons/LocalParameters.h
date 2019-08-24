#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include <Parameters.h>

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

    std::vector<MMseqsParameter *> assembleresults;
    std::vector<MMseqsParameter *> extractorfssubset;
    std::vector<MMseqsParameter *> filternoncoding;
    std::vector<MMseqsParameter *> hybridassembleresults;
    std::vector<MMseqsParameter *> checkcycle;
    std::vector<MMseqsParameter *> assemblerworkflow;
    std::vector<MMseqsParameter *> nuclassemblerworkflow;

    PARAMETER(PARAM_FILTER_PROTEINS)
    PARAMETER(PARAM_PROTEIN_FILTER_THRESHOLD)
    PARAMETER(PARAM_DELETE_TMP_INC)
    PARAMETER(PARAM_CHOP_CYCLE)
    int filterProteins;
    int deleteFilesInc;
    float proteinFilterThreshold;
    bool chopCycle;

private:
    LocalParameters() :
            Parameters(),
            PARAM_FILTER_PROTEINS(PARAM_FILTER_PROTEINS_ID,"--filter-proteins", "Filter Proteins", "filter proteins by a neural network [0,1]",typeid(int), (void *) &filterProteins, "^[0-1]{1}$"),
            PARAM_PROTEIN_FILTER_THRESHOLD(PARAM_PROTEIN_FILTER_THRESHOLD_ID,"--protein-filter-threshold", "Protein Filter Threshold", "filter proteins lower than threshold [0.0,1.0]",typeid(float), (void *) &proteinFilterThreshold, "^0(\\.[0-9]+)?|1(\\.0+)?$"),
            PARAM_DELETE_TMP_INC(PARAM_DELETE_TMP_INC_ID,"--delete-tmp-inc", "Delete temporary files incremental", "delete temporary files incremental [0,1]",typeid(int), (void *) &deleteFilesInc, "^[0-1]{1}$"),
            PARAM_CHOP_CYCLE(PARAM_CHOP_CYCLE_ID,"--chop-cycle", "Chop Cycle", "Remove superflous part of cycle fragment [0,1]",typeid(bool), (void *) &chopCycle, "")


    {
        // assembleresult
        assembleresults.push_back(&PARAM_MIN_SEQ_ID);
        assembleresults.push_back(&PARAM_MAX_SEQ_LEN);
        assembleresults.push_back(&PARAM_THREADS);
        assembleresults.push_back(&PARAM_V);
        assembleresults.push_back(&PARAM_RESCORE_MODE);

        extractorfssubset.push_back(&PARAM_TRANSLATION_TABLE);
        extractorfssubset.push_back(&PARAM_USE_ALL_TABLE_STARTS);
        extractorfssubset.push_back(&PARAM_THREADS);
        extractorfssubset.push_back(&PARAM_V);

        filternoncoding.push_back(&PARAM_PROTEIN_FILTER_THRESHOLD);
        filternoncoding.push_back(&PARAM_THREADS);
        filternoncoding.push_back(&PARAM_V);

        // assembler workflow
        assemblerworkflow = combineList(rescorediagonal, kmermatcher);
        assemblerworkflow = combineList(assemblerworkflow, extractorfs);
        assemblerworkflow = combineList(assemblerworkflow, assembleresults);
        assemblerworkflow = combineList(assemblerworkflow, filternoncoding);

        assemblerworkflow.push_back(&PARAM_FILTER_PROTEINS);
        assemblerworkflow.push_back(&PARAM_NUM_ITERATIONS);
        assemblerworkflow.push_back(&PARAM_DELETE_TMP_INC);
        assemblerworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
        assemblerworkflow.push_back(&PARAM_RUNNER);

        // nucl assembler workflow
        nuclassemblerworkflow = combineList(rescorediagonal, kmermatcher);
        nuclassemblerworkflow = combineList(nuclassemblerworkflow, assembleresults);
        nuclassemblerworkflow.push_back(&PARAM_NUM_ITERATIONS);
        nuclassemblerworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
        nuclassemblerworkflow.push_back(&PARAM_RUNNER);

        // hybridassembleresults
        hybridassembleresults = combineList(rescorediagonal, kmermatcher);
        hybridassembleresults.push_back(&PARAM_NUM_ITERATIONS);
        hybridassembleresults.push_back(&PARAM_REMOVE_TMP_FILES);
        hybridassembleresults.push_back(&PARAM_RUNNER);
        hybridassembleresults.push_back(&PARAM_RESCORE_MODE);

        //checkcycle
        checkcycle.push_back(&PARAM_MIN_SEQ_ID);
        checkcycle.push_back(&PARAM_K);
        checkcycle.push_back(&PARAM_MAX_SEQ_LEN);
        checkcycle.push_back(&PARAM_CHOP_CYCLE);

        filterProteins = 1;
        deleteFilesInc = 1;
        proteinFilterThreshold = 0.2;
        chopCycle = false;

    }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
