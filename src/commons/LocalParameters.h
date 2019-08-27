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
    std::vector<MMseqsParameter *> cyclecheck;
    std::vector<MMseqsParameter *> assemblerworkflow;
    std::vector<MMseqsParameter *> nuclassemblerworkflow;
    std::vector<MMseqsParameter *> hybridassemblerworkflow;

    PARAMETER(PARAM_FILTER_PROTEINS)
    PARAMETER(PARAM_PROTEIN_FILTER_THRESHOLD)
    PARAMETER(PARAM_DELETE_TMP_INC)
    PARAMETER(PARAM_CYCLE_CHECK)
    PARAMETER(PARAM_CHOP_CYCLE)
    int filterProteins;
    int deleteFilesInc;
    float proteinFilterThreshold;
    bool cycleCheck;
    bool chopCycle;

private:
    LocalParameters() :
            Parameters(),
            PARAM_FILTER_PROTEINS(PARAM_FILTER_PROTEINS_ID,"--filter-proteins", "Filter Proteins", "filter proteins by a neural network [0,1]",typeid(int), (void *) &filterProteins, "^[0-1]{1}$"),
            PARAM_PROTEIN_FILTER_THRESHOLD(PARAM_PROTEIN_FILTER_THRESHOLD_ID,"--protein-filter-threshold", "Protein Filter Threshold", "filter proteins lower than threshold [0.0,1.0]",typeid(float), (void *) &proteinFilterThreshold, "^0(\\.[0-9]+)?|1(\\.0+)?$"),
            PARAM_DELETE_TMP_INC(PARAM_DELETE_TMP_INC_ID,"--delete-tmp-inc", "Delete temporary files incremental", "delete temporary files incremental [0,1]",typeid(int), (void *) &deleteFilesInc, "^[0-1]{1}$"),
            PARAM_CYCLE_CHECK(PARAM_CYCLE_CHECK_ID,"--cycyle-check", "Check for circular sequences", "Check for circular sequences",typeid(bool), (void *) &cycleCheck, ""),
            PARAM_CHOP_CYCLE(PARAM_CHOP_CYCLE_ID,"--chop-cycle", "Chop Cycle", "Remove superflous part of circular fragments",typeid(bool), (void *) &chopCycle, "")


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


        //cyclecheck
        cyclecheck.push_back(&PARAM_MAX_SEQ_LEN);
        cyclecheck.push_back(&PARAM_CHOP_CYCLE);
        cyclecheck.push_back(&PARAM_THREADS);

        // nucl assembler workflow
        nuclassemblerworkflow = combineList(rescorediagonal, kmermatcher);
        nuclassemblerworkflow = combineList(nuclassemblerworkflow, assembleresults);
        nuclassemblerworkflow = combineList(nuclassemblerworkflow, cyclecheck);
        nuclassemblerworkflow.push_back(&PARAM_CYCLE_CHECK);
        nuclassemblerworkflow.push_back(&PARAM_NUM_ITERATIONS);
        nuclassemblerworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
        nuclassemblerworkflow.push_back(&PARAM_RUNNER);

        // hybridassembleresults
        hybridassembleresults.push_back(&PARAM_MIN_SEQ_ID);
        hybridassembleresults.push_back(&PARAM_MAX_SEQ_LEN);
        hybridassembleresults.push_back(&PARAM_THREADS);
        hybridassembleresults.push_back(&PARAM_V);
        hybridassembleresults.push_back(&PARAM_RESCORE_MODE); //temporary added until assemble and nuclassemble use same rescoremode


        // hybridassemblerworkflow
        hybridassemblerworkflow = combineList(hybridassemblerworkflow, hybridassembleresults);
        hybridassemblerworkflow = combineList(hybridassemblerworkflow, nuclassemblerworkflow);
        hybridassemblerworkflow.push_back(&PARAM_CYCLE_CHECK);
        hybridassemblerworkflow.push_back(&PARAM_NUM_ITERATIONS);
        hybridassemblerworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
        hybridassemblerworkflow.push_back(&PARAM_RUNNER);


        filterProteins = 1;
        deleteFilesInc = 1;
        proteinFilterThreshold = 0.2;

    }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
