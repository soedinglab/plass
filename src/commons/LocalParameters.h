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

    std::vector<MMseqsParameter *> assembleresults;
    std::vector<MMseqsParameter *> extractorfssubset;
    std::vector<MMseqsParameter *> filternoncoding;
    std::vector<MMseqsParameter *> hybridassembleresults;
    std::vector<MMseqsParameter *> cyclecheck;
    std::vector<MMseqsParameter *> assemblerworkflow;
    std::vector<MMseqsParameter *> nuclassemblerworkflow;
    std::vector<MMseqsParameter *> hybridassemblerworkflow;
    std::vector<MMseqsParameter *> createhdb;

    PARAMETER(PARAM_FILTER_PROTEINS)
    PARAMETER(PARAM_PROTEIN_FILTER_THRESHOLD)
    PARAMETER(PARAM_DELETE_TMP_INC)
    PARAMETER(PARAM_NUM_AA_ITERATIONS)
    PARAMETER(PARAM_NUM_NUCL_ITERATIONS)
    PARAMETER(PARAM_MIN_CONTIG_LEN)
    PARAMETER(PARAM_CYCLE_CHECK)
    PARAMETER(PARAM_CHOP_CYCLE)
    int filterProteins;
    int deleteFilesInc;
    int numAAIterations;
    int numNuclIterations;
    int minContigLen;
    float proteinFilterThreshold;
    bool cycleCheck;
    bool chopCycle;

private:
    LocalParameters() :
            Parameters(),
            PARAM_FILTER_PROTEINS(PARAM_FILTER_PROTEINS_ID,"--filter-proteins", "Filter Proteins", "filter proteins by a neural network [0,1]",typeid(int), (void *) &filterProteins, "^[0-1]{1}$"),
            PARAM_PROTEIN_FILTER_THRESHOLD(PARAM_PROTEIN_FILTER_THRESHOLD_ID,"--protein-filter-threshold", "Protein Filter Threshold", "filter proteins lower than threshold [0.0,1.0]",typeid(float), (void *) &proteinFilterThreshold, "^0(\\.[0-9]+)?|1(\\.0+)?$"),
            PARAM_DELETE_TMP_INC(PARAM_DELETE_TMP_INC_ID,"--delete-tmp-inc", "Delete temporary files incremental", "delete temporary files incremental [0,1]",typeid(int), (void *) &deleteFilesInc, "^[0-1]{1}$"),
            PARAM_NUM_AA_ITERATIONS(PARAM_NUM_AA_ITERATIONS_ID, "--num-aa-iterations", "Number of assembly aa iteration","Number of assembly iterations performed on amino acid level [1, inf]",typeid(int),(void *) &numAAIterations, "^[1-9]{1}[0-9]*$"),
            PARAM_NUM_NUCL_ITERATIONS(PARAM_NUM_NUCL_ITERATIONS_ID, "--num-nucl-iterations", "Number of assembly nucl iteration","Number of assembly iterations performed on nucleotide level [1, inf]",typeid(int),(void *) &numNuclIterations, "^[1-9]{1}[0-9]*$"),
            PARAM_MIN_CONTIG_LEN(PARAM_MIN_CONTIG_LEN_ID, "--min-contig-len", "Minimum contig length", "Minimum length of assembled contig to output", typeid(int), (void *) &minContigLen, "^[1-9]{1}[0-9]*$"),
            PARAM_CYCLE_CHECK(PARAM_CYCLE_CHECK_ID,"--cycle-check", "Check for circular sequences", "Check for circular sequences (avoid infinite extension of circular or long repeated regions) ",typeid(bool), (void *) &cycleCheck, ""),
            PARAM_CHOP_CYCLE(PARAM_CHOP_CYCLE_ID,"--chop-cycle", "Chop Cycle", "Remove superfluous part of circular fragments",typeid(bool), (void *) &chopCycle, "")



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
        cyclecheck.push_back(&PARAM_V);

        //createhdb
        createhdb.push_back(&PARAM_COMPRESSED);
        createhdb.push_back(&PARAM_V);
        
        // nucl assembler workflow
        nuclassemblerworkflow.push_back(&PARAM_CYCLE_CHECK);
        nuclassemblerworkflow.push_back(&PARAM_MIN_CONTIG_LEN);
        nuclassemblerworkflow.push_back(&PARAM_NUM_ITERATIONS);
        nuclassemblerworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
        nuclassemblerworkflow.push_back(&PARAM_RUNNER);
        nuclassemblerworkflow = combineList(nuclassemblerworkflow, kmermatcher);
        nuclassemblerworkflow = combineList(nuclassemblerworkflow, rescorediagonal);
        nuclassemblerworkflow = combineList(nuclassemblerworkflow, assembleresults);
        nuclassemblerworkflow = combineList(nuclassemblerworkflow, cyclecheck);

        // hybridassembleresults
        hybridassembleresults.push_back(&PARAM_MIN_SEQ_ID);
        hybridassembleresults.push_back(&PARAM_MAX_SEQ_LEN);
        hybridassembleresults.push_back(&PARAM_THREADS);
        hybridassembleresults.push_back(&PARAM_V);
        hybridassembleresults.push_back(&PARAM_RESCORE_MODE);


        // hybridassemblerworkflow
        hybridassemblerworkflow.push_back(&PARAM_CYCLE_CHECK);
        hybridassemblerworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
        hybridassemblerworkflow.push_back(&PARAM_RUNNER);
        hybridassemblerworkflow = combineList(hybridassemblerworkflow, hybridassembleresults);
        hybridassemblerworkflow = combineList(hybridassemblerworkflow, nuclassemblerworkflow);

        std::remove(hybridassemblerworkflow.begin(), hybridassemblerworkflow.end(), &PARAM_NUM_ITERATIONS);
        hybridassemblerworkflow.push_back(&PARAM_NUM_AA_ITERATIONS);
        hybridassemblerworkflow.push_back(&PARAM_NUM_NUCL_ITERATIONS);


        filterProteins = 1;
        deleteFilesInc = 1;
        proteinFilterThreshold = 0.2;
        minContigLen = 1000;

    }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
