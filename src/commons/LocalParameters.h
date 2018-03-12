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

    std::vector<MMseqsParameter> assembleresults;
    std::vector<MMseqsParameter> hybridassembleresults;
    std::vector<MMseqsParameter> assemblerworkflow;

private:
    LocalParameters() : Parameters() {
        // assembleresult
        assembleresults.push_back(PARAM_MIN_SEQ_ID);
        assembleresults.push_back(PARAM_V);

        // assembler workflow
        assemblerworkflow = combineList(rescorediagonal, kmermatcher);
        assemblerworkflow = combineList(assemblerworkflow, assembleresults);
        assemblerworkflow.push_back(PARAM_NUM_ITERATIONS);
        assemblerworkflow.push_back(PARAM_REMOVE_TMP_FILES);
        assemblerworkflow.push_back(PARAM_RUNNER);

        //
        hybridassembleresults = combineList(rescorediagonal, kmermatcher);
        hybridassembleresults.push_back(PARAM_NUM_ITERATIONS);
        hybridassembleresults.push_back(PARAM_REMOVE_TMP_FILES);
        hybridassembleresults.push_back(PARAM_RUNNER);
    }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
