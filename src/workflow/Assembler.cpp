#include "Parameters.h"
#include <string>
#include <cassert>
#include <Util.h>
#include "assembler.sh.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"


void setAssemblerWorkflowDefaults(Parameters *p) {
    p->spacedKmer = false;
    p->maskMode = 1;
    p->covThr = 0.0;
    p->evalThr = 0.00001;
    p->seqIdThr = 0.9;
    p->alphabetSize = 21;
    p->kmersPerSequence = 60;
    p->numIterations = 12;
    p->maskMode = 0;
    p->includeOnlyExtendable = true;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
}

int assembler(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setAssemblerWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, 3);
    if(FileUtil::directoryExists(par.db3.c_str())==false){
        Debug(Debug::ERROR) << "Tmp " << par.db3 << " folder does not exist or is not a directory.\n";
        EXIT(EXIT_FAILURE);
    }
    CommandCaller cmd;
    if(par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("NUM_IT", SSTR(par.numIterations).c_str());
    // save some values to restore them later
    size_t alphabetSize = par.alphabetSize;
    size_t kmerSize = par.kmerSize;
    // # 1. Finding exact $k$-mer matches.
    int baseKmerSize = 14;
    bool kmerSizeWasSet = false;
    bool alphabetSizeWasSet = false;
    for (size_t i = 0; i < par.assemblerworkflow.size(); i++) {
        if (par.assemblerworkflow[i].uniqid == par.PARAM_K.uniqid && par.assemblerworkflow[i].wasSet) {
            kmerSizeWasSet = true;
        }
        if (par.assemblerworkflow[i].uniqid == par.PARAM_ALPH_SIZE.uniqid && par.assemblerworkflow[i].wasSet) {
            alphabetSizeWasSet = true;
        }
    }
    if(kmerSizeWasSet==false){
        par.kmerSize = baseKmerSize;
    }
    if(alphabetSizeWasSet == false){
        par.alphabetSize = Parameters::CLUST_LINEAR_DEFAULT_ALPH_SIZE;
    }
    cmd.addVariable("KMER_PER_SEQ", SSTR(par.kmersPerSequence).c_str());
    std::vector<MMseqsParameter> kmerMatcherrWithoutKmerPerseq;
    for(size_t i = 0; i < par.kmermatcher.size(); i++){
        if(par.kmermatcher[i].uniqid != par.PARAM_KMER_PER_SEQ.uniqid ){
            kmerMatcherrWithoutKmerPerseq.push_back(par.kmermatcher[i]);
        }
    }

    for(size_t i = 0; i < par.numIterations; i++){
        std::string key = "KMERMATCHER"+SSTR(i)+"_PAR";
        par.hashShift = i+1;
        cmd.addVariable(key.c_str(), par.createParameterString(kmerMatcherrWithoutKmerPerseq).c_str());
    }
    par.alphabetSize = alphabetSize;
    par.kmerSize = kmerSize;
    // # 2. Hamming distance pre-clustering
    par.filterHits = false;
    par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
    cmd.addVariable("UNGAPPED_ALN_PAR", par.createParameterString(par.rescorediagonal).c_str());
    cmd.addVariable("ASSEMBLE_RESULT_PAR", par.createParameterString(par.assembleresults).c_str());
    FileUtil::writeFile(par.db3 + "/assembler.sh", assembler_sh, assembler_sh_len);
    std::string program(par.db3 + "/assembler.sh");
    cmd.execProgram(program.c_str(), par.filenames);
    return 0;
}
