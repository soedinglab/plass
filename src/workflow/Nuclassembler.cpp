#include <mmseqs/src/commons/DBReader.h>
#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "nuclassembler.sh.h"

void setNuclAssemblerWorkflowDefaults(LocalParameters *p) {
    p->spacedKmer = false;
    p->maskMode = 0;
    p->covThr = 0.0;
    p->evalThr = 0.00001;
    p->seqIdThr = 0.9;
    p->kmersPerSequence = 60;
    p->numIterations = 12;
    p->includeOnlyExtendable = true;
    p->alphabetSize = 5;
    p->kmerSize = 22;
    p->skipNRepeatKmer = 8;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
}

int nuclassembler(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.overrideParameterDescription((Command &)command, par.PARAM_COV_MODE.uniqid, NULL, NULL, par.PARAM_COV_MODE.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_C.uniqid, NULL, NULL, par.PARAM_C.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_MIN_SEQ_ID.uniqid, "overlap sequence identity threshold", NULL,  par.PARAM_MIN_SEQ_ID.category);
    par.overrideParameterDescription((Command &)command, par.PARAM_SEQ_ID_MODE.uniqid, NULL, NULL,  par.PARAM_SEQ_ID_MODE.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_RESCORE_MODE.uniqid, NULL, NULL,  par.PARAM_RESCORE_MODE.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_INCLUDE_ONLY_EXTENDABLE.uniqid, NULL, NULL,  par.PARAM_INCLUDE_ONLY_EXTENDABLE.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_KMER_PER_SEQ.uniqid, NULL, NULL,  par.PARAM_KMER_PER_SEQ.category | MMseqsParameter::COMMAND_EXPERT);

    setNuclAssemblerWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, 2, true, Parameters::PARSE_VARIADIC);

    CommandCaller cmd;
    // paired end reads
    if((par.filenames.size() - 2) % 2  == 0 ){
        std::string readStr = "";
        for(size_t i = 0; i < (par.filenames.size() -2)/2; i++){
            readStr.append(par.filenames[i*2]);
            readStr.append(" ");
            readStr.append(par.filenames[i*2 + 1]);
        }
        cmd.addVariable("PAIRED_END", "1");
        cmd.addVariable("READ_FILES", readStr.c_str());
    }else{
        if(par.filenames.size() == 3){
            cmd.addVariable("READ_FILES", par.filenames[0].c_str());
        }else{
            Debug(Debug::ERROR) << "Read input parameters are wrong. \n";
            Debug(Debug::ERROR) << "For paired end input use READSETA_1.fastq READSETA_2.fastq ... OUTPUT.fasta  \n";
            Debug(Debug::ERROR) << "For single input use READSET.fast(q|a) OUTPUT.fasta  \n";
            EXIT(EXIT_FAILURE);
        }
    }
    cmd.addVariable("OUT_FILE", par.filenames[par.filenames.size() - 2].c_str());

    std::string tmpPath = par.filenames[par.filenames.size() - 1];
    if (FileUtil::directoryExists(tmpPath.c_str()) == false){
        Debug(Debug::INFO) << "Temporary folder " << tmpPath << " does not exist or is not a directory.\n";
        if (FileUtil::makeDir(tmpPath.c_str()) == false){
            Debug(Debug::WARNING) << "Could not crate tmp folder " << tmpPath << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created directory " << tmpPath << "\n";
        }
    }
    size_t hash = par.hashParameter(par.filenames, par.searchworkflow);
    std::string tmpDir(tmpPath + "/" + SSTR(hash));
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::WARNING) << "Could not create sub folder in temporary directory " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);
    FileUtil::symlinkAlias(tmpDir, "latest");
    cmd.addVariable("TMP_PATH", tmpDir.c_str());


    if (par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("NUM_IT", SSTR(par.numIterations).c_str());
    // # 1. Finding exact $k$-mer matches.
    cmd.addVariable("KMERMATCHER_PAR", par.createParameterString(par.kmermatcher).c_str());

    // # 2. Hamming distance pre-clustering
    par.filterHits = false;
    par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
    cmd.addVariable("UNGAPPED_ALN_PAR", par.createParameterString(par.rescorediagonal).c_str());
    cmd.addVariable("ASSEMBLE_RESULT_PAR", par.createParameterString(par.assembleresults).c_str());

    FileUtil::writeFile(tmpDir + "/assembler.sh", nuclassembler_sh, nuclassembler_sh_len);
    std::string program(tmpDir + "/assembler.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
