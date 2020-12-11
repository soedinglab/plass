#include <cassert>

#include "CommandCaller.h"
#include "DBReader.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "Util.h"

namespace nuclassembler {
#include "easyassembler.sh.h"
}

void setEasyNuclAssemblerWorkflowDefaults(LocalParameters *p) {
    p->addBacktrace = false;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
    p->alphabetSize = 5;
    p->chopCycle = true;
    p->covThr = 0.0;
    p->cycleCheck = true;
    p->evalThr = 0.00001;
    p->ignoreMultiKmer = true;
    p->includeOnlyExtendable = true;
    p->kmerSize = 22;
    p->kmersPerSequence = 60;
    p->kmersPerSequenceScale = 0.1;
    p->maskMode = 0;
    p->numIterations = 12;
    p->rescoreMode = Parameters::RESCORE_MODE_GLOBAL_ALIGNMENT;
    p->seqIdThr = 0.99;
    p->spacedKmer = false;
}

void setEasyNuclAssemblerMustPassAlong(LocalParameters *p) {
    p->PARAM_NUM_ITERATIONS.wasSet = true;
    p->PARAM_K.wasSet = true;
    p->PARAM_MIN_SEQ_ID.wasSet = true;
    p->PARAM_ALPH_SIZE.wasSet = true;
    p->PARAM_C.wasSet = true;
    p->PARAM_E.wasSet = true;
    p->PARAM_MASK_RESIDUES.wasSet = true;
    p->PARAM_KMER_PER_SEQ.wasSet = true;
    p->PARAM_KMER_PER_SEQ_SCALE.wasSet = true;
    p->PARAM_SPACED_KMER_MODE.wasSet = true;
    p->PARAM_IGNORE_MULTI_KMER.wasSet = true;
    p->PARAM_INCLUDE_ONLY_EXTENDABLE.wasSet = true;
    p->PARAM_RESCORE_MODE.wasSet = true;
    p->PARAM_ALIGNMENT_MODE.wasSet = true;
    p->PARAM_CYCLE_CHECK.wasSet = true;
    p->PARAM_CHOP_CYCLE.wasSet = true;
}

int easynuclassembler(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();

    par.overrideParameterDescription(par.PARAM_MIN_SEQ_ID, "Overlap sequence identity threshold (range 0.0-1.0)", NULL, 0);
    par.overrideParameterDescription(par.PARAM_NUM_ITERATIONS, "Number of assembly iterations (range 1-inf)", NULL, 0);
    par.overrideParameterDescription(par.PARAM_E, "Extend sequences if the E-value is below (range 0.0-inf)", NULL, 0);
    par.overrideParameterDescription(par.PARAM_COMPRESSED, "Use compressed database format for temporary files", NULL, 0);



    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_COV_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_C.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_SEQ_ID_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_INCLUDE_ONLY_EXTENDABLE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_KMER_PER_SEQ.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_SORT_RESULTS.addCategory(MMseqsParameter::COMMAND_EXPERT);

    setEasyNuclAssemblerWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    setEasyNuclAssemblerMustPassAlong(&par);

    CommandCaller cmd;
    if(par.filenames.size() < 3) {
        Debug(Debug::ERROR) << "Too few input files provided.\n";
        return EXIT_FAILURE;
    }
    else if ((par.filenames.size() - 2) % 2 == 0) {
        cmd.addVariable("PAIRED_END", "1"); // paired end reads
    } else {
        if (par.filenames.size() != 3) {
            Debug(Debug::ERROR) << "Too many input files provided.\n";
            Debug(Debug::ERROR) << "For paired-end input provide READSETA_1.fastq READSETA_2.fastq ... OUTPUT.fasta tmpDir\n";
            Debug(Debug::ERROR) << "For single input use READSET.fast(q|a) OUTPUT.fasta tmpDir\n";
            return EXIT_FAILURE;
        }
        cmd.addVariable("PAIRED_END", NULL); // single end reads
    }

    std::string tmpPath = par.filenames.back();
    if (FileUtil::directoryExists(tmpPath.c_str()) == false) {
        Debug(Debug::INFO) << "Temporary folder " << tmpPath << " does not exist or is not a directory.\n";
        if (FileUtil::makeDir(tmpPath.c_str()) == false) {
            Debug(Debug::ERROR) << "Could not crate tmp folder " << tmpPath << ".\n";
            return EXIT_FAILURE;
        } else {
            Debug(Debug::INFO) << "Created directory " << tmpPath << "\n";
        }
    }
    size_t hash = par.hashParameter(command.databases, par.filenames, *command.params);
    std::string tmpDir = tmpPath + "/" + SSTR(hash);
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Could not create sub folder in temporary directory " << tmpDir << ".\n";
            return EXIT_FAILURE;
        }
    }
    par.filenames.pop_back();
    FileUtil::symlinkAlias(tmpDir, "latest");
    char *p = realpath(tmpDir.c_str(), NULL);
    if (p == NULL) {
        Debug(Debug::ERROR) << "Could not get real path of " << tmpDir << "!\n";
        EXIT(EXIT_FAILURE);
    }
    cmd.addVariable("TMP_PATH", p);
    free(p);

    cmd.addVariable("OUT_FILE", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("RUNNER", par.runner.c_str());

    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createdb).c_str());
    cmd.addVariable("ASSEMBLY_PAR", par.createParameterString(par.nuclassemblerworkflow, true).c_str());
    cmd.addVariable("ASSEMBLY_MODULE", "nuclassembledb");
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    std::string program = tmpDir + "/easyassemble.sh";
    FileUtil::writeFile(program, nuclassembler::easyassembler_sh, nuclassembler::easyassembler_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return 0;
}
