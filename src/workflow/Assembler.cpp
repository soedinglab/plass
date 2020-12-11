#include <cassert>

#include "CommandCaller.h"
#include "DBReader.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "Util.h"

#include "easyassembler.sh.h"

void setEasyAssemblerWorkflowDefaults(LocalParameters *p) {
    p->spacedKmer = false;
    p->maskMode = 0;
    p->covThr = 0.0;
    p->evalThr = 0.00001;
    p->seqIdThr = 0.9;
    p->kmersPerSequence = 60;
    p->numIterations = 12;
    p->alphabetSize = 13;
    p->kmerSize = 14;
    p->orfMinLength = 45;
    p->ignoreMultiKmer = true;
    p->includeOnlyExtendable = true;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
    p->rescoreMode = Parameters::RESCORE_MODE_GLOBAL_ALIGNMENT;
}

void setEasyAssemblerMustPassAlong(LocalParameters *p) {
    p->PARAM_SPACED_KMER_MODE.wasSet = true;
    p->PARAM_MASK_RESIDUES.wasSet = true;
    p->PARAM_C.wasSet = true;
    p->PARAM_E.wasSet = true;
    p->PARAM_MIN_SEQ_ID.wasSet = true;
    p->PARAM_KMER_PER_SEQ.wasSet = true;
    p->PARAM_NUM_ITERATIONS.wasSet = true;
    p->PARAM_ALPH_SIZE.wasSet = true;
    p->PARAM_K.wasSet = true;
    p->PARAM_ORF_MIN_LENGTH.wasSet = true;
    p->PARAM_IGNORE_MULTI_KMER.wasSet = true;
    p->PARAM_INCLUDE_ONLY_EXTENDABLE.wasSet = true;
    p->PARAM_ALIGNMENT_MODE.wasSet = true;
    p->PARAM_RESCORE_MODE.wasSet = true;
}

int easyassembler(int argc, const char **argv, const Command &command) {

    LocalParameters &par = LocalParameters::getLocalInstance();

    par.overrideParameterDescription(par.PARAM_MIN_SEQ_ID, "Overlap sequence identity threshold [0.0, 1.0]", NULL, 0);
    par.overrideParameterDescription(par.PARAM_NUM_ITERATIONS, "Number of assembly iterations [1, inf]", NULL,  0);
    par.overrideParameterDescription(par.PARAM_E, "Extend sequences if the E-value is below [0.0, inf]", NULL,  0);

    par.PARAM_COV_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_C.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ID_OFFSET.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_CONTIG_END_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_CONTIG_START_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ORF_MAX_GAP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ORF_START_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ORF_FORWARD_FRAMES.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ORF_REVERSE_FRAMES.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_SEQ_ID_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_INCLUDE_ONLY_EXTENDABLE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_KMER_PER_SEQ.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_SORT_RESULTS.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_TRANSLATION_TABLE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_USE_ALL_TABLE_STARTS.addCategory(MMseqsParameter::COMMAND_EXPERT);


    setEasyAssemblerWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    setEasyAssemblerMustPassAlong(&par);

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

    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);

    char *p = realpath(tmpDir.c_str(), NULL);
    if (p == NULL) {
        Debug(Debug::ERROR) << "Could not get real path of " << tmpDir << "!\n";
        EXIT(EXIT_FAILURE);
    }
    cmd.addVariable("TMP_PATH", p);
    free(p);
    par.filenames.pop_back();
    cmd.addVariable("OUT_FILE", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("RUNNER", par.runner.c_str());

    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createdb).c_str());
    cmd.addVariable("ASSEMBLY_PAR", par.createParameterString(par.assembleDBworkflow, true).c_str());
    cmd.addVariable("ASSEMBLY_MODULE", "assembledb");
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    std::string program = tmpDir + "/easyassembler.sh";
    FileUtil::writeFile(program, easyassembler_sh, easyassembler_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return 0;
}
