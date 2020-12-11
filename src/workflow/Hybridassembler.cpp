#include <cassert>
#include <typeinfo>

#include "CommandCaller.h"
#include "DBReader.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "MultiParam.h"
#include "Util.h"

namespace hybridassembler {
#include "easyassembler.sh.h"
}

void setEasyHybridAssemblerWorkflowDefaults(LocalParameters *p) {
    //p->createdbMode = Parameters::SEQUENCE_SPLIT_MODE_SOFT;

    p->multiNumIterations = MultiParam<int>(5,5);
    p->multiKmerSize = MultiParam<int>(14,22);
    p->multiSeqIdThr = MultiParam<float>(0.97,0.99);
    p->alphabetSize = MultiParam<int>(13,5);

    p->orfMinLength = 45;
    p->covThr = 0.00;
    p->evalThr = 0.00001;
    p->maskMode = 0;
    p->kmersPerSequence = 60;
    p->kmersPerSequenceScale = 0.1;
    p->spacedKmer = false;
    p->ignoreMultiKmer = true;
    p->includeOnlyExtendable = true;
    p->rescoreMode = Parameters::RESCORE_MODE_GLOBAL_ALIGNMENT;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
    p->maxSeqLen = 200000; //TODO
    p->cycleCheck = true;
    p->chopCycle = true;

    //cluster defaults
    p->clustSeqIdThr = 0.97;
    p->clustCovThr = 0.99;
    p->clusteringMode = 2;
    p->gapOpen = 5;
    p->gapExtend = 2;
    p->zdrop = 200;
}

void setEasyHybridAssemblerMustPassAlong(LocalParameters *p) {
    p->PARAM_MULTI_NUM_ITERATIONS.wasSet = true;
    p->PARAM_MULTI_K.wasSet = true;
    p->PARAM_MULTI_MIN_SEQ_ID.wasSet = true;
    p->PARAM_ALPH_SIZE.wasSet = true;

    p->PARAM_ORF_MIN_LENGTH.wasSet = true;
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
    p->PARAM_MAX_SEQ_LEN.wasSet = true;
    p->PARAM_CYCLE_CHECK.wasSet = true;
    p->PARAM_CHOP_CYCLE.wasSet = true;

    p->PARAM_CLUST_MIN_SEQ_ID_THR.wasSet = true;
    p->PARAM_CLUST_C.wasSet = true;
    p->PARAM_CLUSTER_MODE.wasSet = true;
    p->PARAM_GAP_OPEN.wasSet = true;
    p->PARAM_GAP_EXTEND.wasSet = true;
    p->PARAM_ZDROP.wasSet = true;
}

int easyhybridassembler(int argc, const char **argv, const Command &command) {

    LocalParameters &par = LocalParameters::getLocalInstance();

    par.overrideParameterDescription(par.PARAM_MULTI_MIN_SEQ_ID, "Overlap sequence identity threshold (range 0.0-1.0)", NULL, 0);
    par.overrideParameterDescription(par.PARAM_E, "Extend sequences if the E-value is below (range 0.0-inf)", NULL, 0);
    par.overrideParameterDescription(par.PARAM_GAP_OPEN, "Gap open cost (only for clustering)", NULL, 0);
    par.overrideParameterDescription(par.PARAM_GAP_EXTEND, "Gap extend cost (only for clustering)", NULL, 0);
    par.overrideParameterDescription(par.PARAM_ZDROP, "Maximal allowed difference between score values before alignment is truncated (only for clustering)", NULL, 0);
    par.overrideParameterDescription(par.PARAM_COMPRESSED, "Use compressed database format for temporary files", NULL, 0);

    // reorganize parameters to clarify usage context
    par.PARAM_ALPH_SIZE.replaceCategory(MMseqsParameter::COMMAND_CLUSTLINEAR | MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_INCLUDE_IDENTITY.replaceCategory(MMseqsParameter::COMMAND_CLUSTLINEAR | MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MASK_RESIDUES.replaceCategory(MMseqsParameter::COMMAND_CLUSTLINEAR | MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MASK_LOWER_CASE.replaceCategory(MMseqsParameter::COMMAND_CLUSTLINEAR| MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_SPACED_KMER_MODE.replaceCategory(MMseqsParameter::COMMAND_CLUSTLINEAR | MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_SPLIT_MEMORY_LIMIT.replaceCategory(MMseqsParameter::COMMAND_CLUSTLINEAR | MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_SUB_MAT.replaceCategory(MMseqsParameter::COMMAND_CLUSTLINEAR | MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.replaceCategory(MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT);

    par.PARAM_COV_MODE.replaceCategory(MMseqsParameter::COMMAND_CLUST | MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_GAP_EXTEND.replaceCategory(MMseqsParameter::COMMAND_CLUST | MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_GAP_OPEN.replaceCategory(MMseqsParameter::COMMAND_CLUST | MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ZDROP.replaceCategory(MMseqsParameter::COMMAND_CLUST | MMseqsParameter::COMMAND_EXPERT);

    // hide this parameters (changing them would lead to unexpected behavior)
    par.PARAM_C.replaceCategory(MMseqsParameter::COMMAND_HIDDEN);
    par.PARAM_CONTIG_END_MODE.replaceCategory(MMseqsParameter::COMMAND_HIDDEN);
    par.PARAM_CONTIG_START_MODE.replaceCategory(MMseqsParameter::COMMAND_HIDDEN);
    par.PARAM_ORF_MAX_GAP.replaceCategory(MMseqsParameter::COMMAND_HIDDEN);
    par.PARAM_ORF_START_MODE.replaceCategory(MMseqsParameter::COMMAND_HIDDEN);
    par.PARAM_WRAPPED_SCORING.replaceCategory(MMseqsParameter::COMMAND_HIDDEN);

    // expert
    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_CLUSTER_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_DB_TYPE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ID_OFFSET.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_KMER_PER_SEQ.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_KMER_PER_SEQ_SCALE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_INCLUDE_ONLY_EXTENDABLE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ORF_MAX_LENGTH.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ORF_MIN_LENGTH.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ORF_FORWARD_FRAMES.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ORF_REVERSE_FRAMES.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_SEQ_ID_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_SORT_RESULTS.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_TRANSLATION_TABLE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_TRANSLATE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_USE_ALL_TABLE_STARTS.addCategory(MMseqsParameter::COMMAND_EXPERT);

    setEasyHybridAssemblerWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    setEasyHybridAssemblerMustPassAlong(&par);

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
    cmd.addVariable("ASSEMBLY_PAR", par.createParameterString(par.hybridassembleDBworkflow, true).c_str());
    cmd.addVariable("ASSEMBLY_MODULE", "hybridassembledb");
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    std::string program = tmpDir + "/easyassembler.sh";
    FileUtil::writeFile(program, hybridassembler::easyassembler_sh, hybridassembler::easyassembler_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return 0;
}
