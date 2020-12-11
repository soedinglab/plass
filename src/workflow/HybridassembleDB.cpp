#include "DBReader.h"
#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"

#include "hybridassembledb.sh.h"

void setHybridAssemblerWorkflowDefaults(LocalParameters *p) {

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
    p->addBacktrace = false;
    p->rescoreMode = Parameters::RESCORE_MODE_GLOBAL_ALIGNMENT;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
    p->maxSeqLen = 200000; //TODO
    p->cycleCheck = true;
    p->chopCycle = true;

    //cluster defaults
    p->covThr = 0.99;
    p->clusteringMode = 2;
    p->gapOpen = 5;
    p->gapExtend = 2;
    p->zdrop = 200;
}

int hybridassembledb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setHybridAssemblerWorkflowDefaults(&par);

    par.overrideParameterDescription(par.PARAM_MULTI_MIN_SEQ_ID, "Overlap sequence identity threshold (range 0.0-1.0)", NULL, 0);
    par.overrideParameterDescription(par.PARAM_E, "Extend sequences if the E-value is below (range 0.0-inf)", NULL, 0);
    par.overrideParameterDescription(par.PARAM_GAP_OPEN, "Gap open cost (only for clustering", NULL, 0);
    par.overrideParameterDescription(par.PARAM_GAP_EXTEND, "Gap extend cost (only for clustering)", NULL, 0);
    par.overrideParameterDescription(par.PARAM_ZDROP, "Maximal allowed difference between score values before alignment is truncated (only for clustering)", NULL, 0);

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

    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    CommandCaller cmd;


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
    par.filenames.pop_back();
    free(p);

    cmd.addVariable("OUT_FILE", par.filenames.back().c_str());
    par.filenames.pop_back();

    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("REMOVE_INCREMENTAL_TMP", par.deleteFilesInc ? "TRUE" : NULL);
    cmd.addVariable("RUNNER", par.runner.c_str());

    // set values for protein level assembly
    par.numIterations = par.multiNumIterations.aminoacids;
    par.kmerSize = par.multiKmerSize.aminoacids;
    par.seqIdThr = par.multiSeqIdThr.aminoacids;
    par.alnLenThr = par.multiAlnLenThr.aminoacids;
    cmd.addVariable("NUM_IT", SSTR(par.numIterations).c_str());

    // # 0. Extract ORFs
    // --orf-start-mode 0 --min-length 45 --max-gaps 0
    par.orfStartMode = 0;
    par.orfMaxGaps = 0;
    cmd.addVariable("EXTRACTORFS_LONG_PAR", par.createParameterString(par.extractorfs).c_str());

    // --contig-start-mode 1 --contig-end-mode 0 --orf-start-mode 0 --min-length 30 --max-length 45 --max-gaps 0
    par.contigStartMode = 1;
    par.contigEndMode = 0;
    par.orfStartMode = 0;
    par.orfMaxLength = par.orfMinLength;
    par.orfMinLength = std::min(par.orfMinLength, 20);
    par.orfMaxGaps = 0;
    cmd.addVariable("EXTRACTORFS_START_PAR", par.createParameterString(par.extractorfs).c_str());

    // force parameters for assembly steps
    par.covThr = 0.0;
    par.seqIdMode = 0;
    par.includeOnlyExtendable = true;

    // # 1. Finding exact $k$-mer matches.
    cmd.addVariable("KMERMATCHER_PAR", par.createParameterString(par.kmermatcher).c_str());

    // # 2. Rescore diagonal
    par.filterHits = false;
    bool addBacktrace = par.addBacktrace;
    par.addBacktrace = true;
    cmd.addVariable("UNGAPPED_ALN_PAR", par.createParameterString(par.rescorediagonal).c_str());

    // # 3. Assembly: Extend by left and right extension
    par.seqIdThr = par.multiSeqIdThr.nucleotides;
    cmd.addVariable("ASSEMBLE_RESULT_PAR", par.createParameterString(par.assembleresults).c_str());

    // set mandatory values for nucleotide level assembly step when calling nucleassemble from hybridassemble
    par.numIterations = par.multiNumIterations.nucleotides;
    par.kmerSize = par.multiKmerSize.nucleotides;
    par.seqIdThr = par.multiSeqIdThr.nucleotides;
    par.alnLenThr = par.multiAlnLenThr.nucleotides;
    par.addBacktrace = addBacktrace;
    cmd.addVariable("NUCL_ASM_PAR", par.createParameterString(par.nuclassembleDBworkflow).c_str());

    // set mandatory values for redundancy reduction
    par.seqIdThr = par.clustSeqIdThr;
    par.covMode = 1;
    par.covThr = par.clustCovThr;
    par.wrappedScoring = true;
    par.ignoreMultiKmer = true;
    cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.reduceredundancy).c_str());

    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    FileUtil::writeFile(tmpDir + "/hybridassembledb.sh", hybridassembledb_sh, hybridassembledb_sh_len);
    std::string program(tmpDir + "/hybridassembledb.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
