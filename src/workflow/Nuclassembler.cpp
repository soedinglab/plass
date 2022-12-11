#include "DBReader.h"
#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"

#include "nuclassemble.sh.h"

void setNuclAssemblerWorkflowDefaults(LocalParameters *p) {

    p->numIterations = 8;
    p->kmerSize = 22;
    p->seqIdThr = 0.99;
    p->alphabetSize = 5;

    p->covThr = 0.0;
    p->evalThr = 0.00001;
    p->maskMode = 0;
    p->kmersPerSequence = 60;
    p->kmersPerSequenceScale = 0.1;
    p->spacedKmer = false;
    p->ignoreMultiKmer = true;
    p->includeOnlyExtendable = true;
    p->addBacktrace = false;
    p->rescoreMode = Parameters::RESCORE_MODE_END_TO_END_ALIGNMENT;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
    p->maxSeqLen = 200000;
    p->cycleCheck = true;
    p->chopCycle = true;

}

int nuclassemble(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setNuclAssemblerWorkflowDefaults(&par);

    par.overrideParameterDescription(par.PARAM_MIN_SEQ_ID, "Overlap sequence identity threshold (range 0.0-1.0)", NULL, 0);
    par.overrideParameterDescription(par.PARAM_NUM_ITERATIONS, "Number of assembly iterations (range 1-inf)", NULL, 0);
    par.overrideParameterDescription(par.PARAM_E, "Extend sequences if the E-value is below (range 0.0-inf)", NULL, 0);

    // make parameter visible in short help
    par.overrideParameterDescription( par.PARAM_MAX_SEQ_LEN, NULL, NULL, MMseqsParameter::COMMAND_COMMON);


    par.PARAM_COV_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_C.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_SEQ_ID_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_INCLUDE_ONLY_EXTENDABLE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_KMER_PER_SEQ.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_SORT_RESULTS.addCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

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
    par.filenames.pop_back();
    free(p);

    cmd.addVariable("OUT_FILE", par.filenames.back().c_str());
    par.filenames.pop_back();

    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("REMOVE_INCREMENTAL_TMP", par.deleteFilesInc ? "TRUE" : NULL);

    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("NUM_IT", SSTR(par.numIterations).c_str());
    // # 1. Finding exact $k$-mer matches.
    cmd.addVariable("KMERMATCHER_PAR", par.createParameterString(par.kmermatcher).c_str());

    // # 2. Hamming distance pre-clustering
    par.filterHits = false;
    cmd.addVariable("UNGAPPED_ALN_PAR", par.createParameterString(par.rescorediagonal).c_str());
    cmd.addVariable("ASSEMBLE_RESULT_PAR", par.createParameterString(par.assembleresults).c_str());

    cmd.addVariable("CALL_CYCLE_CHECK", par.cycleCheck ? "TRUE" : NULL);
    cmd.addVariable("CYCLE_CHECK_PAR", par.createParameterString(par.cyclecheck).c_str());

    cmd.addVariable("MIN_CONTIG_LEN", SSTR(par.minContigLen).c_str());

    if (par.contigOutputMode == LocalParameters::OUTPUT_ONLY_EXTENDED_CONTIGS)
        cmd.addVariable("OUTPUT_ONLY_EXTENDED_CONTIGS", "1");
    else
        cmd.addVariable("OUTPUT_ONLY_EXTENDED_CONTIGS", NULL);

    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    cmd.addVariable("DB_MODE", par.dbMode ? "TRUE" : NULL);

    FileUtil::writeFile(tmpDir + "/nuclassemble.sh", nuclassemble_sh, nuclassemble_sh_len);
    std::string program(tmpDir + "/nuclassemble.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
