#include "DBReader.h"
#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"

#include "assembler.sh.h"

void setAssemblerWorkflowDefaults(LocalParameters *p) {
    p->spacedKmer = false;
    p->maskMode = 0;
    p->covThr = 0.0;
    p->evalThr = 0.00001;
    p->seqIdThr = 0.9;
    p->kmersPerSequence = 60;
    p->numIterations = 12;
    p->alphabetSize = 13;
    p->kmerSize = 14;
    p->skipNRepeatKmer = 8;
    p->includeOnlyExtendable = true;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
}

int assembler(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setAssemblerWorkflowDefaults(&par);
    par.overrideParameterDescription((Command &)command, par.PARAM_COV_MODE.uniqid, NULL, NULL, par.PARAM_COV_MODE.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_C.uniqid, NULL, NULL, par.PARAM_C.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_MIN_SEQ_ID.uniqid, "Overlap sequence identity threshold [0.0, 1.0]", NULL,  par.PARAM_MIN_SEQ_ID.category);
    par.overrideParameterDescription((Command &)command, par.PARAM_NUM_ITERATIONS.uniqid, "Number of assembly iterations [1, inf]", NULL,  par.PARAM_NUM_ITERATIONS.category);
    par.overrideParameterDescription((Command &)command, par.PARAM_E.uniqid, "Extend sequences if the E-value is below [0.0, inf]", NULL,  par.PARAM_E.category);

    par.overrideParameterDescription((Command &)command, par.PARAM_SEQ_ID_MODE.uniqid, NULL, NULL,  par.PARAM_SEQ_ID_MODE.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_RESCORE_MODE.uniqid, NULL, NULL,  par.PARAM_RESCORE_MODE.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_INCLUDE_ONLY_EXTENDABLE.uniqid, NULL, NULL,  par.PARAM_INCLUDE_ONLY_EXTENDABLE.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_KMER_PER_SEQ.uniqid, NULL, NULL,  par.PARAM_KMER_PER_SEQ.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_SORT_RESULTS.uniqid, NULL, NULL,  par.PARAM_SORT_RESULTS.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_TRANSLATION_TABLE.uniqid, NULL, NULL, par.PARAM_TRANSLATION_TABLE.category | MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription((Command &)command, par.PARAM_USE_ALL_TABLE_STARTS.uniqid, NULL, NULL, par.PARAM_USE_ALL_TABLE_STARTS.category | MMseqsParameter::COMMAND_EXPERT);
    par.parseParameters(argc, argv, command, 3, true, Parameters::PARSE_VARIADIC);

    CommandCaller cmd;
    if ((par.filenames.size() - 2) % 2 == 0) {
        // paired end reads
        cmd.addVariable("PAIRED_END", "1");
    } else {
        if (par.filenames.size() != 3) {
            Debug(Debug::ERROR) << "Too many input files provided.\n";
            Debug(Debug::ERROR) << "For paired-end input provide READSETA_1.fastq READSETA_2.fastq ... OUTPUT.fasta tmpDir\n";
            Debug(Debug::ERROR) << "For single input use READSET.fast(q|a) OUTPUT.fasta tmpDir\n";
            return EXIT_FAILURE;
        }
        cmd.addVariable("PAIRED_END", NULL);
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
    size_t hash = par.hashParameter(par.filenames, par.assemblerworkflow);
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
    cmd.addVariable("NUM_IT", SSTR(par.numIterations).c_str());
    cmd.addVariable("PROTEIN_FILTER", par.filterProteins == 1 ? "1" : NULL);
    // # 1. Finding exact $k$-mer matches.
    cmd.addVariable("KMERMATCHER_PAR", par.createParameterString(par.kmermatcher).c_str());

    // # 2. Hamming distance pre-clustering
    par.filterHits = false;
    par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
    cmd.addVariable("EXTRACTORFS_SUBSET_PAR", par.createParameterString(par.extractorfssubset).c_str());
    par.addOrfStop = true;
    cmd.addVariable("TRANSLATENUCS_PAR", par.createParameterString(par.translatenucs).c_str());
    cmd.addVariable("UNGAPPED_ALN_PAR", par.createParameterString(par.rescorediagonal).c_str());
    cmd.addVariable("ASSEMBLE_RESULT_PAR", par.createParameterString(par.assembleresults).c_str());

    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    FileUtil::writeFile(tmpDir + "/assembler.sh", assembler_sh, assembler_sh_len);
    std::string program(tmpDir + "/assembler.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
