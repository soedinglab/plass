#include "DBReader.h"
#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "hybridassembler.sh.h"

void setHybridAssemblerWorkflowDefaults(LocalParameters *p) {
    p->spacedKmer = false;
    p->maskMode = 0;
    p->covThr = 0.0;
    p->evalThr = 0.00001;
    p->seqIdThr = 0.97;
//    p->alphabetSize = 21;
    p->kmersPerSequence = 60;
    p->numProtIterations = 12;
    p->numNuclIterations = 12;
    p->includeOnlyExtendable = true;
    p->orfMinLength = 45;
    p->ignoreMultiKmer = true;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
    p->rescoreMode = Parameters::RESCORE_MODE_GLOBAL_ALIGNMENT;

}

int hybridassembler(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setHybridAssemblerWorkflowDefaults(&par);

    par.overrideParameterDescription(par.PARAM_MIN_SEQ_ID, "Overlap sequence identity threshold [0.0, 1.0]", NULL, 0);
//    par.overrideParameterDescription(par.PARAM_ORF_MIN_LENGTH, "Min codons in orf", "minimum codon number in open reading frames", 0);
//    par.overrideParameterDescription(par.PARAM_NUM_ITERATIONS, "Number of assembly iterations [1, inf]", NULL, 0);
    par.overrideParameterDescription(par.PARAM_E, "Extend sequences if the E-value is below [0.0, inf]", NULL, 0);

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

    par.parseParameters(argc, argv, command, true, 0, 0);
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
            Debug(Debug::ERROR) << "Could not create tmp folder " << tmpPath << ".\n";
            return EXIT_FAILURE;
        } else {
            Debug(Debug::INFO) << "Created directory " << tmpPath << "\n";
        }
    }
    size_t hash = par.hashParameter(par.filenames, *command.params);
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
    //cmd.addVariable("NUM_IT", SSTR(par.numNuclIterations).c_str());


    // # 1. Finding exact $k$-mer matches.
    cmd.addVariable("KMERMATCHER_PAR", par.createParameterString(par.kmermatcher).c_str());


    cmd.addVariable("NUM_IT", SSTR(par.numProtIterations).c_str());

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


    // # 2. Hamming distance pre-clustering
    par.filterHits = false;
    par.addBacktrace = true;
    cmd.addVariable("UNGAPPED_ALN_PAR", par.createParameterString(par.rescorediagonal).c_str());
    cmd.addVariable("ASSEMBLE_RESULT_PAR", par.createParameterString(par.assembleresults).c_str());

    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    // set nucleassemble default values when calling nucleassemble from hybridassemble
    par.numIterations = par.numNuclIterations;
    par.kmerSize = 22;
    par.alphabetSize = 5;
    //par.kmersPerSequence = 60;
    par.kmersPerSequenceScale = 0.1;
    par.addBacktrace = false;
    par.cycleCheck = true;
    par.chopCycle = true;
    cmd.addVariable("NUCL_ASM_PAR", par.createParameterString(par.nuclassemblerworkflow).c_str());

    FileUtil::writeFile(tmpDir + "/hybridassembler.sh", hybridassembler_sh, hybridassembler_sh_len);
    std::string program(tmpDir + "/hybridassembler.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
