#include "DBReader.h"
#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"

#include "assemble.sh.h"

void setAssembleDBWorkflowDefaults(LocalParameters *p) {
    p->spacedKmer = false;
    p->maskMode = 0;
    p->covThr = 0.0;
    p->evalThr = 0.00001;
    p->seqIdThr = 0.9;
    p->kmersPerSequence = 60;
//    p->shuffleDatabase = true;
//    p->splitSeqByLen = false;
    p->numIterations = 12;
    p->alphabetSize = 13;
    p->kmerSize = 14;
    p->orfMinLength = 45;
    p->ignoreMultiKmer = true;
    p->includeOnlyExtendable = true;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
    p->rescoreMode = Parameters::RESCORE_MODE_END_TO_END_ALIGNMENT;
}

int assemble(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setAssembleDBWorkflowDefaults(&par);

    par.overrideParameterDescription(par.PARAM_MIN_SEQ_ID, "Overlap sequence identity threshold [0.0, 1.0]", NULL, 0);
    par.overrideParameterDescription(par.PARAM_NUM_ITERATIONS, "Number of assembly iterations [1, inf]", NULL,  0);
    par.overrideParameterDescription(par.PARAM_E, "Extend sequences if the E-value is below [0.0, inf]", NULL,  0);

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
    cmd.addVariable("PROTEIN_FILTER", par.filterProteins == 1 ? "1" : NULL);
    // # 1. Finding exact $k$-mer matches.

    for(int i = 0; i < par.numIterations; i++){
        std::string key = "KMERMATCHER"+SSTR(i)+"_PAR";
        par.hashShift = par.hashShift + i % 2;
        if(par.PARAM_INCLUDE_ONLY_EXTENDABLE.wasSet == false){
            if (i == 0) {
                par.includeOnlyExtendable = false;
            } else {
                par.includeOnlyExtendable = true;
            }
        }
        cmd.addVariable(key.c_str(), par.createParameterString(par.kmermatcher).c_str());
    }

    cmd.addVariable("KMERMATCHER_PAR", par.createParameterString(par.kmermatcher).c_str());

    // # 2. Hamming distance pre-clustering
    par.filterHits = false;

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


    par.addOrfStop = true;
    //cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createdb).c_str());
    cmd.addVariable("TRANSLATENUCS_PAR", par.createParameterString(par.translatenucs).c_str());
    cmd.addVariable("UNGAPPED_ALN_PAR", par.createParameterString(par.rescorediagonal).c_str());
    cmd.addVariable("ASSEMBLE_RESULT_PAR", par.createParameterString(par.assembleresults).c_str());
    cmd.addVariable("FILTERNONCODING_PAR", par.createParameterString(par.filternoncoding).c_str());

    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    FileUtil::writeFile(tmpDir + "/assemble.sh", assemble_sh, assemble_sh_len);
    std::string program(tmpDir + "/assemble.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
