#include <mmseqs/src/commons/DBReader.h>
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
    p->includeOnlyExtendable = true;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
}

int assembler(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    setAssemblerWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, 3);

    const int dbType = DBReader<unsigned int>::parseDbType(par.db1.c_str());
    if(dbType==DBReader<unsigned int>::DBTYPE_AA){
        par.alphabetSize = 21;
        par.kmerSize = 14;
    }
    if (FileUtil::directoryExists(par.db3.c_str()) == false){
        Debug(Debug::INFO) << "Temporary folder " << par.db3 << " does not exist or is not a directory.\n";
        if (FileUtil::makeDir(par.db3.c_str()) == false){
            Debug(Debug::WARNING) << "Could not crate tmp folder " << par.db3 << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created directory " << par.db3 << "\n";
        }
    }
    size_t hash = par.hashParameter(par.filenames, par.searchworkflow);
    std::string tmpDir(par.db3 + "/" + SSTR(hash));
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::WARNING) << "Could not create sub folder in temporary directory " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);
    if (FileUtil::symlinkCreateOrRepleace(par.db3 + "/latest", tmpDir) == false){
        Debug(Debug::WARNING) << "Could not symlink latest folder in temporary directory." << tmpDir << ".\n";
        EXIT(EXIT_FAILURE);
    }

    CommandCaller cmd;
    if (par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("NUM_IT", SSTR(par.numIterations).c_str());
    // nucleotide assembly
    if(dbType == DBReader<unsigned int>::DBTYPE_NUC){
        cmd.addVariable("NUCL", "1");
    }
    // save some values to restore them later
    size_t alphabetSize = par.alphabetSize;
    size_t kmerSize = par.kmerSize;

    // # 1. Finding exact $k$-mer matches.

    cmd.addVariable("KMERMATCHER_PAR", par.createParameterString(par.kmermatcher).c_str());

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

    return EXIT_SUCCESS;
}
