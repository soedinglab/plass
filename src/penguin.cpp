#include "Command.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"

const char* binary_name = "penguin";
const char* tool_name = "PenguiN";
const char* tool_introduction = "protein-guided nucleotide assembler.";
const char* main_author = "Annika Jochheim (annika.jochheim@mpinat.mpg.de)";
const char* show_extended_help = NULL;
const char* show_bash_info = NULL;
bool hide_base_commands = true;
void (*validatorUpdate)(void) = 0;
LocalParameters& localPar = LocalParameters::getLocalInstance();

std::vector<struct Command> commands = {
    // PenguiN workflows
    {"guided_nuclassemble",         guidedNuclAssemble,            &localPar.guidedNuclAssembleworkflow, COMMAND_MAIN,
        "Assemble nucleotide sequences by iterative greedy overlap assembly using protein and nucleotide information",
        NULL,
        "Annika Jochheim <annika.jochheim@mpinat.mpg.de>",
        "<i:fast(a|q)File[.gz]> | <i:fastqFile1_1[.gz] <i:fastqFile1_2[.gz] ... <i:fastqFileN_1[.gz] <i:fastqFileN_2[.gz]> <o:fastaFile> <tmpDir>",
        CITATION_PLASS, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}},
    {"nuclassemble",          nuclassemble,     &localPar.nuclassembleworkflow, COMMAND_MAIN,
        "Assemble nucleotide sequences by iterative greedy overlap assembly using nucleotide information only",
        NULL,
        "Annika Jochheim <annika.jochheim@mpinat.mpg.de>",
        "<i:fast(a|q)File[.gz]> | <i:fastqFile1_1[.gz] <i:fastqFile1_2[.gz] ... <i:fastqFileN_1[.gz] <i:fastqFileN_2[.gz]> <o:fastaFile> <tmpDir>",
        CITATION_PLASS, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}},
    // PenguiN tools
    {"guidedassembleresults", guidedassembleresults, &localPar.guidedassembleresults, COMMAND_HIDDEN,
        "Extending representative sequence to the left and right side using ungapped alignments.",
        NULL,
        "Annika Jochheim <annika.jochheim@mpinat.mpg.de> & Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
        "<i:nuclSequenceDB> <i:aaSequenceDB> <i:nuclAlnResult> <o:nuclAssembly> <o:aaAssembly>",
        CITATION_PLASS, {{"nuclSequenceDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb },
                            {"aaSequenceDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::aaDb },
                            {"nuclAlnResult", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                            {"nuclAssembly", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::nuclDb },
                            {"aaAssembly", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::aaDb }}},
    {"nuclassembleresults",      nuclassembleresult,       &localPar.assembleresults,      COMMAND_HIDDEN,
        "Extending representative sequence to the left and right side using ungapped alignments.",
        NULL,
        "Annika Jochheim <annika.jochheim@mpinat.mpg.de>>",
        "<i:sequenceDB> <i:alnResult> <o:reprSeqDB>",
        CITATION_PLASS, {{"sequenceDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb },
                            {"alnResult", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb  },
                            {"reprSeqDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
    {"mergereads",      mergereads,      &localPar.onlythreads,          COMMAND_HIDDEN,
        "Merge paired-end reads from FASTQ file (powered by FLASH)",
        NULL,
        "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
        "<i:fastaFile1[.gz]> ... <i:fastaFileN[.gz]> <o:sequenceDB>",
        CITATION_PLASS, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}},
    {"cyclecheck",      cyclecheck,      &localPar.cyclecheck,          COMMAND_HIDDEN,
        "Simple cycle detector",
        NULL,
        "Annika Jochheim <annika.jochheim@mpinat.mpg.de>",
        "<i:sequenceDB> <o:sequenceDBcycle>",
        CITATION_PLASS, {{"sequenceDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb },
                            {"cycleResult", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::nuclDb }}},
    {"createhdb",      createhdb,      &localPar.createhdb,          COMMAND_HIDDEN,
        "Generate header db file for given sequence db file",
        NULL,
        "Annika Jochheim <annika.jochheim@mpinat.mpg.de>",
        "<i:sequenceDB> [<i:sequenceDBcycle>] <o:headerDB>",
        CITATION_PLASS, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}}
};
