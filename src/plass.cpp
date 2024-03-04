#include "Command.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"

const char* binary_name = "plass";
const char* tool_name = "Plass";
const char* tool_introduction = "Protein Level Assembler.";
const char* main_author = "Martin Steinegger (martin.steinegger@mpibpc.mpg.de)";
const char* show_extended_help = NULL;
const char* show_bash_info = NULL;
bool hide_base_commands = true;
void (*validatorUpdate)(void) = 0;
LocalParameters& localPar = LocalParameters::getLocalInstance();

std::vector<struct Command> commands = {
        // Plass workflows
        {"assemble",             assemble,            &localPar.assembleworkflow, COMMAND_MAIN,
                "Assemble protein sequences by iterative greedy overlap assembly",
                NULL,
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de> ",
                "<i:fast(a|q)File[.gz]> | <i:fastqFile1_1[.gz] <i:fastqFile1_2[.gz] ... <i:fastqFileN_1[.gz] <i:fastqFileN_2[.gz]> <o:fastaFile> <tmpDir>",
                CITATION_PLASS, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}},
        // Plass tools
        {"assembleresults",      assembleresult,       &localPar.assembleresults,      COMMAND_HIDDEN,
                "Extending representative sequence to the left and right side using ungapped alignments.",
                NULL,
                "Annika Jochheim <annika.jochheim@mpibpc.mpg.de> & Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:sequenceDB> <i:alnResult> <o:reprSeqDB>",
                CITATION_PLASS, {{"sequenceDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                 {"alnResult", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb  },
                                 {"reprSeqDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
        {"findassemblystart",    findassemblystart,    &localPar.onlythreads,          COMMAND_HIDDEN,
                "Compute consensus based new * stop before M amino acid",
                NULL,
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:sequenceDB> <i:alnResult> <o:sequenceDB>",
                CITATION_PLASS, {{"sequenceDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                 {"alnResult", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb  },
                                 {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
        {"filternoncoding",      filternoncoding,      &localPar.filternoncoding,          COMMAND_HIDDEN,
                "Filter non-coding protein sequences",
                NULL,
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:sequenceDB> <o:sequenceDB>",
                CITATION_PLASS, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                 {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
        {"mergereads",      mergereads,      &localPar.onlythreads,          COMMAND_HIDDEN,
                "Merge paired-end reads from FASTQ file (powered by FLASH)",
                NULL,
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:fastaFile1[.gz]> ... <i:fastaFileN[.gz]> <o:sequenceDB>",
                CITATION_PLASS, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}},
        {"createhdb",      createhdb,      &localPar.createhdb,          COMMAND_HIDDEN,
                "Generate header db file for given sequence db file",
                NULL,
                "Annika Jochheim <annika.jochheim@mpibpc.mpg.de>",
                "<i:sequenceDB> [<i:sequenceDBcycle>] <o:headerDB>",
                CITATION_PLASS, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}}
};
