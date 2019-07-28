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
LocalParameters& localPar = LocalParameters::getLocalInstance();

std::vector<struct Command> commands = {
        // Plass tools
        {"assemble",             assembler,            &localPar.assemblerworkflow,    COMMAND_MAIN,
                "Assemble protein sequences by iterative greedy overlap assembly.",
                "Extends sequence to the left and right using ungapped alignments.",
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de> ",
                "<i:fast(a|q)File[.gz]> | <i:fastqFile1_1[.gz] ... <i:fastqFileN_1[.gz] <i:fastqFile1_2[.gz] ... <i:fastqFileN_2[.gz]> <o:fastaFile> <tmpDir>",
                CITATION_PLASS, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}},
        {"nuclassemble",          nuclassembler,     &localPar.nuclassemblerworkflow,    COMMAND_MAIN,
                "Assemble nucleotide sequences by iterative greedy overlap assembly. (experimental)",
                "Extends sequence to the left and right using ungapped alignments.",
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de> ",
                "<i:fast(a|q)File[.gz]> | <i:fastqFile1_1[.gz] ... <i:fastqFileN_1[.gz] <i:fastqFile1_2[.gz] ... <i:fastqFileN_2[.gz]> <o:fastaFile> <tmpDir>",
                CITATION_PLASS, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}},
        {"hybridassemble",         hybridassembler,            &localPar.assemblerworkflow,    COMMAND_HIDDEN,
                "Assemble protein sequences in linear time using protein and nucleotide information.",
                "Assemble protein sequences in linear time using protein and nucleotide information.",
                "Annika Seidel <annika.seidel@mpibpc.mpg.de> & Martin Steinegger <martin.steinegger@mpibpc.mpg.de> ",
                "<i:fast(a|q)File[.gz]> | <i:fastqFile1_1[.gz] ... <i:fastqFileN_1[.gz] <i:fastqFile1_2[.gz] ... <i:fastqFileN_2[.gz]> <o:fastaFile> <tmpDir>",
                CITATION_PLASS, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}},
        {"assembleresults",      assembleresult,       &localPar.assembleresults,      COMMAND_HIDDEN,
                "Extending representative sequence to the left and right side using ungapped alignments.",
                NULL,
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:sequenceDB> <i:alnResult> <o:reprSeqDB>",
                CITATION_PLASS, {{"sequenceDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                 {"alnResult", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb  },
                                 {"reprSeqDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
        {"hybridassembleresults",      hybridassembleresults,       &localPar.assembleresults,      COMMAND_HIDDEN,
                "Extending representative sequence to the left and right side using ungapped alignments.",
                NULL,
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:nuclSequenceDB> <i:aaSequenceDB> <i:nuclAlnResult> <o:nuclAssembly> <o:aaAssembly>",
                CITATION_PLASS, {{"nuclSequenceDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb },
                                 {"aaSequenceDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::aaDb },
                                 {"nuclAlnResult", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                 {"nuclAssembly", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::nuclDb },
                                 {"aaAssembly", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::aaDb }}},
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
        {"correctreads",      correctreads,      &localPar.onlythreads,          COMMAND_HIDDEN,
                "Simple read corrector",
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
                CITATION_PLASS, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}}

};
