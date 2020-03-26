#ifndef LOCALCOMMANDDECLARATIONS_H
#define LOCALCOMMANDDECLARATIONS_H

#include "Command.h"

extern int easyassembler(int argc, const char **argv, const Command &command);
extern int easynuclassembler(int argc, const char **argv, const Command &command);
extern int easyhybridassembler(int argc, const char** argv, const Command &command);
//
extern int assembledb(int argc, const char **argv, const Command &command);
extern int nuclassembledb(int argc, const char** argv, const Command &command);
extern int hybridassembledb(int argc, const char** argv, const Command &command);
extern int assembleresult(int argc, const char** argv, const Command &command);
extern int hybridassembleresults(int argc, const char** argv, const Command &command);
extern int nuclassembleresult(int argc, const char** argv, const Command &command);
extern int filternoncoding(int argc, const char** argv, const Command &command);
extern int mergereads(int argc, const char** argv, const Command &command);
extern int findassemblystart(int argc, const char** argv, const Command &command);
extern int cyclecheck(int argc, const char** argv, const Command &command);
extern int createhdb(int argc, const char** argv, const Command &command);
#endif
