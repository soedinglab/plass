/*
 * Written by Annika Jochheim <annika.jochheim@mpibpc.mpg.de>
 * create header db
 */

#include "DBReader.h"
#include "DBWriter.h"
#include "LocalParameters.h"

#include <algorithm>
#ifdef OPENMP
#include <omp.h>
#endif

#define HEADER_INTERN_SEP  " "

int createhdb(int argc, const char **argv, const Command& command) {

    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> *seqDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str(),  1, DBReader<unsigned int>::USE_INDEX);
    seqDbr->open(DBReader<unsigned int>::NOSORT);

    const bool hasCycleLookup = par.filenames.size() > 2;
    DBReader<unsigned int> *cycleDbr = NULL;
    DBWriter *headerDbw;
    const char *headerData, *headerIndex;
    if (hasCycleLookup) {
        cycleDbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(),  1, DBReader<unsigned int>::USE_INDEX);
        cycleDbr->open(DBReader<unsigned int>::NOSORT);
        headerData = par.hdr3.c_str();
        headerIndex = par.hdr3Index.c_str();

    } else {
        headerData = par.hdr2.c_str();
        headerIndex = par.hdr2Index.c_str();
    }

    if (FileUtil::fileExists(headerData) == true) {
        Debug(Debug::WARNING) << headerData << " exists and will be overwritten.\n";
    }
    headerDbw = new DBWriter (headerData, headerIndex, 1, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    headerDbw->open();

     // header format: 'ID len:<len> cycle:<0|1>', later 'ID len:<len> cycle:<0|1> cov:<cov>'
     for (size_t id = 0; id < seqDbr->getSize(); id++) {


         std::string headerLine =
                 std::to_string(id) + HEADER_INTERN_SEP + "len:" + std::to_string(seqDbr->getSeqLen(id));
         size_t seqKey = seqDbr->getDbKey(id);

         if (cycleDbr != NULL) {
             bool cycle = (cycleDbr->getId(seqKey) != UINT_MAX);
             headerLine = headerLine + HEADER_INTERN_SEP + "cycle:" + std::to_string(cycle);
         }

         // TODO: extend header line by coverage field later
         headerLine += "\n";
         headerDbw->writeData(headerLine.c_str(), headerLine.length(), seqKey, 0);

     }


    headerDbw->close(true);
    seqDbr->close();
    delete headerDbw;
    delete seqDbr;

    if (cycleDbr != NULL) {
        cycleDbr->close();
        delete cycleDbr;
    }

    return EXIT_SUCCESS;
}
