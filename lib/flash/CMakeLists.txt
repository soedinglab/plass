add_library(flash combine_reads.h combine_reads.cpp read.h read.cpp util.h util.cpp)
mmseqs_setup_derived_target(flash)
