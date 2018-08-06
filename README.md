# PLASS - Protein Level ASSembly
[ ![Codeship Status for soedinglab/plass](https://app.codeship.com/projects/fc7c4e70-e188-0135-0db2-569fac09cf96/status?branch=master)](https://app.codeship.com/projects/266646)

Plass (Protein-level-assembler) is a software to assemble short read sequencing data on a protein level. Plass is GPL-licensed open source software that is implemented in C++ and available for Linux and macOS. The software is designed to run on multiple cores. 
 
### Install static Linux version
The following command will download the last Plass version, extract it and set the `PATH` environment variable. This version runs only on Linux. If you want to run it on Mac please compile it from source.

      wget https://mmseqs.com/latest/plass-static_sse41.tar.gz
      tar xvzf plass-static_sse41.tar.gz
      export PATH=$(pwd)/plass/bin/:$PATH

### Compile from source
Compiling PLASS from source has the advantage that it will be optimized to the specific system, which should improve its performance. To compile PLASS `git`, `g++` (4.6 or higher) and `cmake` (3.0 or higher) are required. Afterwards, the PLASS binary will be located in the `build/bin` directory.

      git clone https://github.com/soedinglab/plass.git
      cd plass
      git submodule update --init
      mkdir build && cd build
      cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
      make -j 4 && make install
      export PATH="$(pwd)/bin/:$PATH"
        
:exclamation: If you want to compile PLASS on macOS, please install and use `gcc` from Homebrew. The default macOS `clang` compiler does not support OpenMP and PLASS will not be able to run multithreaded. Use the following cmake call:

      CXX="$(brew --prefix)/bin/g++-8" cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..

## How to assemble

PLASS can assemble both paired-end reads (FASTQ) and single reads (FASTA or FASTQ):

      # assemble paired-end reads 
      plass assemble read_1.fastq read_2.fastq assembly.fas tmp

      # assemble single-end reads 
      plass assemble reads.fasta assembly.fas tmp
      
## Soil Reference Catalog (SRC) and Marine Eukaryotic Reference Catalog (MERC)
SRC was created by assembling 640 soil metagenome samples. MERC was assembled from the the metatranscriptomics datasets created by the TARA ocean expedition. Both catalogues were redundancy reduced to 90% sequence identity at 90% coverage.
Each catalog is a single FASTA file containing the sequences, the header identifiers contain the Sequence Read Archive (SRA) identifiers.
The catalogues can be downloaded [here](http://wwwuser.gwdg.de/~compbiol/plass/current_release/).
