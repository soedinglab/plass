# PLASS protein level assembly
[ ![Codeship Status for soedinglab/plass](https://app.codeship.com/projects/fc7c4e70-e188-0135-0db2-569fac09cf96/status?branch=master)](https://app.codeship.com/projects/266646)

Plass (Protein-level-assembler) is a software to assemble short reads on a protein level. Plass is open source GPL-licensed software implemented in C++ for Linux, MacOS. The software is designed to run on multiple cores.


### Install static Linux version
The following command will download the last Plass version, extract it and set the `PATH` variable. This version runs only on linux. If you want to run it on Mac please compile it or use brew.

       wget https://mmseqs.com/latest/plass-static_sse41.tar.gz
       tar xvzf plass-static_sse41.tar.gz
       export PATH=$(pwd)/plass/bin/:$PATH

### Compile from source
Compiling PLASS from source has the advantage that it will be optimized to the specific system, which should improve its performance. To compile PLASS `git`, `g++` (4.6 or higher) and `cmake` (3.0 or higher) are needed. Afterwards, the PLASS binary will be located in `build/bin/`.

       git clone https://github.com/soedinglab/plass.git
       git submodule init
       git submodule update
       cd plass
       mkdir build
       cd build
       cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
       make
       make install
       export PATH=$(pwd)/bin/:$PATH
      
:exclamation: Please install and use `gcc` from Homebrew, if you want to compile PLASS on MacOS. The default MacOS `clang` compiler does not support OpenMP and PLASS will not be able to run multithreaded. Use the following cmake call:

       CXX="$(brew --prefix)/bin/g++-6" cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
      

## How to assemble

PLASS can assemble paired-end reads (fastq) and single reads

      # assemble paired-end reads
      plass assemble read_1.fastq read_2.fastq assembly.fas tmp
      
      # assemble single-end reads
      plass assemble reads.fasta assembly.fas tmp
    
## Soil Reference Catalog (SRC) and Marine Eukaryotic Reference Catalog (MERC)

SRC was created by assembling 640 soil metagenome samples. MERC contains the metatranscriptomic datasets created by the Tara ocean expedition. Both catalogs were redundancy reduced to 90% sequence identity at 90% coverage.

Each catalog is a single Fasta file containing the sequences, the header identifier is Sequence Read Archive identifier.

The catalogs can be downloaded [here](http://wwwuser.gwdg.de/~compbiol/plass/current_release/)

