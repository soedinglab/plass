# PLASS protein level assembly
[ ![Codeship Status for soedinglab/plass](https://app.codeship.com/projects/fc7c4e70-e188-0135-0db2-569fac09cf96/status?branch=master)](https://app.codeship.com/projects/266646)

Plass (Proten-level-assembler) is a software to assemble short reads on a protein level. Plass is open source GPL-licensed software implemented in C++ for Linux, MacOS. The software is designed to run on multiple cores and servers. 
 
 
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
      