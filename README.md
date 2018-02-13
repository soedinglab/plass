# PLASS protein level assembly
[ ![Codeship Status for soedinglab/plass](https://app.codeship.com/projects/fc7c4e70-e188-0135-0db2-569fac09cf96/status?branch=master)](https://app.codeship.com/projects/266646)

### Compile from source
Compiling PLASS from source has the advantage that it will be optimized to the specific system, which should improve its performance. To compile PLASS `git`, `g++` (4.6 or higher) and `cmake` (3.0 or higher) are needed. Afterwards, the PLASS binary will be located in `build/bin/`.

        git clone https://github.com/soedinglab/plass.git
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

PLASS was evaluated on illumina HiSeq 2500 150x2 paired end reads but should be able to assemble protein fragments of all sizes.

       # combine reads using FLASH (Install flash from https://ccb.jhu.edu/software/FLASH/)
       flash read_1.fastq read_2.fastq
       cat out.extendedFrags.fastq out.notCombined_1.fastq out.notCombined_2.fastq > all_merged_reads_nucl.fastq
       # create internal database
       mmseqs createdb all_merged_reads_nucl.fastq all_merged_reads_nucl
       # assemble
       mmseqs assemble all_merged_reads_nucl all_merged_reads_aa_assembly tmp  --mask 0 --min-seq-id 0.9 --num-iterations 12

