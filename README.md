# PLASS protein level assembly


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
        
:exclamation: Please install and use `gcc` from Homebrew, if you want to compile PLASS on MacOS. The default MacOS `clang` compiler does not support OpenMP and MMseqs2 will not be able to run multithreaded. Use the following cmake call:

        CXX="$(brew --prefix)/bin/g++-6" cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
        
        
                

## How to assemble

       # combine reads using FLASH (Install flash from https://ccb.jhu.edu/software/FLASH/)
       flash read_1.fastq read_2.fastq
       cat out.extendedFrags.fastq out.notCombined_1.fastq out.notCombined_2.fastq > all_merged_reads_nucl.fastq
       mmseqs createdb all_merged_reads_nucl.fastq all_merged_reads_nucl
       mmseqs assemble all_merged_reads_nucl all_merged_reads_aa_assembly tmp  
