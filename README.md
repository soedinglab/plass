# PLASS and PenguiN assembler
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/plass.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/plass)
[![BioContainer Pulls](https://img.shields.io/endpoint?url=https%3A%2F%2Fmmseqs.com%2Fbiocontainer.php%3Fcontainer%3Dplass)](https://biocontainers.pro/#/tools/plass)
[![Build Status](https://travis-ci.org/soedinglab/plass.svg?branch=master)](https://travis-ci.org/soedinglab/plass)
[![DOI](https://zenodo.org/badge/118119513.svg)](https://zenodo.org/badge/latestdoi/118119513)

Plass (Protein-Level ASSembler) is a software to assemble protein sequences from short read sequencing data, while PenguiN (Protein guided nucleotide assembler) assembles DNA/RNA contigs. Both are build to assemble data from complex metagenomic datasets. This software is GPL-licensed open source software that is implemented in C++ and available for Linux and macOS and is designed to run on multiple cores. 

[Plass:](https://github.com/soedinglab/plass/tree/master?tab=readme-ov-file#plass---protein-level-assembler) [Steinegger M, Mirdita M and Soeding J. Protein-level assembly increases protein sequence recovery from metagenomic samples manyfold. Nature Methods, doi: doi.org/10.1038/s41592-019-0437-4 (2019)](https://www.nature.com/articles/s41592-019-0437-4).

[PenguiN:](https://github.com/soedinglab/plass/tree/master?tab=readme-ov-file#penguin---Protein-guided-Nucleotide-Assembler) [Jochheim A, Jochheim FE, Kolodyazhnaya A, Morice E, Steinegger M, Soeding J. Strain-resolved de-novo metagenomic assembly of viral genomes and microbial 16S rRNAs. bioRxiv (2024)](https://www.biorxiv.org/content/10.1101/2024.03.29.587318v1)

<p align="center"><img src="https://raw.githubusercontent.com/soedinglab/plass/master/.github/plass.png" height="256" /></p>

### Soil Reference Catalog (SRC) and Marine Eukaryotic Reference Catalog (MERC)
SRC was created by assembling 640 soil metagenome samples. MERC was assembled from the the metatranscriptomics datasets created by the TARA ocean expedition. Both catalogues were redundancy reduced to 90% sequence identity at 90% coverage.
Each catalog is a single FASTA file containing the sequences, the header identifiers contain the Sequence Read Archive (SRA) identifiers.
The catalogues can be downloaded [here](http://wwwuser.gwdg.de/~compbiol/plass/current_release/).
We provide a [HH-suite3](https://github.com/soedinglab/hh-suite) database called "BFD" containing sequences from the Metaclust, SRC, MERC and Uniport at [here](https://bfd.mmseqs.com/).

# PenguiN - Protein-guided Nucleotide Assembler
PenguiN a software to assemble short read sequencing data on a nucleotide level. In a first step it assembles coding sequences using the information from the translated protein sequences. In a second step it links them across non-coding regions. The main purpose of PenguiN is the assembly of complex metagenomic and metatranscriptomic datasets. It was especially tested for the assembly of viral genomes as well as 16S rRNA gene seqeunces. It assembles 3-40 times more complete viral genomes and 6 times as many 16S rRNA sequences than state of the art assemblers like Megahit and the SPAdes variants. PenguiN is GPL-licensed open source software that is implemented in C++ and available for Linux and macOS. The software is designed to run on multiple cores. 

### Install Plass and PenguiN
Our sofwtare can be install via [conda](https://github.com/conda/conda) or as statically compiled Linux version. It requires a 64-bit Linux/MacOS system (check with `uname -a | grep x86_64`) with at least the SSE2 instruction set.

     # install from bioconda
     conda install -c conda-forge -c bioconda plass 
     # static build with AVX2 (fastest)
     wget https://mmseqs.com/plass/plass-linux-avx2.tar.gz; tar xvfz plass-linux-avx2.tar.gz; export PATH=$(pwd)/plass/bin/:$PATH
     # static build with SSE4.1
     wget https://mmseqs.com/plass/plass-linux-sse41.tar.gz; tar xvfz plass-linux-sse41.tar.gz; export PATH=$(pwd)/plass/bin/:$PATH
     # static build with SSE2 (slowest, for very old systems)
     wget https://mmseqs.com/plass/plass-linux-sse2.tar.gz; tar xvfz plass-linux-sse2.tar.gz; export PATH=$(pwd)/plass/bin/:$PATH
 

## How to assemble
Plass and PenguiN can assemble both paired-end reads (FASTQ) and single reads (FASTA or FASTQ):

      # assemble paired-end reads 
      plass assemble examples/reads_1.fastq.gz examples/reads_2.fastq.gz assembly.fas tmp

      # assemble single-end reads 
      plass assemble examples/reads_1.fastq.gz assembly.fas tmp

      # assemble single-end reads using stdin
      cat examples/reads_1.fastq.gz | plass assemble stdin assembly.fas tmp


Important parameters: 

     --min-seq-id         Adjusts the overlap sequence identity threshold
     --min-length         minimum codon length for ORF prediction (default: 40)
     -e                   E-value threshold for overlaps 
     --num-iterations     Number of iterations of assembly
     --filter-proteins    Switches the neural network protein filter off/on

Plass Modules: 

      plass assemble      Assembles proteins (i:Nucleotides -> o:Proteins)
      
      
PenguiN Modules: 

      penguin guided_nuclassemble  Assembles nucleotides using protein and nucleotide information (i:Nucleotides -> o:Nucleotides)
      penguin nuclassemble         Assembles nucleotides using only nucleotdie information (i:Nucleotides -> o:Nucleotides)

### Assemble using MPI 
Both tools can be distributed over several homogeneous computers. However the TMP folder has to be shared between all nodes (e.g. NFS). The following command assembles several nodes:

    RUNNER="mpirun -np 42" plass assemble examples/reads_1.fastq.gz examples/reads_2.fastq.gz assembly.fas tmp


### Compile from source
Compiling from source has the advantage that it will be optimized to the specific system, which should improve its performance. To compile `git`, `g++` (4.6 or higher) and `cmake` (3.0 or higher) are required. Afterwards, the PLASS and PenguiN binaries will be located in the `build/bin` directory.

      git clone https://github.com/soedinglab/plass.git
      cd plass
      git submodule update --init
      mkdir build && cd build
      cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
      make -j 4 && make install
      export PATH="$(pwd)/bin/:$PATH"
        
:exclamation: If you want to compile PLASS or PenguiN on macOS, please install and use `gcc` from Homebrew. The default macOS `clang` compiler does not support OpenMP and PLASS will not be able to run multithreaded. Use the following cmake call:

      CXX="$(brew --prefix)/bin/g++-8" cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..

#### Dependencies

When compiling from source, our sofwtare requires `zlib` and `bzip`.

### Use the docker image
We also provide a Docker image of Plass. You can mount the current directory containing the reads to be assembled and run plass with the following command:

      docker pull soedinglab/plass
      docker run -ti --rm -v "$(pwd):/app" -w /app plass assemble reads_1.fastq reads_2.fastq assembly.fas tmp

## Hardware requirements
Plass needs roughly 1 byte of memory per residue to work efficiently. Plass will scale its memory consumption based on the available main memory of the machine. Plass needs a CPU with at least the SSE4.1 instruction set to run. 

## Known problems 
* The assembly of Plass includes all ORFs having a start and end codon that includes even very short ORFs < 60 amino acids. Many of these short ORFs are spurious since our neural network cannot distingue them well. We would recommend to use other method to verify the coding potential of these. Assemblies above 100 amino acids are mostly genuine protein sequences. 
* Plass in default searches for ORFs of 40 amino acids or longer. This limits the read length to > 120. To assemble this protein, you need to lower the `--min-length` threshold. Be aware using short reads (< 100 length) might result in lower sensitivity.
