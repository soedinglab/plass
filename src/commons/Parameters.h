// Written by Martin Steinegger martin.steinegger@mpibpc.mpg.de
//
// Represents a parameter of MMseqs
//
#ifndef MMSEQS_PARAMETERS
#define MMSEQS_PARAMETERS
#include <string>
#include <vector>
#include <typeinfo>

#define PARAMETER(x) const static int x##_ID = __COUNTER__; \
    				 MMseqsParameter x;

struct MMseqsParameter {
    const int uniqid;
    const char *name;
    const char *display;
    const char *description;
    const std::type_info &type;
    void * value;
    const char * regex;
    bool wasSet;
    MMseqsParameter(int uid, const char * n, const char *display,
                    const char * d, const std::type_info &hash, void * value, const char * regex):
                    uniqid(uid), name(n), display(display), description(d), type(hash), value(value), regex(regex), wasSet(false){}
};


class Parameters          // Parameters for gap penalties and pseudocounts
{
public:

    static const unsigned int ALIGNMENT_MODE_FAST_AUTO = 0;
    static const unsigned int ALIGNMENT_MODE_SCORE_ONLY = 1;
    static const unsigned int ALIGNMENT_MODE_SCORE_COV = 2;
    static const unsigned int ALIGNMENT_MODE_SCORE_COV_SEQID = 3;
    // prefilter search
    static const int SEARCH_GLOBAL = 0;
    static const int SEARCH_LOCAL = 1;
    static const int SEARCH_LOCAL_FAST = 2;
    // format alignment
    static const int FORMAT_ALIGNMENT_BLAST_TAB = 0;
    static const int FORMAT_ALIGNMENT_PAIRWISE  = 1;
    static const int FORMAT_ALIGNMENT_SAM       = 2;

    static const int PROFILE_MODE_HMM = 0;
    static const int PROFILE_MODE_PSSM = 1;

    static const int SET_COVER = 0;
    static const int CONNECTED_COMPONENT = 1;
    static const int GREEDY = 2;

    static const int APC_ALIGNMENTSCORE=1;
    static const int APC_SEQID=2;
    // split mode
    static const int TARGET_DB_SPLIT = 0;
    static const int QUERY_DB_SPLIT = 1;
    static const int DETECT_BEST_DB_SPLIT = 2;
    // split
    static const int AUTO_SPLIT_DETECTION = 0;
    // includeIdentity
    static const int INCLUDE_HIT_AUTO = 0;
    static const int FORCE_INCLUDE = 1;

    static const int MAX_SEQ_LEN = 32000;

    // COMMON
    const char** argv;            //command line parameters
    char argc;              //dimension of argv
    
    // path to databases
    std::string db1;
    std::string db1Index;
    
    std::string db2;
    std::string db2Index;

    std::string db3;
    std::string db3Index;
    
    std::string db4;
    std::string db4Index;
    
    std::string db5;
    std::string db5Index;

    std::string mmdir;

    std::string scoringMatrixFile;       // path to scoring matrix
    size_t maxSeqLen;                    // sequence length
    size_t maxResListLen;                // Maximal result list length per query
    int    verbosity;                    // log level
    int    querySeqType;                 // Query sequence type (PROFILE, AMINOACIDE, NUCLEOTIDE)
    int    targetSeqType;                // Target sequence type (PROFILE, AMINOACIDE, NUCLEOTIDE)
    int    threads;                      // Amounts of threads
    bool   removeTmpFiles;               // Do not delete temp files
    bool   includeIdentity;              // include identical ids as hit
    // PREFILTER
    float  sensitivity;                  // target sens
    int    kmerSize;                     // kmer size for the prefilter
    int    kmerScore;                    // kmer score for the prefilter
    int    alphabetSize;                 // alphabet size for the prefilter
    int    searchMode;                   // Local search type
    bool   profile;                      // using profile information
    bool   nucl;                         // using nucl informatoin
    int    compBiasCorrection;           // Aminoacid composiont correction
    int    diagonalScoring;              // switch diagonal scoring
    int    minDiagScoreThr;              // min diagonal score
    int    spacedKmer;                   // Spaced Kmers
    int    split;                        // Split database in n equal chunks
    int    splitMode;                    // Split by query or target DB (MPI only)
    bool   splitAA;                      // Split database by amino acid count instead

    // ALIGNMENT
    std::string ffindexPrefDB;           // prefilter database (input for alignment module)
    int alignmentMode;                   // alignment mode 0=fastest on parameters,
                                         // 1=score only, 2=score, cov, start/end pos, 3=score, cov, start/end pos, seq.id,
    float  evalThr;                      // e-value threshold for acceptance
    float  covThr;                       // coverage threshold for acceptance
    int    maxRejected;                  // after n sequences that are above eval stop
    float  seqIdThr;                     // sequence identity threshold for acceptance
    bool   fragmentMerge;                // allow fragments to in the result
    bool   addBacktrace;                 // store backtrace string (M=Match, D=deletion, I=insertion)
    bool   realign;                      // realign hit with more conservative score
	
    // workflow
    std::string runner;

    // CLUSTERING
    std::string ffindexAlnDBBase;
    int    clusteringMode;
    int    validateClustering;
    bool   cascaded;

    // SEARCH WORKFLOW
    int numIterations;
    int startSens;
    int sensStepSize;
    //CLUSTERING
    int maxIteration;                   // Maximum depth of breadth first search in connected component
    int similarityScoreType;            // Type of score to use for reassignment 1=alignment score. 2=coverage 3=sequence identity 4=E-value 5= Score per Column

    //extractorf
    int orfMinLength;
    int orfMaxLength;
    int orfMaxGaps;
    bool orfSkipIncomplete;
    bool orfLongest;
    bool orfExtendMin;
    std::string forwardFrames;
    std::string reverseFrames;

    // createprofiledb
    int profileMode;
    bool useIndex;
    // format alignment
    int formatAlignmentMode;            // BLAST_TAB, PAIRWISE or SAM

    // result2msa
    bool allowDeletion;
    bool addInternalId;
    bool compressMSA;
    bool summarizeHeader;
    std::string summaryPrefix;
    bool onlyRepSeq;

    // result2profile
    float filterMaxSeqId;
    float evalProfile;
    float qsc;
    float qid;
    float cov;
    int Ndiff;
    bool wg;
    float pca;
    float pcb;
    bool noPruning;

    // createdb
    bool useHeader;
    int identifierOffset;
    bool splitSeqByLen;

    // rebuildfasta
    bool useHeaderFile;

    // gff2ffindex
    std::string gffType;

    // translate nucleotide
    int translationTable;

    // addSequences
    int minSequences;

    // filterDb
    int filterColumn;
    std::string filterColumnRegex;
	std::string filteringFile;
	std::string mappingFile;
	bool positiveFilter;
	bool trimToOneColumn;

    // evaluationscores
    bool allVsAll;
    bool randomizedRepresentative;
    bool use_sequenceheader;



    
    void checkSaneEnvironment();
    void setDefaults();
    void parseParameters(int argc, const char* argv[],
                         const std::string &programUsageHeader,
                         std::vector<MMseqsParameter> &par,
                         size_t requiredParameterCount,
                         bool printParameters = true,
                         bool isVariadic = false);
    void printUsageMessage(const std::string &programUsageHeader,
                           const std::vector<MMseqsParameter> &parameters);
    void printParameters(int argc, const char* pargv[],
                         const std::vector<MMseqsParameter> &par);
	
	std::vector<MMseqsParameter> removeParameter(std::vector<MMseqsParameter> par,MMseqsParameter x);
	
    Parameters();
    ~Parameters(){};

    PARAMETER(PARAM_S)
    PARAMETER(PARAM_K)
    PARAMETER(PARAM_THREADS)
    PARAMETER(PARAM_ALPH_SIZE)
    PARAMETER(PARAM_MAX_SEQ_LEN)
    PARAMETER(PARAM_PROFILE)
    PARAMETER(PARAM_NUCL)
    PARAMETER(PARAM_DIAGONAL_SCORING)
    PARAMETER(PARAM_MIN_DIAG_SCORE)
    PARAMETER(PARAM_K_SCORE)
    PARAMETER(PARAM_MAX_SEQS)
    PARAMETER(PARAM_SPLIT)
    PARAMETER(PARAM_SPLIT_MODE)
    PARAMETER(PARAM_SPLIT_AMINOACID)
    PARAMETER(PARAM_SUB_MAT)
    PARAMETER(PARAM_SEARCH_MODE)
    PARAMETER(PARAM_NO_COMP_BIAS_CORR)
    PARAMETER(PARAM_SPACED_KMER_MODE)
    PARAMETER(PARAM_REMOVE_TMP_FILES)
    PARAMETER(PARAM_INCLUDE_IDENTITY)
    std::vector<MMseqsParameter> prefilter;

    // alignment
    PARAMETER(PARAM_ALIGNMENT_MODE)
    PARAMETER(PARAM_E)
    PARAMETER(PARAM_C)
    PARAMETER(PARAM_FRAG_MERGE)
    PARAMETER(PARAM_MAX_REJECTED)
    PARAMETER(PARAM_ADD_BACKTRACE)
    PARAMETER(PARAM_REALIGN)
    PARAMETER(PARAM_MIN_SEQ_ID)

    std::vector<MMseqsParameter> alignment;

    // clustering
    PARAMETER(PARAM_CLUSTER_MODE)

    PARAMETER(PARAM_CASCADED)

    //afinity clustering
    PARAMETER(PARAM_MAXITERATIONS)
    PARAMETER(PARAM_SIMILARITYSCORE)
    // logging
    PARAMETER(PARAM_V)
    std::vector<MMseqsParameter> clustering;

    // create profile (HMM, PSSM)
    PARAMETER(PARAM_PROFILE_TYPE)

    // format alignment
    PARAMETER(PARAM_FORMAT_MODE)

    // result2msa
    PARAMETER(PARAM_ALLOW_DELETION)
    PARAMETER(PARAM_ADD_INTERNAL_ID)
    PARAMETER(PARAM_COMPRESS_MSA)
    PARAMETER(PARAM_SUMMARIZE_HEADER)
    PARAMETER(PARAM_SUMMARY_PREFIX)
    PARAMETER(PARAM_REPSEQ)

    // result2profile
    PARAMETER(PARAM_E_PROFILE)
    PARAMETER(PARAM_FILTER_MAX_SEQ_ID)
    PARAMETER(PARAM_FILTER_QSC)
    PARAMETER(PARAM_FILTER_QID)
    PARAMETER(PARAM_FILTER_COV)
    PARAMETER(PARAM_FILTER_NDIFF)
    PARAMETER(PARAM_WG)
    PARAMETER(PARAM_PCA)
    PARAMETER(PARAM_PCB)
//    PARAMETER(PARAM_NO_PRUNING)


    // workflow
    PARAMETER(PARAM_RUNNER)

    // search workflow
    PARAMETER(PARAM_NUM_ITERATIONS)
    PARAMETER(PARAM_START_SENS)
    PARAMETER(PARAM_SENS_STEP_SIZE)
    PARAMETER(PARAM_USE_INDEX)
    // extractorfs
    PARAMETER(PARAM_ORF_MIN_LENGTH)
    PARAMETER(PARAM_ORF_MAX_LENGTH)
    PARAMETER(PARAM_ORF_MAX_GAP)
    PARAMETER(PARAM_ORF_SKIP_INCOMPLETE)
    PARAMETER(PARAM_ORF_LONGEST)
    PARAMETER(PARAM_ORF_EXTENDMIN)
    PARAMETER(PARAM_ORF_FORWARD_FRAMES)
    PARAMETER(PARAM_ORF_REVERSE_FRAMES)

    // createdb
    PARAMETER(PARAM_USE_HEADER) // also used by extractorf
    PARAMETER(PARAM_ID_OFFSET)  // same
    PARAMETER(PARAM_DONT_SPLIT_SEQ_BY_LEN)

    // rebuildfasta
    PARAMETER(PARAM_USE_HEADER_FILE)

    // gff2ffindex
    PARAMETER(PARAM_GFF_TYPE)

    // translate_nucleotide
    PARAMETER(PARAM_TRANSLATION_TABLE)

    // addsequences
    PARAMETER(PARAM_MIN_SEQUENCES)

    // filterDb
    PARAMETER(PARAM_FILTER_COL)
    PARAMETER(PARAM_FILTER_REGEX)
    PARAMETER(PARAM_FILTER_POS)
    PARAMETER(PARAM_FILTER_FILE)
    PARAMETER(PARAM_MAPPING_FILE)
    PARAMETER(PARAM_TRIM_TO_ONE_COL)
	
    // evaluationScore
    PARAMETER(PARAM_EVALUATION_ALLVSALL)
    PARAMETER(PARAM_EVALUATION_RANDOMIZEDREPRESENTATIVE)
    PARAMETER(PARAM_EVALUATION_USE_SEQUENCEHEADER)

    std::vector<MMseqsParameter> empty;

    std::vector<MMseqsParameter> onlyverbosity;
    std::vector<MMseqsParameter> createFasta;
    std::vector<MMseqsParameter> createprofiledb;
    std::vector<MMseqsParameter> result2profile;
    std::vector<MMseqsParameter> result2msa;
    std::vector<MMseqsParameter> extractorf;
    std::vector<MMseqsParameter> splitffindex;
    std::vector<MMseqsParameter> createindex;
    std::vector<MMseqsParameter> formatalignment;
    std::vector<MMseqsParameter> createdb;
    std::vector<MMseqsParameter> rebuildfasta;
    std::vector<MMseqsParameter> gff2ffindex;
    std::vector<MMseqsParameter> detectredundancy;
    std::vector<MMseqsParameter> searchworkflow;
    std::vector<MMseqsParameter> clusteringWorkflow;
    std::vector<MMseqsParameter> clusterUpdateSearch;
    std::vector<MMseqsParameter> clusterUpdateClust;
    std::vector<MMseqsParameter> clusterUpdate;
    std::vector<MMseqsParameter> translateNucleotide;
    std::vector<MMseqsParameter> addSequences;
    std::vector<MMseqsParameter> filterDb;
    std::vector<MMseqsParameter> swapresults;
    std::vector<MMseqsParameter> substractresult;
    std::vector<MMseqsParameter> result2newick;
    std::vector<MMseqsParameter> diff;
    std::vector<MMseqsParameter> dbconcat;

    std::vector<MMseqsParameter> evaluationscores;

    std::vector<MMseqsParameter> combineList(std::vector<MMseqsParameter> &par1,
                                              std::vector<MMseqsParameter> &par2);

    std::string createParameterString(std::vector < MMseqsParameter > &vector);

};

#endif
