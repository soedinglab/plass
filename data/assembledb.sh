#!/bin/sh -e
# Assembler workflow script
fail() {
    echo "Error: $1"
    exit 1
}

deleteIncremental() {
    if [ -n "$REMOVE_INCREMENTAL_TMP" ] &&  [ -n "$1" ]; then
         "$MMSEQS" rmdb "$1"
    fi
}

notExists() {
	[ ! -f "$1" ]
}

# check input variables
[ -z "${OUT_FILE}" ] && echo "Please provide OUT_FILE" && exit 1
[ -z "${TMP_PATH}" ] && echo "Please provide TMP_PATH" && exit 1

# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[   -f "${OUT_FILE}" ] &&  echo "${OUT_FILE}.dbtype exists already!" && exit 1
[ ! -d "${TMP_PATH}" ] &&  echo "tmp directory ${TMP_PATH} not found!" && mkdir -p "${TMP_PATH}"


INPUT="$1"
if notExists "${TMP_PATH}/nucl_6f_start"; then
    # shellcheck disable=SC2086
    "$MMSEQS" extractorfs "${INPUT}" "${TMP_PATH}/nucl_6f_start" ${EXTRACTORFS_START_PAR} \
        || fail "extractorfs start step died"
fi

if notExists "${TMP_PATH}/aa_6f_start"; then
    # shellcheck disable=SC2086
    "$MMSEQS" translatenucs "${TMP_PATH}/nucl_6f_start" "${TMP_PATH}/aa_6f_start" ${TRANSLATENUCS_PAR} \
        || fail "translatenucs start step died"
fi

if notExists "${TMP_PATH}/nucl_6f_long"; then
    # shellcheck disable=SC2086
    "$MMSEQS" extractorfs "${INPUT}" "${TMP_PATH}/nucl_6f_long" ${EXTRACTORFS_LONG_PAR} \
        || fail "extractorfs longest step died"
fi

if notExists "${TMP_PATH}/aa_6f_long"; then
    # shellcheck disable=SC2086
    "$MMSEQS" translatenucs "${TMP_PATH}/nucl_6f_long" "${TMP_PATH}/aa_6f_long" ${TRANSLATENUCS_PAR} \
        || fail "translatenucs long step died"
fi

if notExists "${TMP_PATH}/aa_6f_start_long"; then
    # shellcheck disable=SC2086
    "$MMSEQS" concatdbs "${TMP_PATH}/aa_6f_long" "${TMP_PATH}/aa_6f_start" "${TMP_PATH}/aa_6f_start_long" ${VERBOSITY_PAR} \
        || fail "concatdbs start long step died"
fi

if notExists "${TMP_PATH}/aa_6f_start_long_h"; then
    #awk 'BEGIN { printf("%c%c%c%c",12,0,0,0); exit; }' > "${TMP_PATH}/nucl_6f_long_h.dbtype"
    #awk 'BEGIN { printf("%c%c%c%c",12,0,0,0); exit; }' > "${TMP_PATH}/nucl_6f_start_h.dbtype"
    # shellcheck disable=SC2086
    "$MMSEQS" concatdbs "${TMP_PATH}/nucl_6f_long_h" "${TMP_PATH}/nucl_6f_start_h" "${TMP_PATH}/aa_6f_start_long_h" ${VERBOSITY_PAR} \
        || fail "concatdbs start long step died"
fi

INPUT="${TMP_PATH}/aa_6f_start_long"
STEP=0
if [ -z "$NUM_IT" ]; then
    NUM_IT=1
fi

while [ "$STEP" -lt "$NUM_IT" ]; do
    echo "STEP: $STEP"
    # 1. Finding exact $k$-mer matches.
    if notExists "${TMP_PATH}/pref_$STEP.done"; then
        PARAM=KMERMATCHER${STEP}_PAR
        eval KMERMATCHER_TMP="\$$PARAM"
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" kmermatcher "$INPUT" "${TMP_PATH}/pref_$STEP" ${KMERMATCHER_TMP} \
            || fail "Kmer matching step died"
        deleteIncremental "$PREV_KMER_PREF"
        touch "${TMP_PATH}/pref_$STEP.done"
        PREV_KMER_PREF="${TMP_PATH}/pref_$STEP"

    fi

    # 2. Ungapped alignment
    if notExists "${TMP_PATH}/aln_$STEP.done"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" rescorediagonal "$INPUT" "$INPUT" "${TMP_PATH}/pref_$STEP" "${TMP_PATH}/aln_$STEP" ${UNGAPPED_ALN_PAR} \
            || fail "Ungapped alignment step died"
        touch "${TMP_PATH}/aln_$STEP.done"
        deleteIncremental "$PREV_ALN"
        PREV_ALN="${TMP_PATH}/aln_$STEP"
    fi

    ALN="${TMP_PATH}/aln_$STEP"
    if [ $STEP -eq 0 ]; then
        if notExists "${TMP_PATH}/corrected_seqs.done"; then
            # shellcheck disable=SC2086
            "$MMSEQS" findassemblystart "$INPUT" "${TMP_PATH}/aln_$STEP" "${TMP_PATH}/corrected_seqs" ${THREADS_PAR} \
                || fail "Findassemblystart alignment step died"
              touch "${TMP_PATH}/corrected_seqs.done"
              # delete at the end of the first iteration
              PREV_ASSEMBLY="${TMP_PATH}/corrected_seqs"
        fi
        INPUT="${TMP_PATH}/corrected_seqs"

        if notExists "${TMP_PATH}/pref_corrected_$STEP.done"; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" kmermatcher "$INPUT" "${TMP_PATH}/pref_corrected_$STEP" ${KMERMATCHER_TMP} \
                || fail "Kmer matching step died"
            deleteIncremental "$PREV_KMER_PREF"
            touch "${TMP_PATH}/pref_corrected_$STEP.done"
            PREV_KMER_PREF="${TMP_PATH}/pref_corrected_$STEP"
        fi

        if notExists "${TMP_PATH}/aln_corrected_$STEP.done"; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" rescorediagonal "$INPUT" "$INPUT" "${TMP_PATH}/pref_corrected_$STEP" "${TMP_PATH}/aln_corrected_$STEP" ${UNGAPPED_ALN_PAR} \
                || fail "Ungapped alignment step died"
           touch "${TMP_PATH}/aln_corrected_$STEP.done"
           deleteIncremental "$PREV_ALN"
           PREV_ALN="${TMP_PATH}/aln_corrected_$STEP"
        fi
        ALN="${TMP_PATH}/aln_corrected_$STEP"
    fi

    # 3. Assemble
    if notExists "${TMP_PATH}/assembly_$STEP.done"; then
        # shellcheck disable=SC2086
        "$MMSEQS" assembleresults "$INPUT" "${ALN}" "${TMP_PATH}/assembly_$STEP" ${ASSEMBLE_RESULT_PAR} \
            || fail "Assembly step died"

        touch "${TMP_PATH}/assembly_$STEP.done"
        deleteIncremental "$PREV_ASSEMBLY"
        PREV_ASSEMBLY="${TMP_PATH}/assembly_$STEP"
    fi

    INPUT="${TMP_PATH}/assembly_$STEP"
    STEP="$((STEP+1))"

done
STEP="$((STEP-1))"

# post processing
RESULT="${TMP_PATH}/assembly_${STEP}"
if [ -n "${PROTEIN_FILTER}" ]; then
    RESULT="${TMP_PATH}/assembly_${STEP}_filtered"
    if notExists "${TMP_PATH}/assembly_${STEP}_filtered"; then
        # shellcheck disable=SC2086
        "$MMSEQS" filternoncoding "${TMP_PATH}/assembly_${STEP}" "${TMP_PATH}/assembly_${STEP}_filtered" ${FILTERNONCODING_PAR} \
            || fail "filternoncoding died"
    fi
fi

# select only assembled sequences
if notExists "${RESULT}_only_assembled.index"; then
    # detect assembled proteins sequences
    awk 'NR == FNR { f[$1] = $0; next } $1 in f { print f[$1], $0 }' "${RESULT}.index" "${TMP_PATH}/aa_6f_start_long.index" > "${RESULT}_tmp.index"
    awk '$3 > $6 { print $1"\t"$2"\t"$3 }' "${RESULT}_tmp.index" > "${RESULT}_only_assembled1.index"
    # detect complete proteins with * at start and end
    awk '/^\x00?\*[A-Z]*\*$/{ f[NR-1]=1; next } $1 in f { print $0 }' "${RESULT}" "${RESULT}.index" > "${RESULT}_only_assembled2.index"
    # keep only non-redundant entries
    cat "${RESULT}_only_assembled1.index" "${RESULT}_only_assembled2.index" | sort | uniq > "${RESULT}_only_assembled.index"
fi

# create db outfile
if notExists "${OUT_FILE}.dbtype"; then
     # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${RESULT}_only_assembled.index" "${RESULT}" "${OUT_FILE}" --subdb-mode 0 \
        || fail "Createsubdb died"
fi

if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files"
    rm -f "${TMP_PATH}/aa_6f_"*
    rm -f "${TMP_PATH}/nucl_6f_"*
    rm -f "${TMP_PATH}/pref_"*
    rm -f "${TMP_PATH}/aln_"*
    rm -f "${TMP_PATH}/assembly_"*
    rm -f "${TMP_PATH}/assembledb.sh"
fi
