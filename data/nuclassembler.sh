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

cyclecheck() {
	if [ -n "$CALL_CYCLE_CHECK" ]; then
        if notExists "${1}_cycle.done"; then
            # shellcheck disable=SC2086
            "$MMSEQS" cyclecheck "$1" "${1}_cycle" ${CYCLE_CHECK_PAR} \
                || fail "Cycle check step died"


            if [ -s "${1}_cycle" ]; then

                if notExists "${1}_noneCycle"; then
                    awk 'NR==FNR { a[$1]=$0; next } !($1 in a) {print $0}' "${1}_cycle.index" \
                    "${1}.index" > "${1}_noneCycle.index"
                    ln -s "$1" "${1}_noneCycle"
                    ln -s "${1}.dbtype" "${1}_noneCycle.dbtype"
                fi

                if [ -z "$PREV_CYCLE_ALL" ]; then
                    # shellcheck disable=SC2086
                    "$MMSEQS" mvdb "${1}_cycle" "${1}_cycle_all"
                else
                    # shellcheck disable=SC2086
                    "$MMSEQS" concatdbs "${PREV_CYCLE_ALL}" "${1}_cycle" "${1}_cycle_all" --preserve-keys
                fi

            else
                ln -s "$1" "${1}_noneCycle"
                ln -s "${1}.index" "${1}_noneCycle.index"
                ln -s "${1}.dbtype" "${1}_noneCycle.dbtype"
            fi
            touch "${1}_cycle.done"
            deleteIncremental "$PREV_CYCLE"
            PREV_CYCLE="${1}_cycle"
        fi

        if [ -s "${1}_cycle_all" ]; then
            deleteIncremental "${PREV_CYCLE_ALL}"
            PREV_CYCLE_ALL="${1}_cycle_all"
        fi

        PREV_ASSEMBLY="${1}_noneCycle"
    fi
}

# check input variables
[ -z "${OUT_FILE}" ] && echo "Please provide OUT_FILE" && exit 1
[ -z "${TMP_PATH}" ] && echo "Please provide TMP_PATH" && exit 1

# check if files exists
[   -f "${OUT_FILE}" ] &&  echo "${OUT_FILE} exists already!" && exit 1
[ ! -d "${TMP_PATH}" ] &&  echo "tmp directory ${TMP_PATH} not found!" && mkdir -p "${TMP_PATH}"

if notExists "${TMP_PATH}/nucl_reads"; then
    if [ -n "${PAIRED_END}" ]; then
        echo "PAIRED END MODE"
        # shellcheck disable=SC2086
        "$MMSEQS" mergereads "$@" "${TMP_PATH}/nucl_reads" ${VERBOSITY_PAR} \
            || fail "mergereads failed"
    else
        # shellcheck disable=SC2086
        "$MMSEQS" createdb "$@" "${TMP_PATH}/nucl_reads" ${VERBOSITY_PAR} \
            || fail "createdb failed"
    fi
fi
INPUT="${TMP_PATH}/nucl_reads"

STEP=0
if [ -z "$NUM_IT" ]; then
    NUM_IT=1;
fi

while [ $STEP -lt $NUM_IT ]; do
    echo "STEP: $STEP"

    # 1. Finding exact $k$-mer matches.
    if notExists "${TMP_PATH}/pref_${STEP}.done"; then
        # shellcheck disable=SC2086
        "$MMSEQS" kmermatcher "$INPUT" "${TMP_PATH}/pref_${STEP}" ${KMERMATCHER_PAR} \
            || fail "Kmer matching step died"
        deleteIncremental "$PREV_KMER_PREF"
        touch "${TMP_PATH}/pref_${STEP}.done"
        PREV_KMER_PREF="${TMP_PATH}/pref_${STEP}"
    fi

    # 2. Ungapped alignment
    if notExists "${TMP_PATH}/aln_${STEP}.done"; then
        # shellcheck disable=SC2086
        "$MMSEQS" rescorediagonal "$INPUT" "$INPUT" "${TMP_PATH}/pref_${STEP}" "${TMP_PATH}/aln_${STEP}" ${UNGAPPED_ALN_PAR} \
            || fail "Ungapped alignment step died"
        touch "${TMP_PATH}/aln_${STEP}.done"
        deleteIncremental "$PREV_ALN"
        PREV_ALN="${TMP_PATH}/aln_${STEP}"
    fi

    # 3. Assemble
    if notExists "${TMP_PATH}/assembly_${STEP}.done"; then
        # shellcheck disable=SC2086
        "$MMSEQS" assembleresults "$INPUT" "${TMP_PATH}/aln_${STEP}" "${TMP_PATH}/assembly_${STEP}" ${ASSEMBLE_RESULT_PAR} \
            || fail "Assembly step died"
        touch "${TMP_PATH}/assembly_${STEP}.done"
        deleteIncremental "$PREV_ASSEMBLY"
        deleteIncremental "$PREV_ASSEMBLY_STEP"
    fi

    PREV_ASSEMBLY="${TMP_PATH}/assembly_${STEP}"
    PREV_ASSEMBLY_STEP="${TMP_PATH}/assembly_${STEP}"
    cyclecheck "${PREV_ASSEMBLY}"

    INPUT="${PREV_ASSEMBLY}"
    STEP="$((STEP+1))"
done
STEP="$((STEP-1))"
RESULT="${TMP_PATH}/assembly_${STEP}"

if [ -n "$PREV_CYCLE_ALL" ]; then

    RESULT="${TMP_PATH}/assembly_merged"
    if notExists "${TMP_PATH}/assembly_merged"; then
        # shellcheck disable=SC2086
        "$MMSEQS" concatdbs "${PREV_ASSEMBLY}" "${PREV_CYCLE_ALL}" "${TMP_PATH}/assembly_merged" --preserve-keys
    fi
fi

# select only assembled sequences
if notExists "${RESULT}_only_assembled.index"; then
    awk 'NR == FNR { f[$1] = $0; next } $1 in f { print f[$1], $0 }' "${RESULT}.index" "${TMP_PATH}/nucl_reads.index" > "${RESULT}_tmp.index"
    awk '$3 > $6 { print }' "${RESULT}_tmp.index" > "${RESULT}_only_assembled.index"
fi

# select only sequences fullfilling a minimum length threshold
if notExists "${RESULT}_filtered.index"; then
    # shellcheck disable=SC208
    awk -v thr="${MIN_CONTIG_LEN}" '$3 > (thr+1) { print }' "${RESULT}_only_assembled.index" > "${TMP_PATH}/assembly_final.index"
fi

# create fasta output
if notExists "${TMP_PATH}/assembly_final"; then
    ln -s "${RESULT}" "${TMP_PATH}/assembly_final"
fi

if notExists "${TMP_PATH}/assembly_final.dbtype"; then
    ln -s "${RESULT}.dbtype" "${TMP_PATH}/assembly_final.dbtype"
fi

# redundancy reduction by using linclust
if notExists "${TMP_PATH}/assembly_final_rep"; then

    CLUST_INPUT="${TMP_PATH}/assembly_final"
    if notExists "${TMP_PATH}/clu.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" linclust "${CLUST_INPUT}" "${TMP_PATH}/clu" "${TMP_PATH}/clu_tmp" ${CLUSTER_PAR} \
            || fail "Search died"
    fi

    # shellcheck disable=SC2086
    "$MMSEQS" result2repseq "${CLUST_INPUT}" "${TMP_PATH}/clu" "${TMP_PATH}/assembly_final_rep" ${THREADS_PAR} \
         || fail "Result2repseq  died"
fi


if notExists "${TMP_PATH}/assembly_final_rep_h"; then
    # shellcheck disable=SC2086
    if notExists "${PREV_CYCLE_ALL}";then
        "$MMSEQS" createhdb "${TMP_PATH}/assembly_final_rep" "${TMP_PATH}/assembly_final_rep" ${VERBOSITY_PAR} \
                || fail "createhdb failed"
    else
        "$MMSEQS" createhdb "${TMP_PATH}/assembly_final_rep" "${PREV_CYCLE_ALL}" "${TMP_PATH}/assembly_final_rep" ${VERBOSITY_PAR} \
                || fail "createhdb failed"
    fi
fi

if notExists "${TMP_PATH}/assembly_final_rep.fasta"; then
    # shellcheck disable=SC2086
    "$MMSEQS" convert2fasta "${TMP_PATH}/assembly_final_rep" "${TMP_PATH}/assembly_final_rep.fasta" ${VERBOSITY_PAR} \
        || fail "convert2fasta died"
fi

mv -f "${TMP_PATH}/assembly_final_rep.fasta" "$OUT_FILE" \
    || fail "Could not move result to $OUT_FILE"

if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files"
    rm -f "${TMP_PATH}/nucl_reads" "${TMP_PATH}/nucl_reads.index"
    rm -f "${TMP_PATH}/pref_"*
    rm -f "${TMP_PATH}/aln_"*
    rm -f "${TMP_PATH}/assembly_"*
fi
