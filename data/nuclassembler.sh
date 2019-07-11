#!/bin/sh -e

# Assembler workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# check input variables
[ ! -n "${OUT_FILE}" ] && echo "Please provide OUT_FILE" && exit 1
[ ! -n "${TMP_PATH}" ] && echo "Please provide TMP_PATH" && exit 1

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
    if notExists "${TMP_PATH}/pref_$STEP"; then
        # shellcheck disable=SC2086
        "$MMSEQS" kmermatcher "$INPUT" "${TMP_PATH}/pref_$STEP" ${KMERMATCHER_PAR} \
            || fail "Kmer matching step died"
    fi

    # 2. Ungapped alignment
    if notExists "${TMP_PATH}/aln_$STEP"; then
        # shellcheck disable=SC2086
        "$MMSEQS" rescorediagonal "$INPUT" "$INPUT" "${TMP_PATH}/pref_$STEP" "${TMP_PATH}/aln_$STEP" ${UNGAPPED_ALN_PAR} \
            || fail "Ungapped alignment step died"
    fi

    # 3. Assemble
    if notExists "${TMP_PATH}/assembly_$STEP"; then
        # shellcheck disable=SC2086
        "$MMSEQS" assembleresults "$INPUT" "${TMP_PATH}/aln_$STEP" "${TMP_PATH}/assembly_$STEP" ${ASSEMBLE_RESULT_PAR} \
            || fail "Assembly step died"
    fi

    INPUT="${TMP_PATH}/assembly_$STEP"
    STEP="$((STEP+1))"
done
STEP="$((STEP-1))"

RESULT="${TMP_PATH}/assembly_${STEP}"

# select only assembled sequences
if notExists "${RESULT}_only_assembled.index"; then
    awk 'NR == FNR { f[$1] = $0; next } $1 in f { print f[$1], $0 }' "${RESULT}.index" "${TMP_PATH}/nucl_reads.index" > "${RESULT}_tmp.index"
    awk '$3 > $6 { print }' "${RESULT}_tmp.index" > "${RESULT}_only_assembled.index"
fi

# create fasta output
if notExists "${RESULT}_only_assembled"; then
    ln -s "${RESULT}" "${RESULT}_only_assembled"
    ln -s "${RESULT},dbtype" "${RESULT}_only_assembled.dbtype"
fi

if notExists "${RESULT}_only_assembled_h"; then
    ln -s "${TMP_PATH}/nucl_reads_h" "${RESULT}_only_assembled_h"
fi

if notExists "${RESULT}_only_assembled_h.index"; then
    ln -s "${TMP_PATH}/nucl_reads_h.index" "${RESULT}_only_assembled_h.index"
fi

if notExists "${RESULT}_only_assembled.fasta"; then
    # shellcheck disable=SC2086
    "$MMSEQS" convert2fasta "${RESULT}_only_assembled" "${RESULT}_only_assembled.fasta" ${VERBOSITY_PAR} \
        || fail "convert2fasta died"
fi

mv -f "${RESULT}_only_assembled.fasta" "$OUT_FILE" \
    || fail "Could not move result to $OUT_FILE"

if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files"
    rm -f "${TMP_PATH}/nucl_reads" "${TMP_PATH}/nucl_reads.index"
    rm -f "${TMP_PATH}/pref_"*
    rm -f "${TMP_PATH}/aln_"*
    rm -f "${TMP_PATH}/assembly_"*
fi
