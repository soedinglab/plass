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
if notExists "${TMP_PATH}/nucl_6f_start"; then
    # shellcheck disable=SC2086
    "$MMSEQS" extractorfs "${INPUT}" "${TMP_PATH}/nucl_6f_start" ${EXTRACTORFS_START_PAR} \
        || fail "extractorfs start step died"
fi

if notExists "${TMP_PATH}/nucl_6f_long"; then
    # shellcheck disable=SC2086
    "$MMSEQS" extractorfs "${INPUT}" "${TMP_PATH}/nucl_6f_long" ${EXTRACTORFS_LONG_PAR} \
        || fail "extractorfs longest step died"
fi

if notExists "${TMP_PATH}/nucl_6f_start_long"; then
    "$MMSEQS" concatdbs "${TMP_PATH}/nucl_6f_long" "${TMP_PATH}/nucl_6f_start" "${TMP_PATH}/nucl_6f_start_long" \
        || fail "concatdbs start long step died"
fi

if notExists "${TMP_PATH}/nucl_6f_start_long_h"; then
    # shellcheck disable=SC2086
    "$MMSEQS" concatdbs "${TMP_PATH}/nucl_6f_long_h" "${TMP_PATH}/nucl_6f_start_h" "${TMP_PATH}/nucl_6f_start_long_h" ${VERBOSITY_PAR} \
        || fail "concatdbs start long step died"
fi

if notExists "${TMP_PATH}/aa_6f_start_long"; then
    "$MMSEQS" translatenucs "${TMP_PATH}/nucl_6f_start_long" "${TMP_PATH}/aa_6f_start_long" --add-orf-stop \
        || fail "translatenucs step died"
fi

INPUT_AA="${TMP_PATH}/aa_6f_start_long"
INPUT_NUCL="${TMP_PATH}/nucl_6f_start_long"
STEP=0
if [ -z "$NUM_IT" ]; then
    NUM_IT=1
fi

while [ $STEP -lt $NUM_IT ]; do
    echo "STEP: $STEP"

    # 1. Finding exact $k$-mer matches.
    if notExists "${TMP_PATH}/pref_$STEP"; then
        # shellcheck disable=SC2086
        "$MMSEQS" kmermatcher "$INPUT_AA" "${TMP_PATH}/pref_$STEP" ${KMERMATCHER_PAR}   \
            || fail "Kmer matching step died"
    fi

    # 2. Ungapped alignment
    if notExists "${TMP_PATH}/aln_$STEP"; then
        # shellcheck disable=SC2086
        "$MMSEQS" rescorediagonal "$INPUT_AA" "$INPUT_AA" "${TMP_PATH}/pref_$STEP" "${TMP_PATH}/aln_$STEP" ${UNGAPPED_ALN_PAR} \
            || fail "Ungapped alignment step died"
    fi

    # 3. Ungapped alignment protein 2 nucl
    if notExists "${TMP_PATH}/aln_nucl_$STEP"; then
        "$MMSEQS" proteinaln2nucl "$INPUT_NUCL" "$INPUT_NUCL" "$INPUT_AA"  "$INPUT_AA"  "${TMP_PATH}/aln_$STEP" "${TMP_PATH}/aln_nucl_$STEP"  \
            || fail "Ungapped alignment 2 nucl step died"
    fi

    # 4. Assemble
    if notExists "${TMP_PATH}/assembly_aa_$STEP" || notExists "${TMP_PATH}/assembly_aa_$STEP"; then
        # shellcheck disable=SC2086
        "$MMSEQS" hybridassembleresults "$INPUT_NUCL" "$INPUT_AA" "${TMP_PATH}/aln_nucl_$STEP" "${TMP_PATH}/assembly_nucl_$STEP" "${TMP_PATH}/assembly_aa_$STEP" ${ASSEMBLE_RESULT_PAR} \
            || fail "Assembly step died"
    fi

    INPUT_AA="${TMP_PATH}/assembly_aa_$STEP"
    INPUT_NUCL="${TMP_PATH}/assembly_nucl_$STEP"
    STEP="$((STEP+1))"
done
STEP="$((STEP-1))"

RESULT_NUCL="${TMP_PATH}/assembly_nucl_$STEP"
#RESULT_AA="${TMP_PATH}/assembly_aa_$STEP"

# select only assembled orfs
if notExists "${RESULT_NUCL}_only_assembled.index"; then
    awk 'NR == FNR { f[$1] = $0; next } $1 in f { print f[$1], $0 }' "${RESULT_NUCL}.index" "${TMP_PATH}/nucl_6f_start_long.index" > "${RESULT_NUCL}_tmp.index"
    awk '$3 > $6 { print }' "${RESULT_NUCL}_tmp.index" > "${RESULT_NUCL}_only_assembled.index"
fi

if notExists "${RESULT_NUCL}_only_assembled"; then
    ln -s "${RESULT_NUCL}" "${RESULT_NUCL}_only_assembled"
fi

if notExists "${RESULT_NUCL}_only_assembled.dbtype"; then
    ln -s "${RESULT_NUCL}.dbtype" "${RESULT_NUCL}_only_assembled.dbtype"
fi

if notExists "${RESULT_NUCL}_h"; then
    ln -s "${TMP_PATH}/nucl_6f_start_long_h" "${RESULT_NUCL}_only_assembled_h"
fi

if notExists "${RESULT_NUCL}_h.index"; then
    ln -s "${TMP_PATH}/nucl_6f_start_long_h.index" "${RESULT_NUCL}_only_assembled_h.index"
fi

if notExists "${RESULT_NUCL}.fasta"; then
    # shellcheck disable=SC2086
    "$MMSEQS" convert2fasta "${RESULT_NUCL}_only_assembled" "${RESULT_NUCL}_only_assembled.fasta" ${VERBOSITY_PAR} \
        || fail "convert2fasta died"
fi

if notExists "${RESULT_NUCL}.merged.fasta"; then
    # shellcheck disable=SC2086
    "$MMSEQS" convert2fasta "${INPUT}" "${INPUT}.fasta"
    cat "${RESULT_NUCL}_only_assembled.fasta" "${INPUT}.fasta" > "${RESULT_NUCL}.merged.fasta"
fi

# shellcheck disable=SC2086
"$MMSEQS" nuclassemble "${RESULT_NUCL}.merged.fasta" "${OUT_FILE}" "${TMP_PATH}/nuclassembly_2" ${NUCL_ASM_PAR}

#mv -f "${TMP_PATH}/assembly_aa_${STEP}" "${2}_aa" || fail "Could not move result to $2"
#mv -f "${TMP_PATH}/assembly_aa_${STEP}.index" "${2}_aa.index" || fail "Could not move result to $2.index"


if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files"
    rm -f "${TMP_PATH}/pref_"*
    rm -f "${TMP_PATH}/aln_"*
    rm -f "${TMP_PATH}/assembly_"*
    rm -f "${TMP_PATH}/nuclassembly_2/latest/"*
fi
