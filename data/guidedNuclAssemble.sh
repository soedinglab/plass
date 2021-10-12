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
        "$MMSEQS" createdb "$@" "${TMP_PATH}/nucl_reads" ${CREATEDB_PAR} \
            || fail "createdb failed"
     fi
fi

INPUT="${TMP_PATH}/nucl_reads"
TMP_PATH_GUIDED_ASSEMBLY="${TMP_PATH}/guidedassembly_tmp"
[ ! -d "${TMP_PATH_GUIDED_ASSEMBLY}" ] &&  echo "tmp directory ${TMP_PATH_GUIDED_ASSEMBLY} not found!" && mkdir -p "${TMP_PATH_GUIDED_ASSEMBLY}"


if notExists "${TMP_PATH_GUIDED_ASSEMBLY}/nucl_6f_start"; then
    # shellcheck disable=SC2086
    "$MMSEQS" extractorfs "${INPUT}" "${TMP_PATH_GUIDED_ASSEMBLY}/nucl_6f_start" ${EXTRACTORFS_START_PAR} \
        || fail "extractorfs start step died"
fi

if notExists "${TMP_PATH_GUIDED_ASSEMBLY}/nucl_6f_long"; then
    # shellcheck disable=SC2086
    "$MMSEQS" extractorfs "${INPUT}" "${TMP_PATH_GUIDED_ASSEMBLY}/nucl_6f_long" ${EXTRACTORFS_LONG_PAR} \
        || fail "extractorfs longest step died"
fi

if notExists "${TMP_PATH_GUIDED_ASSEMBLY}/nucl_6f_start_long"; then
    "$MMSEQS" concatdbs "${TMP_PATH_GUIDED_ASSEMBLY}/nucl_6f_long" "${TMP_PATH_GUIDED_ASSEMBLY}/nucl_6f_start" "${TMP_PATH_GUIDED_ASSEMBLY}/nucl_6f_start_long" \
        || fail "concatdbs start long step died"
fi

if notExists "${TMP_PATH_GUIDED_ASSEMBLY}/nucl_6f_start_long_h"; then
    # shellcheck disable=SC2086
    "$MMSEQS" concatdbs "${TMP_PATH_GUIDED_ASSEMBLY}/nucl_6f_long_h" "${TMP_PATH_GUIDED_ASSEMBLY}/nucl_6f_start_h" "${TMP_PATH_GUIDED_ASSEMBLY}/nucl_6f_start_long_h" ${VERBOSITY_PAR} \
        || fail "concatdbs start long step died"
fi

if notExists "${TMP_PATH_GUIDED_ASSEMBLY}/aa_6f_start_long"; then
    "$MMSEQS" translatenucs "${TMP_PATH_GUIDED_ASSEMBLY}/nucl_6f_start_long" "${TMP_PATH_GUIDED_ASSEMBLY}/aa_6f_start_long" --add-orf-stop \
        || fail "translatenucs step died"
fi

INPUT_AA="${TMP_PATH_GUIDED_ASSEMBLY}/aa_6f_start_long"
INPUT_NUCL="${TMP_PATH_GUIDED_ASSEMBLY}/nucl_6f_start_long"
STEP=0
if [ -z "$NUM_IT" ]; then
    NUM_IT=1
fi

while [ $STEP -lt $NUM_IT ]; do
    echo "STEP: $STEP"

    # 1. Finding exact $k$-mer matches.
    if notExists "${TMP_PATH_GUIDED_ASSEMBLY}/pref_$STEP.done"; then
        # shellcheck disable=SC2086
        "$MMSEQS" kmermatcher "$INPUT_AA" "${TMP_PATH_GUIDED_ASSEMBLY}/pref_$STEP" ${KMERMATCHER_PAR}   \
            || fail "Kmer matching step died"
        deleteIncremental "$PREV_KMER_PREF"
        touch "${TMP_PATH_GUIDED_ASSEMBLY}/pref_${STEP}.done"
        PREV_KMER_PREF="${TMP_PATH_GUIDED_ASSEMBLY}/pref_${STEP}"
    fi

    # 2. Ungapped alignment
    if notExists "${TMP_PATH_GUIDED_ASSEMBLY}/aln_$STEP.done"; then
        # shellcheck disable=SC2086
        "$MMSEQS" rescorediagonal "$INPUT_AA" "$INPUT_AA" "${TMP_PATH_GUIDED_ASSEMBLY}/pref_$STEP" "${TMP_PATH_GUIDED_ASSEMBLY}/aln_$STEP" ${UNGAPPED_ALN_PAR} \
            || fail "Ungapped alignment step died"
        touch "${TMP_PATH_GUIDED_ASSEMBLY}/aln_$STEP.done"
        deleteIncremental "$PREV_ALN"
        PREV_ALN="${TMP_PATH_GUIDED_ASSEMBLY}/aln_$STEP"
    fi

    # 3. Ungapped alignment protein 2 nucl
    if notExists "${TMP_PATH_GUIDED_ASSEMBLY}/aln_nucl_$STEP.done"; then
        "$MMSEQS" proteinaln2nucl "$INPUT_NUCL" "$INPUT_NUCL" "$INPUT_AA"  "$INPUT_AA"  "${TMP_PATH_GUIDED_ASSEMBLY}/aln_$STEP" "${TMP_PATH_GUIDED_ASSEMBLY}/aln_nucl_$STEP"  \
            || fail "Ungapped alignment 2 nucl step died"
        deleteIncremental "$PREV_ALN_NUCL"
        touch "${TMP_PATH_GUIDED_ASSEMBLY}/aln_nucl_${STEP}.done"
        PREV_ALN_NUCL="${TMP_PATH_GUIDED_ASSEMBLY}/aln_nucl_$STEP"
    fi

    # 4. Assemble
    if notExists "${TMP_PATH_GUIDED_ASSEMBLY}/assembly_aa_nucl_$STEP.done"; then
        # shellcheck disable=SC2086
        "$MMSEQS" guidedassembleresults "$INPUT_NUCL" "$INPUT_AA" "${TMP_PATH_GUIDED_ASSEMBLY}/aln_nucl_$STEP" "${TMP_PATH_GUIDED_ASSEMBLY}/assembly_nucl_$STEP" "${TMP_PATH_GUIDED_ASSEMBLY}/assembly_aa_$STEP" ${ASSEMBLE_RESULT_PAR} \
            || fail "Assembly step died"
        touch "${TMP_PATH_GUIDED_ASSEMBLY}/assembly_aa_nucl_$STEP.done"
        deleteIncremental "$PREV_ASSEMBLY_AA"
        deleteIncremental "$PREV_ASSEMBLY_NUCL"
        PREV_ASSEMBLY_AA="${TMP_PATH_GUIDED_ASSEMBLY}/assembly_aa_$STEP"
        PREV_ASSEMBLY_NUCL="${TMP_PATH_GUIDED_ASSEMBLY}/assembly_nucl_$STEP"
    fi

    INPUT_AA="${TMP_PATH_GUIDED_ASSEMBLY}/assembly_aa_$STEP"
    INPUT_NUCL="${TMP_PATH_GUIDED_ASSEMBLY}/assembly_nucl_$STEP"
    STEP="$((STEP+1))"
done
STEP="$((STEP-1))"

RESULT_NUCL="${TMP_PATH_GUIDED_ASSEMBLY}/assembly_nucl_$STEP"
#RESULT_AA="${TMP_PATH}/assembly_aa_$STEP"

# select only assembled orfs
if notExists "${RESULT_NUCL}_only_assembled.index"; then
    awk 'NR == FNR { f[$1] = $0; next } $1 in f { print f[$1], $0 }' "${RESULT_NUCL}.index" "${TMP_PATH_GUIDED_ASSEMBLY}/nucl_6f_start_long.index" > "${RESULT_NUCL}_tmp.index"
    awk '$3 > $6 { print }' "${RESULT_NUCL}_tmp.index" > "${RESULT_NUCL}_only_assembled.index"
fi

if notExists "${RESULT_NUCL}_only_assembled"; then
    ln -s "${RESULT_NUCL}" "${RESULT_NUCL}_only_assembled"
fi

if notExists "${RESULT_NUCL}_only_assembled.dbtype"; then
    ln -s "${RESULT_NUCL}.dbtype" "${RESULT_NUCL}_only_assembled.dbtype"
fi

if notExists "${RESULT_NUCL}_only_assembled_h"; then
    ln -s "${TMP_PATH_GUIDED_ASSEMBLY}/nucl_6f_start_long_h" "${RESULT_NUCL}_only_assembled_h"
fi

if notExists "${RESULT_NUCL}_only_assembled_h.index"; then
    ln -s "${TMP_PATH_GUIDED_ASSEMBLY}/nucl_6f_start_long_h.index" "${RESULT_NUCL}_only_assembled_h.index"
fi

if notExists "${RESULT_NUCL}_only_assembled_h.dbtype"; then
    ln -s "${TMP_PATH_GUIDED_ASSEMBLY}/nucl_6f_start_long_h.dbtype" "${RESULT_NUCL}_only_assembled_h.dbtype"
fi

if notExists "${RESULT_NUCL}.merged.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" concatdbs "${RESULT_NUCL}_only_assembled" "${INPUT}" "${TMP_PATH}/guided_assembly.merged" \
    || fail "Concat hybridassemblies and reads died"
fi

if notExists "${TMP_PATH}/nuclassembly.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" nuclassemble "${TMP_PATH}/guided_assembly.merged" "${TMP_PATH}/nuclassembly" "${TMP_PATH}/nuclassembly_tmp" ${NUCL_ASM_PAR}
fi

# redundancy reduction using linclust
if notExists "${TMP_PATH}/nuclassembly_rep.dbtype"; then

    CLUST_INPUT="${TMP_PATH}/nuclassembly"
    if notExists "${TMP_PATH}/clu.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" linclust "${CLUST_INPUT}" "${TMP_PATH}/clu" "${TMP_PATH}/clu_tmp" ${CLUSTER_PAR} \
            || fail "Redundancy reduction step died"
    fi

    if notExists "${TMP_PATH}/${CLUST_INPUT}_rep"; then
        # shellcheck disable=SC2086
        "$MMSEQS" result2repseq "${CLUST_INPUT}" "${TMP_PATH}/clu" "${CLUST_INPUT}_rep" ${THREADS_PAR} \
            || fail "Result2repseq  died"
    fi
fi

if notExists "${CLUST_INPUT}_rep_cycle.index" && [ -f "${TMP_PATH}/nuclassembly_cycle.index" ]; then
    awk 'NR == FNR { f[$1] = $0; next } $1 in f { print $0 }' "${TMP_PATH}/nuclassembly_cycle.index" "${CLUST_INPUT}_rep.index" > "${CLUST_INPUT}_rep_cycle.index"
fi

if notExists "${TMP_PATH}/nuclassembly_rep_h.dbtype"; then
    # shellcheck disable=SC2086
    if [ -f "${TMP_PATH}/nuclassembly_rep_cycle.index" ]; then
        "$MMSEQS" createhdb "${TMP_PATH}/nuclassembly_rep" "${TMP_PATH}/nuclassembly_rep_cycle" "${TMP_PATH}/nuclassembly_rep" ${VERBOSITY_PAR} \
            || fail "createhdb failed"
    else
        "$MMSEQS" createhdb "${TMP_PATH}/nuclassembly_rep" "${TMP_PATH}/nuclassembly_rep" ${VERBOSITY_PAR} \
            || fail "createhdb failed"
    fi
fi

if notExists "${TMP_PATH}/nuclassembly_rep.fasta"; then
    # shellcheck disable=SC2086
    "$MMSEQS" convert2fasta "${TMP_PATH}/nuclassembly_rep" "${TMP_PATH}/nuclassembly_rep.fasta" ${VERBOSITY_PAR} \
        || fail "convert2fasta died"
fi

mv -f "${TMP_PATH}/nuclassembly_rep.fasta" "$OUT_FILE" \
    || fail "Could not move result to $OUT_FILE"

#mv -f "${TMP_PATH}/assembly_aa_${STEP}" "${2}_aa" || fail "Could not move result to $2"
#mv -f "${TMP_PATH}/assembly_aa_${STEP}.index" "${2}_aa.index" || fail "Could not move result to $2.index"


if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files"
    "$MMSEQS" rmdb "${TMP_PATH}/nucl_reads"
    "$MMSEQS" rmdb "${TMP_PATH}/nucl_reads_h"
    rm -f "${TMP_PATH_GUIDED_ASSEMBLY}/aa_6f_"*
    rm -f "${TMP_PATH_GUIDED_ASSEMBLY}/nucl_6f_"*
    rm -f "${TMP_PATH_GUIDED_ASSEMBLY}/pref_"*
    rm -f "${TMP_PATH_GUIDED_ASSEMBLY}/aln_"*
    rm -f "${TMP_PATH_GUIDED_ASSEMBLY}/assembly_"*
    "$MMSEQS" rmdb "${TMP_PATH}/guided_assembly.merged"
    rm -f "${TMP_PATH}/nuclassembly"*
    "$MMSEQS" rmdb "${TMP_PATH}/clu"
    rm -f "${TMP_PATH}/guidedNuclAssemble.sh"
fi
