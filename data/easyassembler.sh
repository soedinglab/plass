#!/bin/sh -e
# Frame for all assembler workflows
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
        "$MMSEQS" createdb "$@" "${TMP_PATH}/nucl_reads" ${CREATEDB_PAR} \
            || fail "createdb failed"
    fi
fi

INPUT="${TMP_PATH}/nucl_reads"
if notExists "${TMP_PATH}/assembly.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" "${ASSEMBLY_MODULE}" "${INPUT}" "${TMP_PATH}/assembly" "${TMP_PATH}/assembly_tmp" ${ASSEMBLY_PAR} \
        || fail "Assembly died"
fi

if notExists "${TMP_PATH}/assembly_h.dbtype"; then
    # shellcheck disable=SC2086
    if [ -f "${TMP_PATH}/assembly_cycle.index" ]; then
    "$MMSEQS" createhdb "${TMP_PATH}/assembly" "${TMP_PATH}/assembly_cycle" "${TMP_PATH}/assembly" ${VERBOSITY_PAR} \
            || fail "createhdb failed"
    else
        "$MMSEQS" createhdb "${TMP_PATH}/assembly" "${TMP_PATH}/assembly" ${VERBOSITY_PAR} \
            || fail "createhdb failed"
    fi
fi

if notExists "${TMP_PATH}/assembly.fasta"; then
    # shellcheck disable=SC2086
    "$MMSEQS" convert2fasta "${TMP_PATH}/assembly" "${TMP_PATH}/assembly.fasta" ${VERBOSITY_PAR} \
        || fail "convert2fasta died"
fi

mv -f "${TMP_PATH}/assembly.fasta" "$OUT_FILE" \
    || fail "Could not move result to $OUT_FILE"


if [ -n "${REMOVE_TMP}" ]; then
    echo "Removing temporary files"
    "$MMSEQS" rmdb "${TMP_PATH}/assembly"
    "$MMSEQS" rmdb "${TMP_PATH}/assembly_h"
    "$MMSEQS" rmdb "${TMP_PATH}/nucl_reads"
    "$MMSEQS" rmdb "${TMP_PATH}/nucl_reads_h"
    rm -rf "${TMP_PATH}/assembly_tmp"
    rm -f "${TMP_PATH}/easyassembler.sh"
fi