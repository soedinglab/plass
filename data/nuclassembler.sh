#!/bin/bash -e

# Assembler workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [[ $1 == */* ]]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d $(dirname "$1") ]; then
            echo "$(cd $(dirname "$1"); pwd)/$(basename "$1")"
    fi
}


# check amount of input variables
[ ! -n "${READ_FILES}" ] && echo "Please provide ${READ_FILES}" && exit 1;
[ ! -n ${OUT_FILE} ] && echo "Please provide ${OUT_FILE}" && exit 1;
[ ! -n ${TMP_PATH} ] && echo "Please provide ${TMP_PATH}" && exit 1;

# check if files exists
[   -f "${OUT_FILE}" ] &&  echo "${OUT_FILE} exists already!" && exit 1;
[ ! -d "${TMP_PATH}" ] &&  echo "tmp directory ${TMP_PATH} not found!" && mkdir -p ${TMP_PATH};


if notExists "${TMP_PATH}/nucl_reads"; then
    if [ ${PAIRED_END} -eq 1 ]; then
        echo "PAIRED END MODE"
        $MMSEQS mergereads $READ_FILES "${TMP_PATH}/nucl_reads"
    else
        $MMSEQS createdb $READ_FILES "${TMP_PATH}/nucl_reads"
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
        $MMSEQS kmermatcher "$INPUT" "${TMP_PATH}/pref_$STEP" ${KMERMATCHER_PAR}   \
        || fail "Kmer matching step died"
    fi

    # 2. Ungapped alignment
    if notExists "${TMP_PATH}/aln_$STEP"; then
        $MMSEQS rescorediagonal "$INPUT" "$INPUT" "${TMP_PATH}/pref_$STEP" "${TMP_PATH}/aln_$STEP" ${UNGAPPED_ALN_PAR} \
        || fail "Ungapped alignment step died"
    fi

    # 3. Assemble
    if notExists "${TMP_PATH}/assembly_$STEP"; then
        $MMSEQS assembleresults "$INPUT" "${TMP_PATH}/aln_$STEP" "${TMP_PATH}/assembly_$STEP" ${ASSEMBLE_RESULT_PAR} \
    || fail "Assembly step died"
    fi

    INPUT="${TMP_PATH}/assembly_$STEP"
    STEP=$(($STEP+1))
done
STEP=$(($STEP-1))


# select only assembled sequences
RESULT="${TMP_PATH}/assembly_${STEP}"
awk 'NR == FNR { f[$1] = $0; next } $1 in f { print f[$1], $0 }' ${RESULT}.index "${TMP_PATH}/nucl_reads.index" > ${RESULT}_tmp.index
awk '$3 > $6 { print }' ${RESULT}_tmp.index > ${RESULT}_only_assembled.index
mv ${RESULT}.index ${RESULT}_old.index
mv ${RESULT}_only_assembled.index ${RESULT}.index

# create fasta output
if notExists "${RESULT}_h"; then
    ln -s "${TMP_PATH}/nucl_reads_h" "${RESULT}_h"
fi
if notExists "${RESULT}_h.index"; then
    ln -s "${TMP_PATH}/nucl_reads_h.index" "${RESULT}_h.index"
fi
if notExists "${RESULT}.fasta"; then
    $MMSEQS convert2fasta "${RESULT}" "${RESULT}.fasta"
fi

echo "${RESULT}.fasta"

echo "$OUT_FILE"

mv -f "${RESULT}.fasta" "$OUT_FILE" || fail "Could not move result to $OUT_FILE"


if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files"
    rm -f "${TMP_PATH}/pref_"*
    rm -f "${TMP_PATH}/aln_"*
    rm -f "${TMP_PATH}/assembly_"*
fi
