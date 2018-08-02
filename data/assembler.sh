#!/bin/bash -e

# Assembler workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
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
if notExists "${TMP_PATH}/nucl_6f_start"; then
    $MMSEQS extractorfs "${INPUT}" "${TMP_PATH}/nucl_6f_start" --contig-start-mode 1 --contig-end-mode 0 --orf-start-mode 0 --min-length 30 --max-length 45 --max-gaps 0 \
        || fail "extractorfs start step died"
fi


if notExists "${TMP_PATH}/aa_6f_start"; then
    $MMSEQS translatenucs "${TMP_PATH}/nucl_6f_start" "${TMP_PATH}/aa_6f_start" --add-orf-stop \
        || fail "translatenucs start step died"
fi

if notExists "${TMP_PATH}/nucl_6f_long"; then
    $MMSEQS extractorfs ${INPUT} "${TMP_PATH}/nucl_6f_long" --orf-start-mode 0 --min-length 45 --max-gaps 0 \
        || fail "extractorfs longest step died"
fi

if notExists "${TMP_PATH}/aa_6f_long"; then
    $MMSEQS translatenucs "${TMP_PATH}/nucl_6f_long" "${TMP_PATH}/aa_6f_long" --add-orf-stop \
        || fail "translatenucs long step died"
fi

# $MMSEQS concatdbs "${TMP_PATH}/aa_6f_start" "${TMP_PATH}/aa_6f_end" "${TMP_PATH}/aa_6f_start_end"
if notExists "${TMP_PATH}/aa_6f_start_long"; then
    $MMSEQS concatdbs "${TMP_PATH}/aa_6f_long" "${TMP_PATH}/aa_6f_start" "${TMP_PATH}/aa_6f_start_long" \
        || fail "concatdbs start long step died"
fi

# $MMSEQS concatdbs "${TMP_PATH}/aa_6f_start" "${TMP_PATH}/aa_6f_end" "${TMP_PATH}/aa_6f_start_end"
if notExists "${TMP_PATH}/aa_6f_start_long_h"; then
    $MMSEQS concatdbs "${TMP_PATH}/aa_6f_long_h" "${TMP_PATH}/aa_6f_start_h" "${TMP_PATH}/aa_6f_start_long_h" \
        || fail "concatdbs start long step died"
fi

INPUT="${TMP_PATH}/aa_6f_start_long"
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

    if [ $STEP -eq 0 ]; then
        if notExists "${TMP_PATH}/corrected_seqs"; then
            $MMSEQS findassemblystart "$INPUT" "${TMP_PATH}/aln_$STEP" "${TMP_PATH}/corrected_seqs" \
        || fail "Findassemblystart alignment step died"
        fi
        INPUT="${TMP_PATH}/corrected_seqs"
        if notExists "${TMP_PATH}/aln_corrected_$STEP"; then
            $MMSEQS rescorediagonal "$INPUT" "$INPUT" "${TMP_PATH}/pref_$STEP" "${TMP_PATH}/aln_corrected_$STEP" ${UNGAPPED_ALN_PAR} \
        || fail "Ungapped alignment step died"
        fi
        if notExists "${TMP_PATH}/assembly_$STEP"; then
            $MMSEQS assembleresults "$INPUT" "${TMP_PATH}/aln_corrected_$STEP" "${TMP_PATH}/assembly_$STEP" ${ASSEMBLE_RESULT_PAR} \
        || fail "Assembly step died"
        fi
    else
      # 3. Assemble
        if notExists "${TMP_PATH}/assembly_$STEP"; then
            $MMSEQS assembleresults "$INPUT" "${TMP_PATH}/aln_$STEP" "${TMP_PATH}/assembly_$STEP" ${ASSEMBLE_RESULT_PAR} \
        || fail "Assembly step died"
        fi
    fi

    INPUT="${TMP_PATH}/assembly_$STEP"
    STEP=$(($STEP+1))
done
STEP=$(($STEP-1))

    # post processing
RESULT="${TMP_PATH}/assembly_${STEP}"
if [ ${PROTEIN_FILTER} -eq 1 ]; then
    RESULT="${TMP_PATH}/assembly_${STEP}_filtered"
    if notExists "${TMP_PATH}/assembly_${STEP}_filtered"; then
        $MMSEQS filternoncoding "${TMP_PATH}/assembly_${STEP}" "${TMP_PATH}/assembly_${STEP}_filtered" \
            || fail "Filter protein with NN step died"
    fi
fi


# select only assembled sequences
awk 'NR == FNR { f[$1] = $0; next } $1 in f { print f[$1], $0 }' ${RESULT}.index ${TMP_PATH}/aa_6f_start_long.index > ${RESULT}_tmp.index
awk '$3 > $6 { print }' ${RESULT}_tmp.index > ${RESULT}_only_assembled.index
mv ${RESULT}.index ${RESULT}_old.index
mv ${RESULT}_only_assembled.index ${RESULT}.index

# create fasta output
if notExists "${RESULT}_h"; then
    ln -s ${TMP_PATH}/aa_6f_start_long_h "${RESULT}_h"
fi
if notExists "${RESULT}_h.index"; then
    ln -s ${TMP_PATH}/aa_6f_start_long_h.index "${RESULT}_h.index"
fi

if notExists "${RESULT}.fasta"; then
    $MMSEQS convert2fasta "${RESULT}" "${RESULT}.fasta"
fi

mv -f "${RESULT}.fasta" "$OUT_FILE" || fail "Could not move result to $OUT_FILE"

if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files"
    rm -f "${TMP_PATH}/pref_"*
    rm -f "${TMP_PATH}/aln_"*
    rm -f "${TMP_PATH}/assembly_"*
fi
