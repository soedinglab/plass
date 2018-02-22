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
[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[   -f "$2" ] &&  echo "$2 exists already!" && exit 1;
[ ! -d "$3" ] &&  echo "tmp directory $3 not found!" && mkdir -p $3;

INPUT="$(abspath $1)"
TMP_PATH="$(abspath $3)"

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

INPUT="${TMP_PATH}/aa_6f_start_long"
STEP=0
if [ -z "$NUM_IT" ]; then
    NUM_IT=1;
fi

M=$KMER_PER_SEQ
MERGEDBSTR=""
while [ $STEP -lt $NUM_IT ]; do
    echo "STEP: $STEP"

    # 1. Finding exact $k$-mer matches.
    PARAM="KMERMATCHER${STEP}_PAR"
    if notExists "${TMP_PATH}/pref_$STEP"; then
        $MMSEQS kmermatcher "$INPUT" "${TMP_PATH}/pref_$STEP" --kmer-per-seq "$M" ${!PARAM} \
        || fail "Kmer matching step died"
    fi

    # 2. Ungapped alignment
    if notExists "${TMP_PATH}/aln_$STEP"; then
        $MMSEQS rescorediagonal "$INPUT" "$INPUT" "${TMP_PATH}/pref_$STEP" "${TMP_PATH}/aln_$STEP" ${UNGAPPED_ALN_PAR} \
        || fail "Ungapped alignment step died"
    fi

    if [ $STEP -eq 0 ]; then
        if notExists "${TMP_PATH}/corrected_reads"; then
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
    MERGEDBSTR=" ${TMP_PATH}/assembly_$STEP"$MERGEDBSTR
    STEP=$(($STEP+1))
done
STEP=$(($STEP-1))

# post processing
if notExists "${TMP_PATH}/assembly_${STEP}_filtered"; then
    $MMSEQS filternoncoding "${TMP_PATH}/assembly_${STEP}" "${TMP_PATH}/assembly_${STEP}_filtered" \
        || fail "Filter protein with NN step died"
fi

mv -f "${TMP_PATH}/assembly_${STEP}_filtered" "$2" || fail "Could not move result to $2"
mv -f "${TMP_PATH}/assembly_${STEP}_filtered.index" "$2.index" || fail "Could not move result to $2.index"

if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files"
    rm -f "${TMP_PATH}/pref_*"
    rm -f "${TMP_PATH}/aln_*"
    rm -f "${TMP_PATH}/assembly*"
fi
