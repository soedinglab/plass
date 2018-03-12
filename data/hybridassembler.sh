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

INPUT="$(abspath "$1")"
TMP_PATH="$(abspath "$3")"
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

if notExists "${TMP_PATH}/nucl_6f_start_long"; then
    $MMSEQS concatdbs "${TMP_PATH}/nucl_6f_long" "${TMP_PATH}/nucl_6f_start" "${TMP_PATH}/nucl_6f_start_long" \
        || fail "concatdbs start long step died"
fi

INPUT_AA="${TMP_PATH}/aa_6f_start_long"
INPUT_NUCL="${TMP_PATH}/nucl_6f_start_long"
STEP=0
if [ -z "$NUM_IT" ]; then
    NUM_IT=1;
fi

while [ $STEP -lt $NUM_IT ]; do
    echo "STEP: $STEP"

    # 1. Finding exact $k$-mer matches.
    if notExists "${TMP_PATH}/pref_$STEP"; then
        $MMSEQS kmermatcher "$INPUT_AA" "${TMP_PATH}/pref_$STEP" ${KMERMATCHER_PAR}   \
        || fail "Kmer matching step died"
    fi

    # 2. Ungapped alignment
    if notExists "${TMP_PATH}/aln_$STEP"; then
        $MMSEQS rescorediagonal "$INPUT_AA" "$INPUT_AA" "${TMP_PATH}/pref_$STEP" "${TMP_PATH}/aln_$STEP" ${UNGAPPED_ALN_PAR} \
        || fail "Ungapped alignment step died"
    fi

    # 3. Ungapped alignment protein 2 nucl
    if notExists "${TMP_PATH}/aln_nucl_$STEP"; then
        $MMSEQS proteinaln2nucl "$INPUT_NUCL" "$INPUT_NUCL" "${TMP_PATH}/aln_$STEP" "${TMP_PATH}/aln_nucl_$STEP"  \
        || fail "Ungapped alignment step died"
    fi

    # 4. Assemble
    if notExists "${TMP_PATH}/assembly_aa_$STEP" || notExists "${TMP_PATH}/assembly_aa_$STEP"; then
        $MMSEQS hybridassembleresults "$INPUT_NUCL" "$INPUT_AA" "${TMP_PATH}/aln_nucl_$STEP" "${TMP_PATH}/assembly_nucl_$STEP" "${TMP_PATH}/assembly_aa_$STEP"  ${ASSEMBLE_RESULT_PAR} \
    || fail "Assembly step died"
    fi

    INPUT_AA="${TMP_PATH}/assembly_aa_$STEP"
    INPUT_NUCL="${TMP_PATH}/assembly_nucl_$STEP"
    STEP=$(($STEP+1))
done
STEP=$(($STEP-1))

mv -f "${TMP_PATH}/assembly_nucl_${STEP}" "${2}_nucl" || fail "Could not move result to $2"
mv -f "${TMP_PATH}/assembly_nucl_${STEP}.index" "${2}_nucl.index" || fail "Could not move result to $2.index"

mv -f "${TMP_PATH}/assembly_aa_${STEP}" "${2}_aa" || fail "Could not move result to $2"
mv -f "${TMP_PATH}/assembly_aa_${STEP}.index" "${2}_aa.index" || fail "Could not move result to $2.index"


if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files"
    rm -f "${TMP_PATH}/pref_"*
    rm -f "${TMP_PATH}/aln_"*
    rm -f "${TMP_PATH}/assembly_"*
fi
