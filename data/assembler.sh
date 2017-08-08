#!/bin/bash
# Assembler workflow script
calc() { awk "BEGIN{print int($*) }";}

checkReturnCode () { 
	[ $? -ne 0 ] && echo "$1" && exit 1;
}
notExists () { 
	[ ! -f "$1" ] 
}
# check amount of input variables
[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[   -f "$2" ] &&  echo "$2 exists already!" && exit 1;
[ ! -d "$3" ] &&  echo "tmp directory $3 not found!" && mkdir -p $3;

export OMP_PROC_BIND=TRUE

INPUT="$1"
STEP=0
[ -z "$NUM_IT" ] && NUM_IT=1;
TOTAL_INPUT_CNT=$(wc -l ${INPUT}".index"|awk '{print $1}')
M=$KMER_PER_SEQ
while [ $STEP -lt $NUM_IT ]; do
    echo "STEP: "$STEP
    # 1. Finding exact $k$-mer matches.

    PARAM=KMERMATCHER${STEP}_PAR
    notExists "$3/pref_$STEP"          && $MMSEQS kmermatcher "$INPUT" "$3/pref_$STEP" ${!PARAM} --kmer-per-seq $M  && checkReturnCode "Kmer matching step died"
    # 2. Ungapped alignment
    notExists "$3/aln_$STEP" && $MMSEQS rescorediagonal "$INPUT" "$INPUT" "$3/pref_$STEP" "$3/aln_$STEP" ${UNGAPPED_ALN_PAR} && checkReturnCode "Ungapped alignment step died"
    # 3. Assemble
    notExists "$3/assembly_$STEP"         && $MMSEQS assembleresults "$INPUT" "$3/aln_$STEP" "$3/assembly_$STEP" ${ASSEMBLE_RESULT_PAR}  && checkReturnCode "Assembly step died"
    CNT=$(wc -l "$3/assembly_${STEP}.index"|awk '{print $1}')
    M=$(calc $TOTAL_INPUT_CNT/$CNT*$KMER_PER_SEQ)
    echo $M
    echo $TOTAL_INPUT_CNT
    echo $CNT
    echo $KMER_PER_SEQ
    INPUT="$3/assembly_$STEP"
    let STEP=STEP+1
done
let STEP=STEP-1
# post processing
mv -f "$3/assembly_$STEP" "$2"
mv -f "$3/assembly_$STEP.index" "$2.index"
checkReturnCode "Could not move result to $2"

if [ -n "$REMOVE_TMP" ]; then
 echo "Remove temporary files"
 rm -f "$3/pref" "$3/pref.index"
 rm -f "$3/aln" "$3/aln.index"
 rm -f "$3/assembly" "$3/assembly.index"
fi
