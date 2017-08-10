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
MERGEDBSTR=""
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
    INPUT="$3/assembly_$STEP"
    MERGEDBSTR=" $3/assembly_$STEP"$MERGEDBSTR
    let STEP=STEP+1
done

let STEP=STEP-1
echo $MERGEDBSTR 
# merge databases 
notExists "$3/assembly_${STEP}_merge" && mmseqs mergedbs $1 "$3/assembly_${STEP}_merge" $MERGEDBSTR && checkReturnCode "Merge databases step died"
# first line should be the longest assembled sequence
notExists "$3/assembly_${STEP}_filter" && mmseqs filterdb "$3/assembly_${STEP}_merge" "$3/assembly_${STEP}_filter" --extract-lines 1 && checkReturnCode "Filter database step died"
# remove entries with just null bytes
awk '$3 > 3 {print $0}' $3/assembly_${STEP}_filter.index > $3/assembly_${STEP}_filter.nonull.index
mv $3/assembly_${STEP}_filter.nonull.index $3/assembly_${STEP}_filter.index
# post processing
mv -f "$3/assembly_${STEP}_filter" "$2"
checkReturnCode "Could not move result to $2"
mv -f "$3/assembly_${STEP}_filter.index" "$2.index"
checkReturnCode "Could not move result to $2.index"

if [ -n "$REMOVE_TMP" ]; then
 echo "Remove temporary files"
 rm -f "$3/pref_*" 
 rm -f "$3/aln_*" 
 rm -f "$3/assembly*"
fi
