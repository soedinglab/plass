#!/bin/bash
# Assembler workflow script
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

notExists "$3/nucl_6f_start"  && $MMSEQS extractorfs ${INPUT} "$3/nucl_6f_start" --orf-start-state 1 --orf-end-state 0 --longest-orf --min-length 30 --max-length 45 --max-gaps 0 && checkReturnCode "extractorfs start step died"
notExists "$3/aa_6f_start"  && $MMSEQS translatenucs "$3/nucl_6f_start" "$3/aa_6f_start" --add-orf-stop && checkReturnCode "translatenucs start step died"i
notExists "$3/nucl_6f_long"  && $MMSEQS extractorfs ${INPUT} "$3/nucl_6f_long" --longest-orf --min-length 45 --max-gaps 0 && checkReturnCode "extractorfs longest step died"
notExists "$3/aa_6f_long"  && $MMSEQS translatenucs "$3/nucl_6f_long" "$3/aa_6f_long" --add-orf-stop && checkReturnCode "translatenucs long step died"

#notExists "$3/aa_6f_start_end"  && $MMSEQS concatdbs "$3/aa_6f_start" "$3/aa_6f_end" "$3/aa_6f_start_end" && checkReturnCode "concatdbs start end step died"
notExists "$3/aa_6f_start_long"  && $MMSEQS concatdbs "$3/aa_6f_long" "$3/aa_6f_start" "$3/aa_6f_start_long" && checkReturnCode "concatdbs start long step died"
INPUT="$3/aa_6f_start_long"

STEP=0
[ -z "$NUM_IT" ] && NUM_IT=1;
#TOTAL_INPUT_CNT=$(wc -l ${INPUT}".index"|awk '{print $1}')
M=$KMER_PER_SEQ
MERGEDBSTR=""
while [ $STEP -lt $NUM_IT ]; do
    echo "STEP: "$STEP
    # 1. Finding exact $k$-mer matches.

    PARAM=KMERMATCHER${STEP}_PAR
    notExists "$3/pref_$STEP"          && $MMSEQS kmermatcher "$INPUT" "$3/pref_$STEP" ${!PARAM} --kmer-per-seq $M  && checkReturnCode "Kmer matching step died"
    #notExists "$3/pref_$STEP"          && $MMSEQS prefilter "$INPUT" "$INPUT" "$3/pref_$STEP" -s 2  && checkReturnCode "Kmer matching step died"
    #notExists "$3/pref_$STEP"          && $MMSEQS prefilter "$INPUT" "$INPUT" "$3/pref_$STEP"  -s 3 --max-seqs 10000 --min-ungapped-score 60 && checkReturnCode "Kmer matching step died"
    # 2. Ungapped alignment
    notExists "$3/aln_$STEP" && $MMSEQS rescorediagonal "$INPUT" "$INPUT" "$3/pref_$STEP" "$3/aln_$STEP" ${UNGAPPED_ALN_PAR} && checkReturnCode "Ungapped alignment step died"
    if [ $STEP -eq 0 ]; then
       notExists "$3/corrected_reads" && $MMSEQS findassemblystart "$INPUT" "$3/aln_$STEP" "$3/corrected_seqs" && checkReturnCode "Findassemblystart alignment step died"
       INPUT="$3/corrected_seqs"
       notExists "$3/aln_corrected_$STEP" && $MMSEQS rescorediagonal "$INPUT" "$INPUT" "$3/pref_$STEP" "$3/aln_corrected_$STEP" ${UNGAPPED_ALN_PAR} && checkReturnCode "Ungapped alignment step died"
       notExists "$3/assembly_$STEP"      && $MMSEQS assembleresults "$INPUT" "$3/aln_corrected_$STEP" "$3/assembly_$STEP" ${ASSEMBLE_RESULT_PAR}  && checkReturnCode "Assembly step died"
    else
      # 3. Assemble
      notExists "$3/assembly_$STEP"         && $MMSEQS assembleresults "$INPUT" "$3/aln_$STEP" "$3/assembly_$STEP" ${ASSEMBLE_RESULT_PAR}  && checkReturnCode "Assembly step died"
    fi
    #CNT=$(wc -l "$3/assembly_${STEP}.index"|awk '{print $1}')
    INPUT="$3/assembly_$STEP"
    MERGEDBSTR=" $3/assembly_$STEP"$MERGEDBSTR
    let STEP=STEP+1
done

let STEP=STEP-1
#echo $MERGEDBSTR
# merge databases
#$notExists "$3/assembly_${STEP}_merge" && $MMSEQS mergedbs $1 "$3/assembly_${STEP}_merge" $MERGEDBSTR && checkReturnCode "Merge databases step died"
# first line should be the longest assembled sequence
#notExists "$3/assembly_${STEP}_filter" && $MMSEQS filterdb "$3/assembly_${STEP}_merge" "$3/assembly_${STEP}_filter" --extract-lines 1 && checkReturnCode "Filter database step died"
# remove entries with just null bytes
#awk '$3 > 3 {print $0}' $3/assembly_${STEP}_filter.index > $3/assembly_${STEP}_filter.nonull.index
#mv $3/assembly_${STEP}_filter.nonull.index $3/assembly_${STEP}_filter.index
# post processing
notExists "$3/assembly_${STEP}_filtered"      && $MMSEQS filternonecoding "$3/assembly_${STEP}" "$3/assembly_${STEP}_filtered"  && checkReturnCode "Filter protein with NN step died"

mv -f "$3/assembly_${STEP}_filtered" "$2"
#mv -f "$3/assembly_${STEP}" "$2"
checkReturnCode "Could not move result to $2"
mv -f "$3/assembly_${STEP}_filtered.index" "$2.index"
#mv -f "$3/assembly_${STEP}.index" "$2.index"
checkReturnCode "Could not move result to $2.index"

if [ -n "$REMOVE_TMP" ]; then
 echo "Remove temporary files"
 rm -f "$3/pref_*"
 rm -f "$3/aln_*"
 rm -f "$3/assembly*"
fi
