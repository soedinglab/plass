#!/bin/sh -e

PLASS="$1"
MMSEQS="$2"
BASEDIR="$3"

mkdir -p "${BASEDIR}"
wget -qO- http://wwwuser.gwdg.de/~compbiol/plass/plass_regression_data.tar.gz | tar -xzC "${BASEDIR}"

"${PLASS}" assemble "${BASEDIR}/allgenomes_reads_sample_1.fastq" "${BASEDIR}/allgenomes_reads_sample_2.fastq" "${BASEDIR}/final.contigs.aa.fa" "${BASEDIR}/tmp"

"${MMSEQS}" createdb "${BASEDIR}/final.contigs.aa.fa" "${BASEDIR}/final.contigs.aa"
"${MMSEQS}" createdb "${BASEDIR}/prochloroccus_allproteins.fasta" "${BASEDIR}/prochloroccus_allproteins"
"${MMSEQS}" createdb "${BASEDIR}/prochloroccus_allproteins_nr.fasta" "${BASEDIR}/prochloroccus_allproteins_nr"

len_distribution() {
    awk '{ n[$3]++ } END { for (i in n) print i,n[i] }' "$1" | sort -n
}
mapped_distribution() {
    awk -v lencut="$2" \
      'BEGIN{ group = ""; len = 0; cov = "" } $1 != group { if(len >= lencut) { n[cov *len]++; } group = $1; len = $8; cov = (1 + $7 - $6)/$8; } { cov = cov > (1 + $7 - $6)/$8 ? cov : (1 + $7 - $6)/$8} END{ if(len >= lencut) { n[cov*len]++; }; for (i in n) print i,n[i] }' "$1" | sort -n
}

calc() {
    awk "BEGIN { print $* }"
}
mapped_fraction() {
    SUM=$(len_distribution "$1" | awk -v len="$3" 'BEGIN { sum = 0 } $1 > len { sum += $1 * $2 } END { print sum }')
    ALIGNED=$(mapped_distribution "$2" "$3" | awk 'BEGIN { sum = 0 } { sum += $1 * $2 } END { print sum }')
    printf '%.0f\t%.0f\t%.3f\n' "$SUM" "$ALIGNED" "$(calc $ALIGNED/$SUM)"
}

evaluate() {
    ASSEMBLY="$1"
    REFERENCE="$2"
    REFERENCENR="$3"
    RESULT="$4"
    LEN="$5"

    awk -v len="$LEN" '$3 > len { print }' "${ASSEMBLY}.index" > "${ASSEMBLY}.ids"
    "${MMSEQS}" createsubdb "${ASSEMBLY}.ids" "${ASSEMBLY}" "${ASSEMBLY}.${LEN}"
    "${MMSEQS}" createsubdb "${ASSEMBLY}.ids" "${ASSEMBLY}_h" "${ASSEMBLY}.${LEN}_h"

    "${MMSEQS}" search "${ASSEMBLY}.${LEN}" "${REFERENCE}" "${RESULT}/assembly_against_reference" "${RESULT}/tmp" -s 5 --max-seqs 5000 --min-ungapped-score 100 -a --min-seq-id 0.89
    for i in $(seq 90 99 | awk '{ print $1/100 }'); do
        "${MMSEQS}" filterdb "${RESULT}/assembly_against_reference" "${RESULT}/assembly_against_reference_${i}" --filter-column 3 --comparison-value ${i} --comparison-operator ge
        "${MMSEQS}" createtsv "${ASSEMBLY}.${LEN}" "${REFERENCE}" "${RESULT}/assembly_against_reference_${i}" "${RESULT}/assembly_against_reference_${i}.tsv"
        mapped_fraction "${ASSEMBLY}.${LEN}.index" "${RESULT}/assembly_against_reference_${i}.tsv" "$LEN" >> "${RESULT}/precision"
    done

    # sens
    "${MMSEQS}" search "$REFERENCENR" "${ASSEMBLY}.${LEN}" "${RESULT}/reference_against_assembly" "${RESULT}/tmp" --max-seqs 500000 -a --min-seq-id 0.89
    for i in $(seq 90 99 | awk '{ print $1/100 }'); do
        "${MMSEQS}" filterdb "${RESULT}/reference_against_assembly" "${RESULT}/reference_against_assembly_${i}" --filter-column 3 --comparison-value $i --comparison-operator ge
        "${MMSEQS}" createtsv "$REFERENCENR" "${ASSEMBLY}.${LEN}" "${RESULT}/reference_against_assembly_${i}" "${RESULT}/reference_against_assembly_${i}.tsv"
        mapped_fraction "${REFERENCENR}.index" "${RESULT}/reference_against_assembly_${i}.tsv" "$LEN" >> "${RESULT}/sens"
    done
}

evaluate "${BASEDIR}/final.contigs.aa" "${BASEDIR}/prochloroccus_allproteins" "${BASEDIR}/prochloroccus_allproteins_nr" "${BASEDIR}" 100

cat "${BASEDIR}/sens" "${BASEDIR}/precision" > "${BASEDIR}/report"
cat "${BASEDIR}/report"

check() {
    REPORT=$1
    VALS=$2
    if [ ! -s "${REPORT}" ]; then
        exit 1
    fi

    GOOD=$(awk -v check="$VALS" 'BEGIN { split(check,checklist," "); cnt=0 } $3 >= (checklist[NR]-0.005) { cnt = cnt+1 } END { if(cnt == NR) { print "GOOD" } else { print "BAD" } }' "${REPORT}")
    if [ "$GOOD" != "GOOD" ]; then
        >&2 echo "Failed check!"
        exit 1
    fi
}
check "${BASEDIR}/report" "0.495 0.474 0.451 0.422 0.389 0.343 0.295 0.245 0.196 0.133 0.980 0.980 0.979 0.979 0.977 0.974 0.965 0.940 0.864 0.649"

