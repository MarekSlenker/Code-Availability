#!/bin/bash




FASTQ2="${FASTQ1//\_1\./_2.}" # Input, reverse (forward is in ${FASTQ1})
FASTQ="${FASTQ1%_1.f*q*}" # Base name
TRM1="$(basename "${FASTQ}.trm.R1.fq")" # Trimmed, forward
TRM2="$(basename "${FASTQ}.trm.R2.fq")" # Trimmed, reverse
FASTQBASE="$(basename "${FASTQ}")" # Ensure to get name of the directory containing pair of sequences only (not full path)



### fastq-pair: Rewrite paired end fastq files to make sure that all reads have a mate and to separate out singletons.
<PATH>/bin/fastq-pair-1.0/bin/fastq_pair -p $FASTQ1 $FASTQ2

fastp --in1 $FASTQ1.paired.fq --in2 $FASTQ2.paired.fq --out1 "$TRM1" --out2 "$TRM2" \
    --length_required 50 --html "$FASTQ".html --thread $PBS_NUM_PPN --trim_poly_g --qualified_quality_phred 20 --length_required 50

fastqc -o "${FASTQBASE}" -t 18 "${TRM1}" "${TRM2}"


