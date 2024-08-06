#!/bin/bash

#Script to prepare phased haplotype sequences for each for one sample.

prefix=$1
PLOIDY=$2
genelist=$3

#Run bcftools to extract sequences


for((i=1;i<=$PLOIDY;i++)); do
    while read GEN; do
    samtools faidx "$prefix".fasta $GEN | bcftools consensus -H $i ./"$prefix".fasta.snps.whatshap.POLYPHASE.vcf.gz >> ./"$prefix".v"$i".phased.fasta; done < ./"$genelist"
done

exit


