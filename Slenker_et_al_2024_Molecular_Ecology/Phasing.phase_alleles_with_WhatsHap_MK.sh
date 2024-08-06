#!/bin/bash


prefix=$1
PLD=$2


echo "POLYPHASE:"
echo

pip install whatshap --user --upgrade

python -m whatshap polyphase \
$prefix.fasta.snps.vcf \
$prefix.fasta.marked.bam \
--reference=$prefix.fasta \
--output $prefix.fasta.snps.POLYPHASE.whatshap.vcf \
--ploidy $PLD

python -m whatshap stats $prefix.fasta.snps.POLYPHASE.whatshap.vcf \
--gtf $prefix.whatshap.gtf \
--tsv $prefix.whatshap.stats.tsv

bgzip -c $prefix.fasta.snps.POLYPHASE.whatshap.vcf > $prefix.fasta.snps.whatshap.POLYPHASE.vcf.gz
tabix $prefix.fasta.snps.whatshap.POLYPHASE.vcf.gz



exit

