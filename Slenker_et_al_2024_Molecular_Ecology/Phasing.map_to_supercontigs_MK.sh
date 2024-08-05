#!/bin/bash

# This workflow will take the supercontig output of HybPiper and return a supercontig that
# contains heterozygous positions as ambiguity bases. Uses paired reads.

# The script should be run on a FASTA file containing all the supercontigs of interest.



#############COMMAND LINE ARGUMENTS############

prefix=$1
read1fq=$2
read2fq=$3
PLD=$4

#####STEP ZERO: Make Reference Databases

java -jar $PICARD CreateSequenceDictionary REFERENCE=$prefix.fasta
bwa index $prefix.fasta
samtools faidx $prefix.fasta

#####STEP ONE: Map reads

echo "Mapping Reads"

bwa mem -M -t 2 $prefix.fasta $read1fq $read2fq | samtools sort -@ 2 -o $prefix.fasta.sorted.bam

samtools index $prefix.fasta.sorted.bam

java -Xmx2g -jar $PICARD FastqToSam  \
F1=$read1fq \
F2=$read2fq \
O=$prefix.fasta.unmapped.bam \
SM=$prefix.fasta

java -Xmx2g -jar $PICARD MergeBamAlignment \
ALIGNED=$prefix.fasta.sorted.bam \
UNMAPPED=$prefix.fasta.unmapped.bam \
O=$prefix.fasta.merged.bam \
R=$prefix.fasta

#####STEP TWO: Mark duplicates

echo "Marking Duplicates"
java -Xmx2g -jar $PICARD MarkDuplicates \
I=$prefix.fasta.merged.bam \
O=$prefix.fasta.marked.bam \
M=$prefix.fasta.metrics.txt

#######STEP THREE: Identify variants, select only SNPs

echo "Identifying variants"

samtools index $prefix.fasta.marked.bam

gatk --java-options "-Xmx4g" HaplotypeCaller  \
-R "$prefix.fasta" \
-I  $prefix.fasta.marked.bam \
-O $prefix.fasta.vcf \
-ploidy "$PLD" --native-pair-hmm-threads 2


gatk VariantFiltration \
-R "$prefix.fasta" \
-V $prefix.fasta.vcf \
-O $prefix.fasta.filtered.vcf \
--filter-expression 'QD < 2.0' --filter-name QDfilter \
--filter-expression 'DP < 8.0' --filter-name DPfilter \
--filter-expression 'MQ < 40.0' --filter-name MQfilter \
--filter-expression 'FS > 60.0' --filter-name FSfilter


gatk SelectVariants \
-R $prefix.fasta \
--select-type-to-include SNP \
--exclude-filtered \
-V $prefix.fasta.filtered.vcf \
-O $prefix.fasta.snps.vcf



# Preparation for creating masks
# REMOVE 1/1  a 0/0 from vcf, because otherwise it will make 'Ns' if it is located outside the main phased block
# make vcf with "N" as ALT base - for N masking step


case $PLD in
    2)
    NULY='0/0'
    JEDNOTKY='1/1'
    ;;
    3)
    NULY='0/0/0'
    JEDNOTKY='1/1/1'
    ;;
    4)
    NULY='0/0/0/0'
    JEDNOTKY='1/1/1/1'
    ;;
    5)
    NULY='0/0/0/0/0'
    JEDNOTKY='1/1/1/1/1'
    ;;
    6)
    NULY='0/0/0/0/0/0'
    JEDNOTKY='1/1/1/1/1/1'
    ;;
    *)
    echo EEEEEEEEEEEEEEEEEe
    ;;
esac

while read LINE; do
    if [[ $LINE == *"$NULY"* ]]; then
    echo "$LINE" | sed 's/\([0-9A-Za-z_.]\+\t[0-9]\+\t.\t\)\([A-Z]\+\t\)\([A-Z]\+\t\)/\1\2\2/' >> $prefix.fasta.snps.N.vcf
    echo "$LINE" | sed 's/\([0-9A-Za-z_.]\+\t[0-9]\+\t.\t\)\([A-Z]\+\t\)\([A-Z]\+\t\)/\1\2\2/' >> $prefix.fasta.snps.toIupac.vcf
    elif [[ $LINE == *"$JEDNOTKY"* ]]; then
    #echo "$LINE" | sed 's/\([0-9A-Za-z_.]\+\t[0-9]\+\t.\t\)\([A-Z]\+\t\)\([A-Z]\+\t\)/\1\2\3/' >> $prefix.fasta.snps.N.vcf
    #echo "$LINE" | sed 's/\([0-9A-Za-z_.]\+\t[0-9]\+\t.\t\)\([A-Z]\+\t\)\([A-Z]\+\t\)/\1\2\3/' >> $prefix.fasta.snps.toIupac.vcf
    echo "$LINE" >> $prefix.fasta.snps.N.vcf
    echo "$LINE" >> $prefix.fasta.snps.toIupac.vcf
    else # ostatne = N
    echo "$LINE" | sed 's/\([0-9A-Za-z_.]\+\t[0-9]\+\t.\t[A-Z]\+\t\)\([A-Z]\)/\1N/' >> $prefix.fasta.snps.N.vcf
    echo "$LINE" >> $prefix.fasta.snps.toIupac.vcf
    fi 
done < $prefix.fasta.snps.vcf




bgzip -c $prefix.fasta.snps.toIupac.vcf > $prefix.fasta.snps.toIupac.vcf.gz
bgzip -c $prefix.fasta.snps.N.vcf > $prefix.fasta.snps.N.vcf.gz
tabix -p vcf $prefix.fasta.snps.toIupac.vcf.gz
tabix -p vcf $prefix.fasta.snps.N.vcf.gz


######STEP FOUR: Output new supercontig FASTA with ambiguity codes

echo "Generating IUPAC FASTA file"


bcftools consensus --fasta-ref "$prefix".fasta --iupac-codes --output $prefix.fasta.iupac.fasta $prefix.fasta.snps.toIupac.vcf.gz

bcftools consensus --fasta-ref "$prefix".fasta --haplotype A --output $prefix.fasta.Nmasked.fasta $prefix.fasta.snps.N.vcf.gz



echo
echo
echo "Done!"

exit




