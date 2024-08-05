#!/bin/bash



export RETRIVEDSEQDIR   # input sequences
export DEDUPDIR         # fastq reads
export ALLELE_FREQ=0.3  # Minimum allele frequency regarded for assigning ambiguity code. If the alternative allele is less frequent than the chosen value, it is not included in ambiguity coding.
export READ_DEPTH=10    # Required coverage to call variants. Variants are only called when there is a coverage of at least 10.
export ALLELE_COUNT=4   # Minimum count of alleles to be regarded for assigning ambiguity code.


module add samtools-1.9
module add picard-2.9.0
module add bwa-0.7.5a
module add gatk-4.1.6.0
module add parallel-20160622
module add bcftools-1.11


cd "$SCRATCHDIR"

# collect all sequences of a particular sample
for seq in $RETRIVEDSEQDIR/*fasta; do
if grep -q "$SAMPLE" $seq 
then 
echo ">"$(basename $seq) >> "$SAMPLE".ref.fasta
grep "$SAMPLE" $seq -A1 --no-filename | tail -n 1 >> "$SAMPLE".ref.fasta
fi
done

# indexing..
sed -i 's/-//g' "$SAMPLE".ref.fasta
bwa index "$SAMPLE".ref.fasta
samtools faidx "$SAMPLE".ref.fasta
java -jar $PICARD CreateSequenceDictionary R="$SAMPLE".ref.fasta

# copy fastq reads from $DEDUPDIR to current working directory
cp "$DEDUPDIR"/"$SAMPLE".dedup.R1.fq.bz2 .
cp "$DEDUPDIR"/"$SAMPLE".dedup.R2.fq.bz2 .
bunzip2 *.bz2

# map reads
bwa mem -R "@RG\tID:4\tLB:lib1\tPL:ILLUMINA\tSM:$SAMPLE" "$SAMPLE".ref.fasta "$SAMPLE".dedup.R1.fq "$SAMPLE".dedup.R2.fq -t 2 | samtools sort > "$SAMPLE".bam
samtools index "$SAMPLE".bam


# call variants
gatk --java-options "-Xmx4g" HaplotypeCaller  \
    -R "$SAMPLE".ref.fasta \
    -I "$SAMPLE".bam \
    -O "$SAMPLE".raw.vcf \
    -ploidy 2

gatk VariantFiltration \
    -R "$SAMPLE".ref.fasta \
    -V "$SAMPLE".raw.vcf \
    -O "$SAMPLE".filtered.vcf \
    --filter-expression 'QD < 2.0' --filter-name QDfilter \
    --filter-expression 'DP < 8.0' --filter-name DPfilter \
    --filter-expression 'MQRankSum < -12.5' --filter-name MQRS \
    --filter-expression 'ReadPosRankSum < -8.0' --filter-name RPRS \
    --filter-expression 'SOR > 3.0' --filter-name SOR \
    --filter-expression 'MQ < 40.0' --filter-name MQfilter \
    --filter-expression 'FS > 60.0' --filter-name FSfilter 


gatk SelectVariants \
    -R "$SAMPLE".ref.fasta \
    -V "$SAMPLE".raw.vcf \
    -O "$SAMPLE".passed.vcf \
    --exclude-filtered -select-type SNP

bgzip "$SAMPLE".passed.vcf
tabix -f "$SAMPLE".passed.vcf.gz


# and finaly create sequences, with SNPs coded with iupac ambiguity codes
bcftools consensus -I -i "(AD[:1] / FORMAT/DP) >= $ALLELE_FREQ && FORMAT/DP >= $READ_DEPTH && (AD[:1]) >= $ALLELE_COUNT" -f "$SAMPLE".ref.fasta "$SAMPLE".passed.vcf.gz > "$SAMPLE".consensus.fasta


exit
