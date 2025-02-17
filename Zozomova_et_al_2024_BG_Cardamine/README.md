## Variant calling & filtration
Genome of C. amara ([GCA_040955855.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_040955855.1/)) was indexed, 
```ruby
samtools faidx GCA_040955855.1_C_amara_ONT_v2_genomic.fna
java -jar $PICARD CreateSequenceDictionary R="GCA_040955855.1_C_amara_ONT_v2_genomic.fna"
bwa index "GCA_040955855.1_C_amara_ONT_v2_genomic.fna"
```

and reads were mapped to reference using BWA
```ruby
bwa mem  -t $NCUP GCA_040955855.1_C_amara_ONT_v2_genomic.fna "$SAMPLE"*R1*f*q* "$SAMPLE"*R2*f*q* | samtools view -bu | samtools sort -l 9 -o "$SAMPLE".paired.bam
```
and further processed by PICARD.
```ruby
RUNNUMBER="1"
RGLB="$SAMPLE".lib1
RGPU="${RUNNUMBER}".unit1
RGSM="$SAMPLE"

mkdir tmp
java -Xmx16g -XX:+UseSerialGC -Djava.io.tmpdir=tmp -jar "${PICARD}" AddOrReplaceReadGroups INPUT="$SAMPLE".paired.bam OUTPUT="$SAMPLE".rg.bam RGID="${RUNNUMBER}" RGLB="${RGLB}" RGPL="illumina" RGPU="${RGPU}" RGSM="${RGSM}" 

java -Xmx16g -XX:+UseSerialGC -Djava.io.tmpdir=tmp -jar "${PICARD}" BuildBamIndex INPUT="$SAMPLE".rg.bam 

limit=$(echo `ulimit -n` - 50 | bc)
java -Xmx16g -XX:+UseSerialGC -Djava.io.tmpdir=tmp -jar "${PICARD}" MarkDuplicates I="$SAMPLE".rg.bam O="$SAMPLE".dedup.bam MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=$limit M=dup_metrics.log ASSUME_SORTED=true TAGGING_POLICY=All

java -Xmx16g -XX:+UseSerialGC -Djava.io.tmpdir=tmp -jar "${PICARD}" BuildBamIndex INPUT="$SAMPLE".dedup.bam
```

Variant calling was performed for each individual (specifying ploidy level) using modules from the GATK 4.4.0.0.

```ruby
gatk --java-options -Xmx16g HaplotypeCaller -R GCA_040955855.1_C_amara_ONT_v2_genomic.fna -I "$SAMPLE".dedup.bam  -ERC GVCF -ploidy $PLOIDY --min-base-quality-score 20 --max-genotype-count 350 -O "$SAMPLE".gvcf.gz

```

Next, all single-sample GVCFs were imported into GenomicsDB. 

Genomic intervals (-L) were equivalent to contigs. Each interval/contig was processed in a separate job.

```ruby
SAMPLELIST=$(find . -name "*.gvcf.gz" | sed 's/^\.\///' | sed 's/^/-V /' | tr "\n" " ")

mkdir tmp
gatk --java-options "-Xmx80g -XX:+UseSerialGC" GenomicsDBImport $SAMPLELIST \
--genomicsdb-workspace-path db."$INTERVAL" \
-L "$INTERVAL" \
--batch-size 50 --tmp-dir ./tmp --reader-threads 4 >GenomicsDBImport.log.out 2>GenomicsDBImport.log_error.out

gatk --java-options "-Xmx80g -XX:+UseSerialGC" GenotypeGVCFs -R GCA_040955855.1_C_amara_ONT_v2_genomic.fna \
-V gendb://db."$INTERVAL" \
-O "$INTERVAL".vcf.gz \
-L "$INTERVAL" \
--tmp-dir ./tmp  >GenotypeGVCFs.log.out 2>GenotypeGVCFs.log_error.out
```

All 944 VCFs, i.e. genotyped contigs, were concatenated using bcftools.
```ruby
bcftools concat -O z *.vcf.gz > concat.vcf.gz
```

Filtering variant calls was done in the GATK 4.4.0.0.
```ruby
gatk --java-options "-Xmx70g" SelectVariants -V concat.vcf.gz -O concat.bialelic.vcf.gz --restrict-alleles-to BIALLELIC -select-type SNP --exclude-non-variants

gatk --java-options "-Xmx70g" VariantFiltration \
-V concat.bialelic.vcf.gz -O concat.bialelic.filtered.DP8.vcf.gz \
--filter-expression 'QD < 2.0'  --filter-name 'QD' --filter-expression 'FS > 60.0' --filter-name 'FS' \
--filter-expression 'MQ < 40.0' --filter-name 'MQ' --filter-expression 'MQRankSum < -12.5' --filter-name 'MQRS' \
--filter-expression 'ReadPosRankSum < -8.0' --filter-name 'RPRS' --filter-expression 'SOR > 3.0' --filter-name 'SOR' \
--genotype-filter-expression 'DP < 8' --genotype-filter-name 'DP' --set-filtered-genotype-to-no-call

gatk --java-options "-Xmx70g" SelectVariants \
-V concat.bialelic.filtered.DP8.vcf.gz -O concat.bialelic.filtered.DP8.passed.vcf.gz \
--exclude-filtered  --max-nocall-fraction 0.6
```

To determine missing data per sample, VCF file was converted to phylip format using [vcf2phylip.py](https://github.com/edgardomortiz/vcf2phylip/blob/master/vcf2phylip.py) script and the percentage of missing data was calculated as the number of `N` divided by the total number of SNPs. The number of Ns was calculated by the following command:
```ruby
function countchar()
{
    while IFS= read -r i; do printf "%s" "$i" | tr -dc "$1" | wc -m; done
}
countchar 'N' <concat.bialelic.filtered.DP8.passed.vcf.min4.phy
```


## Maximum likelihood (ML) tree
The [vcf2phylip.py](https://github.com/edgardomortiz/vcf2phylip/blob/master/vcf2phylip.py) script was used to transform the data from the VCF file to the PHYLIP, and invariant sites were removed with the script [ascbias.py](https://github.com/btmartin721/raxml_ascbias).
The ML tree was constructed by RAxML-NG v.0.9.0, employing GTR model with Felsenstein’s ascertainment bias correction using script [MLTree.1.bestTree.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/MLTree.1.bestTree.sh).
Bootstrap analyses were performed using 500 replicates [MLTree.2.BS_trees.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/MLTree.2.BS_trees.sh), and the final tree with BS support was inferred using 
```ruby
cat *.raxml.bootstraps > allbootstraps.bootstraps

raxml-ng --support \
--tree concat.bialelic.filtered.DP8.passed.vcf.min4.ascbias_FELS.raxml.bestTree \
--bs-trees allbootstraps.bootstraps \
--prefix concat.bialelic.filtered.DP8.passed.vcf.min4.ascbias_FELS.raxml --threads 1 
```
