This repository stores scripts and commands used for the analysis of RADseq data for the paper [Kantor et al. 2026]().


## Demultiplexing, quality filtering and deduplication
Raw reads were demultiplexed using [rad_1_demultiplexing_1_qsub.sh](https://github.com/V-Z/RAD-Seq-scripts/blob/master/bin/rad_1_demultiplexing_1_qsub.sh) and [rad_1_demultiplexing_2_run.sh](https://github.com/V-Z/RAD-Seq-scripts/blob/master/bin/rad_1_demultiplexing_2_run.sh) scripts.  

Next, [fastq pair](https://github.com/linsalrob/fastq-pair) was used to rewrite paired-end fastq files to ensure all reads have a mate and separate singletons.
```ruby
fastq_pair -p $FASTQ1 $FASTQ2
```

The [fastp](https://github.com/OpenGene/fastp) was used to filter low-quality reads.
```ruby
fastp --in1 $FASTQ1.paired.fq --in2 $FASTQ2.paired.fq --out1 "$FQ".trm.R1.fq --out2 "$FQ".trm.R2.fq \
    --length_required 50 --html "$FASTQ".html --thread $PBS_NUM_PPN --trim_poly_g --qualified_quality_phred 20 --length_required 50
```

and [BBMap](https://sourceforge.net/projects/bbmap/) for deduplication
```ruby
clumpify.sh -Xmx"${JAVAMEM}" in="$FQ".trm.R1.fq in2="$FQ".trm.R2.fq out="$FQ".dedup.R1.fq out2="$FQ".dedup.R2.fq dedupe optical spany adjacent 
```

## Variant calling & filtration
Genome of *Cardamine amara* (GenBank: [GCA_040955855.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_040955855.1/)) was indexed, 

```ruby
REF="GCA_040955855.1_C_amara_ONT_v2_genomic.fna"
samtools faidx "$REF"
java -jar $PICARD CreateSequenceDictionary R="$REF"
bwa index "$REF"
```

and reads were mapped to the reference using BWA
```ruby
bwa mem  -t "$PBS_NUM_PPN" -R "@RG\tID:run1\tLB:"$FQ".lib1\tPL:ILLUMINA\tSM:$FQ" "$REF" "$FQ".dedup*R1*f*q "$FQ".dedup*R2*f*q | samtools view -bu -@ "$PBS_NUM_PPN" | samtools sort -l 9 -@ "$PBS_NUM_PPN" -o "$FQ".paired.bam

java -Xmx"${JAVAMEM}" -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar "${PICARD}" BuildBamIndex INPUT="$FQ".paired.bam
```

Variant calling was performed for each individual (specifying ploidy level) using modules from the GATK 4.4.0.0.

```ruby
gatk --java-options "-Djava.io.tmpdir="${SCRATCHDIR}"/tmp -Xmx"$JAVAMEM" -XX:ParallelGCThreads=2" HaplotypeCaller \
  -R $REF \
  -I "$FQ".paired.bam  \
  -ploidy $PLIODY \
  -ERC GVCF \
  --min-base-quality-score 20 --max-genotype-count 350 \
  -O $FQ.gvcf.gz
```

Next, all single-sample GVCFs were imported into GenomicsDB. 

Genomic intervals (-L) were equivalent to contigs. Each interval/contig was processed in a separate job.

```ruby
SAMPLELIST=$(find . -name "*.gvcf.gz" | sed 's/^\.\///' | sed 's/^/-V /' | tr "\n" " ")

mkdir tmp
gatk --java-options "-Xmx"$MEM"g -XX:+UseSerialGC" GenomicsDBImport \
$SAMPLELIST --genomicsdb-workspace-path db."$INTERVAL" \
-L "$INTERVAL" \
--batch-size 50 --tmp-dir ./tmp --reader-threads 4 >GenomicsDBImport.log.out 2>GenomicsDBImport.log_error.out

gatk --java-options "-Xmx80g -XX:+UseSerialGC" GenotypeGVCFs \
-R $REF \
-V gendb://db."$INTERVAL" \
-O "$INTERVAL".vcf.gz \
-L "$INTERVAL" \
--tmp-dir ./tmp  >GenotypeGVCFs.log.out 2>GenotypeGVCFs.log_error.out
```

All 944 VCFs, i.e. genotyped contigs, were concatenated.
```ruby
SAMPLELIST=$(find . -name "*.vcf.gz" | sed 's/^\.\///' | sed 's/^/-I /' | tr "\n" " ")
gatk MergeVcfs $SAMPLELIST -O ../concat.vcf.gz
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
-V concat.bialelic.filtered.DP8.vcf.gz -O concat.bialelic.filtered.DP8.passed.m02.vcf.gz \
--exclude-filtered  --max-nocall-fraction 0.2
```





## Maximum likelihood (ML) tree
The [vcf2phylip.py](https://github.com/edgardomortiz/vcf2phylip/blob/master/vcf2phylip.py) script was used to transform the data from the VCF file to the PHYLIP, and invariant sites were removed with the script [ascbias.py](https://github.com/btmartin721/raxml_ascbias).
The ML tree was constructed by RAxML-NG v.0.9.0, employing GTR model with Felsenstein’s ascertainment bias correction using script [MLTree.1.bestTree.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Bacak_et_al_2025_Mol_Phylogenet_Evol/MLTree.1.bestTree.sh). Bootstrap analyses were performed using 500 replicates [MLTree.2.BS_trees.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/MLTree.2.BS_trees.sh), and the final tree with BS support was inferred using 


```ruby
cat *.raxml.bootstraps > allbootstraps.bootstraps

raxml-ng --support \
--tree concat.bialelic.filtered.DP8.passed.vcf.min4.ascbias_Felsenstein.raxml.bestTree \
--bs-trees allbootstraps.bootstraps \
--prefix concat.bialelic.filtered.DP8.passed.vcf.min4.ascbias_Felsenstein.raxml --threads 1 
```

## PCA
The principal component analysis (PCA) based on covariance matrix, as implemented in the R package adegenet was calculated using [PCA.R](https://github.com/MarekSlenker/Code-Availability/blob/main/Zozomova_et_al_2025_New_Phytol/PCA.R) script.


## Clonality
To assess clonality, we calculated Nei's genetic distances from a VCF file using [Neis_distances.diploidized.R](https://raw.githubusercontent.com/MarekSlenker/Code-Availability/refs/heads/main/Kantor_et_al_2026_Taxon/Neis_distances.diploidized.R) script. The frequency distributions of Nei's pairwise genetic distances were then examined. Genet assignment threshold following the approach proposed by Tsujimoto et al. (2020). According to this approach, the first peak in the frequency distribution of Nei's pairwise genetic distances indicates somatic mutations and genotypic errors within genets (genetically identical individuals resulting from clonal reproduction), followed by distances representing variation between the genets. 















