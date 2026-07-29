This repository stores scripts and commands used for the analysis of RADseq data for the paper [Majerova et al. 2026]().

---

The ddRADseq reads were demultiplexed using [rad_1_demultiplexing_1_qsub.sh](https://github.com/V-Z/RAD-Seq-scripts/blob/master/bin/rad_1_demultiplexing_1_qsub.sh) and [rad_1_demultiplexing_2_run.sh](https://github.com/V-Z/RAD-Seq-scripts/blob/master/bin/rad_1_demultiplexing_2_run.sh) scripts.  

Next, [fastq pair](https://github.com/linsalrob/fastq-pair) was used to ensure all reads have a mate pare.
```ruby
fastq_pair -p $FASTQ1 $FASTQ2
```

The [fastp](https://github.com/OpenGene/fastp) was used to filter reads with quality scores below Q20, lengths less than 50 bp, and to trim putative polyG tails. 
```ruby
fastp --in1 $FASTQ1.paired.fq --in2 $FASTQ2.paired.fq --out1 "$FQ".trm.R1.fq --out2 "$FQ".trm.R2.fq \
    --length_required 50 --html "$FASTQ".html --thread $PBS_NUM_PPN --trim_poly_g --qualified_quality_phred 20 --length_required 50
```

Sequence duplicates were removed using clumpify.sh from [BBTools](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/). 
```ruby
clumpify.sh -Xmx"${JAVAMEM}" in="$FQ".trm.R1.fq in2="$FQ".trm.R2.fq out="$FQ".dedup.R1.fq out2="$FQ".dedup.R2.fq dedupe optical spany adjacent 
```

The demultiplexed, filtered, and deduplicated reads were mapped to the *Cardamine amara* genome (GenBank accession: [GCA_040955855.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_040955855.1/)) using BWA-MEM algorithm implemented in `BWA` 0.7.5a. 

```ruby
bwa mem  -t "$PBS_NUM_PPN" -R "@RG\tID:run1\tLB:"$FQ".lib1\tPL:ILLUMINA\tSM:$FQ" "$REF" "$FQ".dedup*R1*f*q "$FQ".dedup*R2*f*q | samtools view -bu -@ "$PBS_NUM_PPN" | samtools sort -l 9 -@ "$PBS_NUM_PPN" -o "$FQ".paired.bam

java -Xmx"${JAVAMEM}" -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar "${PICARD}" BuildBamIndex INPUT="$FQ".paired.bam
```

Variants were called for each individual using `HaplotypeCaller` module of `GATK` v4.4.0.0.

```ruby
gatk --java-options "-Djava.io.tmpdir="${SCRATCHDIR}"/tmp -Xmx"$JAVAMEM" -XX:ParallelGCThreads=2" HaplotypeCaller \
  -R $REF \
  -I "$FQ".paired.bam  \
  -ploidy $PLIODY \
  -ERC GVCF \
  --min-base-quality-score 20 --max-genotype-count 350 \
  -O $FQ.gvcf.gz
```

Next, single-sample GVCFs were imported into GenomicsDB. 

Genomic intervals (-L) were equivalent to chromosoms. Each chromosome was processed in a separate job.

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

All partial VCFs, were concatenated.

```ruby
SAMPLELIST=$(find . -name "*.vcf.gz" | sed 's/^\.\///' | sed 's/^/-I /' | tr "\n" " ")
gatk MergeVcfs $SAMPLELIST -O ../concat.vcf.gz
```

Variant calling and biallelic SNP filtering were performed using `SelectVariants` module of `GATK` with a minimum sequencing depth of 8×, following the filter parameters indicated by GATK's Best Practices.  

```ruby
gatk --java-options "-Xmx70g" SelectVariants -V concat.vcf.gz -O concat.bialelic.vcf.gz --restrict-alleles-to BIALLELIC -select-type SNP --exclude-non-variants

gatk --java-options "-Xmx70g" VariantFiltration \
-V concat.bialelic.vcf.gz -O concat.bialelic.filtered.DP8.vcf.gz \
--filter-expression 'QD < 2.0'  --filter-name 'QD' --filter-expression 'FS > 60.0' --filter-name 'FS' \
--filter-expression 'MQ < 40.0' --filter-name 'MQ' --filter-expression 'MQRankSum < -12.5' --filter-name 'MQRS' \
--filter-expression 'ReadPosRankSum < -8.0' --filter-name 'RPRS' --filter-expression 'SOR > 3.0' --filter-name 'SOR' \
--genotype-filter-expression 'DP < 8' --genotype-filter-name 'DP' --set-filtered-genotype-to-no-call

gatk --java-options "-Xmx70g" SelectVariants \
-V concat.bialelic.filtered.DP8.vcf.gz -O concat.bialelic.filtered.DP8.passed.m06.vcf.gz \
--exclude-filtered  --max-nocall-fraction 0.2
```



For analyses based on the adegenet genlight object, including principal component analysis, NeighborNet, and clonality assessment (see below), genotypes were encoded to reflect allele dosage of mixed-ploidy data. For genome-wide admixture computations (see below) requiring diploid genotypes, polyploid genotypes were converted to a pseudo-diploid format while preserving homozygous or heterozygous states of original polyploids. Specifically, genotypes such as AABB or AAAABB were recoded as AB, whereas homozygous genotypes such as AAAA or AAAAAA were recoded as AA.  


**The principal component analysis (PCA)**, based on a covariance matrix, as implemented in the R package adegenet was calculated using [PCA.R](https://github.com/MarekSlenker/Code-Availability/blob/main/Zozomova_et_al_2025_New_Phytol/PCA.R) script.

The **NeighborNet** analysis was based on Nei’s genetic distances, calculated using the R package `StAMPP` ([Neis_distances.R](https://github.com/MarekSlenker/Code-Availability/blob/main/Zozomova_et_al_2025_New_Phytol/Neis_distances.R)), and the resulting network was visualized in `SplitsTree4`. 



For the Bayesian clustering analyses implemented in `STRUCTURE` v.2.3.4, datasets of unlinked SNPs were generated to minimise linkage disequilibrium. First, RADseq loci separated by at least 1,000 bp were identified using the [identiﬁRadLoci.workﬂow](https://github.com/MarekSlenker/vcf_prune/blob/main/identifiRadLoci.workflow). Subsequently, 100 independent datasets were created by selecting a single random SNP per locus from loci harboring at least six SNPs, using the [vcf_prune.py](https://github.com/MarekSlenker/vcf_prune/blob/main/vcf_prune.py) script. For each dataset, we tested a range of clusters from K = 1–10. Each run consisted of a 100,000-generation burn-in followed by 900,000 generations of data collection, assuming an admixture model and correlated allele frequencies, following [STRUCTURE.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/STRUCTURE.sh). Results were summarized using [CLUMPAK](https://tau.evolseq.net/clumpak/), the optimal number of clusters was determined following the Evanno method, and the homogeneity of results was assessed on the graphs produced by [structureSum](https://github.com/MarekSlenker/structureSum). 




Genome-wide admixture in individuals sampled from the four localities showing species co-occurrence and genome size variation (Table 1) was quantified using the hybrid index (Buerkle, 2005) and interspecific heterozygosity implemented in the Introgress R package (Gompert & Buerkle, 2010). Individuals from single-species populations were used as parental references, with diploid C. apennina set as parent 1 and tetraploid C. amporitana set as 0. These analyses were performed using pseudo-diploidized genotypes and a subset of SNPs showing differentiation between the parental species (GST > 0.95; Hedrick, 2005), as implemented in the vcfR R package (Knaus & Grünwald, 2017).
To determine whether the hybrids identified by the analyses described above share a single origin or arose independently within each of the four sampled localities, we searched for diagnostic alleles (both reference and alternate) of the C. apennina and C. amporitana populations and evaluated their presence in hybrids. An allele was considered diagnostic for a given population if it was observed in at least 30% of non-missing samples within that population and was not detected in any other C. apennina and C. amporitana populations. This set of diagnostic alleles was subsequently analyzed in hybrid populations, where an allele was considered present if it occurred again in at least 30% of non-missing samples. Under a scenario of independent origins, we expected hybrids to show increased frequencies of diagnostic alleles corresponding to their sympatric parental populations.
Furthermore, to assess clonality within the populations, particularly for the hybrids, we used the Nei’s genetic distances calculated in the R package StAMPP, as described above. Following the approach proposed by Tsujimoto et al. (2020), we examined the frequency distributions of pairwise genetic distances within each population (and ploidy level) to determine a genet assignment threshold. Under this framework, the first peak in the frequency distribution represents the distribution of somatic mutations and sequencing errors within genets (clones), while subsequent distances represent variation between distinct genets. To empirically verify the position of these clonal peaks, we calculated pairwise distances for two replicated samples (which were excluded from the primary VCF file). These replicate distances were plotted onto the frequency distribution graphs to confirm that they fell below the inferred genet assignment threshold.
To infer the direction of interspecific crosses, we constructed a maximum-likelihood (ML) tree based on plastome SNP data. The demultiplexed, quality-filtered, and deduplicated ddRADseq reads were mapped to the Cardamine amara plastid reference sequence (KY562580.1) and processed as described above. The resulting VCF file was converted to PHYLIP format using vcf2phylip.py (Ortiz, 2019), and invariant sites were removed with ascbias.py following the recommendations of Leaché et al. (2015). The maximum-likelihood tree was inferred in RAxML-NG v0.9.0 (Kozlov et al., 2019) under the GTR model, incorporating Felsenstein’s ascertainment bias correction. 
