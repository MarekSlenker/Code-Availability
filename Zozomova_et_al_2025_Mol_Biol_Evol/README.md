This repository stores scripts and commands used for the analysis of RADseq data for the paper Zozomova et al. 2025

## Demultiplexing, quality filtering and deduplication
Raw reads were demultiplexed using [rad_1_demultiplexing_1_qsub.sh](https://github.com/V-Z/RAD-Seq-scripts/blob/master/bin/rad_1_demultiplexing_1_qsub.sh) and [rad_1_demultiplexing_2_run.sh](https://github.com/V-Z/RAD-Seq-scripts/blob/master/bin/rad_1_demultiplexing_2_run.sh) scripts.  

next, [fastq pair](https://github.com/linsalrob/fastq-pair) was used to rewrite paired-end fastq files to ensure all reads have a mate and separate singletons.
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
Genome of C. amara ([GCA_040955855.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_040955855.1/)) was indexed, 
```ruby
REF="GCA_040955855.1_C_amara_ONT_v2_genomic.fna"
samtools faidx "$REF"
java -jar $PICARD CreateSequenceDictionary R="$REF"
bwa index "$REF"
```

and reads were mapped to reference using BWA
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

All 944 VCFs, i.e. genotyped contigs, were concatenated using bcftools.
```ruby
bcftools concat -O z *.vcf.gz > concat.vcf.gz
tabix concat.vcf.gz
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
--exclude-filtered  --max-nocall-fraction 0.2 -L regions1000apart.bed
```

The last `SelectVariants` command filtered out any SNPs outside identified RAD loci (regions1000apart.bed). Loci were identified following [identifiRadLoci.workflow](https://github.com/MarekSlenker/vcf_prune/blob/main/identifiRadLoci.workflow).  

To determine missing data per sample, VCF file was converted to phylip format using [vcf2phylip.py](https://github.com/edgardomortiz/vcf2phylip/blob/master/vcf2phylip.py) script and the percentage of missing data was calculated as the number of `N` divided by the total number of SNPs. The number of Ns was calculated by the following command:
```ruby
function countchar()
{
    while IFS= read -r i; do printf "%s" "$i" | tr -dc "$1" | wc -m; done
}
countchar 'N' <concat.bialelic.filtered.DP8.passed.vcf.min4.phy
```




## STRUCTURE
The STRUCTURE analysis was conducted on 100 datasets produced by selecting a single random SNP from each RADseq locus containing at least six SNPs, using [vcf_prune.py](https://github.com/MarekSlenker/vcf_prune/blob/main/vcf_prune.py) script. The STRUCTURE analysis itself was run like in the [STRUCTURE.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/STRUCTURE.sh) script. The results were visualised by [CLUMPAK](https://tau.evolseq.net/clumpak/) and the homogeneity of results was assessed on the graphs produced by [structureSum](https://github.com/MarekSlenker/structureSum).



## A neighbor-net network
A neighbor-net network was created using the NeighborNet algorithm in SplitsTree4 based on Nei’s genetic distances calculated in the StAMPP R package [Neis_distances.R](https://github.com/MarekSlenker/Code-Availability/blob/main/Zozomova_et_al_2025/Neis_distances.R).



## PCA
The principal component analysis (PCA) based on covariance matrix, as implemented in the R package adegenet was calculated using [PCA.R](https://github.com/MarekSlenker/Code-Availability/blob/main/Zozomova_et_al_2025/PCA.R) script.




## Maximum likelihood (ML) tree
The [vcf2phylip.py](https://github.com/edgardomortiz/vcf2phylip/blob/master/vcf2phylip.py) script was used to transform the data from the VCF file to the PHYLIP, and invariant sites were removed with the script [ascbias.py](https://github.com/btmartin721/raxml_ascbias).
The ML tree was constructed by RAxML-NG v.0.9.0, employing GTR model with Lewis’s ascertainment bias correction using script [MLTree.1.bestTree.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Zozomova_et_al_2025/MLTree.1.bestTree.sh). Bootstrap analyses were performed using 500 replicates [MLTree.2.BS_trees.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/MLTree.2.BS_trees.sh), and the final tree with BS support was inferred using 
```ruby
cat *.raxml.bootstraps > allbootstraps.bootstraps

raxml-ng --support \
--tree concat.bialelic.filtered.DP8.passed.vcf.min4.ascbias_Lewis.raxml.bestTree \
--bs-trees allbootstraps.bootstraps \
--prefix concat.bialelic.filtered.DP8.passed.vcf.min4.ascbias_Lewis.raxml --threads 1 
```


## PolyRelatedness
The relatedness coefficients were estimated in PolyRelatedness. VCF file was formatted to PolyRelatedness input file, and coefficients were calculated [PolyRelatedness.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/PolyRelatedness.sh). Heatmap and violin plots were made in R [PolyRelatedness.plots.R](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/PolyRelatedness.plots.R). 

## NewHybrids
This analysis was performed on subsets of 100 loci, with each subset differentiating parental groups, namely matthioli-rivularis, amara-rivularis and acris-rivularis. We will demonstrate it using matthioli-rivularis.  

We used `bcftools` to select samples from pure populations of matthioli and rivularis. The resulting vcf file was processed in R.
```ruby
library(adegenet)
library(vcfR)

vcf <- read.vcfR("BGhybridy.bialelic.filtered.DP8.passed.m02.inRegs.bezOutgroups.noInvs.MatRiv.vcf.gz")

# Converting VCF data to a genlight object

genlight <- vcfR2genlight.tri.MK(vcf)  # you can find this function in the script Neis_distances.R
locNames(genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")   # add real SNP.names

pop(aa.genlight) = c("riv","riv","riv", ... "mat","mat","mat")
pops <- as.factor(c("riv","riv","riv", ... "mat","mat","mat"))

diffs <- genetic_diff(vcf, pops = pops, method = 'nei')

write.table(diffs, file = "BGhybridy.bialelic.filtered.DP8.passed.m02.inRegs.bezOutgroups.noInvs.MatRiv.diff", quote = F, row.names = F)
```

Next, we selected only SNPs with Gprimest == 1 (see BGhybridy.bialelic.filtered.DP8.passed.m02.inRegs.bezOutgroups.noInvs.MatRiv.diff file). To minimize the effects of linkage disequilibrium, one SNP per scaffold was selected, prioritizing those with the least missing data. The dataset was reduced to 100 SNPs by random selection. Next, we created several data sets, each containing samples from pure parental populations combined with samples from a single hybrid zone, and the subset of just identified SNPs. The SNP subsets were converted to the NewHybrids input format using [vcf_to_newhybrids_format.py](https://github.com/mscharmann/tools/blob/master/vcf_to_newhybrids_format.py). The `NewHybrids` was run using the following command.

```ruby
newhybrids-no-gui-linux.exe -d BGhybridy.....newhybrids.txt --burn-in 1000000 --num-sweeps 9000000 --no-gui > std_err_output.txt 2>&1
```



## Hybrid index 

The hybrid index was calculated in `GenoDive v. 3.06`, following GenoDive manual. Another estimation of hybrid index and interspecific heterozygosity was done using `Introgress` R package.  

  
Unlike GenoDive (which accepts the full data matrix and polyploids), Introgress was run on a subset of SNPs that differed between the parental populations (those identified for `NewHybrids`, but not reduced to 100 SNPs). The triploid genotypes were converted into the diploid ones (preserving homozygote and heterozygote genotypes; i.e. 0/0/0 -> 0/0, 0/0/1 -> 0/1, 0/1/1 -> 0/1, 1/1/1 -> 1/1). Indices were calculated using [introgress.R](https://github.com/MarekSlenker/Code-Availability/blob/main/Zozomova_et_al_2025/introgress.R) script.


























