This repository stores scripts and commands used for the analysis of RADseq data for the paper Zozomova et al. 2025



# RADseq data processing

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
Genome of Erysimum cheiranthoides (GenBank: [GCA_011420285.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_011420285.1/)) was indexed, 

```ruby
REF="GCA_011420285.1_BTI_Eche1.2_genomic.fasta"
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

All 223 VCFs, i.e. genotyped contigs, were concatenated.
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




## Maximum likelihood (ML) tree
The [vcf2phylip.py](https://github.com/edgardomortiz/vcf2phylip/blob/master/vcf2phylip.py) script was used to transform the data from the VCF file to the PHYLIP, and invariant sites were removed with the script [ascbias.py](https://github.com/btmartin721/raxml_ascbias).
The ML tree was constructed by RAxML-NG v.0.9.0, employing GTR model with Felsenstein’s ascertainment bias correction using script [MLTree.1.bestTree.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Bacak_et_al_2025/MLTree.1.bestTree.sh). Bootstrap analyses were performed using 500 replicates [MLTree.2.BS_trees.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/MLTree.2.BS_trees.sh), and the final tree with BS support was inferred using 


```ruby
cat *.raxml.bootstraps > allbootstraps.bootstraps

raxml-ng --support \
--tree concat.bialelic.filtered.DP8.passed.vcf.min4.ascbias_Lewis.raxml.bestTree \
--bs-trees allbootstraps.bootstraps \
--prefix concat.bialelic.filtered.DP8.passed.vcf.min4.ascbias_Lewis.raxml --threads 1 
```

### quartet sampling

Branch support and the amount of discordance were further assessed by the quartet sampling method (Pease et al., 2018), enabling to distinguish between lack of support and conflicting support in the phylogenetic tree.
Branch support and the degree of discordance were further assessed using the quartet sampling method (Pease et al., 2018), which enables distinction between a lack of support and conflicting support within the phylogenetic tree.

```ruby

python3 <PATH>/quartetsampling/bin/pysrc/quartet_sampling.py \
--tree Erysimum.odoratum_karpatske.filtered.DP8.passed.m02.inRegs.min4.ascbias__FELS.raxml.bestTree \
--align Erysimum.odoratum_karpatske.filtered.DP8.passed.m02.inRegs.min4.ascbias.phy.fasta.phylip \
--reps 1000 \
--threads 16 --lnlike 2
```
The results were visualised using [quartetsampling.r](https://github.com/tomas-fer/scripts/blob/master/quartetsampling.r) script from Tomas Fer's repo.




Further insight into the overall genetic structure was obtained using a Bayesian clustering approach implemented in STRUCTURE v. 2.3.4 (Pritchard et al. 2000) and a neighbor-net network in SplitsTree4 (Huson and Bryant, 2006). STRUCTURE analysed 100 datasets, each containing a single randomly selected SNP from RADseq loci with at least six SNPs, using the vcf_prune.py script (Šlenker, 2024). Calculations were performed and summarized as described in Šlenker et al. (2021). For the neighbor-net (NN) analysis, Nei’s genetic distances (Nei 1972) were calculated in the R package StAMPP (Pembleton et al. 2013) using R 4.4.0 (R Core Team 2024).  
Furthermore, a Bayes factor species delimitation analysis (BFD*, Leaché et al., 2014; Leaché and Bouckaert, 2018) was performed to statistically validate the genetic clusters within the Carpathian diploids (corresponding to E. witmannii s.l. and E. vagicum). Marginal likelihoods of species trees were derived using the Path Sampling approach with SNAPP v.1.4.2 (Bryant et al., 2012) and BEAST v. 2.5.0 (Bouckaert et al., 2014). The dataset of unlinked SNPs was used, which was reduced to three samples per genetic cluster, except for the smallest cluster Vihorlat, where only 2 samples were selected. Analyses were run in eight steps for each model, with 1,000,000 MCMC iterations, sampling every 1,000th, and a burn-in cutoff of 10%. Competing species delimitation models were ranked by comparing their marginal likelihood estimates and their support was assessed by calculating the Bayes factor (Kass and Raftery, 1995), as suggested by Leaché and Bouckaert (2018). Seven alternative species models were explored, either keeping the Carpathian diploids as one unit or splitting them into two to four entities, taking into account the ML tree, NN and STRUCTURE clustering results. Moreover, TreeAnnotator (Drummond and Rambaut, 2007) was used to summarise the posterior distribution of species trees and to identify the topology with the best posterior support, using the species model with the highest support in BFD* described above.
Polyploid origins and potential reticulation and introgression events in diploids were also examined using SNaQ (Solís-Lemus and Ané, 2016; Solís-Lemus et al., 2017) and Dsuite (Malinsky et al., 2021). For the SNaQ analysis, concordance factors were calculated from unlinked SNPs using the R function SNPs2CF (Olave and Meyer, 2020), and a starting tree was inferred using the Quartet MaxCut algorithm (Snir and Rao, 2012). To reduce computational demands,  the Carpathian diploids (E. witmannii s.l. and E. vagicum) were kept as a single entity. To test complex patterns of heterogeneous introgression along the genome using the ABBA–BABA and related statistics (Durand et al., 2011), Dsuite (Malinsky et al., 2021) was used, which calculates the D, f4-ratio, and f-branch statistics. For this analysis, Carpathian diploids were sorted into four entities following the results of BFD*.
The RADseq reads were also utilized to obtain data from plastomes. The reads were mapped to the plastome of E. cheiranthoides (GenBank accession number MN207123.1), processed, and the ML tree was constructed in RAxML-NG as described above.

Hyb-Seq data processing
The Hyb-Seq reads were processed using HybPiper v. 2.2.0 (Johnson et al., 2016) to extract consensus sequences of targeted exons. Highly variable sequences (indicative of potential paralogs), where the proportion of SNPs exceeded 5% were excluded from subsequent processing (identified using HybPhaser; https://github.com/LarsNauheimer/HybPhaser), resulting in 964 sequences (exons/supercontigs). Consensus sequences were aligned using MAFFT v. 7.450 (Katoh and Standley, 2013), and flanking regions and sites with gaps in more than 25% of sequences were removed using the R package ips (Heibl, 2008 onward) in R 4.4.0 (R Core Team, 2024). ML trees were inferred in RAxML-NG using the best-fitting substitution models as determined by the ModelFinder function of IQ-TREE v.1.6.12 (Chernomor et al., 2016; Kalyaanamoorthy et al., 2017) based on the Bayesian information criterion. Bootstrap analyses were performed using 500 replicates. For the species tree reconstruction, internal branches with bootstrap support ≤20% were collapsed using Newick-Utilities v. 1.6 (Junier and Zdobnov, 2010). The species tree was constructed employing a multispecies coalescent model implemented in ASTRAL-III (Zhang et al., 2018), including computation of local posterior probabilities to evaluate branch support (Sayyari and Mirarab, 2016).
The supercontig sequences of polyploid accessions of E. odoratum were further processed for read-backed phasing to infer allele sequences, as described in detail in Šlenker et al. (2021). We applied four different methods to identify homeologous diploid subgenomes and the most likely parental species or lineages: PhyloSD (Sancho et al., 2022), EPA-ng (Barbera et al., 2019), AlleleSorting (Šlenker et al., 2021), and GRAMPA (Thomas et al., 2017). In the PhyloSD approach we followed the pipeline by Sancho et al. (2022) with some modifications. Due to the pipeline's requirement for a single representative of diploid genomes, we calculated species tree of each gene using ASTRAL- III. Due to unacceptable loss of data, we did not discard incongruent diploid skeletons (unlike in Sancho et al. 2022), but rather applied stricter criteria in the Bootstrapping Refinement step,  keeping only the homeologs that were confirmed by at least 20% of bootstrap replicates. Only the major homeolog-types (those with at least 12-15% representation in the polyploid genome) were further processed with the “Subgenome Assignment” algorithm and used for the subgenomic tree constructions in RAxML-NG. The subgenomic ML trees were finally summarized in ASTRAL-III. In the AlleleSorting approach (applicable to the tetraploids only), alleles were sorted into two homeologs based on sequence divergence, labelled to attribute them to different subgenomes (Šlenker et al. 2021), and treated as independent accessions in the coalescent based species tree inference in ASTRAL-III. EPA-ng (Barbera et al., 2019), a reimplementation of the evolutionary placement algorithm (EPA), performs maximum likelihood-based placement of allelic sequences onto a reference phylogenetic tree. Only placements with a likelihood weight ratio greater than 0.9 were considered significant and subsequently used for species tree inference in ASTRAL-III. Finally, GRAMPA (Gene-tree Reconciliation Algorithm with MUL-trees for Polyploid Analysis) uses an algorithm for counting gene duplications and losses to identify polyploidy events, distinguishing between allo- and autopolyploid, and place them on a phylogeny (Thomas et al., 2017). 
