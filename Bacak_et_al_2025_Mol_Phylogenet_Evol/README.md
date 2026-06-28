This repository stores scripts and commands used for the analysis of RADseq data for the paper [Bacak et al. 2025](https://doi.org/10.1016/j.ppees.2025.125921).

### Table of Contents
**[RADseq data processing](#radseq-data-processing)**<br>
&nbsp;&nbsp;&nbsp;&nbsp;[Demultiplexing, quality filtering and deduplication](#demultiplexing-quality-filtering-and-deduplication)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[Variant calling & filtration](#variant-calling--filtration)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[Maximum likelihood (ML) tree](#maximum-likelihood-ml-tree)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[Species tree](#species-tree)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[STRUCTURE](#structure)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[A neighbor-net network](#a-neighbor-net-network)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[Bayes factor species delimitation analysis (BFD*)](#bayes-factor-species-delimitation-analysis-bfd)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[Dsuite](#dsuite)<br>
<br>

**[Hyb-Seq data processing](#hyb-seq-data-processing)**<br>
&nbsp;&nbsp;&nbsp;&nbsp;[HybPiper](#hybpiper)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[HybPhaser](#hybphaser)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[Maximum likelihood (ML) trees](#maximum-likelihood-ml-trees)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[ASTRAL](#astral)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[Read-backed phasing](#read-backed-phasing)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[PhyloSD](#phylosd)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[EPA-ng](#epa-ng)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[AlleleSorting](#allelesorting)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[GRAMPA](#grampa)<br>




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
The ML tree was constructed by RAxML-NG v.0.9.0, employing GTR model with Felsenstein’s ascertainment bias correction using script [MLTree.1.bestTree.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Bacak_et_al_2025_Mol_Phylogenet_Evol/MLTree.1.bestTree.sh). Bootstrap analyses were performed using 500 replicates [MLTree.2.BS_trees.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/MLTree.2.BS_trees.sh), and the final tree with BS support was inferred using 


```ruby
cat *.raxml.bootstraps > allbootstraps.bootstraps

raxml-ng --support \
--tree concat.bialelic.filtered.DP8.passed.vcf.min4.ascbias_Lewis.raxml.bestTree \
--bs-trees allbootstraps.bootstraps \
--prefix concat.bialelic.filtered.DP8.passed.vcf.min4.ascbias_Lewis.raxml --threads 1 
```


## Species tree
The polyploid genotypes were converted into the diploid ones, preserving homozygote and heterozygote genotypes using sed commands.

```ruby
sed -i 's/\.\/\.\/\./.\/./g' $VCF
sed -i 's/1\/1\/1/1\/1/g' $VCF
sed -i 's/0\/1\/1/0\/1/g' $VCF
sed .... 
```


The fasta sequences for each RAD locus (`$REGS`) of each sample were generated using GATK's FastaAlternateReferenceMaker tool, with heterozygous SNP sites outputed using IUPAC ambiguity codes.

```ruby
gatk --java-options "-Xmx120g -XX:+UseSerialGC" SelectVariants \
-V $VCF \
-O "$SAMPLE".vcf.gz \
--sample-name $SAMPLE


gatk --java-options "-Xmx120g -XX:+UseSerialGC" FastaAlternateReferenceMaker \
-R /auto/pruhonice1-ibot/nfs4/home/mslenker/Projects/ref/Erysimum_cheiranthoides/GCA_011420285.1_BTI_Eche1.2_genomic.fasta \
-O "$SAMPLE".fasta \
-L $(basename $REGS) \
-V "$SAMPLE".vcf.gz \
--use-iupac-sample $SAMPLE
```


The FASTA sequences for each RAD locus of each sample were collected using the following commands.

```ruby
# how many regions do we have?
echo $(wc -l $REGS)
# 19356


# Removes line breaks from a fasta files
for file in *.fasta; do
    echo "$file"
    awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $file  > temp_file
    mv temp_file $file
done

# FastaAlternateReferenceMaker uses "." if SNP is missing
sed -i 's/\./-/g' *fasta

parallel --jobs 32 '
  seq={};
  while read SAMPLE; do
      echo ">$SAMPLE" >>res/"$seq".fasta
      grep -A 1 ">$seq " "$SAMPLE".fasta | grep -v ">" >>res/"$seq".fasta
  done <"$SAMPLELIST"
' ::: {1..19356}

```

Phylogenetic trees were inferred in RAxML-NG, using the best-fitting substitution models (see below), from loci that were at least 150 bp in length and contained more than 10 phylogenetically informative sites (calculated by AMAS). A species tree was then estimated with ASTRAL-III, using individual gene trees in which branches with low support (≤20%; based on 500 bootstrap replicates) were collapsed. Additional details are provided below.


## quartet sampling



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



## STRUCTURE
The STRUCTURE analysis was conducted on 100 datasets produced by selecting a single random SNP from each RADseq locus containing at least six SNPs, using [vcf_prune.py](https://github.com/MarekSlenker/vcf_prune/blob/main/vcf_prune.py) script. The STRUCTURE analysis itself was run as in the [STRUCTURE.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/STRUCTURE.sh) script. The results were visualised by [CLUMPAK](https://tau.evolseq.net/clumpak/) and the homogeneity of results was assessed on the graphs produced by [structureSum](https://github.com/MarekSlenker/structureSum).


## A neighbor-net network
A neighbor-net network was created using the NeighborNet algorithm in SplitsTree4 based on Nei’s genetic distances calculated in the StAMPP R package [Neis_distances.R](https://github.com/MarekSlenker/Code-Availability/blob/main/Bacak_et_al_2025_Mol_Phylogenet_Evol/Neis_distances.R).





## Bayes factor species delimitation analysis (BFD*)
VCF file was subsampled, invariant SNPs were removed, and one SNP per RAD locus was selected. The VCF file was converted to phylip using `vcf2phylip.py` and the XML input file for the BEAST was created in BEAUTi using SNAPP template. The BEAST was run from the command line. 

```ruby
beast -threads 16 $INFILE > "$INFILE".log
```
Different delimitation models were compared according to the value of the "marginal L estimate", and the Bayes factor (BF) was computed as $`BF = {2* (MLE1 − MLE0)}`$.  
The SNAPP package was further employed to estimate a coalescent-based species tree directly from SNP data. The XML input file was created in BEAUTi, and the BEAST was run from the command line as above. Topology with the best posterior support was inferred in the TreeAnnotator.


## Dsuite
The [Dsuite](https://github.com/millanek/Dsuite) was employed to test patterns of heterogeneous introgression along the genome using the ABBA–BABA statistics.
The D, f4-ratio, and f-branch statistics were calculated following the Dsuite manual.
```ruby
Dsuite Dtrios -c \
-n Erysimum.odoratum_karpatske.Dsuite \
-t Erysimum.odoratum_karpatske.filtered.DP8.passed.m02.inRegs.ascbias__FELS.bezDiffusumCroaticum.bezOdoratum.astralTree \
Erysimum.odoratum_karpatske.filtered.DP8.passed.m02.inRegs.bezDiffusumCroaticum.bezOdoratum.vcf.gz \
Erysimum.odoratum_karpatske.bezDiffusumCroaticum.bezOdoratum.map

Dsuite Fbranch -p 0.01 \
Erysimum.odoratum_karpatske.filtered.DP8.passed.m02.inRegs.ascbias__FELS.bezDiffusumCroaticum.bezOdoratum.astralTree \
Erysimum.odoratum_karpatske.bezDiffusumCroaticum_Erysimum.odoratum_karpatske.bezDiffusumCroaticum_tree.txt \
> Erysimum.odoratum_karpatske.bezDiffusumCroaticum_Erysimum.odoratum_karpatske.bezDiffusumCroaticum_Fbranch.01.txt

python /home/mint/bin/Dsuite/utils/dtools.py \
Erysimum.odoratum_karpatske.bezDiffusumCroaticum_Erysimum.odoratum_karpatske.bezDiffusumCroaticum_Fbranch.01.txt \
Erysimum.odoratum_karpatske.filtered.DP8.passed.m02.inRegs.ascbias__FELS.bezDiffusumCroaticum.bezOdoratum.astralTree
```

# Hyb-Seq data processing

## HybPiper
The Hyb-Seq reads were processed using [HybPiper v. 2.2.0](https://github.com/mossmatters/HybPiper/releases/tag/v2.2.0)  


```ruby
hybpiper assemble  \
--readfiles "$SAMPLE".trm.R1.fq --readfiles "$SAMPLE".trm.R2.fq \
--targetfile_dna "$(basename "$BAITFILE")" --bwa \
--cpu 4 --prefix "$SAMPLE" --hybpiper_output "$SAMPLE"

hybpiper retrieve_sequences supercontig --targetfile_dna $BAITFILE --sample_names namelist
```

Consensus sequences were aligned using MAFFT v. 7.450, flanking regions and sites with gaps in more than 25% of sequences were removed using the R package ips in R 4.0.0. 
```R
library("ape")
library("ips")

seq = read.dna(file=args[1], format="fasta")
aln = mafft(x=seq, method="auto", maxiterate=100, exec="/software/mafft/7.313/bin/mafft")

aln.trm = deleteEmptyCells(DNAbin=aln)
aln.trm = trimEnds(aln.trm, min.n.seq = nrow(aln.trm)*0.98)
aln.trm = deleteGaps(x=aln.trm, gap.max=round(nrow(aln.trm)/4))
aln.trm = del.colgapsonly(x=aln.trm, threshold=0.1, freq.only=FALSE)
aln.trm = deleteEmptyCells(DNAbin=aln.trm)

write.dna(x=aln.trm, file=args[2], format="fasta", append=FALSE, nbcol=-1, colsep="", colw=80) 
```


## HybPhaser
Next, we used HybPhaser to identify highly variable sequences (indicative of potential paralogs; employing 2x samples only), following [1. SNP assessment](https://github.com/LarsNauheimer/HybPhaser?tab=readme-ov-file#1-snp-assessment), although with some modifications, since we did not use the HybPiper results folder as input data, but only cleaned sequences of 2x samples.  
Sequences, with SNPs coded with iupac ambiguity codes were created by script 
[HybPhaser.1.Consensus_sequence_generation.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/HybPhaser.1.Consensus_sequence_generation.sh), 
table with the proportions of SNPs in each locus and sample was created by [HybPhaser.2a.count_snps.R](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/HybPhaser.2a.count_snps.R), and tables and graphs to assess the variability of sequences were created by the script 
[HybPhaser.2b.assess_dataset.R](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/HybPhaser.2b.assess_dataset.R).

Highly variable sequences, where the proportion of SNPs exceeded 5% were excluded from subsequent processing.


## Maximum likelihood (ML) trees
The best-fitting substitution model was determined for each alignment using the ModelFinder function of `IQ-TREE`. 
Inferred model was parsed to RAxML format using [ML_Trees.PhyhlogenyModelParser.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/ML_Trees.PhyhlogenyModelParser.sh), and the best-scoring ML tree with bootstrap support was inferred by `raxml-ng` (see [ML_Trees.RAxML-NG.IQ_TREE.1seq.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/ML_Trees.RAxML-NG.IQ_TREE.1seq.sh) script).


## ASTRAL
For the species tree reconstruction, internal branches with bootstrap support ≤20% were collapsed using Newick-Utilities v. 1.6. 
```ruby
parallel "nw_ed  {} 'i & b<20' o > ./{}.BS20" ::: *support
```
The species tree was constructed employing a multispecies coalescent model implemented in ASTRAL-III, including the computation of local posterior probabilities to evaluate branch support.
```ruby
cat *.BS20 > bs20_trees
java -jar ~/bin/astral.5.7.8/astral.5.7.8.jar -i bs20_trees -o astralTree.tree --namemapfile namemapfile -t 4 -r 10000
```



## Read-backed phasing
The code (roughly) follows the procedures outlined in the [alleles_workflow](https://github.com/mossmatters/phyloscripts/tree/master/alleles_workflow) GitHub repository. The following scripts implement these procedures.  
The sequences of each sample are phased by [Phasing.1.phasing_oneSample.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/Phasing.1.phasing_oneSample.sh). This script takes consensus sequences and fastq reads at the input and produces phased sequences `"$SAMPLE".v1.phased.fasta, "$SAMPLE".v2.phased.fasta, ...` and unphased `"$SAMPLE".unPhased.fasta` sequences. You need to run this script for all samples.  
Phased sequences are sometimes represented by multiple mutually unphased blocks (take a look at `"$SAMPLE".whatshap.gtf`). Selection of the longest phased block and masking of the remaining variant is the responsibility of [Phasing.2.masking.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/Phasing.2.masking.sh) script. Phased sequences are written to RESDIR directory. This script works with all samples simultaneously, using files produced by the previous script.   

The SAMPLEPLOIDYLIST file needed for <ins>Phasing.2.masking.sh</ins> has the following structure (sample1 \n ploidy of sample1 \n sample2 \n ploidy of sample2 \n ....)  
```
acraBAB6
2
acraC003_104
2
acraC018_101
4
acraC095_109
3
acraC149_8
2
```

## PhyloSD
Due to the pipeline's requirement for a single representative of diploid genomes, the species tree from each gene tree was calculated by ASTRAL-III following [PhyloSD.1.ASTRALGeneTrees.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Bacak_et_al_2025_Mol_Phylogenet_Evol/PhyloSD.1.ASTRALGeneTrees.sh) script, merging the respective 2x samples to a single taxon, but keeping polyploids with separate (phased) alleles.  

#### 1) NEAREST DIPLOID SPECIES NODE algorithm
1.5) Root and sort nodes in trees, losing bootstrap and aLRT supports on the way. (see also [Root and sort nodes in trees...](https://github.com/eead-csic-compbio/allopolyploids?tab=readme-ov-file#15-root-and-sort-nodes-in-trees-loosing-bootstrap-and-alrt-supports-on-the-way))
```ruby
mkdir 1.5_RootAndSort
parallel -j 8 "echo {}; perl5.38.2 <PATH>/PhyloSD/bin/PhyloSD/_reroot_tree.pl {} > 1.5_RootAndSort/{.}.root.ph" ::: *astralTree
```
1.6) Check diploid skeleton (topology) for each tree (see also [Check diploid skeleton (topology) for each tree)](https://github.com/eead-csic-compbio/allopolyploids?tab=readme-ov-file#16-check-diploid-skeleton-topology-for-each-tree))
```ruby
cd 1.5_RootAndSort
for FILE in *root.ph; do
   perl5.38.2 <PATH>/PhyloSD/bin/PhyloSD/_check_diploids.pl $FILE;
done > ../diploids.log
```
Due to an unacceptable loss of data, incongruent diploid skeletons were not discarded (unlike in the original pipeline), and all sequences and trees were moved to `1.7_congruent_and_labelled_files` folder.
```ruby
mkdir 1.7_congruent_and_labelled_files
cp 1.5_RootAndSort/*ph 1.7_congruent_and_labelled_files # trees
cp inputSequences/*fna 1.7_congruent_and_labelled_files # sequences
```

1.8) Labelling polyploid homeologs (see also [Labelling polyploid homeologs](https://github.com/eead-csic-compbio/allopolyploids?tab=readme-ov-file#18-labelling-polyploid-homeologs))
```ruby
cd 1.7_congruent_and_labelled_files
parallel -j 8 "echo {}; perl5.38.2 <PATH>/PhyloSD/bin/PhyloSD/_check_lineages_polyploids.pl -v -f {} -t {.}.raxml.bestTree.root.ph > {}.log" ::: *.fna
```
Homeologs of polyploids, that fit the criteria in [polyconfig.Erysimum.pm](https://github.com/MarekSlenker/Code-Availability/blob/main/Bacak_et_al_2025_Mol_Phylogenet_Evol/polyconfig.Erysimum.pm) config file (defined according to the phylogenetic tree), and thus can be attributed to one of the diploid parents, were written to `label.reduced.fna` files.

#### 2) BOOTSTRAPPING REFINEMENT algorithm

2.1) Each polyploid homeolog from `label.reduced.fna` files was merged with 2x samples+outgroup (supercontigs_aln_gt31_skontrolovane_consens2Fazovane_aln_bezOutoci_cons2x), and the correctness of attribution to 2x parental taxa was tested by bootstrapping (see also [Set the pruned FASTA alignments (diploids + outgroups + one polyploid homeolog)](https://github.com/eead-csic-compbio/allopolyploids?tab=readme-ov-file#21-set-the-pruned-fasta-alignments-diploids--outgroups--one-polyploid-homeolog)).  

```ruby

for f in *label.reduced.fna; do
  bn=${f%.label.reduced.fna}
  for r in $(grep ">" $f); do
    echo $r
    grep -A 1 "$r" "$f" > ../2.1.one_allopolyploid_plus_diploids/"$bn"."$r".fna
    cat ../supercontigs_aln_gt31_skontrolovane_consens2Fazovane_aln_bezOutoci_cons2x/"$bn".fasta >> ../2.1.one_allopolyploid_plus_diploids/"$bn"."$r".fna
  done
done
```

2.2) Run 500 non-parametric bootstrapping replicates & Labelling polyploid homeologs. [PhyloSD.2.LabelBSTrees.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Bacak_et_al_2025_Mol_Phylogenet_Evol/PhyloSD.2.LabelBSTrees.sh). The results are in `counts` files. Those files summarize the results of re-labelling polyploid homeologs. We required confirmation by at least 20% of bootstrap replicates. That means if some homeolog was originally labelled as "witmannii" (step 1.8), we keep that particular sequence only if more than 100 BS trees (20%) were re-labelled as "witmannii".  

2.18) "Homeologs' ML consensus tree" (see also [2.18) Phylogenomic analysis of concatenated labelled, filtered and corrected genes/MSAs](https://github.com/eead-csic-compbio/allopolyploids?tab=readme-ov-file#218-phylogenomic-analysis-of-concatenated-labelled-filtered--and-corrected-genesmsas-homeologs-ml-consensus-tree)). We concatenated sequences of 2x samples and the labelled homeologs of each polyploid (those with at least 12-15% representation in the polyploid genome). If more than one homeolog of the gene was labelled with the same 2x label, a homeolog with higher BS support was chosen. The phylogenetic tree was computed in RAxML-NG from the concatenated alignment, as described above.


#### 3) SUBGENOME ASSIGNMENT algorithm
In this step, homeologs referring to the same subgenome are amalgamated, based on PCoA-MST (Principal coordinates analysis-minimum spanning tree). Using the following code, we computed patristic distances and PCoA-MST in R.

```R
library(adephylo)
library(ape)
library(stats)

tree = read.tree("Ery.Karpatske.fazovane.presli.BS.concatenated.nad12percent.raxml.bestTree")

patristicDists = distTips(tree, method = "patristic")

distMatrix= as.matrix(patristicDists)

poylploids = distMatrix[grep(
"odoratumRetezat_witmannii|odoratumRetezat_saxosum|odoratum22chrom_saxosum|odoratum22chrom_witmannii|odoratum22chrom_crassistylum|odoratum6x_saxosum|odoratum6x_witmannii|odoratum6x_crassistylum|odoratum6x_cuspidatum|odoratum4x_witmannii|odoratum4x_cuspidatum",
colnames(distMatrix)),
grep(
"odoratumRetezat_witmannii|odoratumRetezat_saxosum|odoratum22chrom_saxosum|odoratum22chrom_witmannii|odoratum22chrom_crassistylum|odoratum6x_saxosum|odoratum6x_witmannii|odoratum6x_crassistylum|odoratum6x_cuspidatum|odoratum4x_witmannii|odoratum4x_cuspidatum",
colnames(distMatrix))]

maf.coa <- dudi.pco(as.dist(poylploids), scannf = FALSE, nf = 3)

maf.mst <- ade4::mstree(dist.dudi(maf.coa), 1)

s.label(maf.coa$li, label = row.names(maf.coa$li),clab = 0.8,  cpoi = 2, neig = maf.mst, cnei = 1)

plot(maf.coa$li[,1],maf.coa$li[,2], asp=1)
```

3.5) Amalgamate homeologs: however, the homeologs were too differentiated to be merged. Therefore, we used ASTRAL-III to infer the final subgenomic tree.


## EPA-ng
The EPA-ng performs maximum likelihood-based placement of allelic sequences onto a reference phylogenetic tree. Only placements with a likelihood weight ratio greater than 0.9 were considered significant and subsequently used for species tree inference in ASTRAL-III. Since the samples were not evenly represented, we used usearch to cluster the sequences for each taxon, and the centroid sequence was used as a representative sequence for the entire taxon.


```ruby

PP="Kavlar2" # polypolod sample
# ALN is a fasta file, containing sequences of 2x samples and one polyploid ($PP)

for ALN in *fasta; do
    echo $ALN
    grep -A 1 "$PP" $ALN > querySeq.fasta 
    sed -i "/$PP/,+1 d" $ALN

    # USEARCH
    for SAMP in samples*; do
        SMPL=${SAMP#*.}
        grep --no-group-separator -A 1 -f $SAMP $ALN > "$SAMP".toClust

        sed -i 's/a/A/g' "$SAMP".toClust
        sed -i 's/c/C/g' "$SAMP".toClust
        sed -i 's/t/T/g' "$SAMP".toClust
        sed -i 's/g/G/g' "$SAMP".toClust

        <PATH>/usearch11.0.667_i86linux32 -cluster_fast "$SAMP".toClust \
        -id 0.1 \
        -centroids "$ALN"."$SAMP".centroids.fasta \
        -consout "$ALN"."$SAMP".consout.fasta \
        -uc "$ALN"."$SAMP".clusters.uc 

        echo ">$SMPL" >> ${ALN%.*}.centroids.fasta
        grep -v ">" "$ALN"."$SAMP".centroids.fasta >> ${ALN%.*}.centroids.fasta

        rm "$SAMP".toClust
    done

    cat querySeq.fasta >>${ALN%.*}.centroids.fasta

    mafft ${ALN%.*}.centroids.fasta > ${ALN%.*}.centroids.aln.fasta

    awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${ALN%.*}.centroids.aln.fasta  > temp_file
    mv temp_file ${ALN%.*}.centroids.aln.fasta


    grep -A 1 "$PP" ${ALN%.*}.centroids.aln.fasta > querySeq.fasta 
    sed -i "/$PP/,+1 d" ${ALN%.*}.centroids.aln.fasta



    ./iqtree -s ${ALN%.*}.centroids.aln.fasta -m TESTONLY -st DNA -nt $PBS_NUM_PPN -pre iqtreeoutput -rcluster 10 -redo -quiet

    # find the best model
    MODEL=$(grep "Best-fit model according to BIC:" iqtreeoutput.iqtree | cut -d ' ' -f 6)

    # https://raw.githubusercontent.com/MarekSlenker/Code-Availability/refs/heads/main/Slenker_et_al_2024_Molecular_Ecology/ML_Trees.PhyhlogenyModelParser.sh
    NEWMODEL=$( ./ML_Trees.PhyhlogenyModelParser.sh $MODEL )
    ./raxml-ng --msa ${ALN%.*}.centroids.aln.fasta --model "$NEWMODEL" --tree pars{5},rand{5} --blopt nr_safe --redo  --threads $PBS_NUM_PPN --force --prefix ${ALN%.*}.centroids >> ${ALN%.*}.centroids.log || ./raxml-ng --msa ${ALN%.*}.centroids.aln.fasta --model "$NEWMODEL" --tree pars{5},rand{5} --blopt nr_safe --threads 1 --prefix ${ALN%.*}.centroids >> ${ALN%.*}.centroids.log

    epa-ng --ref-msa ${ALN%.*}.centroids.aln.fasta --tree  ${ALN%.*}.centroids.raxml.bestTree --query querySeq.fasta --model ${ALN%.*}.centroids.raxml.bestModel --redo
    mv epa_result.jplace ${ALN%.*}.centroids.RAxMLTree.jplace

done
```

Results were evaluated as follow:  

**1)** jplace files were edited:
```ruby
sed -i 's/:[0-9.]*{/{/g' *jplace
```

**2)** The best placement for each edge_num was associated with the respective taxon or left unassociated, if it points to an internal node (internal nodes couldn't be compared among trees, as input trees had different topology). R code follows.
```ruby

for (F in c("BREITs1","GLAVs1", "Maliscak9", "SCO3", "CIUC6R", "Kavlar2", "RET1s", "VGrad1")) {
  wd=paste("<PATH>/epa/vyhodnotenie/",F, sep = "")
  setwd(wd)
  
  
      for (jplace_path in list.files(".", "jplace")) {
        jplace <- fromJSON(jplace_path, simplifyVector = FALSE)
        
        # Extract the tree string and read it as phylo object
        tree_newick <- jplace$tree
        tree <- read.tree(text = tree_newick)
        
        # tree$tip.label
        tipLabelTable=matrix(unlist(strsplit(tree$tip.label, "\\{|\\}")), ncol=2, byrow=T)
        
        placements <- jplace$placements
        
        # write.table(paste(jplace_path), 
        #            file = paste("../",F,".vyhodnotene.txt", sep = ""), append = T, quote = F, row.names = F, col.names = F)
        
        for (i in 1:length(placements)) {
          p <- placements[[i]]
          p$n[[1]]
          edge_num=p[1]$p[[1]][[1]]
          like_weight_ratio=p[1]$p[[1]][[3]]
          write.table(paste(jplace_path, ":",edge_num, ":",like_weight_ratio, ":",tipLabelTable[which(tipLabelTable[,2] == edge_num), 1],":",p$n[[1]], sep = ""), 
                      file = paste("../",F,".vyhodnotene.txt", sep = ""), append = T, quote = FALSE, row.names = FALSE, col.names = FALSE)
          
        }
        
      }

}
```

**3)** The final step was to count allelic sequences placed with confidence above 95%.

```ruby
setwd("<PATH>/epa/vyhodnotenie/")



for (F in c("BREITs1","GLAVs1", "Maliscak9", "SCO3", "CIUC6R", "Kavlar2", "RET1s", "VGrad1")) {
  dd = read.delim(paste(F, ".vyhodnotene.txt", sep = ""), sep = ":", header = FALSE)
  dd = dd[- which(dd$V4 == ""),]
  aa = table(dd$V4[which(dd$V3 > 0.95)])
  
  
  write.table(t(aa), file=paste(F,".95.res", sep = ""), col.names = T, row.names = FALSE)
  write.table(t(aa/sum(aa)), file=paste(F,".95.res", sep = ""), col.names = FALSE, row.names = FALSE, append = T)
  
  write.table(dd[which(dd$V3 > 0.95),], file=paste(F,".above095.res", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE)
   
}

```



## AlleleSorting
In the [AlleleSorting](https://github.com/MarekSlenker/AlleleSorting) approach (applicable to the tetraploids only), alleles were sorted into two homeologs based on sequence divergence, following the pipeline proposed on that repository.  


## GRAMPA
The GRAMPA uses an algorithm for counting gene duplications and losses to identify polyploidy events, distinguishing between allo- and autopolyploid, and place them on a phylogeny. 

```ruby
SPECIESTREE="BREITs1_trees.Root.astralTree"
GENETREE="BREITs1_trees.Root.nwk"
H1="BREITs1"

grampa \
-s $SPECIESTREE \
-g $GENETREE \
-h1 $H1 \
-p 8 \
--overwrite

```





