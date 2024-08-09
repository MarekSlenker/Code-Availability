<table><tr><td>https://github.com/lifeparticle/Markdown-Cheatsheet</td></tr></table>

### Table of Contents
**[Hyb-Seq data processing](#hyb-seq-data-processing)**<br>
&nbsp;&nbsp;&nbsp;&nbsp;[HybPiper](#hybpiper)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[HybPhaser](#hybphaser)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[Read-backed phasing](#read-backed-phasing)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[Maximum likelihood (ML) trees](#maximum-likelihood-ml-trees)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[ASTRAL](#astral)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[STRUCTURE](#structure)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[SNaQ (PhyloNetworks)](#snaq-phylonetworks)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[PhyloSD](#phylosd)<be>

**[RADseq data processing](#radseq-data-processing)**<br>
&nbsp;&nbsp;&nbsp;&nbsp;[Variant calling & filtration](#variant-calling--filtration)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[Maximum likelihood (ML) tree](#maximum-likelihood-ml-tree)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[Bayes factor species delimitation analysis (BFD*)](#bayes-factor-species-delimitation-analysis-bfd)<be>

&nbsp;&nbsp;&nbsp;&nbsp;[STRUCTURE](#structure)<be>
&nbsp;&nbsp;&nbsp;&nbsp;[STRUCTURE](#structure)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[STRUCTURE](#structure)<be>
&nbsp;&nbsp;&nbsp;&nbsp;[STRUCTURE](#structure)<br>
&nbsp;&nbsp;&nbsp;&nbsp;[STRUCTURE](#structure)<be>
&nbsp;&nbsp;&nbsp;&nbsp;[STRUCTURE](#structure)<br>


# Hyb-Seq data processing
## HybPiper
The Hyb-Seq reads were processed using [HybPiper v. 1.3](https://github.com/mossmatters/HybPiper/releases/tag/v1.3.1_final)  
```ruby
python2 HybPiper/reads_first.py --bwa -r "$SAMPLE".R{1,2}.fastq -b "$BAITFILE" --cpu "$NCPU" --prefix "$SAMPLE"
python2 HybPiper/intronerate.py --prefix "$SAMPLE" --addN 
python2 HybPiper/cleanup.py "$SAMPLE" 
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

Sequences, where the proportion of single nucleotide polymorphisms exceeded 2% were excluded from subsequent processing.

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

## Maximum likelihood (ML) trees
Maximum likelihood (ML) trees were inferred:  
* for each gene separately using RAxML-NG. Best-fitting substitution model was determined through the ModelFinder function of IQ-TREE. [ML_Trees.RAxML-NG.IQ_TREE.1seq.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/ML_Trees.RAxML-NG.IQ_TREE.1seq.sh)
* from concatenated genes (concatenated by [AMAS](https://github.com/marekborowiec/AMAS)).
  * Each gene was treated as a separate partition with its best-fitting substitution model determined through the ModelFinder function of IQ-TREE based on the BIC. [ML_Trees.concat.1.bestFittingModels.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/ML_Trees.concat.1.bestFittingModels.sh). <em>(partitions have to start with 'DNA, ')</em>
  * Inferred models have to be parsed to RAxML format using [ML_Trees.PhyhlogenyModelParser.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/ML_Trees.PhyhlogenyModelParser.sh), and overall structure as well. The resulting file (partitions.best_scheme) should look like this:  
    ```ruby
    GTR+F+I+G4, p1_Assembly_000000000043_Assembly_000000004815_Assembly_000000007919 = 1-560, 101083-102022, 117256-117802;
    GTR+F+I+G4, p2_Assembly_000000000043__Assembly_000000002965 = 561-907, 58272-59025;
    TVM+F+I+G4, p3_Assembly_000000000055 = 908-1636;
    ```    

  * The best-scoring ML tree was inferred by [ML_Trees.concat.2.bestTree.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/ML_Trees.concat.2.bestTree.sh)
  * Bootstrap analyses were performed using 500 replicates. [ML_Trees.concat.3.BS_trees.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/ML_Trees.concat.3.BS_trees.sh)
  * The final ML tree with support values was inferred based on the BS replicates and best-scoring ML tree.
    
    ```ruby
    # merge all partial BS
    cat *.raxml.bootstraps > allbootstraps.bootstraps
    
    raxml-ng --support --tree "${ALIGNMENT%.*}".raxml.bestTree --bs-trees allbootstraps.bootstraps --threads 1 --prefix "${ALIGNMENT%.*}"  >> "${ALIGNMENT%.*}".log
    ```

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

## STRUCTURE
Homogeneous genetic clusters were inferred using Bayesian clustering algorithm implemented in STRUCTURE v. 2.3.4. The [snipStrup pipeline](https://github.com/MarekSlenker/snipStrup) was used to map the Hyb-Seq reads to the consensus sequences, call variants and convert the VCF file. The STRUCTURE analysis itself was run like in the [STRUCTURE.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/STRUCTURE.sh) script. The results were visualised by [CLUMPAK](https://tau.evolseq.net/clumpak/) and the homogeneity of results was assessed on the graphs produced by [structureSum](https://github.com/MarekSlenker/structureSum).


## SNaQ (PhyloNetworks)
Reticulation events were inferred with 500 gene trees, computed from the longest gene alignments. The gene trees were summarized using quartet concordance factors. The species tree reconstructed in ASTRAL-III served as the initial tree, and the snaq method inferred the best phylogenetic network, testing 0 to 4 reticulation events, each optimized with 50 independent runs. Determination of the optimal number of hybrid edges followed Solís-Lemus & Ané (2016). Finally, support for all branches and reticulations in the network was estimated by 500 bootstrap replicates. More details are provided below:
1. fasta to nexus
   ```
   for f in *fasta; do
   echo $f;
   seqret -sequence "$f" -outseq "$f".nex -osformat nexus
   done
   
   
   # move each nexus to its respective folder
   ls *nex | sed 's/.fasta.nex//' > list
   while read L; do
     mkdir $L; mv "$L".fasta.nex "$L";
   done < list
   rm *fasta
   ```
2. RAxML and ASTRAL trees for individual genes. [PhyloNetworks.1.RAxML.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/PhyloNetworks.1.RAxML.sh)
3. Overall ASTRAL tree from individual gene trees. This will be the H0 starting tree. [PhyloNetworks.2.ASTRAL.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/PhyloNetworks.2.ASTRAL.sh)
4. Calculate the quartet concordance factors (CF) observed in the gene trees, and replace samples by the species names (based on the mapping file `allele-species-map.csv`). [PhyloNetworks.3.CFtable.jl](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/PhyloNetworks.3.CFtable.jl) (julia script). The `allele-species-map.csv` has following structure:
   ```
   allele,species
   acraC015_107,acris
   acraC019_103,acris
   ambC014_101,amara
   amoC046_107,amara
   ```
5. The `tableCF.astral.speciesNames.csv` contains uninformative 4-taxon sets (eg., acris,acris,acris,amara), which have to be removed.
   ```
   sed -i '/.*EBalkan.*EBalkan.*EBalkan/d' tableCF.astral.speciesNames.csv
   sed -i '/.*Dinaric.*Dinaric.*Dinaric/d' tableCF.astral.speciesNames.csv
   sed -i '/.*SharGramos.*SharGramos.*SharGramos/d' tableCF.astral.speciesNames.csv
   sed -i '/.*Pindicola.*Pindicola.*Pindicola/d' tableCF.astral.speciesNames.csv
   sed -i '/.*Vardousiae.*Vardousiae.*Vardousiae/d' tableCF.astral.speciesNames.csv
   sed -i '/.*lazica.*lazica.*lazica/d' tableCF.astral.speciesNames.csv
   sed -i '/.*matthioli.*matthioli.*matthioli/d' tableCF.astral.speciesNames.csv
   sed -i '/.*rivularis.*rivularis.*rivularis/d' tableCF.astral.speciesNames.csv
   sed -i '/.*anatolica.*anatolica.*anatolica/d' tableCF.astral.speciesNames.csv
   ```
6. Run snaq. The first run uses astral.tre, next runs start with the best hmax-1 network. [PhyloNetworks.4.snaq.jl](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/PhyloNetworks.4.snaq.jl). Estimation of the optimal number of hybridizations followed [Choice of number of hybridizations](https://github.com/JuliaPhylo/PhyloNetworks.jl/wiki/Choice-of-number-of-hybridizations).
7. BS support was determined from ASTRAL bootstrap replicates (from raxml's bootstrap gene trees). ASTRAL BS trees were collected by
   ```ruby
   for F in astral.Assembly*; do
     echo $F
     cd $F
     head -n 100 astral.tre > astral.BStrees
     echo "$F"/astral.BStrees >> ../astral.BSlistfiles
     cd ..
   done
   ```
   and BS support was inferred by script [PhyloNetworks.5.BS.jl](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/PhyloNetworks.5.BS.jl).

## PhyloSD
Due to the pipeline's requirement for a single representative of diploid genomes, the species tree from each gene tree was calculated by ASTRAL-III using [PhyloSD.1.ASTRALGeneTrees.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/PhyloSD.1.ASTRALGeneTrees.sh) script, merging 2x samples to the single representative, but keeping polyploids with separate (phased) alleles.  
#### 1) NEAREST DIPLOID SPECIES NODE algorithm
1.5) Root and sort nodes in trees, loosing bootstrap and aLRT supports on the way. (see also [Root and sort nodes in trees...](https://github.com/eead-csic-compbio/allopolyploids?tab=readme-ov-file#15-root-and-sort-nodes-in-trees-loosing-bootstrap-and-alrt-supports-on-the-way))
```ruby
mkdir 1.5_RootAndSort
parallel -j 8 "echo {}; perl5.38.2 ./PhyloSD/bin/PhyloSD/_reroot_tree.pl {} > 1.5_RootAndSort/{.}.root.ph" ::: *astralTree
```
1.6) Check diploid skeleton (topology) for each tree (see also [Check diploid skeleton (topology) for each tree)](https://github.com/eead-csic-compbio/allopolyploids?tab=readme-ov-file#16-check-diploid-skeleton-topology-for-each-tree))
```ruby
cd 1.5_RootAndSort
for FILE in *root.ph; do
   perl5.38.2 ./PhyloSD/bin/PhyloSD/_check_diploids.pl $FILE;
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
parallel -j 8 "echo {}; perl5.38.2 ./PhyloSD/bin/PhyloSD/_check_lineages_polyploids.pl -v -f {} -t {.}.raxml.bestTree.root.ph > {}.log" ::: *.fna
```
Homeologs of polyploids, that fit the criteria in `PhyloSD.polyconfig.pm` config file (defined according to the phylogenetic tree), and thus can be attributed to one of the diploid parents, were written to `label.reduced.fna` files.

#### 2) BOOTSTRAPPING REFINEMENT algorithm

2.1) Set the pruned FASTA alignments (diploids + outgroups + one polyploid homeolog) (see also [Set the pruned FASTA alignments (diploids + outgroups + one polyploid homeolog)](https://github.com/eead-csic-compbio/allopolyploids?tab=readme-ov-file#21-set-the-pruned-fasta-alignments-diploids--outgroups--one-polyploid-homeolog))  
Each sequence from `label.reduced.fna` files has to be merged with 2x samples, and the correctness of attribution to 2x parental taxa will be tested by bootstrapping.
```ruby
# keep sequences of 2x only
cp -r inputSequences inputSequences2x && cd inputSequences2x && sed -i '/acrisPP$/,+1 d' *fna

cd 1.7_congruent_and_labelled_files
mkdir ../2.1.one_allopolyploid_plus_diploids
for f in *label.reduced.fna; do
  bn=${f%.label.reduced.fna}
  for r in $(grep ">" $f | cut -f 1 -d ' '); do
    S=${r#>}
    S=${S%%.*}
    mkdir -p ../2.1.one_allopolyploid_plus_diploids/$S
    grep -A 1 "$r" "$f" > ../2.1.one_allopolyploid_plus_diploids/$S/"$bn""$r".fna
    cat ../inputSequences2x/"$bn".fna >> ../2.1.one_allopolyploid_plus_diploids/$S/"$bn""$r".fna
  done
done
```

2.2) Run 500 non-parametric bootstrapping replicates & Labelling polyploid homeologs. [PhyloSD.2.LabelBSTrees.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/PhyloSD.2.LabelBSTrees.sh). The results are in `counts` files. Those files summarize the results of re-labelling polyploid homeologs. We required confirmation by at least 40% of bootstrap replicates. That means if some homeolog was originally labelled as "SharGramos" (step 1.8), we keep that particular sequence only if more than 200 BS trees (40%) were re-labelled as "SharGramos".

2.18) "Homeologs' ML consensus tree" (see also [2.18) Phylogenomic analysis of concatenated labelled, filtered and corrected genes/MSAs](https://github.com/eead-csic-compbio/allopolyploids?tab=readme-ov-file#218-phylogenomic-analysis-of-concatenated-labelled-filtered--and-corrected-genesmsas-homeologs-ml-consensus-tree)). We concatenated sequences of 2x samples and the labelled homeologs of each polyploid (with at least 15% representation in the polyploid species). If more than one homeolog of the gene was labelled with the same 2x label, a homeolog with higher BS support was chosen. The phylogenetic tree `PhyloSD.Cacris.raxml.bestTree` was computed in RAxML-NG from the concatenated alignment, as described above.


#### 3) SUBGENOME ASSIGNMENT algorithm
Sample acraC095.109 contains 2 close homeologs, EBalkan and Dinaric. This was done to see if those homeologs refer to the same subgenome. Using the following code, we computed patristic distances, PCoA-MST (Principal coordinates analysis-minimum spanning tree), in R.

```R
library(adephylo)
library(ape)
library(stats)

tree = read.tree("PhyloSD.Cacris.raxml.bestTree")

patristicDists = distTips(tree,   method = "patristic")

distMatrix= as.matrix(patristicDists)

poylploids = distMatrix[grep("acrisPP", colnames(distMatrix)),grep("acrisPP", colnames(distMatrix))]

maf.coa <- dudi.pco(as.dist(poylploids), scannf = FALSE, nf = 3)

maf.mst <- ade4::mstree(dist.dudi(maf.coa), 1)

s.label(maf.coa$li, label = row.names(maf.coa$li),clab = 0.8,  cpoi = 2, neig = maf.mst, cnei = 1)

plot(maf.coa$li[,1],maf.coa$li[,2], asp=1)
```

3.5) Amalgamate homeologs. Two homeologs of acraC095.109 were amalgamated. If both homeologs were present for the same gene, that with higher BS support was kept.

3.6) Compute the Subgenomic ML consensus tree. The tree was computed in RAxML-NG.



***



# RADseq data processing
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
## Bayes factor species delimitation analysis (BFD*)
VCF file was subsampled, invariant SNPs were removed, and one SNP per RAD locus was selected. The VCF file was converted to phylip using `vcf2phylip.py` and the XML input file for the BEAST was created in BEAUTi using SNAPP template. The BEAST was rum from the command line. 
```ruby
beast -threads 16 $INFILE > "$INFILE".log
```
Different delimitation models were compared according to the value of the "marginal L estimate", and the Bayes factor (BF) was computed as $`BF = {2* (MLE1 − MLE0)}`$.  
The SNAPP package was further employed to estimate a coalescent-based species tree directly from SNP data. The XML input file was created in BEAUTi, and the BEAST was rum from the command line as above. Topology with the best posterior support was inferred in the TreeAnnotator.



The STRUCTURE analysis was conducted as stated above for the Hyb-Seq data, based on 100 datasets produced by selecting a single random SNP from each RADseq locus containing at least six SNPs, using vcf_prune.py script (Šlenker, 2024).
RADseq: Assessment of reticulation events and introgression
Potential admixture and reticulation events were inferred using SNaQ (Solís-Lemus & Ané, 2016; Solís-Lemus et al., 2017), neighbor-net network, STRUCTURE (Pritchard et al., 2000), and Dsuite (Malinsky et al., 2021). The SNaQ analysis was calculated as described above for the Hyb-Seq data, with the following modifications: Concordance Factors were computed from unlinked SNP by R function SNPs2CF (Olave & Meyer, 2020), and a starting tree was inferred using Quartet MaxCut algorithm (Snir & Rao, 2012). To reduce computational demands, only a single individual per population was used, except for C. acris subsp. vardousiae where two samples were included. A neighbor-net network was created using the NeighborNet algorithm in SplitsTree4 (Huson & Bryant, 2006) based on Nei’s genetic distances (Nei, 1972) calculated in the StAMPP R package (Pembleton et al., 2013). To test complex patterns of heterogeneous introgression along the genome using the ABBA–BABA and related statistics (Durand et al., 2011), Dsuite (Malinsky et al., 2021) was employed. The D, f4-ratio, and f-branch statistics was calculated, omitting the putative hybrid populations. 
RADseq: The origin of polyploid populations
To explore the origin of polyploid populations, the relatedness coefficient between the diploid lineages and each polyploid population was estimated, using the method-of-moment estimator implemented in PolyRelatedness (Huang et al., 2014). The distribution of relatedness was presented as violin plots and heatmap (Adler & Kelly, 2021; R Core Team, 2020)
RADseq: Demographic modeling, patterns of genetic diversity and rarity
The demographic history underlying the observed patterns of divergence within the C. acris complex was investigated using the diffusion approximations to the allele frequency spectrum implemented in the dadi python package (Gutenkunst et al., 2009), utilizing the routine proposed by Portik et al. (2017). The site frequency spectra were analysed and projected using easySFS (https://github.com/isaacovercast/easySFS). Each model was inferred in 50 independent runs, and the best-fitting model was selected according to AIC (Akaike information criterion) and ∆AIC scores. The 2D analysis pipeline was applied to pairwise comparisons of the genetic lineages resolved within the C. acris complex, examining whether the observed divergence patterns resulted from vicariance with ancient or more recent gene flow (secondary contact), or due to past range expansion (founder event).
Summary statistics of genetic diversity were calculated in each diploid population, comprising the nucleotide diversity (π), expected heterozygosity (He), observed heterozygosity (Ho), and private allele number (Ap), using the population program in Stacks v. 2.62 (Catchen et al., 2013). To avoid unequal sample sizes, six individuals with the lowest proportion of missing genotypes were selected per population.




