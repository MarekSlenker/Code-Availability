<table><tr><td>https://github.com/lifeparticle/Markdown-Cheatsheet</td></tr></table>

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



# Read-backed phasing

The code (roughly) follows the procedures outlined in the [alleles_workflow](https://github.com/mossmatters/phyloscripts/tree/master/alleles_workflow) GitHub repository. The following scripts implement these procedures.  
The sequences of each sample are phased by [Phasing.1.phasing_oneSample.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/Phasing.1.phasing_oneSample.sh). This script takes consensus sequences and fastq reads at the input and produces phased sequences `"$SAMPLE".v1.phased.fasta, "$SAMPLE".v2.phased.fasta, ...` and unphased `"$SAMPLE".unPhased.fasta` sequences. Phased sequences are sometimes represented by multiple mutually unphased blocks (take a look at `"$SAMPLE".whatshap.gtf`). Selection of the longest phased block and masking of remaining variant is responsibility of [Phasing.2.masking.sh](https://github.com/MarekSlenker/Code-Availability/blob/main/Slenker_et_al_2024_Molecular_Ecology/Phasing.2.masking.sh) script.  

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









# don't read any further, I'm still working on it 
# don't read any further, I'm still working on it
# don't read any further, I'm still working on it
# don't read any further, I'm still working on it
# don't read any further, I'm still working on it



The resulting sequences were forwarded for read-backed phasing to infer allele sequences, as described in detail in Šlenker et al. (2021). Consensus exon sequences were concatenated to genes using AMAS (Borowiec, 2016), while phased exon sequences were not, as they represent unphasable blocks. Both the phased exons and consensus genes were used in downstream analyses.
Hyb-Seq: Phylogenetic analyses and Bayesian clustering
Maximum likelihood (ML) trees were inferred both from concatenated genes and for each gene separately using RAxML-NG v.0.9.0 (Kozlov et al., 2019). Each gene was treated as a separate partition with its best-fitting substitution model determined through the ModelFinder function of IQ-TREE v.1.6.12 (Chernomor et al., 2016; Kalyaanamoorthy et al., 2017) based on the Bayesian information criterion. Bootstrap analyses were performed using 500 replicates. For the species tree reconstruction, internal branches with bootstrap support ≤20% were collapsed using Newick-Utilities v. 1.6 (Junier & Zdobnov, 2010). The species tree was constructed employing a multispecies coalescent model implemented in ASTRAL-III (Zhang et al., 2018), including computation of local posterior probabilities to evaluate branch support (Sayyari & Mirarab, 2016).
In addition, homogeneous genetic clusters were inferred using Bayesian clustering algorithm implemented in STRUCTURE v. 2.3.4 (Pritchard et al., 2000). The snipStrup pipeline (https://github.com/MarekSlenker/snipStrup) was used to map the Hyb-Seq reads to the consensus sequences, call variants and convert the VCF file. The STRUCTURE computations were performed and summarised as described in Šlenker et al. (2021).

Hyb-Seq: Assessment of reticulation events
Reticulation events were inferred using SNaQ (Solís-Lemus & Ané, 2016; implemented in PhyloNetworks v. 0.16.2, Solís-Lemus et al., 2017) with 500 gene trees, computed from the longest gene alignments, where internal branches with bootstrap support ≤20% were collapsed. The gene trees were summarized using quartet concordance factors. The species tree reconstructed in ASTRAL-III served as the initial tree, and the snaq method inferred the best phylogenetic network, testing 0 to 4 reticulation events, each optimized with 50 independent runs. Determination of the optimal number of hybrid edges followed Solís-Lemus & Ané (2016). Finally, support for all branches and reticulations in the network was estimated by 500 bootstrap replicates.
Hyb-Seq: The origin of polyploid populations
Phased exon sequences of polyploids were processed by the PhyloSD pipeline (Sancho et al., 2022) to identify the homeologous diploid genomes. Due to the pipeline's requirement for a single representative of diploid genomes, species tree of each gene was calculated using ASTRAL-III. Due to unacceptable loss of data, incongruent diploid skeletons were not discarded (unlike in Sancho et al., 2022), but stricter criteria were applied in the Bootstrapping Refinement step, requiring confirmation by at least 40% of bootstrap replicates. Homeologs’ ML consensus tree was constructed from the homeologs with at least 15% representation in the polyploid species. The Subgenome Assignment algorithm collapsed homeologs referring to the same subgenomes according to the principal coordinate analysis and superimposed minimum spanning tree generated from patristic distances (following Sancho et al., 2022), and the final Subgenomic ML consensus tree was constructed. 

RADseq data processing
The RADseq reads were mapped on the genome of C. amara (Šlenker et al., in prep.) using BWA 0.7.5a (Li, 2013) and the resulting BAM files were processed with Picard Toolkit 2.22.1. (https://broadinstitute.github.io/picard/). Variant calling was performed for each individual using the HaplotypeCaller module from the GATK 4.4.0.0 (McKenna et al., 2010). As next, variants were aggregated and genotyping across all individuals was performed using the GenomicsDBImport and GenotypeGVCFs modules. Biallelic sites with a minimum sequencing depth of 8x, passing the filter parameters indicated by GATK’s best practices (Van der Auwera et al., 2013), and with no more than 30% of missing genotypes were captured using VariantFiltration and SelectVariants modules. Finally, samples with more than 60 % missing genotypes were excluded. Certain analyses required unlinked SNPs (see below), and this was achieved by selecting a single random SNP from each RADseq locus. The loci were identified following identifiRadLoci.workflow (Šlenker, 2024), requiring a minimum sequencing depth of 8x observed in at least 70% of the samples, and a minimum distance of 1,000 bp, collapsing regions less than 1,000 bp apart into a single RAD locus. 

RADseq: Phylogenetic analyses, species delimitation and Bayesian clustering
Phylogenetic relationships were inferred using concatenation and species tree methods. The ML tree was constructed by RAxML-NG v.0.9.0 (Kozlov et al., 2019), employing GTR model with Felsenstein’s ascertainment bias correction. The vcf2phylip.py script (Ortiz, 2019) was used to transform the data from the VCF file to the PHYLIP, and invariant sites were removed with the script ascbias.py (https://github.com/btmartin721/raxml_ascbias) as recommended by Leaché et al. (2015). 
To provide statistical support for the delimitation of inferred genetic clusters in C. acris, Bayes factor species delimitation analysis (BFD*, Leaché et al., 2014) was performed following Leaché & Bouckaert (2018). Marginal likelihoods of species trees were inferred using the Path Sampling approach with the SNAPP v.1.4.2 (Bryant et al., 2012) and BEAST v. 2.5.0 (Bouckaert et al., 2014). The dataset was reduced to unlinked SNPs and three samples per each genetic cluster, except for C. acris subsp. vardousiae with two representative samples. Analyses were run in eight steps for each model, with 1,000,000 MCMC iterations, sampling every 1,000th, and a burn-in cutoff of 10%. Competing species delimitation models were ranked by comparing their marginal likelihood estimates and their support was assessed by calculating the Bayes factor (Kass & Raftery, 1995), as suggested by Leaché & Bouckaert (2018). The current taxonomy model, encompassing four taxa (C. anatolica and C. acris with three subspecies), was compared with four alternative species models based on the ML tree and Bayesian clustering outputs, where C. acris subsp. acris was split into several entities. The SNAPP package was further employed to estimate a coalescent-based species tree directly from SNP data (reduced as in BFD*), using the species model with the highest support in BFD* described above. Unlinked SNPs in binary nexus format were processed in the BEAUti, and 5,000,000 MCMC iterations were performed, logging every 1,000th tree. The consensus tree topology with the best posterior support was identified by TreeAnnotator (Drummond & Rambaut, 2007) with 10% burn-in. Samples from two populations (C005, C241) were excluded from the tree reconstruction methods and species delimitation analysis, as they were revealed to be of hybrid (allopolyploid) origins (see below).
The STRUCTURE analysis was conducted as stated above for the Hyb-Seq data, based on 100 datasets produced by selecting a single random SNP from each RADseq locus containing at least six SNPs, using vcf_prune.py script (Šlenker, 2024).
RADseq: Assessment of reticulation events and introgression
Potential admixture and reticulation events were inferred using SNaQ (Solís-Lemus & Ané, 2016; Solís-Lemus et al., 2017), neighbor-net network, STRUCTURE (Pritchard et al., 2000), and Dsuite (Malinsky et al., 2021). The SNaQ analysis was calculated as described above for the Hyb-Seq data, with the following modifications: Concordance Factors were computed from unlinked SNP by R function SNPs2CF (Olave & Meyer, 2020), and a starting tree was inferred using Quartet MaxCut algorithm (Snir & Rao, 2012). To reduce computational demands, only a single individual per population was used, except for C. acris subsp. vardousiae where two samples were included. A neighbor-net network was created using the NeighborNet algorithm in SplitsTree4 (Huson & Bryant, 2006) based on Nei’s genetic distances (Nei, 1972) calculated in the StAMPP R package (Pembleton et al., 2013). To test complex patterns of heterogeneous introgression along the genome using the ABBA–BABA and related statistics (Durand et al., 2011), Dsuite (Malinsky et al., 2021) was employed. The D, f4-ratio, and f-branch statistics was calculated, omitting the putative hybrid populations. 
RADseq: The origin of polyploid populations
To explore the origin of polyploid populations, the relatedness coefficient between the diploid lineages and each polyploid population was estimated, using the method-of-moment estimator implemented in PolyRelatedness (Huang et al., 2014). The distribution of relatedness was presented as violin plots and heatmap (Adler & Kelly, 2021; R Core Team, 2020)
RADseq: Demographic modeling, patterns of genetic diversity and rarity
The demographic history underlying the observed patterns of divergence within the C. acris complex was investigated using the diffusion approximations to the allele frequency spectrum implemented in the dadi python package (Gutenkunst et al., 2009), utilizing the routine proposed by Portik et al. (2017). The site frequency spectra were analysed and projected using easySFS (https://github.com/isaacovercast/easySFS). Each model was inferred in 50 independent runs, and the best-fitting model was selected according to AIC (Akaike information criterion) and ∆AIC scores. The 2D analysis pipeline was applied to pairwise comparisons of the genetic lineages resolved within the C. acris complex, examining whether the observed divergence patterns resulted from vicariance with ancient or more recent gene flow (secondary contact), or due to past range expansion (founder event).
Summary statistics of genetic diversity were calculated in each diploid population, comprising the nucleotide diversity (π), expected heterozygosity (He), observed heterozygosity (Ho), and private allele number (Ap), using the population program in Stacks v. 2.62 (Catchen et al., 2013). To avoid unequal sample sizes, six individuals with the lowest proportion of missing genotypes were selected per population.




