#!/bin/bash



####
##   Step 0: concatenate all partial R E libs       >merged-library.fa

# mitefinder
cat <PATH>/mitefinder_output/C087_203_hap1_mitefinder.fasta >merged-library.fa

# HelitronScanner
cat <PATH>/HelitronScanner/C087_203_hap1.fa_helitron.hel.fa >>merged-library.fa

# SineScan
cat <PATH>/sinescan_output/C087_203_hap1.sine.fa >>merged-library.fa

# RepeatModeler2
cat <PATH>/RepeatModeler2/C087_203_hap1-families.fa >>merged-library.fa

# LTR_FINDER + LTR_retriever
cat <PATH>/LTR_retriever/C087_203_hap1.fa.LTRlib.fa >>merged-library.fa

# Dfam and RepBase
cat <PATH>/dfam-tetools/Libraries/Brassicaceae.famdbLib.fasta >>merged-library.fa



####
## Step 2: Second VSEARCH Command

vsearch-2.30.0-linux-x86_64/bin/vsearch \
--cluster_fast merged-library.fa \
--id 0.95 \
--centroids merged-library.centroids.fa \
--uc result.uc \
-consout consensus.fa \
-msaout aligned.fasta \
--log vsearch2.log


# ####
# ## Step 3 - Renaming and classification of final consensus
sed -i 's/centroid=*//' final.nr.consensus.fa
sed -i 's/;seqs=[0-9]*$//' final.nr.consensus.fa




