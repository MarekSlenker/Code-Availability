#!/bin/bash




gatk --java-options "-Xmx120g" SelectVariants \
-V BGhybridy.vcf.gz.bialelic.filtered.DP8.passed.m02.inRegs.vcf.gz \
-O ABBA_BABA/BGhybridy.vcf.gz.bialelic.filtered.DP8.passed.m02.inRegs.ABABABA.noInvs.vcf.gz \
--exclude-non-variants \
--sample-name toTake.args


cd ABBA_BABA

gatk --java-options "-Xmx120g" VariantsToTable \
-V BGhybridy.vcf.gz.bialelic.filtered.DP8.passed.m02.inRegs.ABABABA.noInvs.vcf.gz \
-O BGhybridy.geno \
-F CHROM -F POS -GF GT
sed -i 's/\.\/\.\/\./N\/N\/N/g' BGhybridy.geno
sed -i 's/\.\/\./N\/N/g' BGhybridy.geno

sed -i 's/\.|\.|./N|N|N/g' BGhybridy.geno
sed -i 's/\.|\./N|N/g' BGhybridy.geno

sed -i 's/.GT//g' BGhybridy.geno


# https://github.com/simonhmartin/genomics_general
python <PATH>/ABBA_BABA/genomics_general/freq.py \
-g BGhybridy.geno \
-p BEG_amara -p BEG_riv -p BEG_rivAm -p BEG_rivMat -p BIT_amara -p BIT_mat -p BIT_riv -p BIT_rivAm -p JEB_amara -p JEB_mat -p JEB_riv -p JEB_rivAm -p JEB_rivMat -p KAT_mat -p KAT_riv -p KAT_rivAm -p KAT_rivMat -p POB_amara -p POB_mat -p POB_riv -p POB_rivAm -p YUN_mat -p YUN_riv -p YUN_rivMat -p ZEL_mat -p ZEL_rivAc -p ZEL_rivMatAc -p mat -p riv -p acris -p amara -p outgrp \
--popsFile BGhybridy.pop \
--ploidyFile BGhybridy.pld \
--target derived \
-o BGhybridy.tsv.gz

# -p PRD_mat  -p JEB_mat -p JEB_rivMat -p outgrp \



# MOVE TO R
module add r/4.4.0-gcc-10.2.1-oxdi5pz 
R


D.stat <- function(p1, p2, p3) {
    ABBA <- (1 - p1) * p2 * p3
    BABA <- p1 * (1 - p2) * p3
    (sum(ABBA, na.rm=T) - sum(BABA, na.rm=T)) / (sum(ABBA, na.rm=T) + sum(BABA, na.rm=T))
    }

freq_table = read.table("BGhybridy.tsv.gz", header=T, as.is=T)

nrow(freq_table)

head(freq_table)



########################################
# The D statistic
########################################




P1 = "riv"
P2 <- "BEG_riv"
P3 <- "BEG_rivMat"



D <- D.stat(freq_table[,P1], freq_table[,P2], freq_table[,P3])
print(paste("D =", round(D,4)))


########################################
# Block Jackknife
########################################
# source("<PATH>/ABBA_BABA/genomics_general/jackknife.R")


block_indices <- get.block.indices(block_size=1e6, positions=freq_table$position, chromosomes=freq_table$scaffold)
n_blocks <- length(block_indices)
print(paste("Genome divided into", n_blocks, "blocks."))

D_jackknife <- block.jackknife(block_indices=block_indices,
                               FUN=D.stat,
                               freq_table[,P1], freq_table[,P2], freq_table[,P3])
print(paste("D jackknife mean =", round(D_jackknife$mean,4)))


D_Z <- D_jackknife$mean / D_jackknife$standard_error
print(paste("D Z score = ", round(D_Z,3)))
write.table(D_jackknife$pseudovalues, paste(P3,".",P1,".",P2,".",".D_jackknife", sep=""), col.names=F, row.names=F)




