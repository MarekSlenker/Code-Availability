# adopted from https://github.com/LarsNauheimer/HybPhaser/blob/main/1b_assess_dataset.R

#####################################################################################
### Generating files for assessment of missing data,. paralogs and heterozygosity ###
#####################################################################################


# load packages
library(ape)


tab_snps <- readRDS(file="Table_SNPs.Rds")
tab_length <- readRDS(file="Table_consensus_length.Rds")

tab_snps <- as.matrix(tab_snps)
loci <- t(tab_snps)


##########################################################

nloci <- dim(loci)[2]
nsamples <- dim(loci)[1]

failed_loci <- which(colSums(is.na(loci))==nrow(loci))
failed_samples <- which(colSums(is.na(tab_snps))==nrow(tab_snps))

# per locus

seq_per_locus <- vector()
for(i in 1:nloci){
  seq_per_locus[i] <- length(which(!(is.na(loci[,i]))))
}
names(seq_per_locus) <- colnames(loci)
seq_per_locus_prop <- seq_per_locus/nsamples


# per sample

seq_per_sample <- vector()
for(i in 1:length(colnames(tab_snps))){
  seq_per_sample[i] <- length(which(!(is.na(tab_snps[,i]))))
}
names(seq_per_sample) <- colnames(tab_snps)
seq_per_sample_prop <- seq_per_sample/nloci



tab_snps_cl1 <- tab_snps
loci_cl1 <- t(tab_snps_cl1)




############################################################################################
### Dataset optimization step 2, removing paralogs for a) all samples and b) each sample ### 
############################################################################################


### 2a) Paralogs across multiple samples (removing loci with unusually high proportions of SNPs across all samples)
###################################################################################################################

loci_cl1_colmeans <- colMeans(as.matrix(loci_cl1), na.rm = T)
nloci_cl1 <- length(colnames(loci_cl1))
nsamples_cl1 <- length(colnames(tab_snps_cl1))

loci_cl1_colmeans_mean <- round(mean(loci_cl1_colmeans),4)
loci_cl1_colmeans_median <- round(median(loci_cl1_colmeans),4)


# applying chosen threshold

# "outliers"
  threshold_value <- 1.5*IQR(loci_cl1_colmeans, na.rm = TRUE )+quantile(loci_cl1_colmeans, na.rm = TRUE )[4]
  outloci_para_all_values <- loci_cl1_colmeans[which(loci_cl1_colmeans > threshold_value)]
  outloci_para_all <- names(outloci_para_all_values)
  outloci_para_all
  write.table(outloci_para_all, file="outloci_IQR", row.names=F, col.names=F, quote=F)
  

# threshold_value = NUMERIC
  threshold_value = 0.02  # 2%
  outloci_para_all_values <- loci_cl1_colmeans[which(loci_cl1_colmeans > threshold_value)]
  outloci_para_all <- names(outloci_para_all_values)
  outloci_para_all
  write.table(outloci_para_all, file="outloci_2%", row.names=F, col.names=F, quote=F)
  
  threshold_value = 0.05  # 5%
  outloci_para_all_values <- loci_cl1_colmeans[which(loci_cl1_colmeans > threshold_value)]
  outloci_para_all <- names(outloci_para_all_values)
  outloci_para_all
  write.table(outloci_para_all, file="outloci_5%", row.names=F, col.names=F, quote=F)

# color outliers red
colour_outparaall <- rep("black",nloci_cl1)
colour_outparaall[which(colnames(loci_cl1[,order(loci_cl1_colmeans)]) %in% outloci_para_all)] <- "red"
loci_cl1_order_means <- loci_cl1[,order(loci_cl1_colmeans)]


# generate bar graph
for(i in 1:2){
  
  if(i==1){
    pdf(file="2a_Paralogs_for_all_samples.pdf", width = 11, height=7)
  } else {
    png(file="2a_Paralogs_for_all_samples.png", width = 1400, height=1000)
    par(cex.axis=2, cex.lab=2, cex.main=2)
  }
  
  layout(matrix(c(1,2),2,2, byrow=TRUE), widths=c(5,1))
  
  barplot(sort(loci_cl1_colmeans), col=colour_outparaall, border = NA, las=2,
          main=paste("Mean % SNPs across samples (n=",nsamples_cl1,") for each locus (n=", nloci_cl1,")", sep=""))
  if(length(threshold_value)>0){abline(h=threshold_value, col="red", lty=2)}
  boxplot(loci_cl1_colmeans, las=2)
  if(length(threshold_value)>0){abline(h=threshold_value, col="red", lty=2)}
  
  dev.off()
}
