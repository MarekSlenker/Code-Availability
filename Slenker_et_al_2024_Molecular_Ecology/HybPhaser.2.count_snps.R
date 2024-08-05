# adapted from https://github.com/LarsNauheimer/HybPhaser/blob/main/1a_count_snps.R



##################
### SNPs count ###
##################


# load packages
library(ape)
library(seqinr)
library(stringr)



#####################################################################################
### Counting polymorphic sites (masked as ambiguity codes in conseneus sequences) ###
#####################################################################################

# getting targets and sample names
targets_name = read.delim("seqList", header=F)
samples <- readLines("2xSanmpes")


head(targets_name)
length(targets_name)


tab_snps <- NULL
tab_length <- NULL

# fill tables with information on snps and sequence length for each sample and locus 


for(sample in samples){
  print(paste(sample))
  
  tab_sample <- data.frame(targets=targets_name, seq_length=NA, ambis=NA, ambi_prop=NA)
  # tab_sample <- data.frame(targets=targets_name$V1, seq_length=NA, ambis=NA, ambi_prop=NA)

  fasta = read.fasta(paste(sample, ".consensus.fasta", sep = ""), as.string=TRUE, set.attributes = FALSE )

  for (seq in labels(fasta)) {
    cat(sample, ": ",seq, "\n")
    fasta[seq]
    seq2 <- gsub("N|[?]|-","",fasta[seq])
    
    tab_sample$seq_length[tab_sample$targets == names(fasta[seq])] = round(str_length(seq2),0)
    tab_sample$ambis[tab_sample$targets == names(fasta[seq])] = (str_count(seq2,"Y|K|R|S|M|y|k|r|s|m|w") + str_count(seq2,"W|D|H|B|V|w|d|h|b|v")*2)
  }

  tab_sample$ambi_prop <- tab_sample$ambis/tab_sample$seq_length  
  
  tab_snps = cbind(tab_snps, tab_sample$ambi_prop)
  tab_length = cbind(tab_length, tab_sample$seq_length)

  # tab_snps[match(sample,samples)] <- tab_sample$ambi_prop
  # tab_length[match(sample,samples)] <- tab_sample$seq_length

}


colnames(tab_snps) <- samples
rownames(tab_snps) <- targets_name
#rownames(tab_snps) <- targets_name$V1

colnames(tab_length) <- samples
rownames(tab_length) <- targets_name
# rownames(tab_length) <- targets_name$V1

### Generate output tables and save data in Robjects


saveRDS(tab_snps,file="Table_SNPs.Rds")
saveRDS(tab_length,file="Table_consensus_length.Rds")

write.table(tab_snps,file="Table_SNPs.txt")
write.table(tab_length,file="Table_consensus_length.txt")
