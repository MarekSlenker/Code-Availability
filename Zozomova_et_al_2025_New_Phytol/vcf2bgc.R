

library(vcfR)

inVCF="RivMatHybr.unlinked.vcf.gz"
prefix="RivMatHybr"
p0="riv.samples"
p1="mat.samples"
adm="MatRivHybs.samples"



vcf=read.vcfR(inVCF)


# AdmixData

AdmixData<-extract.gt(vcf)

P0=read.table(p0, header = F)
P1=read.table(p1, header = F)
ADM=read.table(adm, header = F)
cols=colnames(AdmixData)
loci=rownames(AdmixData)

dt= AdmixData[,which(colnames(AdmixData) %in% P0$V1)]
for (i in 1:length(loci)) {
  ll = as.matrix(dt[i,]  )
  ll[which(ll == "0/0")] = "2  0"
  ll[which(ll == "0|0")] = "2  0"
  ll[which(ll == "1/1")] = "0  2"
  ll[which(ll == "1|1")] = "0  2"
  ll[which(ll == "0/1")] = "1  1"
  ll[which(ll == "0|1")] = "1  1"
  ll[which(is.na(ll))] = "0  0"
  rbind(paste("locus", i),ll)
  write.table(rbind(paste("locus", i),ll), paste(prefix,"_p0in.txt", sep = ""),
              append = T, quote = F, row.names = F, col.names = F)
}



dt= AdmixData[,which(colnames(AdmixData) %in% P1$V1)]
for (i in 1:length(loci)) {
  ll = as.matrix(dt[i,]  )
  ll[which(ll == "0/0")] = "2  0"
  ll[which(ll == "0|0")] = "2  0"
  ll[which(ll == "1/1")] = "0  2"
  ll[which(ll == "1|1")] = "0  2"
  ll[which(ll == "0/1")] = "1  1"
  ll[which(ll == "0|1")] = "1  1"
  ll[which(is.na(ll))] = "0  0"
  rbind(paste("locus", i),ll)
  write.table(rbind(paste("locus", i),ll), paste(prefix,"_p1in.txt", sep = ""),
              append = T, quote = F, row.names = F, col.names = F)
}



dt= AdmixData[,which(colnames(AdmixData) %in% ADM$V1)]
for (i in 1:length(loci)) {
  ll = as.matrix(dt[i,]  )
  ll[which(ll == "0/0")] = "2  0"
  ll[which(ll == "0|0")] = "2  0"
  ll[which(ll == "1/1")] = "0  2"
  ll[which(ll == "1|1")] = "0  2"
  ll[which(ll == "0/1")] = "1  1"
  ll[which(ll == "0|1")] = "1  1"
  ll[which(is.na(ll))] = "-9  -9"
  rbind(paste("locus", i),ll)
  write.table(rbind(paste("locus", i),"pop 0",ll), paste(prefix,"_admixedin.txt", sep = ""),
              append = T, quote = F, row.names = F, col.names = F)
}
