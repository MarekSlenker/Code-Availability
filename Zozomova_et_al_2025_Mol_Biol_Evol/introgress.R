
library(introgress)
library(vcfR)



# AdmixData
vcf=read.vcfR("BGhybridy.bialelic.filtered.DP8.passed.m02.inRegs.bezOutgroups.noInvs.MatRiv+HZ.najDifSNPS.Ako2x.vcf")


AdmixData<-extract.gt(vcf)

riv=read.table("riv", header = F)  # names of samples of C. rivularis.
cols=colnames(AdmixData)

# sometimes rivularis is 0/0 and sometimes 1/1. I need to check it and set P1/P1 for rivularis, regardless of 0/0 or 1/1
for (i in 1:dim(AdmixData)[1]) {
  if ( # RIV is 0/0
    length(which((na.omit(unlist(strsplit(AdmixData[i, which(cols %in%  riv$V1)], split = "/"))))=="0")) > length(which((na.omit(unlist(strsplit(AdmixData[i, which(cols %in%  riv$V1)], split = "/"))))=="1"))
  ) {
    AdmixData[i, which(AdmixData[i,] == "0/0")] = "P1/P1"
    AdmixData[i, which(AdmixData[i,] == "1/1")] = "P2/P2"
  } else { #  RIV is 1/1
    AdmixData[i, which(AdmixData[i,] == "1/1")] = "P1/P1"
    AdmixData[i, which(AdmixData[i,] == "0/0")] = "P2/P2"
  }
}

AdmixData[which(AdmixData[,] == "0/1")] = "P1/P2"
AdmixData[which(AdmixData[,] == "NA")] = "NA/NA"



# LociData
LociData = cbind(locus=paste(vcf@fix[,1], vcf@fix[,2], sep = "."),
      type="C",
      lg=as.numeric(as.factor(vcf@fix[,1])),
      marker.pos=paste(as.numeric(as.factor(vcf@fix[,1])),vcf@fix[,2], sep = "."))

# prepare data
count.matrix=prepare.data(admix.gen = AdmixData, loci.data = LociData,
                          parental1 = "P1", parental2 = "P2", pop.id = F,ind.id = F, fixed = T)



# estimate HIBRID INDEX
hi.index = est.h(introgress.data = count.matrix, loci.data = LociData, 
                     fixed = T, p1.allele = "P1", p2.allele = "P2")
write.table(hi.index, "h.index.MatRiv.txt", quote = F, sep = ",")

mk.image(count.matrix, LociData, marker.order = NULL,
         hi.index = hi.index, ylab.image = "indivs",
         xlab.h = "pop 2 ancesty", pdf = T, out.file = "img2.pdf")



# Calculate Interspecific Heterozygosity
int.het = calc.intersp.het(introgress.data=count.matrix)
write.table(int.het, "int.het.MatRiv.txt", quote = F, sep = ",")



