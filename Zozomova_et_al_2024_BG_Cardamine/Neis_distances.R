
library(adegenet)
library(vcfR)
library(StAMPP)

# IMPORT SNP dat from VCF - necessary to change path using setwd()
vcf <- read.vcfR("BGhybridy.bialelic.filtered.DP8.subs.passed.m02.regs.vcf.gz")   #read in all data

genlight <- vcfR2genlight.tri.MK(vcf) # run function at the end of the script

pop(genlight)<-c("243OSA","243OSA",..............,"ulULU","ulULU")



### Calculate Nei's distances between individuals/pops
D.ind <- stamppNeisD(aa.genlight, pop = FALSE)  # Nei's 1972 distance between indivs
stamppPhylip(D.ind, file="Cacris.indiv_Neis_distance.dst") # export matrix - for SplitsTree

aa.D.pop <- stamppNeisD(aa.genlight, pop = TRUE)   # Nei's 1972 distance between pops
stamppPhylip(aa.D.pop, file="Cacris.pops_Neis_distance.dst") # export matrix - for SplitsTree




vcfR2genlight.tri.MK <- function (x, n.cores = 1) 
  {
    bi <- is.biallelic(x)
    if (sum(!bi) > 0) {
      msg <- paste("Found", sum(!bi), "loci with more than two alleles.")
      msg <- c(msg, "\n", paste("Objects of class genlight only support loci with two alleles."))
      msg <- c(msg, "\n", paste(sum(!bi), "loci will be omitted from the genlight object."))
      warning(msg)
      x <- x[bi, ]
    }
    x <- addID(x)
    CHROM <- x@fix[, "CHROM"]
    POS <- x@fix[, "POS"]
    ID <- x@fix[, "ID"]
    x <- extract.gt(x)

  x[x == "0/0"] = 0      
  x[x == "0/0/0"] = 0    
  
  
  x[x == "0/0/1"] = 1    
  x[x == "0/1/1"] = 1    
  x[x == "0/1"] = 1      
  
  x[x == "1/1/1"] = 2    
  x[x == "1/1"] = 2      

  x[x == "0|0"] = 0      
  x[x == "0|0|0"] = 0    
  
  x[x == "0|0|1"] = 1    
  x[x == "0|1"] = 1      
  x[x == "0|1|1"] = 1    
  
  x[x == "1|1|1"] = 2    
  x[x == "1|1"] = 2      
  
    if (requireNamespace("adegenet")) {
      x <- new("genlight", t(x), n.cores = n.cores)
    }
    else {
      warning("adegenet not installed")
    }
    adegenet::chromosome(x) <- CHROM
    adegenet::position(x) <- POS
    adegenet::locNames(x) <- ID
    return(x)
}


