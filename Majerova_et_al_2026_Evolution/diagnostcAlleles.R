#!/usr/bin/env Rscript

library(vcfR)

vcfAmpApe= "Cardamine_AmpApe123.filtered.DP8.bilallelic.passed.m06.bezBlbych.bezOpakovani.AmpApe.1408.vcf"
mapAmpApe= "ampApe.map"
vcfHybr="Cardamine_AmpApe123.filtered.DP8.bilallelic.passed.m06.bezBlbych.bezOpakovani.HYBR.1408.vcf"
mapHybr="Hybr.map"
requiredSampleFractionWithAllele= 0.3
allowedMissingDataPerSNP=0.6


# DIAGNOSTIC ALLELES

vcf = read.vcfR(vcfAmpApe)
AmpApeMap = read.delim(mapAmpApe, header = F)
colnames(AmpApeMap) = c("sample", "pop")

gt = extract.gt(vcf, element = "GT", as.numeric = FALSE)
fix = as.data.frame(vcf@fix, stringsAsFactors = FALSE)
  
diagnosticAlleles <- data.frame()
  
for (snp in seq(1,nrow(gt))) {
    chrom <- fix$CHROM[snp]
    pos <- fix$POS[snp]
    ref <- fix$REF[snp]
    alt <- fix$ALT[snp]
    
    # alleles <- c(ref, unlist(strsplit(alt, ",", fixed = TRUE)))
    
    # Determine each allele
    for (allele in c(0,1)) {  # Expecting only biallelic SNPs
      #counts <- carrier_counts_by_pop(gt[i, ], AmpApeMap, allele_index)
      
      countsByPopulation = data.frame()
      
      for (pop in unique(AmpApeMap$pop)) {
        popSamples <- AmpApeMap$sample[AmpApeMap$pop == pop]
        
        alleleCount=0
        presentSamples = 0
        missingSamples = 0
        
          for (s in popSamples) {
            # gt[snp,popSamples]
            
            if (is.na(gt[snp,s])) {
              missingSamples = missingSamples+1
              next
            }
            alleles <- unlist(strsplit(gt[snp,s], "[/|]"))
            presentSamples <- presentSamples + 1
           
            if (allele %in% alleles) {
              alleleCount=alleleCount + 1
            }
            
          }
        
        sampleFractionWithAllele = ifelse(presentSamples > 0, alleleCount / presentSamples,0)
        
        countsByPopulation = rbind(countsByPopulation, data.frame(
                              population = pop,
                              presentSamples = presentSamples,
                              missingSamples = missingSamples,
                              totalSamples = length(popSamples),
                              samplesWithAllele = alleleCount,
                              sampleFractionWithAllele = sampleFractionWithAllele,
                              missingFraction=missingSamples/length(popSamples),
                              stringsAsFactors = FALSE
          )
        )
      }
      
      # Allele counts and frequencies for each population
      countsByPopulation
      
      # Find diagnostic SNPs
      for (j in seq(1, nrow(countsByPopulation))) {
        if (countsByPopulation$missingFraction[j] > allowedMissingDataPerSNP) {
          next
        }
        if (countsByPopulation$sampleFractionWithAllele[j] < requiredSampleFractionWithAllele) {
          next
        }
        focalPop <- countsByPopulation$population[j]
        
        otherPops <- countsByPopulation[countsByPopulation$population != focalPop, ]
        
        absentInOtherPops = all(otherPops$samplesWithAllele == 0)
        
        if (absentInOtherPops) {
          diagnosticAlleles = rbind(diagnosticAlleles, data.frame(
                                                          CHROM = chrom,
                                                          POS = pos,
                                                          REF = ref,
                                                          ALT = alt,
                                                          allele = allele,
                                                          diagnosticForPopulation = focalPop,
                                                          presentSamples = countsByPopulation$presentSamples[j],
                                                          missingSamples = countsByPopulation$missingSamples[j],
                                                          totalSamples = countsByPopulation$totalSamples[j],
                                                          samplesWithAllele = countsByPopulation$samplesWithAllele[j],
                                                          sampleFractionWithAllele = countsByPopulation$sampleFractionWithAllele[j],
                                                          missingFraction=countsByPopulation$missingFraction[j],
                                                          stringsAsFactors = FALSE
            )
          )
        }
      }
      
      
      
    }
  
}
  
diagnosticAlleles
rm(list=setdiff(ls(), c("diagnosticAlleles", "vcfHybr", "mapHybr" )))



################################################
# COUNT ALLELES
################################################


vcf = read.vcfR(vcfHybr)
HybrMapFile= read.delim(mapHybr, header = F)
colnames(HybrMapFile) = c("sample", "pop")


gt = extract.gt(vcf, element = "GT", as.numeric = FALSE)
fix = as.data.frame(vcf@fix, stringsAsFactors = FALSE)



fix$ID=paste(fix$CHROM, fix$POS, fix$REF, fix$ALT, sep = "_")
diagnosticAlleles$ID=paste(diagnosticAlleles$CHROM,diagnosticAlleles$POS,diagnosticAlleles$REF,diagnosticAlleles$ALT,sep = "_")


countsInHybridPopulations = data.frame()

for (i in seq(1,nrow(diagnosticAlleles))) {
  
      hybridVcfRow = which(fix$ID == diagnosticAlleles[i, ]$ID)
      
      
      for (pop in unique(HybrMapFile$pop)) {
        popSamples = HybrMapFile$sample[HybrMapFile$pop == pop]
        
        alleleCount=0
        presentSamples = 0
        missingSamples = 0
        
        for (s in popSamples) {
          # gt[hybridVcfRow,popSamples]
          
          if (is.na(gt[hybridVcfRow,s])) {
            missingSamples = missingSamples+1
            next
          }
          alleles <- unlist(strsplit(gt[hybridVcfRow,s], "[/|]"))
          presentSamples <- presentSamples + 1
          
          if (diagnosticAlleles[i, ]$allele %in% alleles) {
            alleleCount=alleleCount + 1
          }
          
        }
        
        sampleFractionWithAllele = ifelse(presentSamples > 0, alleleCount / presentSamples,0)
        if (sampleFractionWithAllele > 0) {
                  countsInHybridPopulations = rbind(countsInHybridPopulations, data.frame(
                    CHROM=diagnosticAlleles$CHROM[i],
                    POS=diagnosticAlleles$POS[i],
                    REF=diagnosticAlleles$REF[i],
                    ALT=diagnosticAlleles$ALT[i],
                    population = pop,
                    presentSamples = presentSamples,
                    missingSamples = missingSamples,
                    totalSamples = length(popSamples),
                    samplesWithAllele = alleleCount,
                    sampleFractionWithAllele = sampleFractionWithAllele,
                    missingFraction=missingSamples/length(popSamples),
                    stringsAsFactors = FALSE)
                  )
          }
        }
    
}
  
# pocty alely v kazdej pop, a aj pocetnosti, ...
countsInHybridPopulations
  
  
write.table(countsInHybridPopulations, 
            file = "diagnosticAllelesInHybrids.txt",
            sep = "\t", 
            quote = FALSE, row.names = FALSE)

