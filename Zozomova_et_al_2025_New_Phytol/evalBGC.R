



library("ClineHelpR")

genesDIR="<PATH>/BGC/RivMat.unlinked/BGCOutputFiles"
plotDIR="<PATH>/BGC/RivMat.unlinked"
prefix="RivMatHybr"
admixPop="RivMatHybr"
popmap="<PATH>/BGC/RivMat.unlinked/BGC.popmat"
genes.loci.file="<PATH>/BGC/RivMat.unlinked/BGC.loci.txt"

###############################################################################
# Run BGC plotting functions
###############################################################################

##### Transcripomic alignment ######

# Combine multiple independent BGC runs together
bgc.genes <-
  combine_bgc_output(results.dir = genesDIR,
                     prefix = prefix)

# # Plot the likelihood and parameter traces
# plot_traces(df.list = bgc.genes,
#             prefix = prefix,
#             plotDIR = plotDIR)

# Detect BGC outliers for known genes
gene.outliers <-
  get_bgc_outliers(
    df.list = bgc.genes,
    admix.pop = admixPop,
    popmap = popmap,
    loci.file = genes.loci.file
  )

write.table(gene.outliers[1], "gene.outliers.1.txt", row.names=F )
write.table(gene.outliers[2], "gene.outliers.2.txt", row.names=F)
write.table(gene.outliers[3], "gene.outliers.3.txt", row.names=F)

##### Scaffold alignment #####

# Aggregate runs together

# Make the phi plot for the transcriptomic alignment
# Any of these parameters can be adjusted as needed.
# Here, both.outlier.tests is FALSE
# This means that outliers are flagged if they are significant in either method
phiPlot(
  outlier.list = gene.outliers,
  popname = paste0(admixPop, " Genes"),
  line.size = 0.35,
  saveToFile = paste0(prefix, "_genes"),
  plotDIR = plotDIR,
  both.outlier.tests = FALSE,
  neutral.color = "gray60",
  alpha.color = "cornflowerblue",
  beta.color = "firebrick",
  both.color = "purple",
  hist.y.origin = 1.2,
  hist.height = 1.8,
  margins = c(160.0, 5.5, 5.5, 5.5),
  hist.binwidth = 0.05
)


# ABPLOT LOKALNE

gene.outliers.1=read.delim("gene.outliers.1.txt", sep = "", header = T)
gene.outliers.2=read.delim("gene.outliers.2.txt", sep = "", header = T)
gene.outliers.3=read.delim("gene.outliers.3.txt", sep = "", header = T)

gene.outliers=(list(gene.outliers.1, gene.outliers.2, gene.outliers.3))

library("ClineHelpR")


prefix="RivMatHybr"

alphaBetaPlot(
  gene.outliers,
  alpha.color = "cornflowerblue",
  beta.color = "orange",
  neutral.color = "gray60",
  saveToFile = prefix,
  plotDIR = ".",
  padding = 0.2,
)



###############xx

cat gene.outliers.1.txt |cut -f 5 -d " " |sort |uniq -c
cat gene.outliers.1.txt |cut -f 6 -d " " |sort |uniq -c
cat gene.outliers.1.txt |cut -f 7 -d " " |sort |uniq -c
cat gene.outliers.1.txt |cut -f 8 -d " " |sort |uniq -c



