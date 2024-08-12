

library(gplots)


mm = as.matrix(read.table("Cacris_sl.genotypes.t.polyrelatedness.OUT.txt"))


heatmap(mm, Rowv = NA, Colv = NA, revC = T)   
# heatmap(mm, revC = T, scale="column")   


ss = read.table("C243.txt", header = F)
head(ss)
ss2 = ss[!(ss$V1 %in% c("C018","C095", "C192", "C241", "C242", "C243")),] # removing polyploids

vioplot::vioplot(ss2$V2~ss2$V1)

library(ggplot2)
ggplot(ss2, aes(x=V1, y=V2, fill=V1)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin()
