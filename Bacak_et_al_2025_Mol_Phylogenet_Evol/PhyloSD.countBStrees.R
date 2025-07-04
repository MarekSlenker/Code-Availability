c = read.delim("vybrane.counts")

colSums(c[,-1])


write.table(t(colSums(c[,-1])), "t", quote = F, row.names = F)
