library(tissuesGeneExpression)
data(tissuesGeneExpression)

library(RColorBrewer) 
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

library(genefilter)
rv <- rowVars(e)
idx <- order(-rv)[1:50]

library(gplots) 
cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(tissue)]
head(cbind(colnames(e),cols))
heatmap.2(e[idx,], labCol=tissue,
          trace="none", 
          ColSideColors=cols, 
          col=hmcol)

