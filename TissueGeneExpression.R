library(tissuesGeneExpression)
data(tissuesGeneExpression)

library(RColorBrewer) 
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

library(genefilter)
rv <- rowVars(e)
idx <- order(-rv)[1:100]

library(gplots) 
cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(tissue)]
head(cbind(colnames(e),cols))
heatmap.2(e[idx,], labCol=tissue,
          trace="none", 
          ColSideColors=cols, 
          col=hmcol)

library(GSE5859Subset)
data(GSE5859Subset)
install.packages("matrixStats")
library(matrixStats)

library(gplots)

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
month = format( sampleInfo$date, "%m")
rv <- rowVars(geneExpression)
idx <- order(-rv)[1:25]
cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(as.character(sampleInfo$group))]

heatmap.2(geneExpression[idx,], 
          trace = 'none', labRow = geneAnnotation[idx,]$CHR,
          col = hmcol, labCol = month,
          ColSideColors = cols)

set.seed(17)
m = nrow(geneExpression)
n = ncol(geneExpression)
x = matrix(rnorm(m*n),m,n)
g = factor(sampleInfo$g )

pvals <- rowttests(x, g)$p.value
idx <- order(pvals)[1:50]
cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(as.character(sampleInfo$g))]
heatmap.2(x[idx,], 
          trace = 'none', labRow = geneAnnotation[idx,]$CHR,
          col = hmcol, labCol = month,
          ColSideColors = cols)

sds <- genefilter::rowSds(x,g)
idx <- order(-sds)[1:50]
cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(as.character(sampleInfo$g))]
heatmap.2(x[idx,], 
          trace = 'none', labRow = geneAnnotation[idx,]$CHR,
          col = hmcol, labCol = month,
          ColSideColors = cols)

