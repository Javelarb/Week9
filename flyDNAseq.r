library("DESeq2")

sampleInfo = read.table("RNAseq384_SampleCoding.txt", header = T)
sampleInfo$FullSampleName = as.character(sampleInfo$FullSampleName)
sampleInfo$PlateColumn = as.character(sampleInfo$PlateColumn)

## I am assuming feature counts finished
countdata = read.table("fly_counts.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
countdata = countdata[ ,6:ncol(countdata)]

# Remove crap from colnames in countdata
temp = colnames(countdata)
temp = gsub("RNASeq_out.","",temp)
temp = gsub("\\_.*","",temp)
colnames(countdata) = temp

##  does everything match up...
sampleInfo <- sampleInfo[sampleInfo$SampleNumber %in% colnames(countdata),]

# create DEseq2 object & run DEseq
dds = DESeqDataSetFromMatrix(countData=countdata, colData=sampleInfo, design=~PlateColumn)
dds <- DESeq(dds)
res <- results( dds )
res

plotMA( res, ylim = c(-1, 1) )
plotDispEsts( dds )
hist( res$pvalue, breaks=20, col="grey" )

###  throw out lowly expressed genes?? ... I leave this as an exercise
###  add external annotation to "gene ids"
# log transform
rld = rlog( dds )
## this is where you could just extract the log transformed normalized data for each sample ID, and then roll your own analysis
head( assay(rld) )
mydata = assay(rld)

sampleDists = dist( t( assay(rld) ) )
# heat map
sampleDistMatrix = as.matrix( sampleDists )
rownames(sampleDistMatrix) = rld$PlateColumn
colnames(sampleDistMatrix) = NULL

library( "gplots" )
library( "RColorBrewer" )
#dev.off()

#colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none")

# PCs
# wow you can sure tell tissue apart
print( plotPCA( rld, intgroup = "PlateColumn") )
# heat map with gene clustering
library( "genefilter" )
# these are the top genes (that tell tissue apart no doubt)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", trace="none", dendrogram="column", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
# volcano plot this is an exercise
