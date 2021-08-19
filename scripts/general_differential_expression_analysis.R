#####################################################################
############################ load required packages #################
#####################################################################
library(DESeq2)
library(RUVSeq)
library(ggplot2)
library(vsn)


##please take corresponding vignettes as reference for more details of the functions used

rm(list = ls())

#####################################################################
message("differential expression analysis related to miR-146a-5p data")
#####################################################################

###load required data
load("./data/datMeta_mir146.RData")


#datMeta contains data related to miR-146a-5p overexpression in 
#immortalized microglial (iMG) culture. 
#dataExpr contains raw microRNA counts in each sample
#condition contains information about the groups
#samples from both groups were equally distributed in different lanes for sequencing. 


### data preprocessing and differential expression analysis
zfGenes=dataExpr
nSamples=ncol(zfGenes)
filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=nSamples/2)#removing genes with low counts
filtered <- zfGenes[filter,]
genes <- rownames(filtered)
set <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(condition,  row.names=colnames(filtered)))

###differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = filtered, colData = pData(set), design=~condition)
dds <- DESeq(dds)
rld <- rlog(dds)
meanSdPlot(assay(rld)) #plot mean SD plot



#####################################################################
#PCA plot
#####################################################################
z=plotPCA(rld, intgroup = "condition", ntop = 500,  returnData = FALSE)
nudge <- position_nudge(y = 1)
print(z + theme_bw() )


### find the number of differentially expressed genes
res<- results(dds, contrast=c("condition","mimic","nc")) 
res <- res[order(res$padj), ]
normalized_counts=as.data.frame(counts(dds, normalized=TRUE))
diffexpr_mir146 <- merge(as.data.frame(res), normalized_counts, by="row.names", sort=FALSE)
names(diffexpr_mir146)[1] <- "Gene"

###number of signficantly deregulated genes
table(diffexpr_mir146$padj<0.05)[[2]]


#print(head(resdata))

#####################################################################
#volcano plot
#####################################################################

with(res, plot(log2FoldChange, -log10(padj), pch=20,  xlim=c( (min(res$log2FoldChange, na.rm=T)-1), (max(res$log2FoldChange, na.rm=T)+1) ), ylim=c( (min(-log10(res$padj), na.rm=T)-1), (max(-log10(res$padj), na.rm=T)+1) )))
with(subset(res, padj<0.05 & log2FoldChange >0), points(log2FoldChange, -log10(padj), pch=20, col="red"))
with(subset(res, padj<0.05 & log2FoldChange < 0), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
abline(h=-log10(0.05), lwd = 2, lty = 2)




