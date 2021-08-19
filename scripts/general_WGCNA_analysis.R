

#this script runs WGCNA analysis and find modules in aging mice model treated with 3-miR inhibitor
#and scramble control. Same analysis is also used to analyze the data from APP/PS1 mice only with a change in softpower as described in Methods and Materials. 

options(stringsAsFactors  =  FALSE)
#load packages
library(WGCNA)
library(flashClust)
library(AnnotationDbi)
library("Mus.musculus")
library("org.Mm.eg.db")
library("clusterProfiler")


load("./data/datMeta_agingWGCNA.Rdata") #load data #datNorm contains normalized counts, detail contains phenotypic information


datExpr = as.data.frame(t(datNorm)) #transpose data
datExpr[datExpr == 0] <- 1
datExpr=na.omit(datExpr)
datTraits = detail
names(datTraits) = "condition_short"
datTraits$group = c(rep("aging_young_nc", 6), rep("aging_old_nc", 6), rep("aging_old_inhibitor", 7))
datTraits$group_numeric = as.numeric(as.factor(datTraits$group))
datatrait_numeric=sapply(datTraits, function (x) as.numeric(x))
datatrait_numeric=datatrait_numeric[,-c(1,2)]
rownames(datExpr) = rownames(datTraits)
table(rownames(datTraits)==rownames(datExpr)) 


A = adjacency(t(datExpr),type="signed") 
k = as.numeric(apply(A,2,sum))-1 
Z.k = scale(k)
thresholdZ.k = -2.5 
outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
traitColors = data.frame(numbers2colors(datatrait_numeric,signed=FALSE))
datColors = data.frame(outlier = outlierColor,traitColors)

plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample Dendrogram Heatmap")



allowWGCNAThreads() 


# Pick soft threshold
powers = c(c(1:10), seq(from =10, to=30, by=1)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(datExpr,networkType = "signed", corFnc = "bicor",verbose = 5,powerVector = powers)
sizeGrWindow(9,5)
par(mfrow= c(1,2))
cex1=1.0
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")


plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

softPower = 17  ## Choose based on fit to scale-free topology at R2 > 0.8
adjacency = adjacency(datExpr, corFnc = "bicor", type = "signed", power = softPower)
TOM = TOMsimilarity(adjacency,TOMType = "signed", verbose = 0)
dissTOM = 1-TOM
geneTree = flashClust(as.dist(dissTOM), method = "average")


###Run WGCNA with selected parameters

minModuleSize =100; ds = 2;cutHeight = 0.99999
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method = "hybrid",
                            deepSplit = ds, pamRespectsDendro = T,pamStage = T,
                            minClusterSize = minModuleSize, cutHeight = cutHeight)
dynamicColors = labels2colors(dynamicMods)


#calculate eigevalue
MEList = moduleEigengenes(datExpr, colors = dynamicColors,softPower = softPower)
MEs = MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")

#merge similar clusters
MEDissThres = 0.15
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
moduleColors = merge$colors
MEs = merge$newMEs


# Calculate module membership for each module

modNames=substring(names(MEs),3)
geneModulecor=corAndPvalue(datExpr, MEs, use = "p")
geneModuleMembership = as.data.frame(geneModulecor$cor)
MMPvalue = as.data.frame(geneModulecor$p)
names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")

geneInfo = data.frame(geneID = colnames(datExpr), moduleColor = moduleColors, geneModuleMembership, MMPvalue)


##preparing files to plot 7B
group_mice = data.frame(c(rep("3m_nc", 6), rep("16m_nc", 6), rep("16m_in",7)))
data_eigengene = data.frame(cbind(datTraits$group,group_mice[,1],  MEs))
names(data_eigengene)[2] = "plot_pheno.group"

comparison = list(c("3m_nc","16m_nc"),   c("3m_nc", "16m_in"))


#plot 7B

ggboxplot(data_eigengene, x = "plot_pheno.group", y = "MEblue",
ylab = "eigenvalue",
          title = "MEblue") + stat_compare_means(comparisons = comparison) 


