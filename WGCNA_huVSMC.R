#Weighted Gene Coexpression Network Analysis

setwd("~/")

# Load Packages
library(Biobase)
library(GEOquery)
library(limma)
library(dplyr)
library(WGCNA)
library(sva)
library(colorspace)
library(qvalue)
library(cluster)
library(flashClust)

options(stringsAsFactors = FALSE);
allowWGCNAThreads()

#Read in Vascular Smooth Muscle Gene Expression File (Proliferative or Quiescent, log2 normalized)
annot = (read.table("Gene_Expression.txt",sep = "\t", header=T))

#Format according to below:
#Rows = Samples
#Row Names = Sample IDs
#Columns = ENS ID
#Column Headers = ENS ID
#Class = "numeric"
input = annot

#Keep the genes for which 80% of all samples show expression:

step1 = colSums(input !=0) #counts how many samples express the gene
step2 = which(step1 %in% c(0:120)) #creates a col vector of the gene indexes that don't have enough sample expression
testtable = input[,-c(step2)] #gets rid of the genes that don't have enough sample expression
input = testtable
#run the new input through the gsg code; there should be no additional change from gsg
gsg <- goodSamplesGenes(input, verbose=3)
dim0 = dim(input)
#If not,
# we remove the offending genes and samples from the data set
if(!gsg$allOK) {  
  if(sum(!gsg$goodGenes) > 0)
    printFlush(paste('Removing genes:', paste(names(input)[!gsg$goodGenes], collapse = ', ')))
  if(sum(!gsg$goodSamples) > 0)
    printFlush(paste('Removing samples:', paste(rownames(input)[!gsg$goodSamples], collapse = ', ')))
  # Remove offending genes and samples from the data
  input <- input[gsg$goodSamples, gsg$goodGenes]
}
if(dim0[1]!=nrow(input)){paste0("Null samples were removed")}


# Cluster the Null samples to ID outliers (Euclidian distance)
sampleTree = flashClust(dist(input), method = "average")
plot(sampleTree, main = "sample clustering to detect outliers", sub="", xlab="", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# apply a cut after reviewing this plot (where?)
cut = 100 #height of the line to cut out sample outliers, change according to Gene Expression file, no sample outliers were detected
# Plot a line to show the cut
pdf(file= "prePro_Null_sample_clustering_to_detect_outliers.pdf", width=24, height=9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "sample clustering to detect outliers", sub="", xlab="", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = cut, col = "red");
dev.off()
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cut, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
input = input[keepSamples, ]
length(which(keepSamples==F))

###############################################################
###### initial analysis: identify exponential power for sparsification

# range of powers for consideration
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# select power(s)
sft = pickSoftThreshold(input, powerVector = powers, verbose = 5)
pdf('wgcnaBasic_powerAssessment.pdf', width=10, height=5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# specify power for further network analysis
pow = 6
###############################################################
###### initial analysis: Run WGCNA on the datasets

# First work on the Null data for scale free topology
# Calculate the overall connectivity for each gene
# Plot a histogram of k and a scale free topology plot
k=softConnectivity(datE=input, corFnc="bicor", type="unsigned", power=pow, blockSize=25000, verbose=2)
pdf('wgcnaBasic_connectivity_histogram_and_scale_free_topology.pdf', width=10, height=5)
par(mfrow=c(1,2))
hist(k)
k_omit=na.omit(k)
scaleFreePlot(k_omit, main="Check scale free topology\n")
dev.off()


#If there are issues with function blockwise, try to restart the R session by going to Session and hitting restart
#then reuploading library(WGCNA); WGCNA may mask the cor function from the stats package which is necessary for function blockwise

#blockwise

bwnet = blockwiseModules(input, maxBlockSize = 2500,
                         power = 6, TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.1, numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "file_name_TOM-blockwise", verbose = 3)

bwLables = bwnet$colors
bwModuleColors = labels2colors(bwLables)


#Visualization of Modules
# open a graphics window
sizeGrWindow(6,6)
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 2
plotDendroAndColors(bwnet$dendrograms[[2]], bwModuleColors[bwnet$blockGenes[[2]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 2",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

# Plot the dendrogram and the module colors underneath for block 3
plotDendroAndColors(bwnet$dendrograms[[3]], bwModuleColors[bwnet$blockGenes[[3]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 3",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 4
plotDendroAndColors(bwnet$dendrograms[[4]], bwModuleColors[bwnet$blockGenes[[4]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 4",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 5
plotDendroAndColors(bwnet$dendrograms[[5]], bwModuleColors[bwnet$blockGenes[[5]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 5",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()


#Create Adjacency matrix and topological overlap
ADJ = adjacency(input, power=pow, corFnc='bicor')
dissTOM = TOMdist(ADJ, TOMType="unsigned")
collectGarbage()

#Call the hierarchical clustering function
hierTOM = hclust(as.dist(dissTOM),method="average")
pdf("wgcnaBasic_dendrogram_mphage.pdf",height=6,width=16)
par(mfrow=c(1,2))
plot(hierTOM,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity", labels=FALSE,hang=0.04);
dev.off()


###############################################################
###### principle components for eigengene visualization
###############################################################

modules = bwModuleColors
## Null model
# PCA, MDS
PCs = moduleEigengenes(input,  colors=modules)
ME = PCs$eigengenes
distPC = 1-abs(cor(ME,use="p"))
distPC = ifelse(is.na(distPC), 0, distPC)
pcTree = hclust(as.dist(distPC),method="a")
MDS  = cmdscale(as.dist(distPC),2)
colors = names(table(modules))

# PCA, MDS plots
pdf("ModuleEigengeneVisualizations_mphage.pdf",height=6,width=6)
par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)
plot(pcTree, xlab="",ylab="",main="",sub="")
plot(MDS, col= colors,  main="MDS plot", cex=2, pch=19)
ordergenes = hierTOM$order
plotMat(scale(input[,ordergenes]) , rlabels= modules[ordergenes], clabels= colnames(input), rcols=modules[ordergenes])
for (which.module in names(table(modules))){
  ME2 = ME[, paste("ME",which.module, sep="")]
  barplot(ME2, col=which.module, main="", cex.main=2,
          ylab="eigengene expression",xlab="array sample")
}
dev.off()

# Seperately plot the eigengene tree to see of there are any modules to be merged
# ColorDynamicHybrid, MEDissThres
# Plot the cut line into the dendogram
pdf('Module_Eigengene_Clustering_with_merge_line_mphage.pdf')
plot(pcTree, main = "Module_Eigengene_Clustering",xlab = "", sub = "")
MEDissThres=0.1
abline(h=MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
# Merge module colors
merged = mergeCloseModules(input, modules, cutHeight = MEDissThres, verbose = 3)
mergedModules = merged$colors


#########

# Eigengenes of the new merged modules:
# In order to obtain the variance explained by eigengene metric, calculate the eigengenes with the new colors
mergedMEs = merged$newMEs;
PCs_merged = moduleEigengenes(input,  colors=merged$colors)
mergedModules.varExplained=rbind(colnames(PCs_merged$eigengenes), PCs_merged$varExplained)
write.table(mergedModules.varExplained, file="Variance_Explained_by_Eigengene_Modules.txt", quote=F, row.names=F, col.names=F, sep="\t")

# See what the merging did to module colors
pdf("Module_Dendogram_before_and_after_module_merging_mphage.pdf", height=12, width=15)
plotDendroAndColors(hierTOM, cbind(modules, mergedModules), c("Before Merge", "After Merge"), dendroLabels = FALSE, hang = 0.05, addGuide = TRUE, guideHang = 0.05)
dev.off()

# Plot only the merged colors and the dendogram for Null
pdf("Module_Dendogram_after_module_merging_mphage.pdf", height=12, width=15)
plotDendroAndColors(hierTOM, mergedModules, "DynamicHybrid", dendroLabels = FALSE, hang = 0.05, addGuide = TRUE, guideHang = 0.05)
dev.off()


#######################################################################
#Calculate intramodular connectivity and write it to a file
#for colorDynamicHybrid
#######################################################################

colr = merged$colors
Alldegrees=intramodularConnectivity(ADJ, colr)
head(Alldegrees)
Alldegrees<-cbind(Alldegrees, colr)
head(Alldegrees)
Alldegrees$probesetID<-rownames(Alldegrees)

#Order them based on the color and within connectivity
o=order(Alldegrees$colr, -Alldegrees$kWithin)
Alldegrees.ordered = Alldegrees[o,]
write.table(Alldegrees.ordered, file="intramodular_connectivity_for_all_genes_mphage.txt", quote=FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

#Generalize intramodular connectivity for all genes on the array by calculating
#a correlation between a gene and each module eigengene
#for colorDynamicHybrids
datKMEcD=signedKME(input, ME, outputColumnName="MM.")
# Display the first few rows of the data frame
head(datKMEcD)
write.table(datKMEcD, file="intramodular_connectivity_for_all_genes_with_all_modules_colorDynamicHybridTOM_mphage.txt",quote=FALSE,
            sep = "\t",row.names = FALSE,col.names = TRUE)

#save workspace before working with phenotypes
#save.image("iterative_without_shared_uptoTRAITS.RData")

##########################################################################################
#Identify genes in each module, converts ENS to Gene Symbol

genenames = annot[,2]
ENSnames = annot[,1]
modGenesALL = c()
colortableALL = c()

for (color in colors)
{
  module = color
  probes = colnames(input)
  inModule = is.finite(match(colr, module));
  modProbes = as.matrix(probes[inModule])
  modGenes = as.matrix(genenames[match(modProbes, ENSnames)]);
  len = length(modGenes)
  colortable = matrix(color, len,1)
  colortableALL = c(colortableALL, colortable )
  modGenesALL = c(modGenesALL, modGenes)
  
}

geneswithcolor = rbind(colortableALL,modGenesALL)
geneswithcolor_trans = t(geneswithcolor)


#write file with module assignments
modassignments = t(geneswithcolor)
write.table(modassignments, file = "Module_Assignments.txt", row.names = FALSE, col.names = TRUE);

#Pathway Enrichment
library(AnnotationDbi)
library(GO.db)
library(org.Hs.eg.db)
library(anRichmentMethods)
library(anRichment)
library(topGO)

options(stringsAsFactors = FALSE);

#GO Terms
entrez = convert2entrez(organism = "human", symbol = modGenesALL); 
# How many conversions were successful? 
table(is.finite(entrez))


GOcollection = buildGOcollection(organism = "human") 

GOenrichment = enrichmentAnalysis(
  classLabels = colortableALL, identifiers = entrez, refCollection = GOcollection,
  useBackground = "given",
  threshold = 0.25,
  thresholdType = "FDR",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "grey");


collectGarbage();

names(GOenrichment)

names(GOenrichment$enrichmentTable)

table.display = GOenrichment$enrichmentTable;
table.display$overlapGenes = shortenStrings(table.display$overlapGenes, maxLength = 70, split = "|");
head(table.display);

write.csv(GOenrichment$enrichmentTable, file = "GOenrichment.csv", row.names = FALSE);

#KEGG Pathways
biosysCollection = BioSystemsCollection("human")
#msdbColl = MSigDBCollection(file = "download_file.jsp.xml", organism = "human")

entrez = convert2entrez(organism = "human", symbol = modGenesALL); 
# How many conversions were successful? 
table(is.finite(entrez))
keggEnrichment =  enrichmentAnalysis(classLabels = colortableALL, identifiers = entrez,
                                     refCollection = biosysCollection,
                                     useBackground = "given", threshold = 0.25, thresholdType = "FDR");

write.csv(keggEnrichment$enrichmentTable, file = "KEGGenrichment.csv", row.names = FALSE)

#Hallmark Pathways

#msdbColl = MSigDBCollection(file = "download_file.jsp.xml", organism = "human")

hallmarkEnrichment =  enrichmentAnalysis(classLabels = colortableALL, identifiers = entrez,
                                         refCollection = msdbColl,
                                         useBackground = "given", threshold = 0.25, thresholdType = "FDR");
write.csv(hallmarkEnrichment$enrichmentTable, file = "HALLMARKenrichment.csv", row.names = FALSE)


GOtermnames = GOenrichment$enrichmentTable$dataSetName
names(GOenrichment$dataSetDetails)
names(GOenrichment$dataSetDetails[[1]][[1]])
GOenrichment$dataSetDetails$black[[3]]$commonGeneEntrez
knownGroups(GOcollection)
GO.BPcollection = subsetCollection(GOcollection, tags = "GO.BP")
internalColl = internalCollection(organism = "human")
knownGroups(internalColl, sortBy = "size")

#######################


#Identify genes in CAD GWAS loci
GWASloci = read.table("Insert Gene List.txt", stringsAsFactors = FALSE)
GWASloci = data.frame(GWASloci)

matchedCADgenes = match(t(GWASloci),modGenesALL)
index = t(geneswithcolor[(1:2),matchedCADgenes])
CADgeneswithmodules = index[!is.na(index[,1]) & !is.na(index[,2]),]

#hub Genes
hubinput = input
ensvalues = colnames(input)
ensannotvalues = ENSnames
geneindex = match(ensvalues, ensannotvalues)
genesforhub = genenames[geneindex]
colnames(hubinput) = genesforhub
hubgenes = chooseTopHubInEachModule(hubinput, colr, omitColors = "grey")


#Fisher Exact Test
#overlap table
#Overlap counts and Fisher exact tests for two sets of module labels
Fishergenesall = colnames(hubinput)
FisherCADGenes = CADgeneswithmodules[,2]
########
FisherMatch = match(FisherCADGenes,Fishergenesall)
FisherCADexpressed = Fishergenesall
FisherCADexpressed[FisherMatch] = 1
FisherCADexpressed[-c(FisherMatch)] = 0
fisherindex = which(colr %in% "grey")
FisherColors = colr[-c(fisherindex)]
FisherCADexpressed = FisherCADexpressed[-c(fisherindex)]
FisherResults = overlapTable(FisherCADexpressed,FisherColors)
View(FisherResults$pTable)

#save workspace

save.image("WGCNA.RData")

