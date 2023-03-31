##############################################
#Preservation Analysis Smooth Muscle Cells
##############################################
library(WGCNA)
library(org.Hs.eg.db)
library(GO.db)
library(anRichmentMethods)



#run module creation method prior to using this script (WGCNA)
#Load in expression files and module assignment files

#load in Proliferative (withFBS) expression file
withExpRAW= read.table("Proliferative_expression.txt",sep = "\t", header=T)
withExp = withExpRAW[,-c(1,2)]
withExp = convertNumericColumnsToNumeric(withExp)
withnewrownames = withExpRAW[,2]
rownames(withExp) = withnewrownames
#load in Quiescent (withoutFBS) expression file
withoutExpRAW = read.table("Quiescent_expression.txt",sep = "\t", header=T)
withoutExp = withoutExpRAW[,-c(1,2)]
withoutExp = convertNumericColumnsToNumeric(withoutExp)
newrownames = withoutExpRAW[,2]
rownames(withoutExp) = newrownames

withoutExp = t(withoutExp)
withExp = t(withExp)
#load in module assignment with_FBS (Proliferative) file
withMOD = read.table("Module_Assignments_withFBS.txt", header = TRUE)
#load in module assignment without_FBS (Quiescent) file
withoutMOD_iter = read.table("Module_Assignments_withoutFBS.txt", header = TRUE)



#sort module column based off of gene name, gene order must be matching between the two datasets

withMOD = withMOD[match(withoutMOD_iter[,2],withMOD[,2]),]
withoutMOD = withoutMOD_iter


#Data formatting

#the control/healthy group is the data without FBS

mergedModulesControl=as.character(withoutMOD[,1]) #Quiescent, withoutFBS becomes the "Control"
mergedModulesTest=as.character(withMOD[,1]) #Proliferative, withFBS, becomes the "Test"

#####Preservation Analysis


#Create a contingency (or overlap) table of Control and Test modules 
#together with Fisher's exact test
#p-values for observing the observed overlaps by chance. 
#Use the overlapTable function that does the necessary calculations.
# Calculate the contingency table and p-values
overlap = overlapTable(mergedModulesControl, mergedModulesTest);
# The numMat will encode color. We use -log of the p value.
numMat = -log10(overlap$pTable);
numMat[numMat >50] = 50;
# Prepare for generating a color-coded plot of the overlap table. The text of the table will consist of
# counts and corresponding p-values.
textMat = paste(overlap$countTable, "\n", signif(overlap$pTable, 2));
dim(textMat) = dim(numMat)
# Additional information for the plot. These will be used shortly.
xLabels = paste("M", sort(unique(mergedModulesTest)));
yLabels = paste("M", sort(unique(mergedModulesControl)));
xSymbols = paste(sort(unique(mergedModulesTest)), ": ", table(mergedModulesTest), sep = "")
ySymbols = paste(sort(unique(mergedModulesControl)), ": ", table(mergedModulesControl), sep = "")
# Plot the overlap table
pdf(file = "Quiescent_Proliferative_Module_Membership_Overlap_Table.pdf", w = 14, h = 14); fp = TRUE
fcex = 0.55;
pcex = 0.55;
fcexl = 0.55;
pcexl = 0.55;
par(mar = c(10, 10, 2, 1.0));
labeledHeatmap(Matrix = numMat, xLabels = xLabels, xSymbols = xSymbols, yLabels = yLabels, ySymbols = ySymbols, colorLabels = TRUE,
               colors =  blueWhiteRed(100)[50:100], textMatrix = textMat, cex.text = if (fp) fcex else pcex, setStdMargins = FALSE, cex.lab = if (fp) fcexl else pcexl,
               xColorWidth = 0.08, main = "Quiescent modules (rows) vs. Proliferative modules (columns)", cex.main = 1.2);
dev.off()

save.image(file="OVERLABTABLE.Rdata")

#To quantify this result, we take advantage of the modulePreservation function.
#This function assesses how well a module in one study is preserved in another study using a number of strategies
#and outputs a single Z-score summary.
# Number of data sets that we work with
nSets = 2;
setLabels = c("Quiescent", "Proliferative");

multiExpr  = list(Control=list(data=withoutExp),Test=list(data=withExp))
multiColor = list(Control = mergedModulesControl, Test= mergedModulesTest)
lapply(multiExpr, lapply, dim)

mp=modulePreservation(multiExpr,multiColor,referenceNetworks=c(1,2),verbose=3,corFnc="bicor",networkType="signed",
                      savePermutedStatistics = TRUE, permutedStatisticsFile= 'QuiescentProliferativePreservation', calculateQvalue= TRUE,
                      nPermutations=500,maxGoldModuleSize=1000,maxModuleSize=1000)
stats = mp$preservation$Z$ref.Control$inColumnsAlsoPresentIn.OxPAPC
stats[order(-stats[,2]),c(1:2)]

save.image(file="permutations")


###################################
#Preservation Interpretation

#### load Adjacency matrices for both datasets one at a time
withoutADJ = ADJ #Quiescent WGCNA ADJ
withADJ = ADJ #Proliferative WGCNA ADJ

#Look for Control module preservation in Test data
#Isolate the observed statistics and their Z scores:
ref = 1 #Control
test = 2 #Test
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

#Look at the main output: the preservation medianRank and Zsummary statistics.
# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
withoutscores = cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
                      signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))
#Plot the statistics
# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
pdf(file="Quiescent_Module_Preservation_in_Proliferative_Zsummary_and_medianRank.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
dev.off();

#We now plot the density and connectivity statistics all in one plot. We include the module quality measures for
#comparison:
# Re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];
# Exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));
# Create numeric labels for each module
labs = match(modColors[plotMods], standardColors(50));
# Start the plot: open a suitably sized graphical window and set sectioning and margins.
sizeGrWindow(12, 9);
par(mfrow = c(3,5))
par(mar = c(3,3,2,1))
par(mgp = c(1.6, 0.4, 0));
# Plot each Z statistic in a separate plot.
pdf(file="Quiescent_Module_Preservation_in_Proliferative_All_Conservation_Statistics.pdf", wi=5, h=5)
for (s in 1:ncol(statsZ))
{
  min = min(statsZ[plotMods, s], na.rm = TRUE);
  max = max(statsZ[plotMods, s], na.rm = TRUE);
  if (min > -max/5) min = -max/5
  plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
       main = colnames(statsZ)[s],
       cex = 1.7,
       ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
       ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
       xlim = c(20, 1000))
  labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 0.7, offs = 0.04);
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}
dev.off()

#Look for Test module preservation in Control data
#Isolate the observed statistics and their Z scores:
ref = 2 #Test
test = 1 #Control
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

#Look at the main output: the preservation medianRank and Zsummary statistics.
# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
withscores = cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
                   signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))
#Plot the statistics
# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
pdf(file="Proliferative_Module_Preservation_in_Quiescent_Zsummary_and_medianRank.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
dev.off();

#We now plot the density and connectivity statistics all in one plot. We include the module quality measures for
#comparison:
# Re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];
# Exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));
# Create numeric labels for each module
labs = match(modColors[plotMods], standardColors(50));
# Start the plot: open a suitably sized graphical window and set sectioning and margins.
sizeGrWindow(12, 9);
par(mfrow = c(3,5))
par(mar = c(3,3,2,1))
par(mgp = c(1.6, 0.4, 0));
# Plot each Z statistic in a separate plot.
pdf(file="Proliferative_Module_Preservation_in_Quiescent_All_Conservation_Statistics.pdf", wi=5, h=5)
for (s in 1:ncol(statsZ))
{
  min = min(statsZ[plotMods, s], na.rm = TRUE);
  max = max(statsZ[plotMods, s], na.rm = TRUE);
  if (min > -max/5) min = -max/5
  plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
       main = colnames(statsZ)[s],
       cex = 1.7,
       ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
       ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
       xlim = c(20, 1000))
  labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 0.7, offs = 0.04);
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}
dev.off()

#We now plot the density and connectivity statistics all in one plot. We include the module quality measures for
#comparison:
# Re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];
# Exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));
# Create numeric labels for each module
labs = match(modColors[plotMods], standardColors(50));
# Start the plot: open a suitably sized graphical window and set sectioning and margins.
sizeGrWindow(12, 9);
par(mfrow = c(3,5))
par(mar = c(3,3,2,1))
par(mgp = c(1.6, 0.4, 0));
# Plot each Z statistic in a separate plot.
pdf(file="Proliferative_Module_Preservation_in_Quiescent_All_Conservation_Statistics.pdf", wi=5, h=5)
for (s in 1:ncol(statsZ))
{
  min = min(statsZ[plotMods, s], na.rm = TRUE);
  max = max(statsZ[plotMods, s], na.rm = TRUE);
  if (min > -max/5) min = -max/5
  plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
       main = colnames(statsZ)[s],
       cex = 1.7,
       ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
       ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
       xlim = c(20, 1000))
  labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 0.7, offs = 0.04);
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}
dev.off()

#Obtain the translation table between module colors and the numeric labels
data.frame(color = modColors[plotMods], label = labs)

# We now output the complete results
# This variable will contain the summary table
summaryTable = NULL
# Loop over all combinations of reference and tests sets
for (ref in 1:nSets) for (test in 1:nSets) if (ref!=test)
{
  modules = rownames(mp$preservation$Z[[ref]][[test]]);
  nMods = length(modules);
  sizes = mp$preservation$Z[[ref]][[test]][, 1];
  acc = matrix(NA, nMods, 3);
  if (test!=4)
  {
    acc[match(rownames(mp$accuracy$observed[[ref]][[test]]), modules), ] =
      mp$accuracy$observed[[ref]][[test]][, -1, drop = FALSE];
    colnames(acc) = colnames(mp$accuracy$observed[[ref]][[test]])[-1];
    accZ = mp$accuracy$Z[[ref]][[test]][, -1, drop = FALSE];
    acc.log.p = mp$accuracy$log.p[[ref]][[test]][, -1, drop = FALSE];
    acc.log.pBonf = mp$accuracy$log.pBonf[[ref]][[test]][, -1, drop = FALSE];
  } else {
    accZ = matrix(NA, nMods, 3);
    acc.log.p = matrix(NA, nMods, 3);
    acc.log.pBonf = matrix(NA, nMods, 3);
    colnames(acc) = colnames(mp$accuracy$observed[[1]][[2]])[-1];
    colnames(accZ) = colnames(mp$accuracy$Z[[1]][[2]])[-1];
    colnames(acc.log.p) = colnames(mp$accuracy$log.p[[1]][[2]])[-1];
    colnames(acc.log.pBonf) = colnames(mp$accuracy$log.pBonf[[1]][[2]])[-1];
  }
  # Table of results for this reference-test combination
  tab = cbind(referenceSet = rep(setLabels[ref], nMods),
              testSet = rep(setLabels[test], nMods),
              moduleLabel = modules,
              moduleSize = sizes,
              mp$quality$observed[[ref]][[test]][, -1, drop = FALSE],
              mp$preservation$observed[[ref]][[test]][, -1, drop = FALSE],
              acc,
              mp$referenceSeparability$observed[[ref]][[test]][, -1, drop = FALSE],
              mp$testSeparability$observed[[ref]][[test]][, -1, drop = FALSE],
              mp$quality$Z[[ref]][[test]][, -1, drop = FALSE],
              mp$quality$log.p[[ref]][[test]][, -1, drop = FALSE],
              mp$quality$log.pBonf[[ref]][[test]][, -1, drop = FALSE],
              mp$preservation$Z[[ref]][[test]][, -1, drop = FALSE],
              mp$preservation$log.p[[ref]][[test]][, -1, drop = FALSE],
              mp$preservation$log.pBonf[[ref]][[test]][, -1, drop = FALSE],
              accZ,
              acc.log.p,
              acc.log.pBonf,
              mp$referenceSeparability$Z[[ref]][[test]][, -1, drop = FALSE],
              mp$referenceSeparability$log.p[[ref]][[test]][, -1, drop = FALSE],
              mp$referenceSeparability$log.pBonf[[ref]][[test]][, -1, drop = FALSE],
              mp$testSeparability$Z[[ref]][[test]][, -1, drop = FALSE],
              mp$testSeparability$log.p[[ref]][[test]][, -1, drop = FALSE],
              mp$testSeparability$log.pBonf[[ref]][[test]][, -1, drop = FALSE]
  )
  # Add the table to the main table.
  if (is.null(summaryTable)) summaryTable = tab else summaryTable = rbind(summaryTable, tab);
}
# Save the table in txt format.
write.table(summaryTable, file = "Quiescent_Proliferative_ModulePreservation_completeResults.txt", row.names = FALSE, sep = "\t", quote = FALSE);

save.image("Permutations_all.Rdata")