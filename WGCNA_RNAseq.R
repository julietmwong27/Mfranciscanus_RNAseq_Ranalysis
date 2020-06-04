# Juliet M Wong
# Weighted Gene Co-Expression Analysis (WGCNA) (Langfelder and Horvath 2008)
# Using limma (Ritchie et al. 2015) - voom-corrected data to determine clustering of similarly expressed genes and correlating them to treatments
# Temperature and pCO2 effects on red sea urchin embryo gene expression
# RNA-seq

# Download WGCNA
# Only run the following commands once to install WGCNA and flashClust on your computer
install.packages("BiocManager")
BiocManager::install("WGCNA")
install.packages("flashClust")

# Load WGCNA and flashClust libraries every time you open R
library(WGCNA)
library(flashClust)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# formatting data for WGCNA
# Take an object that contains the normalized counts file 
head(v2$E)
datExpr <- v2$E
head(datExpr)
dim(datExpr)

# Manipulate file so it matches the format WGCNA needs
datExpr = as.data.frame(t(datExpr))
dim(datExpr)
names(datExpr)


# Run this to check if there are gene outliers
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg#allOK
# The last statement returned TRUE, so all genes have passed the cuts

# Cluster the samples to see if there are any obvious outliers
sampleTree = hclust(dist(datExpr), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main =2)
# Plot a line to show the cut
abline(h = 15, col = "red")

# No cluster is an outlier (below the line)

#Create an object called "datTraits" that contains your trait data
datTraits = read.csv("Traits.csv")
head(datTraits)
row.names(datExpr) = datTraits$Sample
#form a data frame analogous to expression data that will hold the clinical traits.
rownames(datTraits) = datTraits$Sample
datTraits$Sample = NULL
table(rownames(datTraits)==rownames(datExpr)) #should return TRUE if datasets align correctly, otherwise your names are out of order
collectGarbage()

# now expression data is datExpr, and the clinical traits is datTraits
# visualize how the clinical traits relate to the sample dendrogram
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)
#Plot the sample dendrogram and the colors underneath
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")

save(datExpr, datTraits, file = "SamplesAndTraits.RData")


###-START Step-by-step module construction----------
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
#f you use RStudio or other environments not supplied by R Core team, please disable parallel execution
disableWGCNAThreads()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")
# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# The soft power value needs to be at least 7 (above where the scale-free topology fit index reaches 0.90)
# 16 was selected based off of recommendations by the package authors (see https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html)
# for sample size number for signed networks

#load data
load("SamplesAndTraits.RData");

#softpower selected from above
softPower = 16
allowWGCNAThreads()
adjacency = adjacency(datExpr, power = softPower, type="signed"); #specify network type

# Construct Networks: This may need to be run on a cluster (generated RData files can be transferred back to Rstudio)
# translate the adjacency into topological overlap matrix and caluclate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="signed")
dissTOM = 1-TOM
save(dissTOM, file= "dissTOM.RData")

# Generate Modules
# Generate a clustered gene tree; call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method="average")
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main="Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
#This sets the minimum number of genes to cluster into a module
minModuleSize = 30
# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Plot the module assignment under the gene dendrogram
# Convert numeric labels into colors
dynamicColors= labels2colors(dynamicMods)
table(dynamicColors)
#black        blue       brown       green greenyellow        grey 
#211       16164        1879         305         155          32 
#magenta        pink      purple         red         tan   turquoise 
#191         204         161         239         129       17316 
#yellow 
#591 

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHand = 0.05, main = "Gene dendrogram and module colors")
# Merging of modules whose expression profiles are very similar
# To quantify co-expression similarity of entire modules, we calculate their eigengenes and cluster them on their correlation
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors= dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method= "average")
# Plot the result; how the eigengenes cluster together
sizeGrWindow(7,6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_allSamples_signed_RLDfiltered.RData")

#load("Network_allSamples_signed_RLDfiltered.RData")
# set a threshold for merging modules
dynamicMergeCut(24) # calculates the threshold for module merging using the number of samples (in this case, 24)
MEDissThres = 0.2559876
# Plot the cut line into the dendrogram; adjust threshold if necessary
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight=MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
table(mergedColors)
#mergedColors
#black        blue       brown       green greenyellow        grey 
#372       16164        1879         305         155          32 
#magenta        pink         red   turquoise      yellow 
#191         204         368       17316         591 

# Eigengenes of the new merged modules
mergedMEs = merge$newMEs
# To see what the merging did to our module colors, we plat the gene dendrogram again, with both original and merged modules
sizeGrWindow(12,9)
# pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang=0.05)
dev.off()
# use the merged module colors in mergedColors; save relevant variables for next analysis
# Rename to moduleColors
moduleColors = mergedColors
table(mergedColors)

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method= "average")
# Plot the result; how the eigengenes cluster together
sizeGrWindow(7,6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
# Save module colors and labels for use in subsequent part
save(MEs, moduleLabels, moduleColors, geneTree, file= "NetworkConstruction-stepBystep_2.RData")

### END STEP-BY-STEP module construction

# Relate gene expression modules to traits
# Correlate traits
# Define number of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Print correlation heatmap between modules and traits
textMatrix = paste(signif(moduleTraitCor, 2), " (",
                   signif(moduleTraitPvalue, 1), ")", sep="")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

# Display the correlation values with a heatmap plot
labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(datTraits),
               yLabels= names(MEs),
               ySymbols= names(MEs),
               colorLabels= FALSE,
               colors= blueWhiteRed(50),
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 0.8,
               zlim= c(-1,1),
               main= paste("Module-Treatment relationships"))
dev.off()


#plotting massive table of all information for GO_MWU analysis
annot=read.table("SeqName_GO_annotations.txt",header=F)
names(datExpr)<-gsub("-tr","",names(datExpr))
names(datExpr)<-sub("transcript:","", names(datExpr))

probes = names(datExpr)
probes2annot = match(probes,annot$V1)
datGS.Traits=data.frame(cor(datExpr,datTraits,use="p"))
names(datGS.Traits)=paste("cor",names(datGS.Traits),sep=".")
datME=moduleEigengenes(datExpr,moduleColors)$eigengenes
datKME=signedKME(datExpr, datME, outputColumnName="kme.")
datOutput=data.frame(ProbeID=names(datExpr), annot[probes2annot,],moduleColors,datKME,datGS.Traits)
dim(datOutput) #37577    21

write.table(datOutput,"SeqNames_WGCNA_modulevalues.txt",row.names=F,sep="\t", col.names=T, quote=F)


