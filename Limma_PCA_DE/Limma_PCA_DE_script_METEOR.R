# Juliet M Wong
# Temperature and pCO2 effects on red sea urchin embryo gene expression
# RNA-seq
# PCA and Differential Expression analysis 

# R v 3.6.2
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install(version = "3.10")
#source("http://bioconductor.org/biocLite.R") # Older versions of R

# Download limma (Ritchie et al. 2015) and load required libraries
#biocLite("edgeR") # analysis performed in limma, but uses a few of the functions in edgeR
#biocLite("limma")
#biocLite("statmod")
#BiocManager::install(c("edgeR", "limma","statmod","ggplot2"))
library("edgeR")
library("limma")
library("statmod")
library("ggplot2")
#BiocManager::install('PCAtools')
library("PCAtools")
#install.packages("vegan")
library("vegan")
#BiocManager::install('EnhancedVolcano')
library("EnhancedVolcano")
library(readr)

#ALL SAMPLES----------------
# Read in counts data (.matrix files)
GeneCounts <- read.delim("results_gene.matrix", header=TRUE, row.names=1)

# All samples
All<-GeneCounts[ ,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)] 

# Filtering for more that 0.5 cpm across at least 3 samples
keep <- rowSums(cpm(All)>0.5) >=3
y<- All[keep,] 
# Check how many sequences remain after filtering
dim(y)
#37577    24

# Read in the model (contains stage, temperature treatment, pCO2 treatment, Biological replicate #s)
model<-read.delim("model.txt", header=TRUE)
# Create a treatment factor (combined treatment)
f_all <- paste(model$Stage,model$Temperature,model$pCO2,sep=".")
f_all <- factor(f_all)

# Create a model matrix of samples x treatment
design <-model.matrix(~0+f_all)
colnames(design) <- levels(f_all)
design

# Apply TMM normalization to RNA-seq read counts
# Create a DGEList object
dge <- DGEList(counts=y, group=f_all)
dge <- calcNormFactors(dge)
dge$samples

# Create MDS plot 
colors <- rep(c("pink","yellow","yellowgreen","lightblue","red","orange","green","blue"))
plotMDS(dge, labels=1:24,col=colors[f])
labelMDS <- rep(c("HH_Ga","HL_Ga","LH_Ga","LL_Ga","HH_Pr","HL_Pr","LH_Pr","LL_Pr"))
plotMDS(dge, labels=labelMDS[f_all],col=colors[f_all])

# Log transform the data for the PCA
logCPM <- cpm(dge, log=TRUE, prior.count = 3)

# Create PCA plot across log transformed data using PCAtools (https://github.com/kevinblighe/PCAtools)
# Read in the metadata for the PCA
metadata_all <- read_csv("Metadata_all.csv", 
                        col_types = cols(Temp = col_number(), 
                                         pCO2 = col_number(),Stage = col_character()))
View(metadata_all)
rownames(metadata_all) <- metadata_all$Sample

# Check the sample names match between the expression data and metadata
all(colnames(logCPM) == rownames(metadata_all))

# PCA for all stages
pca_all <- pca(logCPM, metadata = metadata_all)

# Biplot of first 2 PCs
all_colors <- c("black",rep("#D55E00",6),rep("#993399",6),rep("#009E73",6),rep("#0072B2",6),"black")
all_shapes <- rep(c(16, 17), times = 13)
pcaplot_all <- biplot(pca_all, pointSize = 4.0, labSize = 0)
pcaplot_all <- pcaplot_all + aes(shape = metadata_all$Stage, color = metadata_all$Stage)
pcaplot_all <- pcaplot_all + scale_colour_manual(values=all_colors)
pcaplot_all <- pcaplot_all + scale_shape_manual(values=all_shapes)
pcaplot_all <- pcaplot_all + theme_bw()
pcaplot_all <- pcaplot_all + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=14),legend.position = 'none',
                                   panel.background = element_blank(), axis.line = element_line(color = "black"), axis.text.y = element_text(angle = 90))
pcaplot_all

# Calculate the number of PCs that contribute >80% of the explained variation
which(cumsum(pca_all$variance) > 80) [1]
#7 PCs

# Correlate the top PCs back to the metadata (Temp, pCO2)
eigencorplot(pca_all, components = getComponents(pca_all, 1:7),
             col = c("mediumpurple3","mediumpurple2","mediumpurple1","white","goldenrod1","goldenrod2","goldenrod3"),
             metavars = c('pCO2','Temp'))
# Significant PC correlates exist for PC2, PC3, and PC4

# Write out loading results for PC2, PC3, and PC4
all_load <- pca_all$loadings[2:4]
write.table(all_load, "all_load.txt")


#multivariate analysis 
library(vegan)
colData <- data.frame(row.names = colnames(dge$counts), bucket=c("HH1_Ga","HH1_Pr","HH2_Ga","HH2_Pr","HH3_Ga","HH3_Pr","HL1_Ga","HL1_Pr","HL2_Ga","HL2_Pr","HL3_Ga","HL3_Pr","LH1_Ga","LH1_Pr","LH2_Ga","LH2_Pr","LH3_Ga","LH3_Pr","LL1_Ga","LL1_Pr","LL2_Ga","LL2_Pr","LL3_Ga","LL3_Pr"))
colData$treat = substr(colData$bucket, 1,2)
colData$treat_temp = substr(colData$bucket, 1,1)
colData$treat_pCO2 = substr(colData$bucket, 2,2)
colData$stage = substr(colData$bucket, 5,6)

conditions=colData
conditions$stage=as.factor(conditions$stage)
conditions$treat=as.factor(conditions$treat)
conditions$treat_temp=as.factor(conditions$treat_temp)
conditions$treat_pCO2=as.factor(conditions$treat_pCO2)

#analysis using euclidean distance method
adonis(t(na.omit(logCPM))~stage,data=conditions,method="eu")       
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#stage      1    246847  246847  43.543 0.66434  0.001 ***
#  Residuals 22    124720    5669         0.33566           
#Total     23    371567                 1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(t(na.omit(logCPM))~stage*treat_temp*treat_pCO2,data=conditions,method="eu")
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#stage                        1    246847  246847  52.853 0.66434  0.001 ***
#  treat_temp                   1     15010   15010   3.214 0.04040  0.041 *  
#  treat_pCO2                   1      8620    8620   1.846 0.02320  0.207    
#stage:treat_temp             1     11914   11914   2.551 0.03207  0.071 .  
#stage:treat_pCO2             1      5242    5242   1.122 0.01411  0.239    
#treat_temp:treat_pCO2        1      4617    4617   0.988 0.01242  0.278    
#stage:treat_temp:treat_pCO2  1      4590    4590   0.983 0.01235  0.272    
#Residuals                   16     74727    4670         0.20111           
#Total                       23    371567                 1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
#pie chart
slices <- c(66.4, 4.1, 2.3, 7.1, 20.1)
labs <- c("Stage","Temperature","pCO2","Interactions","Residuals")
labs <- paste(labs, slices) #add percents to labels
labs <- paste(labs,"%",sep="") #add % to labels
piecol <- c("grey46","grey70","grey10","grey28","white")
pie(slices,labels=labs, col=piecol,main="All")

# Apply limma pipelines for differential expression irrespective of stage
model
# Set treatment factors 
TempTreat <- factor(model$Temperature, levels=c("high","low"))
pCO2Treat <- factor(model$pCO2, levels=c("high","low"))

# Create a model matrix of samples x temperature treatment
designTemp <-model.matrix(~0+TempTreat)
colnames(designTemp) <- c("TempHigh","TempLow")
designTemp

fit <- lmFit(logCPM,designTemp) 
fit <- eBayes(fit, trend = TRUE) 
topTable(fit, coef=ncol(designTemp))

all_temp_data_output <- topTable(fit,coef=ncol (designTemp), adjust.method="BH", number=2000, p.value=0.05) 
write.table(all_temp_data_output, "all_temp_data_output.txt")

# Define a contrast matrix to compare high versus low temperature
cont.matrix <- makeContrasts(TempHigh-TempLow, levels=designTemp)
cont.matrix

# Extract the linear model fit for the temperature contrast
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
plotSA(fit2)

# Assess differential expression (high versus low temp)
colnames(fit2)
alltempcomp <- topTable(fit2,adjust="BH", number=Inf)

# High vs Low temp (all stages)
all.temp.highvslow <- topTable(fit2,coef=1,adjust="BH", number=Inf, p.value=0.05)
write.table(all.temp.highvslow, "all.temp.highvslow_sig.txt")

# write out for GO_MWU
all.temp.highvslow_all <- topTable(fit2,coef=1,adjust="BH", number=Inf)
write.table(all.temp.highvslow_all, "all.temp.highvslow_all.txt")

# Store temperature results and summarize
results <- decideTests(fit2)
summary(results)
#TempHigh - TempLow
#Down                  307
#NotSig              37011
#Up                    259

# Volcano plot for High vs. Low temp (all stages)
# Set custom colors for volcano
keyvals_alltemp <- ifelse(all.temp.highvslow_all$adj.P.Val > 0.05, 'grey30',
                          ifelse(all.temp.highvslow_all$logFC > 1, 'firebrick1',
                                 ifelse(all.temp.highvslow_all$logFC < 1 & all.temp.highvslow_all$logFC > 0,'lightpink2',
                                        ifelse(all.temp.highvslow_all$logFC > -1 & all.temp.highvslow_all$logFC < 0, 'lightblue3','dodgerblue2'))))
names(keyvals_alltemp) [keyvals_alltemp=='grey30'] <- 'NS'
names(keyvals_alltemp) [keyvals_alltemp=='firebrick1'] <- 'Up + logFC > 1'
names(keyvals_alltemp) [keyvals_alltemp=='lightpink2'] <- 'Up'
names(keyvals_alltemp) [keyvals_alltemp=='dodgerblue2'] <- 'Down + logFC < -1'
names(keyvals_alltemp) [keyvals_alltemp=='lightblue3'] <- 'Down'

# Plot volcano with labels
temp_all_volcano <- EnhancedVolcano(all.temp.highvslow_all, lab = rownames(all.temp.highvslow_all), x = 'logFC', y = 'adj.P.Val', 
                                    colCustom = keyvals_alltemp,pCutoff = 0.05, colAlpha = 0.8,
                                    labSize = 2.0, gridlines.major = FALSE, gridlines.minor = FALSE)
temp_all_volcano <- temp_all_volcano + ggplot2::coord_cartesian(xlim=c(-3,3),ylim=c(0,10))
temp_all_volcano

# Plot volcano without labels
temp_all_volcano <- EnhancedVolcano(all.temp.highvslow_all, lab = NA, x = 'logFC', y = 'adj.P.Val', 
                                    colCustom = keyvals_alltemp,pCutoff = 0.05, colAlpha = 0.8,
                                    labSize = 2.0, gridlines.major = FALSE, gridlines.minor = FALSE)
temp_all_volcano <- temp_all_volcano + ggplot2::coord_cartesian(xlim=c(-3,3),ylim=c(0,10))
temp_all_volcano


# Create a model matrix of samples x pCO2 treatment
designpCO2 <-model.matrix(~0+pCO2Treat)
colnames(designpCO2) <- c("pCO2High","pCO2Low")
designpCO2

fit <- lmFit(logCPM,designpCO2) 
fit <- eBayes(fit, trend = TRUE) 
topTable(fit, coef=ncol(designpCO2))

all_pCO2_data_output <- topTable(fit,coef=ncol (designpCO2), adjust.method="BH", number=2000, p.value=0.05) 
write.table(gastrula_pCO2_data_output, "all_pCO2_data_output.txt")

# Define a contrast matrix to compare pCO2 treatments (high vs low)
cont.matrix <- makeContrasts(pCO2High-pCO2Low, levels=designpCO2)
cont.matrix

# Extract the linear model fit for contrast of high vs low pCO2 
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
plotSA(fit2)

# Assess differential expression
colnames(fit2)
allpCO2comp <- topTable(fit2,adjust="BH", number=Inf)

# gastrula High vs Low pCO2
all.pCO2.highvslow <- topTable(fit2,coef=1,adjust="BH", number=Inf, p.value=0.05)
write.table(all.pCO2.highvslow, "all.pCO2.highvslow_sig.txt")

# write out for GO_MWU
all.pCO2.highvslow_all <- topTable(fit2,coef=1,adjust="BH", number=Inf)
write.table(all.pCO2.highvslow_all, "all.pCO2.highvslow_all.txt")

# Store results and summarize
results <- decideTests(fit2)
summary(results)
#pCO2High - pCO2Low
#Down                  119
#NotSig              37458
#Up                      0

# Volcano plot for High vs. Low pCO2 (all stages)
# Set custom colors for volcano
keyvals_allpCO2 <- ifelse(all.pCO2.highvslow_all$adj.P.Val > 0.05, 'grey30',
                          ifelse(all.pCO2.highvslow_all$logFC > 1, 'firebrick1',
                                 ifelse(all.pCO2.highvslow_all$logFC < 1 & all.pCO2.highvslow_all$logFC > 0,'lightpink2',
                                        ifelse(all.pCO2.highvslow_all$logFC > -1 & all.pCO2.highvslow_all$logFC < 0, 'lightblue3','dodgerblue2'))))
names(keyvals_allpCO2) [keyvals_allpCO2=='grey30'] <- 'NS'
names(keyvals_allpCO2) [keyvals_allpCO2=='firebrick1'] <- 'Up + logFC > 1'
names(keyvals_allpCO2) [keyvals_allpCO2=='lightpink2'] <- 'Up'
names(keyvals_allpCO2) [keyvals_allpCO2=='dodgerblue2'] <- 'Down + logFC < -1'
names(keyvals_allpCO2) [keyvals_allpCO2=='lightblue3'] <- 'Down'

# Plot volcano with labels
pCO2_all_volcano <- EnhancedVolcano(all.pCO2.highvslow_all, lab = rownames(all.pCO2.highvslow_all), x = 'logFC', y = 'adj.P.Val', 
                                    colCustom = keyvals_allpCO2,pCutoff = 0.05, colAlpha = 0.8,
                                    labSize = 2.0, gridlines.major = FALSE, gridlines.minor = FALSE)
pCO2_all_volcano <- pCO2_all_volcano + ggplot2::coord_cartesian(xlim=c(-3,3),ylim=c(0,10))
pCO2_all_volcano

# Plot volcano without labels
pCO2_all_volcano <- EnhancedVolcano(all.pCO2.highvslow_all, lab = NA, x = 'logFC', y = 'adj.P.Val', 
                                    colCustom = keyvals_allpCO2,pCutoff = 0.05, colAlpha = 0.8,
                                    labSize = 2.0, gridlines.major = FALSE, gridlines.minor = FALSE)
pCO2_all_volcano <- pCO2_all_volcano + ggplot2::coord_cartesian(xlim=c(-3,3),ylim=c(0,10))
pCO2_all_volcano

# To make a Grouped Bar Plot of count data
DEcounts <- structure(list(Temp= c(259, 307), pCO2 = c(0,119)),
                      .Names = c("Temperature (17 vs. 13C)", "pCO2 (1050 vs. 475uatm)"), class = "data.frame", row.names = c(NA,-2L))
barcolors <- c("firebrick1","dodgerblue2")
barplot(as.matrix(DEcounts), ylab = "Gene counts", ylim = c(0,400),cex.lab = 1.2, beside=TRUE, col=barcolors)
legend(7,3300, c("Down-regulated","Up-regulated"), cex=1, bty="n", fill=barcolors)

#GASTRULA ONLY---------------------------------

# Read in only gastrula samples only
gastrula<-GeneCounts[ ,c(1,3,5,7,9,11,13,15,17,19,21,23)] 

# Filtering for more that 0.5 cpm across at least 6 samples
keep <- rowSums(cpm(gastrula)>0.5) >=3
y<- gastrula[keep,] 
# Check how many sequences remain after filtering
dim(y)
# 33654    12

# Read in the model (contains stage, temperature treatment, pCO2 treatment, Biological replicate #s)
model<-read.delim("gastrula_model.txt", header=TRUE)
# Create a treatment factor (combined treatment)
f_gas <- paste(model$Stage,model$Temperature,model$pCO2,sep=".")
f_gas <- factor(f_gas)

# Create a model matrix of samples x treatment
design <-model.matrix(~0+f_gas)
colnames(design) <- levels(f_gas)
design

# Apply TMM normalization to RNA-seq read counts
# Create a DGEList object
dge <- DGEList(counts=y, group=f_gas)
dge <- calcNormFactors(dge)
dge$samples

# Create MDS plot 
colors <- rep(c("red","orange","green","blue"))
plotMDS(dge, labels=1:12,col=colors[f_gas])
par(xpd=T, mar=par()$mar+c(0,0,0,4))
mdsplot <- plotMDS(dge,pch=16,col=colors[f_gas])
plot(mdsplot, pch=16, col=colors[f_gas],xlab="Leading logFC dim 1", ylab="Leading logFC dim 2")
legend(2, 0 ,legend = c("HH_Ga","HL_Ga","LH_Ga","LL_Ga"), col=c("red","orange","green","blue"), pch=16, bty="n")
dev.off()

# Log transform the data for the PCA
logCPM <- cpm(dge, log=TRUE, prior.count = 3)

# Create PCA plot across log transformed data using PCAtools (https://github.com/kevinblighe/PCAtools)
# Read in the metadata for the PCA
metadata_ga <- read_csv("Metadata_ga.csv", 
                        col_types = cols(Temp = col_number(), 
                                         length = col_number(), pCO2 = col_number()))
View(metadata_ga)
rownames(metadata_ga) <- metadata_ga$Sample

# Check the sample names match between the expression data and metadata
all(colnames(logCPM) == rownames(metadata_ga))

# PCA for gastrula
pca_ga <- pca(logCPM, metadata = metadata_ga)

# Biplot of first 2 PCs
ga_colors <- c(rep("#D55E00",3),rep("#993399",3),rep("#009E73",3),rep("#0072B2",3))
pcaplot_gas <- biplot(pca_ga, pointSize = 4.0, labSize = 0)
pcaplot_gas <- pcaplot_gas + scale_colour_manual(values=ga_colors, name=NULL)
pcaplot_gas <- pcaplot_gas + theme_bw()
pcaplot_gas <- pcaplot_gas + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=14),legend.position="none",
                                   panel.background = element_blank(), axis.line = element_line(color = "black"), axis.text.y = element_text(angle = 90))
pcaplot_gas

# Calculate the number of PCs that contribute >80% of the explained variation
which(cumsum(pca_ga$variance) > 80) [1]
#8 PCs

# Correlate the top PCs back to the metadata (Temp, pCO2, other phys data)
eigen_ga <- eigencorplot(pca_ga, components = getComponents(pca_ga, 1:8),
             col = c("mediumpurple3","mediumpurple2","mediumpurple1","white","goldenrod1","goldenrod2","goldenrod3"),
             metavars = c('length','pCO2','Temp'))
eigen_ga
# Significant PC correlates exist for PC1 and PC2

# Write out loading results for PC1 and PC2
ga_load <- pca_ga$loadings[1:2]
write.table(ga_load, "ga_load.txt")

#multivariate analysis 
colData <- data.frame(row.names = colnames(dge$counts), bucket=c("HH1_Ga","HH2_Ga","HH3_Ga","HL1_Ga","HL2_Ga","HL3_Ga","LH1_Ga","LH2_Ga","LH3_Ga","LL1_Ga","LL2_Ga","LL3_Ga"))
colData$treat = substr(colData$bucket, 1,2)
colData$treat_temp = substr(colData$bucket, 1,1)
colData$treat_pCO2 = substr(colData$bucket, 2,2)

conditions=colData
conditions$treat=as.factor(conditions$treat)
conditions$treat_temp=as.factor(conditions$treat_temp)
conditions$treat_pCO2=as.factor(conditions$treat_pCO2)

#analysis using euclidean distance method
adonis(t(na.omit(logCPM))~treat_temp*treat_pCO2,data=conditions,method="eu")
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#treat_temp             1     10949 10948.7 2.72933 0.20258  0.001 ***
#  treat_pCO2             1      7164  7164.4 1.78596 0.13256  0.027 *  
#  treat_temp:treat_pCO2  1      3841  3841.4 0.95759 0.07108  0.440    
#Residuals              8     32092  4011.5         0.59379           
#Total                 11     54047                 1.00000           
---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#pie chart
slices <- c(20.3, 13.2, 7.1, 59.4)
labs <- c("Temperature","pCO2","Interaction","Residuals")
labs <- paste(labs, slices) #add percents to labels
labs <- paste(labs,"%",sep="") #add % to labels
piecol <- c("grey70","grey10","grey28","white")
pie(slices,labels=labs, col=piecol,main="Gastrula")


# Apply limma pipelines for differential expression at gastrula stage
# Set treatment factors 
TempTreat <- factor(model$Temperature, levels=c("high","low"))
pCO2Treat <- factor(model$pCO2, levels=c("high","low"))

# Create a model matrix of samples x temperature treatment
designTemp <-model.matrix(~0+TempTreat)
colnames(designTemp) <- c("TempHigh","TempLow")
designTemp

# Differential expression: limma-trend
fit <- lmFit(logCPM,designTemp) 
fit <- eBayes(fit, trend = TRUE) 
topTable(fit, coef=ncol(designTemp))

gastrula_temp_data_output <- topTable(fit,coef=ncol (designTemp), adjust.method="BH", number=2000, p.value=0.05) 
write.table(gastrula_temp_data_output, "gastrula_temp_data_output.txt")

# Define a contrast matrix to compare high versus low temperature
cont.matrix <- makeContrasts(TempHigh-TempLow, levels=designTemp)
cont.matrix

# Extract the linear model fit for the temperature contrast
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
plotSA(fit2)

# Assess differential expression (gastrula high versus low temp)
colnames(fit2)
gastrulatempcomp <- topTable(fit2,adjust="BH", number=Inf)

# gastrula High vs Low temp
ga.temp.highvslow <- topTable(fit2,coef=1,adjust="BH", number=Inf, p.value=0.05)
write.table(ga.temp.highvslow, "ga.temp.highvslow_sig.txt")

# write out for GO_MWU
ga.temp.highvslow_all <- topTable(fit2,coef=1,adjust="BH", number=Inf)
write.table(ga.temp.highvslow_all, "ga.temp.highvslow_all.txt")

# Store gastrula temperature results and summarize
results <- decideTests(fit2)
summary(results)
#TempHigh - TempLow
#Down                 1955
#NotSig              29650
#Up                   2049

# Volcano plot for High vs. Low temp (gastrula)
# Set custom colors for volcano
keyvals_gatemp <- ifelse(ga.temp.highvslow_all$adj.P.Val > 0.05, 'grey30',
                          ifelse(ga.temp.highvslow_all$logFC > 1, 'firebrick1',
                                 ifelse(ga.temp.highvslow_all$logFC < 1 & ga.temp.highvslow_all$logFC > 0,'lightpink2',
                                        ifelse(ga.temp.highvslow_all$logFC > -1 & ga.temp.highvslow_all$logFC < 0, 'lightblue3','dodgerblue2'))))
names(keyvals_gatemp) [keyvals_gatemp=='grey30'] <- 'NS'
names(keyvals_gatemp) [keyvals_gatemp=='firebrick1'] <- 'Up + logFC > 1'
names(keyvals_gatemp) [keyvals_gatemp=='lightpink2'] <- 'Up'
names(keyvals_gatemp) [keyvals_gatemp=='dodgerblue2'] <- 'Down + logFC < -1'
names(keyvals_gatemp) [keyvals_gatemp=='lightblue3'] <- 'Down'

# Plot volcano with labels
temp_ga_volcano <- EnhancedVolcano(ga.temp.highvslow_all, lab = rownames(ga.temp.highvslow_all), x = 'logFC', y = 'adj.P.Val', 
                                    colCustom = keyvals_gatemp,pCutoff = 0.05, colAlpha = 0.8,
                                    labSize = 2.0, gridlines.major = FALSE, gridlines.minor = FALSE)
temp_ga_volcano <- temp_ga_volcano + ggplot2::coord_cartesian(xlim=c(-5,5),ylim=c(0,8))
temp_ga_volcano

# Plot volcano without labels
temp_ga_volcano <- EnhancedVolcano(ga.temp.highvslow_all, lab = NA, x = 'logFC', y = 'adj.P.Val', 
                                    colCustom = keyvals_gatemp,pCutoff = 0.05, colAlpha = 0.8,
                                    axisLabSize = 10, legendIconSize = 1, legendLabSize = 10, gridlines.major = FALSE, gridlines.minor = FALSE)
temp_ga_volcano <- temp_ga_volcano + ggplot2::coord_cartesian(xlim=c(-5,5),ylim=c(0,8))
temp_ga_volcano


# Create a model matrix of samples x pCO2 treatment
designpCO2 <-model.matrix(~0+pCO2Treat)
colnames(designpCO2) <- c("pCO2High","pCO2Low")
designpCO2

fit <- lmFit(logCPM,designpCO2) 
fit <- eBayes(fit, trend = TRUE) 
topTable(fit, coef=ncol(designpCO2))
gastrula_pCO2_data_output <- topTable(fit,coef=ncol (designpCO2), adjust.method="BH", number=2000, p.value=0.05) 
write.table(gastrula_pCO2_data_output, "gastrula_pCO2_data_output.txt")

# Define a contrast matrix to compare pCO2 treatments (high vs low)
cont.matrix <- makeContrasts(pCO2High-pCO2Low, levels=designpCO2)
cont.matrix

# Extract the linear model fit for contrast of high vs low pCO2 
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
plotSA(fit2)

# Assess differential expression
colnames(fit2)
gastrulapCO2comp <- topTable(fit2,adjust="BH", number=Inf)

# gastrula High vs Low pCO2
ga.pCO2.highvslow <- topTable(fit2,coef=1,adjust="BH", number=Inf, p.value=0.05)
write.table(ga.pCO2.highvslow, "ga.pCO2.highvslow_sig.txt")

# write out for GO_MWU
ga.pCO2.highvslow_all <- topTable(fit2,coef=1,adjust="BH", number=Inf)
write.table(ga.pCO2.highvslow_all, "ga.pCO2.highvslow_all.txt")

# Store results and summarize
results <- decideTests(fit2)
summary(results)
#pCO2High - pCO2Low
#Down                  166
#NotSig              33479
#Up                      9

# Volcano plot for High vs. Low pCO2 (gastrula)
# Set custom colors for volcano
keyvals_gapCO2 <- ifelse(ga.pCO2.highvslow_all$adj.P.Val > 0.05, 'grey30',
                         ifelse(ga.pCO2.highvslow_all$logFC > 1, 'firebrick1',
                                ifelse(ga.pCO2.highvslow_all$logFC < 1 & ga.pCO2.highvslow_all$logFC > 0,'lightpink2',
                                       ifelse(ga.pCO2.highvslow_all$logFC > -1 & ga.pCO2.highvslow_all$logFC < 0, 'lightblue3','dodgerblue2'))))
names(keyvals_gapCO2) [keyvals_gapCO2=='grey30'] <- 'NS'
names(keyvals_gapCO2) [keyvals_gapCO2=='firebrick1'] <- 'Up + logFC > 1'
names(keyvals_gapCO2) [keyvals_gapCO2=='lightpink2'] <- 'Up'
names(keyvals_gapCO2) [keyvals_gapCO2=='dodgerblue2'] <- 'Down + logFC < -1'
names(keyvals_gapCO2) [keyvals_gapCO2=='lightblue3'] <- 'Down'

# Plot volcano with labels
pCO2_ga_volcano <- EnhancedVolcano(ga.pCO2.highvslow_all, lab = rownames(ga.pCO2.highvslow_all), x = 'logFC', y = 'adj.P.Val', 
                                   colCustom = keyvals_gapCO2,pCutoff = 0.05, colAlpha = 0.8,
                                   labSize = 2.0, gridlines.major = FALSE, gridlines.minor = FALSE)
pCO2_ga_volcano <- pCO2_ga_volcano + ggplot2::coord_cartesian(xlim=c(-5,5),ylim=c(0,8))
pCO2_ga_volcano

# Plot volcano without labels
pCO2_ga_volcano <- EnhancedVolcano(ga.pCO2.highvslow_all, lab = NA, x = 'logFC', y = 'adj.P.Val', 
                                   colCustom = keyvals_gapCO2,pCutoff = 0.05, colAlpha = 0.8,
                                   axisLabSize = 10, legendIconSize = 1, legendLabSize = 10, gridlines.major = FALSE, gridlines.minor = FALSE)
pCO2_ga_volcano <- pCO2_ga_volcano + ggplot2::coord_cartesian(xlim=c(-5,5),ylim=c(0,8))
pCO2_ga_volcano


# To make a Grouped Bar Plot of count data
DEcounts <- structure(list(Temp= c(2049, 1955), pCO2 = c(9,166)),
                      .Names = c("Temperature (17 vs. 13C)", "pCO2 (1050 vs. 475uatm)"), class = "data.frame", row.names = c(NA,-2L))
barcolors <- c("firebrick1","dodgerblue2")
barplot(as.matrix(DEcounts), ylab = "Gene counts", ylim = c(0,4000),cex.lab = 1.2, beside=TRUE, col=barcolors)
legend(7,3300, c("Down-regulated","Up-regulated"), cex=1, bty="n", fill=barcolors)


#PRISM ONLY---------------------------------

# prism samples only
prism <-GeneCounts[ ,c(2,4,6,8,10,12,14,16,18,20,22,24)] 

# Filtering for more that 0.5 cpm across at least 6 samples
keep <- rowSums(cpm(prism)>0.5) >=3
y<- prism[keep,] 
# Check how many sequences remain after filtering
dim(y)
#33471    12

# Read in the model (contains stage, temperature treatment, pCO2 treatment, Biological replicate #s)
model<-read.delim("prism_model.txt", header=TRUE)
# Create a treatment factor (combined treatment)
f_pri <- paste(model$Stage,model$Temperature,model$pCO2,sep=".")
f_pri <- factor(f_pri)

# Create a model matrix of samples x treatment
design <-model.matrix(~0+f_pri)
colnames(design) <- levels(f_pri)
design

# Apply TMM normalization to RNA-seq read counts
# Create a DGEList object
dge <- DGEList(counts=y, group=f_pri)
dge <- calcNormFactors(dge)
dge$samples

# Create MDS plot 
colors <- rep(c("red","orange","green","blue"))
plotMDS(dge, labels=1:12,col=colors[f_pri])

par(xpd=T, mar=par()$mar+c(0,0,0,4))
mdsplot <- plotMDS(dge,pch=16,col=colors[f_pri])
plot(mdsplot, pch=16, col=colors[f_pri],xlab="Leading logFC dim 1", ylab="Leading logFC dim 2")
legend(2, 0 ,legend = c("HH_Pr","HL_Pr","LH_Pr","LL_Pr"), col=c("red","orange","green","blue"), pch=16, bty="n")
dev.off()

# Log transform the data for the PCA
logCPM <- cpm(dge, log=TRUE, prior.count = 3)

# Create PCA plot across log transformed data using PCAtools (https://github.com/kevinblighe/PCAtools)
# Read in the metadata for the PCA
metadata_pr <- read_csv("Metadata_pr.csv", 
                        col_types = cols(LT50 = col_number(),Temp = col_number(), 
                                         length = col_number(), pCO2 = col_number(),Stage = col_character()))
View(metadata_pr)
rownames(metadata_pr) <- metadata_pr$Sample

# Check the sample names match between the expression data and metadata
all(colnames(logCPM) == rownames(metadata_pr))

# PCA for gastrula
pca_pr <- pca(logCPM, metadata = metadata_pr)

# Biplot of first 2 PCs
pr_colors <- c(rep("#D55E00",3),rep("#993399",3),rep("#009E73",3),rep("#0072B2",3))
pr_shapes <- rep(17,12)
pcaplot_pri <- biplot(pca_pr, pointSize = 4.0, labSize = 0)
pcaplot_pri <- pcaplot_pri + aes(shape = metadata_pr$Stage)
pcaplot_pri <- pcaplot_pri + scale_colour_manual(values=pr_colors)
pcaplot_pri <- pcaplot_pri + scale_shape_manual(values=pr_shapes)
pcaplot_pri <- pcaplot_pri + theme_bw()
pcaplot_pri <- pcaplot_pri + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=14),legend.position="none",
                                   panel.background = element_blank(), axis.line = element_line(color = "black"), axis.text.y = element_text(angle = 90))
pcaplot_pri

# Calculate the number of PCs that contribute >80% of the explained variation
which(cumsum(pca_ga$variance) > 80) [1]
#8 PCs

# Correlate the top PCs back to the metadata (Temp, pCO2, other phys data)
eigen_pr <- eigencorplot(pca_pr, components = getComponents(pca_pr, 1:8),
             col = c("mediumpurple3","mediumpurple2","mediumpurple1","white","goldenrod1","goldenrod2","goldenrod3"),
             metavars = c('LT50','length','pCO2','Temp'))
eigen_pr
# Significant PC correlates exist for PC1 and PC2

# Write out loading results for PC1 and PC2
pr_load <- pca_pr$loadings[1:2]
write.table(pr_load, "pr_load.txt")


#multivariate analysis
colData <- data.frame(row.names = colnames(dge$counts), bucket=c("HH1_Pr","HH2_Pr","HH3_Pr","HL1_Pr","HL2_Pr","HL3_Pr","LH1_Pr","LH2_Pr","LH3_Pr","LL1_Pr","LL2_Pr","LL3_Pr"))
colData$treat = substr(colData$bucket, 1,2)
colData$treat_temp = substr(colData$bucket, 1,1)
colData$treat_pCO2 = substr(colData$bucket, 2,2)

conditions=colData
conditions$treat=as.factor(conditions$treat)
conditions$treat_temp=as.factor(conditions$treat_temp)
conditions$treat_pCO2=as.factor(conditions$treat_pCO2)

#analysis using euclidean distance method
adonis(t(na.omit(logCPM))~treat_temp*treat_pCO2,data=conditions,method="eu")
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#treat_temp             1     13474 13473.9  3.8678 0.27196  0.001 ***
#  treat_pCO2             1      4614  4613.7  1.3244 0.09312  0.091 .  
#treat_temp:treat_pCO2  1      3587  3587.0  1.0297 0.07240  0.362    
#Residuals              8     27869  3483.6         0.56252           
#Total                 11     49544                 1.00000           
---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#pie chart
slices <- c(27.2, 9.3, 7.2, 56.3)
labs <- c("Temperature","pCO2","Interaction","Residuals")
labs <- paste(labs, slices) #add percents to labels
labs <- paste(labs,"%",sep="") #add % to labels
piecol <- c("grey70","grey10","grey28","white")
pie(slices,labels=labs, col=piecol,main="Prism")

# Differential expression: limma-trend
# Apply limma pipelines for differential expression at prism stage
# Set treatment factors 
TempTreat <- factor(model$Temperature, levels=c("high","low"))
pCO2Treat <- factor(model$pCO2, levels=c("high","low"))

# Create a model matrix of samples x temperature treatment
designTemp <-model.matrix(~0+TempTreat)
colnames(designTemp) <- c("TempHigh","TempLow")
designTemp

fit <- lmFit(logCPM,designTemp) 
fit <- eBayes(fit, trend = TRUE) 
topTable(fit, coef=ncol(designTemp))
prism_temp_data_output <- topTable(fit,coef=ncol (designTemp), adjust.method="BH", number=2000, p.value=0.05) 
write.table(prism_temp_data_output, "prism_temp_data_output.txt")

# Define a contrast matrix to compare high versus low temperature
cont.matrix <- makeContrasts(TempHigh-TempLow, levels=designTemp)
cont.matrix

# Extract the linear model fit for the temperature contrast
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
plotSA(fit2)

# Assess differential expression (prism high versus low temp)
colnames(fit2)
prismtempcomp <- topTable(fit2,adjust="BH", number=Inf)

# prism High vs Low temp
pr.temp.highvslow <- topTable(fit2,coef=1,adjust="BH", number=Inf, p.value=0.05)
write.table(pr.temp.highvslow, "pr.temp.highvslow_sig.txt")

# write out for GO_MWu
pr.temp.highvslow_all <- topTable(fit2,coef=1,adjust="BH", number=Inf)
write.table(pr.temp.highvslow_all, "pr.temp.highvslow_all.txt")

# Store gastrula temperature results and summarize
results <- decideTests(fit2)
summary(results)
#TempHigh - TempLow
#Down                 3434
#NotSig              26195
#Up                   3842

# Volcano plot for High vs. Low temp (prism)
# Set custom colors for volcano
keyvals_prtemp <- ifelse(pr.temp.highvslow_all$adj.P.Val > 0.05, 'grey30',
                         ifelse(pr.temp.highvslow_all$logFC > 1, 'firebrick1',
                                ifelse(pr.temp.highvslow_all$logFC < 1 & pr.temp.highvslow_all$logFC > 0,'lightpink2',
                                       ifelse(pr.temp.highvslow_all$logFC > -1 & pr.temp.highvslow_all$logFC < 0, 'lightblue3','dodgerblue2'))))
names(keyvals_prtemp) [keyvals_prtemp=='grey30'] <- 'NS'
names(keyvals_prtemp) [keyvals_prtemp=='firebrick1'] <- 'Up + logFC > 1'
names(keyvals_prtemp) [keyvals_prtemp=='lightpink2'] <- 'Up'
names(keyvals_prtemp) [keyvals_prtemp=='dodgerblue2'] <- 'Down + logFC < -1'
names(keyvals_prtemp) [keyvals_prtemp=='lightblue3'] <- 'Down'

# Plot volcano with labels
temp_pr_volcano <- EnhancedVolcano(pr.temp.highvslow_all, lab = rownames(pr.temp.highvslow_all), x = 'logFC', y = 'adj.P.Val', 
                                   colCustom = keyvals_prtemp,pCutoff = 0.05, colAlpha = 0.8,
                                   labSize = 2.0, gridlines.major = FALSE, gridlines.minor = FALSE)
temp_pr_volcano <- temp_pr_volcano + ggplot2::coord_cartesian(xlim=c(-5,5),ylim=c(0,8))
temp_pr_volcano

# Plot volcano without labels
temp_pr_volcano <- EnhancedVolcano(pr.temp.highvslow_all, lab = NA, x = 'logFC', y = 'adj.P.Val', 
                                   colCustom = keyvals_prtemp,pCutoff = 0.05, colAlpha = 0.8, 
                                   axisLabSize = 10, legendIconSize = 1, legendLabSize = 10, gridlines.major = FALSE, gridlines.minor = FALSE)
temp_pr_volcano <- temp_pr_volcano + ggplot2::coord_cartesian(xlim=c(-5,5),ylim=c(0,8))
temp_pr_volcano


# Create a model matrix of samples x pCO2 treatment
designpCO2 <-model.matrix(~0+pCO2Treat)
colnames(designpCO2) <- c("pCO2High","pCO2Low")
designpCO2

fit <- lmFit(logCPM,designpCO2) 
fit <- eBayes(fit, trend = TRUE) 
topTable(fit, coef = ncol(designpCO2))
prism_pCO2_data_output <- topTable(fit,coef=ncol (designpCO2), adjust.method="BH", number=2000, p.value=0.05) 
write.table(prism_pCO2_data_output, "prism_pCO2_data_output.txt")

# Define a contrast matrix to compare pCO2 treatments (high vs low)
cont.matrix <- makeContrasts(pCO2High-pCO2Low, levels=designpCO2)
cont.matrix

# Extract the linear model fit for contrast of high vs low pCO2 
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
plotSA(fit2)

# Assess differential expression
colnames(fit2)
prismpCO2comp <- topTable(fit2,adjust="BH", number=Inf)

# prism High vs Low pCO2
pr.pCO2.highvslow <- topTable(fit2,coef=1,adjust="BH", number=Inf, p.value=0.05)
write.table(pr.pCO2.highvslow, "pr.pCO2.highvslow_sig.txt")

# write out for GO_MWU
pr.pCO2.highvslow_all <- topTable(fit2,coef=1,adjust="BH", number=Inf)
write.table(pr.pCO2.highvslow_all, "pr.pCO2.highvslow_all.txt")

# Store results and summarize
results <- decideTests(fit2)
summary(results)
#pCO2High - pCO2Low
#Down                   64
#NotSig              33403
#Up                      4

# Volcano plot for High vs. Low pCO2 (prism)
# Set custom colors for volcano
keyvals_prpCO2 <- ifelse(pr.pCO2.highvslow_all$adj.P.Val > 0.05, 'grey30',
                         ifelse(pr.pCO2.highvslow_all$logFC > 1, 'firebrick1',
                                ifelse(pr.pCO2.highvslow_all$logFC < 1 & pr.pCO2.highvslow_all$logFC > 0,'lightpink2',
                                       ifelse(pr.pCO2.highvslow_all$logFC > -1 & pr.pCO2.highvslow_all$logFC < 0, 'lightblue3','dodgerblue2'))))
names(keyvals_prpCO2) [keyvals_prpCO2=='grey30'] <- 'NS'
names(keyvals_prpCO2) [keyvals_prpCO2=='firebrick1'] <- 'Up + logFC > 1'
names(keyvals_prpCO2) [keyvals_prpCO2=='lightpink2'] <- 'Up'
names(keyvals_prpCO2) [keyvals_prpCO2=='dodgerblue2'] <- 'Down + logFC < -1'
names(keyvals_prpCO2) [keyvals_prpCO2=='lightblue3'] <- 'Down'

# Plot volcano with labels
pCO2_pr_volcano <- EnhancedVolcano(pr.pCO2.highvslow_all, lab = rownames(pr.pCO2.highvslow_all), x = 'logFC', y = 'adj.P.Val', 
                                   colCustom = keyvals_prpCO2,pCutoff = 0.05, colAlpha = 0.8, borderWidth = 0.5,
                                   labSize = 2.0, gridlines.major = FALSE, gridlines.minor = FALSE)
pCO2_pr_volcano <- pCO2_pr_volcano + ggplot2::coord_cartesian(xlim=c(-5,5),ylim=c(0,8))
pCO2_pr_volcano

# Plot volcano without labels
pCO2_pr_volcano <- EnhancedVolcano(pr.pCO2.highvslow_all, lab = NA, x = 'logFC', y = 'adj.P.Val', 
                                   colCustom = keyvals_prpCO2,pCutoff = 0.05, colAlpha = 0.8,
                                   axisLabSize = 10, legendIconSize = 1, legendLabSize = 10, gridlines.major = FALSE, gridlines.minor = FALSE)
pCO2_pr_volcano <- pCO2_pr_volcano + ggplot2::coord_cartesian(xlim=c(-5,5),ylim=c(0,8))
pCO2_pr_volcano


# To make a Grouped Bar Plot of count data
DEcounts <- structure(list(Temp= c(3842,3434), pCO2 = c(4,64)),
                      .Names = c("Temperature (17 vs. 13C)", "pCO2 (1050 vs. 475uatm)"), class = "data.frame", row.names = c(NA,-2L))
barcolors <- c("firebrick1", "dodgerblue2")
barplot(as.matrix(DEcounts), ylab = "Gene counts", ylim = c(0,4000),cex.lab = 1.2, beside=TRUE, col=barcolors)
legend(7,3300, c("Down-regulated","Up-regulated"), cex=1, bty="n", fill=barcolors)


#FIGURES--------------------------
library("ggpubr")
ggarrange(pcaplot_all,pcaplot_all,pcaplot_gas,pcaplot_gas,pcaplot_pri,pcaplot_pri, labels = c("A","B","C","D","E","F"), ncol = 2, nrow =3, align = "hv")
ggarrange(eigen_ga,eigen_pr, labels = c("A","B"), ncol = 1, nrow =2, align = "hv")
ggarrange(temp_all_volcano,pCO2_all_volcano, labels = c("A","B"), ncol = 2, nrow =1, align = "hv")
ggarrange(temp_ga_volcano,pCO2_ga_volcano,temp_pr_volcano,pCO2_pr_volcano, labels = c("A","B","C","D"), ncol = 2, nrow =2, align = "hv")

