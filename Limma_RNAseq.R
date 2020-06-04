# Juliet M Wong
# limma (Ritchie et al. 2015)
# Temperature and pCO2 effects on red sea urchin embryo gene expression
# RNA-seq

# Done using R v 3.6.2
# Download limma and load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")
#source("http://bioconductor.org/biocLite.R") # Older versions of R
#biocLite("edgeR") # analysis performed in limma, but uses a few of the functions in edgeR
#biocLite("limma")
#biocLite("statmod")
BiocManager::install(c("edgeR", "limma","statmod","ggplot2"))
library("edgeR")
library("limma")
library("statmod")
library("ggplot2")
install.packages("vegan")
library("vegan")

# Read in counts data (.matrix files)
GeneCounts <- read.delim("results_gene.matrix", header=TRUE, row.names=1)

#ALL SAMPLES-------------------------------
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
f <- paste(model$Stage,model$Temperature,model$pCO2,sep=".")
f <- factor(f)

# Create a model matrix of samples x treatment
design <-model.matrix(~0+f)
colnames(design) <- levels(f)
design

# Apply TMM normalization to RNA-seq read counts
# Create a DGEList object
dge <- DGEList(counts=y, group=f)
dge <- calcNormFactors(dge)
dge$samples
# Create MDS plot (not yet voom-transformed)
colors <- rep(c("pink","yellow","yellowgreen","lightblue","red","orange","green","blue"))
plotMDS(dge, labels=1:24,col=colors[f])

# voom transformation using voomWithQualityWeights function
# combine observational-level weighting from voom with sample-specific quality weights
v1 <- voomWithQualityWeights(dge, design=design, lib.size=dge$samples$lib.size, plot = TRUE)
labelMDS <- rep(c("HH_Ga","HL_Ga","LH_Ga","LL_Ga","HH_Pr","HL_Pr","LH_Pr","LL_Pr"))
plotMDS(v1, labels=labelMDS[f],col=colors[f])

# Run voom again, this time adding a blocking variable and estimating correlation
# Samples are blocked by biological replicate
corfit <- duplicateCorrelation(v1,design,block=model$BioRep) 
corfit$consensus
v2 <- voomWithQualityWeights(dge,design,plot=TRUE,lib.size=dge$samples$lib.size,block=model$BioRep,correlation=corfit$consensus)

# MDS plot of voom-corrected data
par(xpd=T, mar=par()$mar+c(0,0,0,4))
mdsplot <- plotMDS(v2,pch=16,col=colors[f])
plot(mdsplot, pch=16, col=colors[f],xlab="Leading logFC dim 1", ylab="Leading logFC dim 2")
legend(2, 0 ,legend = c("HH_Ga","HL_Ga","LH_Ga","LL_Ga","HH_Pr","HL_Pr","LH_Pr","LL_Pr"), col=c("pink","yellow","yellowgreen","lightblue","red","orange","green","blue"), pch=16, bty="n")
dev.off()

# PCA across voom-corrected data
pcavoom<- prcomp(t(na.omit(v2$E)))
summary(pcavoom)
pca_all <- as.data.frame(prcomp(t(na.omit(v2$E)))$x)
pca_all$f<-f
pcaplot_all <- qplot(x=PC1, y=PC2, data=pca_all, colour=pca_all$f, shape=pca_all$f,size=I(4))
pcaplot_all <- pcaplot_all + scale_colour_manual(values=c("#D55E00", "#993399", "#009E73","#0072B2","#D55E00", "#993399", "#009E73","#0072B2"), name=NULL)
pcaplot_all <- pcaplot_all + scale_shape_manual(values=c(16,16,16,16,17,17,17,17))
pcaplot_all <- pcaplot_all + theme_bw()
pcaplot_all <- pcaplot_all + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=14),legend.text=element_text(size=14),
                panel.background = element_blank(), axis.line = element_line(color = "black"), axis.text.y = element_text(angle = 90), legend.key = element_blank())
pcaplot_all

#multivariate analysis 
colData <- data.frame(row.names = colnames(v2$E), bucket=c("HH1_Ga","HH1_Pr","HH2_Ga","HH2_Pr","HH3_Ga","HH3_Pr","HL1_Ga","HL1_Pr","HL2_Ga","HL2_Pr","HL3_Ga","HL3_Pr","LH1_Ga","LH1_Pr","LH2_Ga","LH2_Pr","LH3_Ga","LH3_Pr","LL1_Ga","LL1_Pr","LL2_Ga","LL2_Pr","LL3_Ga","LL3_Pr"))
colData$treat = substr(colData$bucket, 1,2)
colData$treat_temp = substr(colData$bucket, 1,1)
colData$treat_pCO2 = substr(colData$bucket, 2,2)
colData$stage = substr(colData$bucket, 5,6)

conditions=colData
conditions$stage=as.factor(conditions$stage)
conditions$treat=as.factor(conditions$treat)
conditions$treat_temp=as.factor(conditions$treat_temp)
conditions$treat_pCO2=as.factor(conditions$treat_pCO2)

adonis(t(na.omit(v2$E))~stage*treat_temp*treat_pCO2,data=conditions,method="eu")
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#stage                        1    383414  383414  34.650 0.57785  0.001 ***
#  treat_temp                   1     27437   27437   2.480 0.04135  0.051 .  
#treat_pCO2                   1     18525   18525   1.674 0.02792  0.184    
#stage:treat_temp             1     23311   23311   2.107 0.03513  0.085 .  
#stage:treat_pCO2             1     11990   11990   1.084 0.01807  0.225    
#treat_temp:treat_pCO2        1     10963   10963   0.991 0.01652  0.282    
#stage:treat_temp:treat_pCO2  1     10833   10833   0.979 0.01633  0.272    
#Residuals                   16    177047   11065         0.26683           
#Total                       23    663519                 1.00000           
---
  #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
#pie chart
slices <-c(57.8, 4.1, 2.8, 8.6, 26.7)
labs <- c("Stage","Temperature","pCO2","Interactions","Residuals")
labs <- paste(labs, slices) #add percents to labels
labs <- paste(labs,"%",sep="") #add % to labels
piecol <- c("grey46","grey70","grey10","grey28","white")
pie(slices,labels=labs, col=piecol,main="All")



#GASTRULA ONLY---------------------------------

# gastrula samples only
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
f <- paste(model$Stage,model$Temperature,model$pCO2,sep=".")
f <- factor(f)

# Create a model matrix of samples x treatment
design <-model.matrix(~0+f)
colnames(design) <- levels(f)
design

# Apply TMM normalization to RNA-seq read counts
# Create a DGEList object
dge <- DGEList(counts=y, group=f)
dge <- calcNormFactors(dge)
dge$samples
# Create MDS plot (not yet voom-transformed)
colors <- rep(c("red","orange","green","blue"))
plotMDS(dge, labels=1:12,col=colors[f])

# voom transformation using voomWithQualityWeights function
# combine observational-level weighting from voom with sample-specific quality weights
v1 <- voomWithQualityWeights(dge, design=design, lib.size=dge$samples$lib.size, plot = TRUE)
labelMDS <- rep(c("HH_Ga","HL_Ga","LH_Ga","LL_Ga"))
plotMDS(v1, labels=labelMDS[f],col=colors[f])

# Run voom again, this time adding a blocking variable and estimating correlation
# Samples are blocked by biological replicate
corfit <- duplicateCorrelation(v1,design,block=model$BioRep) 
corfit$consensus
v2 <- voomWithQualityWeights(dge,design,plot=TRUE,lib.size=dge$samples$lib.size,block=model$BioRep,correlation=corfit$consensus)

# MDS plot of voom-corrected data
par(xpd=T, mar=par()$mar+c(0,0,0,4))
mdsplot <- plotMDS(v2,pch=16,col=colors[f])
plot(mdsplot, pch=16, col=colors[f],xlab="Leading logFC dim 1", ylab="Leading logFC dim 2")
legend(2, 0 ,legend = c("HH_Ga","HL_Ga","LH_Ga","LL_Ga"), col=c("pink","yellow","yellowgreen","lightblue"), pch=16, bty="n")
dev.off()

# PCA across voom-corrected data
pcavoom<- prcomp(t(na.omit(v2$E)))
summary(pcavoom)
pca_gas <- as.data.frame(prcomp(t(na.omit(v2$E)))$x)
pca_gas$f<-f
pcaplot_gas <- qplot(x=PC1, y=PC2, data=pca_gas,colour=pca_gas$f, size=I(4))
pcaplot_gas <- pcaplot_gas + scale_colour_manual(values=c("#D55E00", "#993399", "#009E73","#0072B2"), name=NULL)
pcaplot_gas <- pcaplot_gas + theme_bw()
pcaplot_gas <- pcaplot_gas + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=14),legend.text=element_text(size=14),
                panel.background = element_blank(), axis.line = element_line(color = "black"), axis.text.y = element_text(angle = 90), legend.key = element_blank())

#multivariate analysis
colData <- data.frame(row.names = colnames(v2$E), bucket=c("HH1_Ga","HH2_Ga","HH3_Ga","HL1_Ga","HL2_Ga","HL3_Ga","LH1_Ga","LH2_Ga","LH3_Ga","LL1_Ga","LL2_Ga","LL3_Ga"))
colData$treat = substr(colData$bucket, 1,2)
colData$treat_temp = substr(colData$bucket, 1,1)
colData$treat_pCO2 = substr(colData$bucket, 2,2)

conditions=colData
conditions$treat=as.factor(conditions$treat)
conditions$treat_temp=as.factor(conditions$treat_temp)
conditions$treat_pCO2=as.factor(conditions$treat_pCO2)

#analysis using euclidean distance method
adonis(t(na.omit(v2$E))~treat_temp*treat_pCO2,data=conditions,method="eu")
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#treat_temp             1     19704 19703.5 2.52562 0.19216  0.001 ***
#  treat_pCO2             1     12837 12837.3 1.64550 0.12519  0.034 *  
#  treat_temp:treat_pCO2  1      7586  7586.4 0.97243 0.07399  0.419    
#Residuals              8     62412  7801.5         0.60866           
#Total                 11    102539                 1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#pie chart
slices <-c(19.2, 12.5, 7.4, 60.9)
labs <- c("Temperature","pCO2","Interaction","Residuals")
labs <- paste(labs, slices) #add percents to labels
labs <- paste(labs,"%",sep="") #add % to labels
piecol <- c("grey70","grey10","grey28","white")
pie(slices,labels=labs, col=piecol,main="Gastrula")


# Apply limma pipelines for differential expression at gastrula stage
model<-read.delim("gastrula_model.txt", header=TRUE)
# Set treatment factors 
TempTreat <- factor(model$Temperature, levels=c("high","low"))
pCO2Treat <- factor(model$pCO2, levels=c("high","low"))

# Create a model matrix of samples x temperature treatment
designTemp <-model.matrix(~0+TempTreat)
colnames(designTemp) <- c("TempHigh","TempLow")
designTemp

fit <- lmFit(v2,designTemp,correlation=corfit$consensus) 
fit <- eBayes(fit) 
topTable(fit)
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
write.table(ga.temp.highvslow, "ga.temp.highvslow.txt")

# write out for GO_MWU
ga.temp.highvslow_all <- topTable(fit2,coef=1,adjust="BH", number=Inf)
write.table(ga.temp.highvslow_all, "ga.temp.highvslow_all.txt")

# Store gastrula temperature results and summarize
results <- decideTests(fit2)
summary(results)
#TempHigh - TempLow
#Down                 1735
#NotSig              27283
#Up                   4636

# Create a model matrix of samples x pCO2 treatment
designpCO2 <-model.matrix(~0+pCO2Treat)
colnames(designpCO2) <- c("pCO2High","pCO2Low")
designpCO2

fit <- lmFit(v2,designpCO2,correlation=corfit$consensus) 
fit <- eBayes(fit) 
topTable(fit)
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
ga.pCO2.highvslow_all <- topTable(fit2,coef=1,adjust="BH", number=Inf)

# gastrula High vs Low pCO2
ga.pCO2.highvslow <- topTable(fit2,coef=1,adjust="BH", number=Inf, p.value=0.05)
write.table(ga.pCO2.highvslow, "ga.pCO2.highvslow.txt")

# write out for GO_MWU
ga.pCO2.highvslow_all <- topTable(fit2,coef=1,adjust="BH", number=Inf)
write.table(ga.pCO2.highvslow_all, "ga.pCO2.highvslow_all.txt")

# Store results and summarize
results <- decideTests(fit2)
summary(results)
#pCO2High - pCO2Low
#Down                  202
#NotSig              33403
#Up                     49

# To make a Grouped Bar Plot of count data
DEcounts <- structure(list(Temp= c(4636,1735), pCO2 = c(49,202)),
                      .Names = c("Temperature (17 vs. 13C)", "pCO2 (1050 vs. 475uatm)"), class = "data.frame", row.names = c(NA,-2L))
attach(data)
print(data)
barcolors <- c("firebrick1", "dodgerblue2")
barplot(as.matrix(DEcounts), ylab = "Gene counts", ylim = c(0,5000),cex.lab = 1.2, beside=TRUE, col=barcolors)
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
f <- paste(model$Stage,model$Temperature,model$pCO2,sep=".")
f <- factor(f)

# Create a model matrix of samples x treatment
design <-model.matrix(~0+f)
colnames(design) <- levels(f)
design

# Apply TMM normalization to RNA-seq read counts
# Create a DGEList object
dge <- DGEList(counts=y, group=f)
dge <- calcNormFactors(dge)
dge$samples
# Create MDS plot (not yet voom-transformed)
colors <- rep(c("red","orange","green","blue"))
plotMDS(dge, labels=1:12,col=colors[f])

# voom transformation using voomWithQualityWeights function
# combine observational-level weighting from voom with sample-specific quality weights
v1 <- voomWithQualityWeights(dge, design=design, lib.size=dge$samples$lib.size, plot = TRUE)
labelMDS <- rep(c("HH_Pr","HL_Pr","LH_Pr","LL_Pr"))
plotMDS(v1, labels=labelMDS[f],col=colors[f])

# Run voom again, this time adding a blocking variable and estimating correlation
# Samples are blocked by biological replicate
corfit <- duplicateCorrelation(v1,design,block=model$BioRep) 
corfit$consensus
v2 <- voomWithQualityWeights(dge,design,plot=TRUE,lib.size=dge$samples$lib.size,block=model$BioRep,correlation=corfit$consensus)

# MDS plot of voom-corrected data
par(xpd=T, mar=par()$mar+c(0,0,0,4))
mdsplot <- plotMDS(v2,pch=16,col=colors[f])
plot(mdsplot, pch=16, col=colors[f],xlab="Leading logFC dim 1", ylab="Leading logFC dim 2")
legend(2, 0 ,legend = c("HH_Pr","HL_Pr","LH_Pr","LL_Pr"), col=c("red","orange","green","blue"), pch=16, bty="n")
dev.off()

# PCA across voom-corrected data
pcavoom<- prcomp(t(na.omit(v2$E)))
summary(pcavoom)
pca_pri <- as.data.frame(prcomp(t(na.omit(v2$E)))$x)
pca_pri$f<-f
pcaplot_pri <- qplot(x=PC1, y=PC2, data=pca_pri,colour=pca_pri$f, shape=pca_pri$f, size=I(4))
pcaplot_pri <- pcaplot_pri + scale_colour_manual(values=c("#D55E00", "#993399", "#009E73","#0072B2"), name=NULL)
pcaplot_pri <- pcaplot_pri + scale_shape_manual(values=c(17,17,17,17))
pcaplot_pri <- pcaplot_pri + theme_bw()
pcaplot_pri <- pcaplot_pri + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=14),legend.text=element_text(size=14),
                panel.background = element_blank(), axis.line = element_line(color = "black"), axis.text.y = element_text(angle = 90), legend.key = element_blank())
pcaplot_pri

#multivariate analysis
colData <- data.frame(row.names = colnames(v2$E), bucket=c("HH1_Pr","HH2_Pr","HH3_Pr","HL1_Pr","HL2_Pr","HL3_Pr","LH1_Pr","LH2_Pr","LH3_Pr","LL1_Pr","LL2_Pr","LL3_Pr"))
colData$treat = substr(colData$bucket, 1,2)
colData$treat_temp = substr(colData$bucket, 1,1)
colData$treat_pCO2 = substr(colData$bucket, 2,2)

conditions=colData
conditions$treat=as.factor(conditions$treat)
conditions$treat_temp=as.factor(conditions$treat_temp)
conditions$treat_pCO2=as.factor(conditions$treat_pCO2)

#analysis using euclidean distance method
adonis(t(na.omit(v2$E))~treat_temp*treat_pCO2,data=conditions,method="eu")
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#treat_temp             1     20201 20200.7  3.0438 0.22790  0.001 ***
#  treat_pCO2             1      8475  8474.8  1.2770 0.09561  0.084 .  
#treat_temp:treat_pCO2  1      6870  6869.5  1.0351 0.07750  0.296    
#Residuals              8     53094  6636.8         0.59899           
#Total                 11     88639                 1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#pie chart
slices <-c(22.8, 9.6, 7.7, 59.9)
labs <- c("Temperature","pCO2","Interaction","Residuals")
labs <- paste(labs, slices) #add percents to labels
labs <- paste(labs,"%",sep="") #add % to labels
piecol <- c("grey70","grey10","grey28","white")
pie(slices,labels=labs, col=piecol,main="Prism")

# Apply limma pipelines for differential expression at prism stage
model<-read.delim("prism_model.txt", header=TRUE)
# Set treatment factors 
TempTreat <- factor(model$Temperature, levels=c("high","low"))
pCO2Treat <- factor(model$pCO2, levels=c("high","low"))

# Create a model matrix of samples x temperature treatment
designTemp <-model.matrix(~0+TempTreat)
colnames(designTemp) <- c("TempHigh","TempLow")
designTemp

fit <- lmFit(v2,designTemp,correlation=corfit$consensus) 
fit <- eBayes(fit) 
topTable(fit)
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
write.table(pr.temp.highvslow, "pr.temp.highvslow.txt")

# write out for GO_MWu
pr.temp.highvslow_all <- topTable(fit2,coef=1,adjust="BH", number=Inf)
write.table(pr.temp.highvslow_all, "pr.temp.highvslow_all.txt")

# Store gastrula temperature results and summarize
results <- decideTests(fit2)
summary(results)
#TempHigh - TempLow
#Down                 4286
#NotSig              25053
#Up                   4132

# Create a model matrix of samples x pCO2 treatment
designpCO2 <-model.matrix(~0+pCO2Treat)
colnames(designpCO2) <- c("pCO2High","pCO2Low")
designpCO2

fit <- lmFit(v2,designpCO2,correlation=corfit$consensus) 
fit <- eBayes(fit) 
topTable(fit)
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
write.table(pr.pCO2.highvslow, "pr.pCO2.highvslow.txt")

# write out for GO_MWU
pr.pCO2.highvslow_all <- topTable(fit2,coef=1,adjust="BH", number=Inf)
write.table(pr.pCO2.highvslow_all, "pr.pCO2.highvslow_all.txt")

# Store results and summarize
results <- decideTests(fit2)
summary(results)
#pCO2High - pCO2Low
#Down                   60
#NotSig              33408
#Up                      3

# To make a Grouped Bar Plot of count data
DEcounts <- structure(list(Temp= c(4132,4286), pCO2 = c(3,60)),
                      .Names = c("Temperature (17 vs. 13C)", "pCO2 (1050 vs. 475uatm)"), class = "data.frame", row.names = c(NA,-2L))
attach(data)
print(data)
barcolors <- c("firebrick1", "dodgerblue2")
barplot(as.matrix(DEcounts), ylab = "Gene counts", ylim = c(0,5000),cex.lab = 1.2, beside=TRUE, col=barcolors)
legend(7,3300, c("Down-regulated","Up-regulated"), cex=1, bty="n", fill=barcolors)

