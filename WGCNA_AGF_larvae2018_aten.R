
library(WGCNA)
library(stringr)
library(flashClust)
options(stringsAsFactors=FALSE)

lnames= load(file="AGF_larval_2018_WaldGroupContrasts.Rdata")

colData=read.csv("J19188meta.csv", header=T)
colData=colData[!duplicated(colData$sample), ]
rownames(colData)<-colData$sample
dim(colData) #96  4
colData$cross=as.factor(colData$cross)
colData$treatment=as.factor(colData$treatment)
colData$dam=as.factor(str_sub(colData$cross, 1,2))
colData$sire =as.factor(str_sub(colData$cross,4,5))

colData->allTraits

Pre<-allTraits$treatment=="Pre"
Pre[Pre==T]<- 1
Pre[Pre==F]<- 0
allTraits$Pre=Pre

Hot<-allTraits$treatment=="Hot"
Hot[Hot ==T]<- 1
Hot[Hot ==F]<- 0
allTraits$Hot = Hot

Ambient<-allTraits$treatment=="Ambient"
Ambient[Ambient ==T]<- 1
Ambient[Ambient ==F]<- 0
allTraits$Ambient = Ambient

Ndam<-(allTraits$dam=="SB"| allTraits$dam=="CU" | allTraits$dam=="LS")
Ndam[Ndam ==T]<- 1
Ndam[Ndam ==F]<- 0
allTraits$Ndam = Ndam

Cdam<-(allTraits$dam=="DR" | allTraits$dam=="BK")
Cdam[Cdam ==T]<- 1
Cdam[Cdam ==F]<- 0
allTraits$Cdam = Cdam

Nsire<-(allTraits$sire=="SB"| allTraits$sire=="CU" | allTraits$sire=="LS")
Nsire[Nsire ==T]<- 1
Nsire[Nsire ==F]<- 0
allTraits$Nsire = Nsire

Csire<-(allTraits$sire=="DR" | allTraits$sire=="BK")
Csire[Csire ==T]<- 1
Csire[Csire ==F]<- 0
allTraits$Csire = Csire

SiteXSite<-(allTraits$cross=="DRxDR"| allTraits$cross=="CUxCU" | allTraits$cross=="LSXLS" | allTraits$cross=="BKxBK"| allTraits$cross=="SBxSB")
SiteXSite[SiteXSite ==T]<- 1
SiteXSite[SiteXSite ==F]<- 0
allTraits$SiteXSite = SiteXSite

NXNcrosses<-(allTraits$cross=="SBxSB"| allTraits$cross=="CUxCU" | allTraits$cross=="LSXLS" | allTraits$cross=="CUxSB"| allTraits$cross=="LSxCU" | allTraits$cross=="LSxSB")
NXNcrosses[NXNcrosses ==T]<- 1
NXNcrosses[NXNcrosses ==F]<- 0
allTraits$NXNcrosses = NXNcrosses

CrossSiteCrosses<-(allTraits$cross=="SBxBK"| allTraits$cross=="CUxBK" | allTraits$cross=="DRxSB" )
CrossSiteCrosses[CrossSiteCrosses ==T]<- 1
CrossSiteCrosses[CrossSiteCrosses ==F]<- 0
allTraits$CrossSiteCrosses = CrossSiteCrosses

NdamXCsire<-(allTraits$cross=="SBxBK"| allTraits$cross=="CUxBK" )
NdamXCsire[NdamXCsire ==T]<- 1
NdamXCsire[NdamXCsire ==F]<- 0
allTraits$NdamXCsire = NdamXCsire

CdamXNsire<-(allTraits$cross=="DRxSB" )
CdamXNsire[CdamXNsire ==T]<- 1
CdamXNsire[CdamXNsire ==F]<- 0
allTraits$CdamXNsire = CdamXNsire


allTraits=as.data.frame(allTraits[-c(2:3,6:7)])
datTraits=allTraits
######################################
dim(vstG.df) #11192    96
datExpr= as.data.frame(t(vstG.df[,]))
dim(datExpr) #96 11192

################################# LOAD INTO LS5, REST OF THIS PERFORMED THERE
library(WGCNA)
library(flashClust)

gsg=goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK #If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:

#if (!gsg$allOK)
#	{if (sum(!gsg$goodGenes)>0)
#		printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse= ", ")));
#		if (sum(!gsg$goodSamples)>0)
#			printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse=", ")))
#		datExpr0= datExpr0[gsg$goodSamples, gsg$goodGenes]
#		}

#Now we cluster the samples to find outliers
sampleTree= flashClust(dist(datExpr), method="average")
sizeGrWindow(12,9)
par(cex=0.6)
par(mar= c(0,4,2,0))
plot(sampleTree, main= "Sample Clustering to Detect Outliers", sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)
abline(h=30, col="red")

#form a data frame analogous to expression data that will hold the traits.
rownames(datTraits) <- datTraits$sample
datTraits$sample <- NULL
datTraits= as.data.frame(datTraits)
table(rownames(datTraits)==rownames(datExpr)) #should return TRUE if datasets align correctly, otherwise your names are out of order

#expression data is in datExpr, corresponding clinical traits are datTraits
sampleTree2=flashClust(dist(datExpr), method="average")
traitColors= numbers2colors(datTraits, signed= FALSE)
plotDendroAndColors(sampleTree2, traitColors, groupLabels= names(datTraits), main="Sample Dendrogram and Trait heatmap")

save(datExpr, datTraits, file="AGF_larvae2018_bm10_SamplesAndTraits_aten.RData")

########################## Run this as a script "SoftThresh_TACC.R", ran for <5 mins
library(WGCNA)
library(flashClust)
ll=load("AGF_larvae2018_bm10_SamplesAndTraits_aten.RData")
options(stringsAsFactors = FALSE);
###Step by step network construction and module detection

powers= c(c(1:10), seq(from =12, to=20, by=2)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(datExpr, powerVector=powers, verbose =5, networkType="signed") #call network topology analysis function

sizeGrWindow(9,5)
par(mfrow= c(1,2))
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
############################ 

#from this plot, we would choose a power of 18 becuase it's the lowest power for which the scale free topology index reaches 0.90

########################### Run this as script "TOM.R" or "TOM_merged.R"
library(WGCNA)
library(flashClust)
ll=load("AGF_larvae2018_bm10_SamplesAndTraits_aten.RData")
options(stringsAsFactors = FALSE);

softPower = 18;
adjacency = adjacency(datExpr, power = softPower, type='signed');

TOM = TOMsimilarity(adjacency, TOMType='signed');
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
#pdf(file=Dendrogram_signed_BM10.pdf, width=20, height=20)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
#dev.off()

minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = 2, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize);
table(dynamicMods)
#dynamicMods
#   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 711 1748 1737 1679 1040  766  607  605  421  394  314  243  201  166  162  122 
#  16   17   18 
# 114   85   77

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
#dynamicColors
#       black         blue        brown         cyan        green  greenyellow 
#         605         1737         1679          162          766          243 
#        grey       grey60    lightcyan   lightgreen      magenta midnightblue 
#         711           85          114           77          394          122 
#        pink       purple          red       salmon          tan    turquoise 
#         421          314          607          166          201         1748 
#      yellow 
#        1040 

# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
#pdf(file=Dendrogram_signed_BM10_colors.pdf, width=20, height=20)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
#dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
#sizeGrWindow(7, 6)
#pdf(file=ClusteringEigengenes.pdf, width=20, height=20)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
#dev.off()

MEDissThres = 0.1 #merging threshold
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "AGF_larvae2018_networkConstruct_signedsft14_bm10_merged_aten.RData")

############################################# finished in less than 5 mins
library(WGCNA)
library(flashClust)
options(stringsAsFactors = FALSE);
ll=load("AGF_larvae2018_bm10_SamplesAndTraits_aten.RData")
ll2= load(file="AGF_larvae2018_networkConstruct_signedsft14_bm10_merged_aten.RData")

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors,softPower=18)$eigengenes #change softPower
MEs = orderMEs(MEs0)

#correlations of traits with eigengenes
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))

#sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

#############################################           
#plotting massive table of all information - module membership, genes, gene names, etc.
#setwd("/Users/mariestrader/Dropbox/MEDIP/R_AnalysisMarch2018")
#annot=read.table("gene_info_table_header.txt",header=T,sep="	\t",quote=NULL)
#names(datExpr)<-gsub("-tr","",names(datExpr))
#names(datExpr)<-sub("transcript:","", names(datExpr))

probes = names(datExpr)
probes2annot = match(probes,annot$spu_id)
datGS.Traits=data.frame(cor(datExpr,datTraits,use="p"))
names(datGS.Traits)=paste("cor",names(datGS.Traits),sep=".")
datME=moduleEigengenes(datExpr,moduleColors)$eigengenes
datKME=signedKME(datExpr, datME, outputColumnName="MM.")
datOutput=data.frame(ProbeID=names(datExpr),moduleColors,datKME,datGS.Traits)
dim(datOutput) #11192    29

write.table(datOutput,"AGF_larvae2018_signedsft18_bm10_merged0.1_datOutput_aten.txt",row.names=F,sep="\t", col.names=T, quote=F)

###Making tables for GO analysis, categorical, interesting modules vs entire expression set

table(moduleColors)
#moduleColors
#      black        blue       brown       green greenyellow        grey      grey60 
#        605        1737        1841         766         243         711        2069 
# lightgreen     magenta        pink      purple         red      salmon         tan 
#         77         394         421         314         607         166         201 
#     yellow 
#       1040 

##############Go categorical
datOutput <- read.table("AGF_larvae2018_signedsft18_bm10_merged0.1_datOutput_aten.txt", header=T, sep="\t")

col="magenta"

#Categorical
tab=datOutput[,c(1,2)]

#tab2=rbind(tab,rest)
tab$moduleColors=as.character(tab$moduleColors)

tab$moduleColors[tab$moduleColors!=col]<-0
tab$moduleColors[tab$moduleColors==col]<-1 
tab$moduleColors=as.factor(tab$moduleColors) 
summary(tab) #do counts match table of module colors?
#tab$ProbeID<- sub("^", "gene:", tab$ProbeID )

write.csv(tab,file="GO/GO_magenta_categorical.csv",quote=F,row.names=F)


#################

##########To output ME by sample and plot boxplots
datOutput <- read.table("AGF_larvae2018_signedsft18_bm10_merged0.1_datOutput_aten.txt", header=T, sep="\t")
ll=load("AGF_larvae2018_bm10_SamplesAndTraits_aten.RData")
ll2= load(file="AGF_larvae2018_networkConstruct_signedsft18_bm10_merged_aten.RData")

colData=read.csv("J19188meta.csv", header=T)
colData=colData[!duplicated(colData$sample), ]
rownames(colData)<-colData$sample
dim(colData) #96  3

all(rownames(colData) %in% row.names(datExpr)) #TRUE

meout<-data.frame(cbind(row.names(datExpr),MEs,colData))
meout$treatment=as.factor(meout$treatment)

boxplot(MEturquoise~survivalheat, data=meout, ylab="Brown Eigengene Expression")
boxplot(MEblack~survivalheat, data=meout, ylab="Black Eigengene Expression")
boxplot(MEtan~survivalheat, data=meout, ylab="Tan Eigengene Expression")
boxplot(MEsalmon~survivalheat, data=meout, ylab="Salmon Eigengene Expression")
boxplot(MEmidnightblue~survivalheat, data=meout, ylab="Midnightblue Eigengene Expression")
boxplot(MElightcyan~cross, data=meout_ambo, ylab="LightCyan Eigengene Expression")

meout $treat_f = factor(meout$treatment, levels=c("Pre","Ambient","Hot"))
meout $cross_f = factor(meout$cross, levels=c('LSxSB','LSxLS','DRxDR','SBxSB','LSxCU','DRxSB','BKxBK','CUxCU','CUxSB','CUxBK','SBxBK'))
ggplot(aes(y = MEblue, x = cross_f, fill = treat_f), data = meout) + 
	geom_boxplot()+
	scale_fill_manual(values=c("Pre"="#EBCC2A","Ambient"="#3B9AB2", "Hot"="#F21A00"))+
	geom_hline(yintercept = 0)+
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
			theme_bw()

write.csv(meout,"MEbySample_forboxplots.csv",quote=F,row.names=F)





