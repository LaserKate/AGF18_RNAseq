
library('DESeq2')
library('arrayQualityMetrics')
library('vegan')
library('rgl')
library('ape')
library('pheatmap')
library('dplyr')
library('adegenet')
library('VennDiagram')
library('stringr')
library(tidyr)
library(readr)
library(Rmisc)
library(ggplot2)
library(wesanderson)
library(limma)
library(KOGMWU)
library(ggridges)
library(cowplot)

gcountsALL=read.table("AGF_juvenile_2018_geneCounts_aten.txt", header=T) 

# clean up file names
colnames(gcountsALL)<-sub("X","", colnames(gcountsALL))
colnames(gcountsALL) <- sub("_[^_]+$", "", colnames(gcountsALL))
colnames(gcountsALL)<-gsub("\\.","_", colnames(gcountsALL))

length(gcountsALL[,1]) #25349
dim(gcountsALL) #24952   146

#removing SS samples that are 1/2 D fracs
#gcountsALL <- gcountsALL[ -c(19,42,105,121,145) ]
#sample 30 (column 105), is SBxSB, SS, Hot
#sample 57 (column 121), is BKxBK, SS, Hot
#sample 94 (column 145), is SBxSB, SS, Hot
#sample 125 (column 19) is DRxSB, SS, Hot
#sample 157 (column 42) is DRxSB, SS, Hot
#sample 93 is BKxBK, SS, Hot, also has fair bit of D transcripts but not as much as the others. 

summary(colSums(gcountsALL))
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 263294  491815  607431  641193  736221 1650922


### REMOVE GENES WITH LOW MEAN COUNTS ###

mns = apply(gcountsALL, 1, mean)

gcounts=gcountsALL[mns>10,] #get rid of genes that show little or no expression


table(mns > 10)
#FALSE  TRUE 
#15612  9340  
dim(gcounts) #9340  146


### BUILD A DATAFRAME ASSOCIATING SAMPLE NAMESWITH TREATMENT CONDITIONS ###
colData=read.csv("J19234meta.csv", header=T)

colData=colData[!duplicated(colData$sample), ]
rownames(colData)<-colData$sample
dim(colData) #146  10
colData$group1=factor(paste0(colData$Cross,colData$Zoox))
colData$group2=factor(paste0(colData$Zoox,colData$Treatment.1))
colData$group3=factor(paste0(colData$Cross,colData$Treatment.1))

colDataH=colData[colData$Treatment.1=="Hot",] #27 13

purebred<-(colData$Cross=="SB"| colData$dam=="CU" | colData$dam=="LS")
#Ndam[Ndam ==T]<- 1
#Ndam[Ndam ==F]<- 0
#geneL$Ndam = as.factor(Ndam)

all(colnames(gcounts) %in% (rownames(colData))) #TRUE
colData <- colData[colnames(gcounts),]
all(rownames(colData) == colnames(gcounts))


### OUTLIERS ### - not going to remove any 

dds<-DESeqDataSetFromMatrix(countData=gcounts, colData=colData, design= ~ Treatment.1+Cross+Zoox)
vsd=varianceStabilizingTransformation(dds, blind=T)
e=ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
arrayQualityMetrics(e, intgroup=c("Zoox"), force=T, outdir= "report_for_genes_treat_2022")
#arrayQualityMetrics(e, intgroup=c("Treatment.1"), force=T, outdir= "report_for_genes_cross")

################################ WALD TEST - FULL MODEL, no interactions ###

dds<-DESeqDataSetFromMatrix(gcounts,
	colData = colData, 
	design = ~Cross+Zoox+Treatment.1)
vst=vst(dds)
vstCZT.df = assay(vst)

######## test for overall treatment effect in juveniles
dds<-DESeq(dds, minReplicatesForReplace=Inf) 
resHT<-results(dds, contrast=c('Treatment.1', 'Hot', 'Ambient')) #here is where the two contrasting conditions get defined
mcols(resHT,use.names=TRUE)
summary(resHT)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2311, 25%
#LFC < 0 (down)     : 2584, 28%
#outliers [1]       : 1, 0.011%
#low counts [2]     : 0, 0%

table(resHT$padj < 0.05)
#FALSE  TRUE 
#  5108  4231 
res=data.frame(cbind("gene"=row.names(resHT),"stat"= resHT$stat))
head(res)
write.csv(res,file="GO_symb/resJuves.heat_stat_2022.csv",quote=F, row.names=F)

save(vstCZT.df, resHT, colData,file="AGF_juvenile_2018_WaldCrossTreatZoox.Rdata")


################################## WALD TEST - FULL MODEL - INTERACTIONS ###

dds<-DESeqDataSetFromMatrix(gcounts,
	colData = colData, 
	design = ~group3)

vstG3=vst(dds)
vstG3.df = assay(vstG3)
dds<-DESeq(dds, minReplicatesForReplace=Inf) 

### cross contrasts ###
# population by population 
# BKBK - Backnumbers = Central

resBKBK.HA<-results(dds, contrast=c('group3', 'BKxBKHot', 'BKxBKAmbient')) #here is where the two contrasting conditions get defined
mcols(resBKBK.HA,use.names=TRUE)
summary(resBKBK.HA)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 318, 3.4%
#LFC < 0 (down)     : 292, 3.1%
#outliers [1]       : 1, 0.011%
#low counts [2]     : 543, 5.8%

table(resBKBK.HA $padj < 0.05)
#FALSE  TRUE 
# 8409   387   
res=data.frame(cbind("gene"=row.names(resBKBK.HA),"stat"= resBKBK.HA$stat))
head(res)
write.csv(res,file="resBKBK.HA_stat.csv",quote=F, row.names=F)

# DR - Davies = Central
resDRDR.HA<-results(dds, contrast=c('group3', 'DRxDRHot', 'DRxDRAmbient')) #here is where the two contrasting conditions get defined
mcols(resDRDR.HA,use.names=TRUE)
summary(resDRDR.HA)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 635, 6.8%
#LFC < 0 (down)     : 635, 6.8%
#outliers [1]       : 1, 0.011%

table(resDRDR.HA $padj < 0.05)
#FALSE  TRUE 
# 8601   738
res=data.frame(cbind("gene"=row.names(resDRDR.HA),"stat"= resDRDR.HA$stat))
head(res)
write.csv(res,file="resDRDR.HA_stat.csv",quote=F, row.names=F)

# SB - Sand Bank 7 - North 
resSBSB.HA<-results(dds, contrast=c('group3', 'SBxSBHot', 'SBxSBAmbient')) #here is where the two contrasting conditions get defined
mcols(resSBSB.HA,use.names=TRUE)
summary(resSBSB.HA)
#out of 9311 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 214, 2.3%
#LFC < 0 (down)     : 290, 3.1%
#outliers [1]       : 1, 0.011%

table(resSBSB.HA $padj < 0.05)
#FALSE  TRUE 
# 8823   516 
res=data.frame(cbind("gene"=row.names(resSBSB.HA),"stat"= resSBSB.HA $stat))
head(res)
write.csv(res,file="resSBSB.HA_stat.csv",quote=F, row.names=F)

# CU - Curd Reef - North 
# no hot juves

# LS - Long Sandy - North 
resLSLS.HA<-results(dds, contrast=c('group3', 'LSxLSHot', 'LSxLSAmbient')) #here is where the two contrasting conditions get defined
mcols(resLSLS.HA,use.names=TRUE)
summary(resLSLS.HA)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 286, 3.1%
#LFC < 0 (down)     : 308, 3.3%
#outliers [1]       : 1, 0.011%
#low counts [2]     : 362, 3.9%

table(resLSLS.HA $padj < 0.05)
#FALSE  TRUE 
# 8689   288   
res=data.frame(cbind("gene"=row.names(resLSLS.HA),"stat"= resLSLS.HA $stat))
head(res)
write.csv(res,file="resLSLS.HA_stat.csv",quote=F, row.names=F)

## North by North population crosses ##
## CUXSB
resCUSB.HA<-results(dds, contrast=c('group3', 'CUxSBHot', 'CUxSBAmbient')) #here is where the two contrasting conditions get defined
mcols(resCUSB.HA,use.names=TRUE)
summary(resCUSB.HA)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 4, 0.043%
#LFC < 0 (down)     : 37, 0.4%
#outliers [1]       : 1, 0.011%
#low counts [2]     : 1448, 16%
 
table(resCUSB.HA $padj < 0.05)
#FALSE  TRUE 
#  7873    18   
res=data.frame(cbind("gene"=row.names(resCUSB.HA),"stat"= resCUSB.HA $stat))
head(res)
write.csv(res,file="resCUSB.HA_stat.csv",quote=F, row.names=F)

## LSXCU
resLSCU.HA<-results(dds, contrast=c('group3', 'LSxCUHot', 'LSxCUAmbient')) #here is where the two contrasting conditions get defined
mcols(resLSCU.HA,use.names=TRUE)
summary(resLSCU.HA)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 117, 1.3%
#LFC < 0 (down)     : 131, 1.4%
#outliers [1]       : 1, 0.011%
#low counts [2]     : 2715, 29%

table(resLSCU.HA $padj < 0.05)
#FALSE  TRUE 
#6504   120      
res=data.frame(cbind("gene"=row.names(resLSCU.HA),"stat"= resLSCU.HA $stat))
head(res)
write.csv(res,file="resLSCU.HA_stat.csv",quote=F, row.names=F)


## Central by Central population crosses ## - There are none! 

## Flip Flop Crosses- North Mom X Central Dad ##
##SBXBK
# no juves for this cross

##CUXDR
resCUDR.HA<-results(dds, contrast=c('group3', 'CUxDRHot', 'CUxDRAmbient')) #here is where the two contrasting conditions get defined
mcols(resCUDR.HA,use.names=TRUE)
summary(resCUDR.HA)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 8, 0.086%
#LFC < 0 (down)     : 33, 0.35%
#outliers [1]       : 1, 0.011%
#low counts [2]     : 0, 0%

table(resCUDR.HA$padj < 0.05)
#FALSE  TRUE 
# 9326    13    
res=data.frame(cbind("gene"=row.names(resCUDR.HA),"stat"= resCUDR.HA$stat))
head(res)
write.csv(res,file="resCUDR.HA_stat.csv",quote=F, row.names=F)


## Flip Flop Crosses- Central Mom X North Dad ##
##DRXSB
resDRSB.HA<-results(dds, contrast=c('group3', 'DRxSBHot', 'DRxSBAmbient')) #here is where the two contrasting conditions get defined
mcols(resDRSB.HA,use.names=TRUE)
summary(resDRSB.HA)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 471, 5%
#LFC < 0 (down)     : 652, 7%
#outliers [1]       : 1, 0.011%
#low counts [2]     : 0, 0%

table(resDRSB.HA $padj < 0.05)
#FALSE  TRUE 
#  8707   632      
res=data.frame(cbind("gene"=row.names(resDRSB.HA),"stat"= resDRSB.HA $stat))
head(res)
write.csv(res,file="resDRSB.HA_stat.csv",quote=F, row.names=F)

save(vstG3,vstG3.df,colData,resBKBK.HA,resCUDR.HA,resCUSB.HA,resDRDR.HA,resDRSB.HA,resLSCU.HA,resLSLS.HA,resSBSB.HA,file="AGF_juvenile_2018_WaldGroup3Contrasts_crossTemp_aten.Rdata")

#################################### WALD TEST - FULL MODEL - INTERACTIONS ###

dds<-DESeqDataSetFromMatrix(gcounts,
	colData = colData, 
	design = ~group2)

vstG2=vst(dds)
vstG2.df = assay(vstG2)
dds<-DESeq(dds, minReplicatesForReplace=Inf) 

#####SED
resSED.HA<-results(dds, contrast=c('group2', 'SEDHot', 'SEDAmbient')) #here is where the two contrasting conditions get defined
mcols(resSED.HA,use.names=TRUE)
summary(resSED.HA)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1793, 19%
#LFC < 0 (down)     : 2019, 22%
#outliers [1]       : 3, 0.032%
#low counts [2]     : 0, 0%

table(resSED.HA$padj < 0.05)
#FALSE  TRUE 
# 6148  3189     
res=data.frame(cbind("gene"=row.names(resSED.HA),"stat"= resSED.HA$stat))
head(res)
write.csv(res,file="resSED.HA_stat.csv",quote=F, row.names=F)

#####SS
resSS.HA<-results(dds, contrast=c('group2', 'SSHot', 'SSAmbient')) #here is where the two contrasting conditions get defined
mcols(resSS.HA,use.names=TRUE)
summary(resSS.HA)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 959, 10%
#LFC < 0 (down)     : 1118, 12%
#outliers [1]       : 3, 0.032%
#low counts [2]     : 0, 0%

table(resSS.HA$padj < 0.05)
#FALSE  TRUE 
# 7771  1566    
res=data.frame(cbind("gene"=row.names(resSS.HA),"stat"= resSS.HA$stat))
head(res)
write.csv(res,file="resSS.HA_stat.csv",quote=F, row.names=F)

#####D1
resD1.HA<-results(dds, contrast=c('group2', 'D1Hot', 'D1Ambient')) #here is where the two contrasting conditions get defined
mcols(resD1.HA,use.names=TRUE)
summary(resD1.HA)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 997, 11%
#LFC < 0 (down)     : 1179, 13%
#outliers [1]       : 3, 0.032%
#low counts [2]     : 0, 0%

table(resD1.HA$padj < 0.05)
#FALSE  TRUE 
#7780  1557   
res=data.frame(cbind("gene"=row.names(resD1.HA),"stat"= resD1.HA$stat))
head(res)
write.csv(res,file="resD1.HA_stat.csv",quote=F, row.names=F)

save(vstG2,vstG2.df,colData,resSED.HA,resSS.HA,resD1.HA,file="AGF_juvenile_2018_WaldGroup2Contrasts_zooxTemp_aten.Rdata")

#################################### WALD TEST - Symbiosis only at ambient ###

colDataA=colData[colData$Treatment.1=="Ambient",]

C<-colDataA$Zoox=="C1"
C[C ==T]<- 1
C[C ==F]<- 0
colDataA $C1 = C

D<-colDataA$Zoox=="D1"
D[D ==T]<- 1
D[D ==F]<- 0
colDataA $D1 = D

SS<-colDataA$Zoox=="SS"
SS[SS ==T]<- 1
SS[SS ==F]<- 0
colDataA $SS = SS

SD<-colDataA$Zoox=="SED"
SD[SD ==T]<- 1
SD[SD ==F]<- 0
colDataA $SD = SD

colDataA$C1=as.factor(colDataA$C1)
colDataA$D1=as.factor(colDataA$D1)
colDataA$SS=as.factor(colDataA$SS)
colDataA$SD=as.factor(colDataA$SD)

all(rownames(colDataA) %in% colnames(gcounts)) #TRUE
gcountsA <- gcounts[, rownames(colDataA)]
all(rownames(colDataA) == colnames(gcountsA))


dds<-DESeqDataSetFromMatrix(gcountsA,
	colData = colDataA, 
	design = ~Zoox)

vstZ=vst(dds)
vstZ.df = assay(vstZ)
dds<-DESeq(dds, minReplicatesForReplace=Inf) 

resCD<-results(dds, contrast=c('Zoox', 'C1', 'D1')) #here is where the two contrasting conditions get defined
mcols(resCD,use.names=TRUE)
summary(resCD)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 591, 6.3%
#LFC < 0 (down)     : 554, 5.9%
#outliers [1]       : 10, 0.11%
#low counts [2]     : 0, 0%

table(resCD $padj < 0.05)
#FALSE  TRUE 
# 8621   709      
res=data.frame(cbind("gene"=row.names(resCD),"stat"= resCD$stat))
head(res)
write.csv(res,file="resCD_stat_2022.csv",quote=F, row.names=F)

resCSED<-results(dds, contrast=c('Zoox', 'C1', 'SED')) #here is where the two contrasting conditions get defined
mcols(resCSED,use.names=TRUE)
summary(resCSED)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1872, 20%
#LFC < 0 (down)     : 1918, 21%
#outliers [1]       : 10, 0.11%
#low counts [2]     : 0, 0%

table(resCSED$padj < 0.05)
#FALSE  TRUE 
# 6169  3161 
res=data.frame(cbind("gene"=row.names(resCSED),"stat"= resCSED$stat))
head(res)
write.csv(res,file="resCSED_stat_2022.csv",quote=F, row.names=F)

resCSS<-results(dds, contrast=c('Zoox', 'C1', 'SS')) #here is where the two contrasting conditions get defined
mcols(resCSS,use.names=TRUE)
summary(resCSS)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 645, 6.9%
#LFC < 0 (down)     : 700, 7.5%
#outliers [1]       : 10, 0.11%

table(resCSS$padj < 0.05)
#FALSE  TRUE 
#8364   966  
res=data.frame(cbind("gene"=row.names(resCSS),"stat"= resCSS$stat))
head(res)
write.csv(res,file="resCSS_stat_2022.csv",quote=F, row.names=F)

resDSED<-results(dds, contrast=c('Zoox', 'D1', 'SED')) #here is where the two contrasting conditions get defined
mcols(resDSED,use.names=TRUE)
summary(resDSED)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1025, 11%
#LFC < 0 (down)     : 888, 9.5%
#outliers [1]       : 10, 0.11%
#low counts [2]     : 0, 0%

table(resDSED$padj < 0.05)
#FALSE  TRUE 
#7868  1462 
res=data.frame(cbind("gene"=row.names(resDSED),"stat"= resDSED$stat))
head(res)
write.csv(res,file="resDSED_stat_2022.csv",quote=F, row.names=F) 

resDSS<-results(dds, contrast=c('Zoox', 'D1', 'SS')) #here is where the two contrasting conditions get defined
mcols(resDSS,use.names=TRUE)
summary(resDSS)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 166, 1.8%
#LFC < 0 (down)     : 163, 1.7%
#outliers [1]       : 10, 0.11%
#low counts [2]     : 1082, 12%

table(resDSS$padj < 0.05)
#FALSE  TRUE 
#  8041   207
res=data.frame(cbind("gene"=row.names(resDSS),"stat"= resDSS$stat))
head(res)
write.csv(res,file="resDSS_stat_2022.csv",quote=F, row.names=F)

resSSSED<-results(dds, contrast=c('Zoox', 'SS', 'SED')) #here is where the two contrasting conditions get defined
mcols(resSSSED,use.names=TRUE)
summary(resSSSED)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1260, 13%
#LFC < 0 (down)     : 1208, 13%
#outliers [1]       : 10, 0.11%
#low counts [2]     : 0, 0%

table(resSSSED$padj < 0.05)
#FALSE  TRUE 
# 7459  1871 
res=data.frame(cbind("gene"=row.names(resSSSED),"stat"= resSSSED$stat))
head(res)
write.csv(res,file="resSSSED_stat_2022.csv",quote=F, row.names=F)

############# run ambo juve samples for each sym type separately for KOG_MWU
dds<-DESeqDataSetFromMatrix(gcountsA,
	colData = colDataA, 
	design = ~C1)

dds<-DESeq(dds, minReplicatesForReplace=Inf) 

resC1all<-results(dds, contrast=c('C1', '0', '1')) #here is where the two contrasting conditions get defined
mcols(resC1all,use.names=TRUE)
summary(resC1all)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1320, 14%
#LFC < 0 (down)     : 1266, 14%
#outliers [1]       : 10, 0.11%
#low counts [2]     : 0, 0%

### D1
dds<-DESeqDataSetFromMatrix(gcountsA,
	colData = colDataA, 
	design = ~D1)

dds<-DESeq(dds, minReplicatesForReplace=Inf) 

resD1all<-results(dds, contrast=c('D1', '0', '1')) #here is where the two contrasting conditions get defined
mcols(resD1all,use.names=TRUE)
summary(resD1all)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 86, 0.92%
#LFC < 0 (down)     : 79, 0.85%
#outliers [1]       : 17, 0.18%
#low counts [2]     : 0, 0%

### SS
dds<-DESeqDataSetFromMatrix(gcountsA,
	colData = colDataA, 
	design = ~SS)

dds<-DESeq(dds, minReplicatesForReplace=Inf) 

resSSall<-results(dds, contrast=c('SS', '0', '1')) #here is where the two contrasting conditions get defined
mcols(resSSall,use.names=TRUE)
summary(resSSall)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 244, 2.6%
#LFC < 0 (down)     : 324, 3.5%
#outliers [1]       : 15, 0.16%
#low counts [2]     : 181, 1.9%

### SD
dds<-DESeqDataSetFromMatrix(gcountsA,
	colData = colDataA, 
	design = ~SD)

dds<-DESeq(dds, minReplicatesForReplace=Inf) 

resSDall<-results(dds, contrast=c('SD', '0', '1')) #here is where the two contrasting conditions get defined
mcols(resSDall,use.names=TRUE)
summary(resSDall)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1669, 18%
#LFC < 0 (down)     : 1633, 17%
#outliers [1]       : 14, 0.15%

save(vstZ.df, resCD, resCSED, resCSS, resDSED, resDSS, resSSSED, resSSall, resSDall, resD1all, resC1all, colData,file="AGF_juvenile_2018_WaldZooxAmbient_aten.Rdata")


#################################### WALD TEST - Symbiosis only at heat ###
colDataH=colData[colData$Treatment.1=="Hot",]

all(rownames(colDataH) %in% colnames(gcounts)) #TRUE
gcountsH <- gcounts[, rownames(colDataH)]
all(rownames(colDataH) == colnames(gcountsH))


dds<-DESeqDataSetFromMatrix(gcountsH,
	colData = colDataH, 
	design = ~Zoox)
vstZ=vst(dds)
vstZ.df = assay(vstZ)
dds<-DESeq(dds, minReplicatesForReplace=Inf) 

resSSDhot<-results(dds, contrast=c('Zoox', 'SS', 'D1')) #here is where the two contrasting conditions get defined
mcols(resSSDhot,use.names=TRUE)
summary(resSSDhot)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 56, 0.6%
#LFC < 0 (down)     : 135, 1.4%
#outliers [1]       : 10, 0.11%

table(resSSDhot$padj < 0.05)
#FALSE  TRUE 
# 9198   132 
res=data.frame(cbind("gene"=row.names(resSSDhot),"stat"= resSSDhot$stat))
head(res)
write.csv(res,file="resSSDhot_stat_2022.csv",quote=F, row.names=F)

resSSSEDhot<-results(dds, contrast=c('Zoox', 'SS', 'SED')) #here is where the two contrasting conditions get defined
mcols(resSSSEDhot,use.names=TRUE)
summary(resSSSEDhot)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 21, 0.22%
#LFC < 0 (down)     : 38, 0.41%
#outliers [1]       : 10, 0.11%
#low counts [2]     : 0, 0%

table(resSSSEDhot$padj < 0.05)
#FALSE  TRUE 
# 9294    36     
res=data.frame(cbind("gene"=row.names(resSSSEDhot),"stat"= resSSSEDhot$stat))
head(res)
write.csv(res,file="resSSSEDhot_stat_2022.csv",quote=F, row.names=F)

resD1SEDhot<-results(dds, contrast=c('Zoox', 'D1', 'SED')) #here is where the two contrasting conditions get defined
mcols(resD1SEDhot,use.names=TRUE)
summary(resD1SEDhot)
#out of 9340 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 201, 2.2%
#LFC < 0 (down)     : 131, 1.4%
#outliers [1]       : 10, 0.11%
#low counts [2]     : 0, 0%

table(resD1SEDhot $padj < 0.05)
#FALSE  TRUE 
#9097   233   
res=data.frame(cbind("gene"=row.names(resD1SEDhot),"stat"= resD1SEDhot $stat))
head(res)
write.csv(res,file="resD1SEDhot_stat_2022.csv",quote=F, row.names=F)


################ KOG_MWU #######################
library(KOGMWU)

###### look at effects of symbiosis on juveniles in ambient conditions
ll=load("AGF_juvenile_2018_WaldZooxAmbient_aten.Rdata") 
gene2kog=read.table("Aten_gene2kogClass.tab",header=T,sep="\t",quote=NULL,fill=T)

C_LFC=as.data.frame((cbind(row.names(resC1all), resC1all $log2FoldChange))) 
C_LFC$V1=sub("-RA","", C_LFC$V1)
C=kog.mwu(C_LFC,gene2kog,Alternative="t")

D_LFC=as.data.frame((cbind(row.names(resD1all), resD1all $log2FoldChange))) 
D_LFC$V1=sub("-RA","", D_LFC$V1)
D=kog.mwu(D_LFC,gene2kog,Alternative="t")

SED_LFC=as.data.frame((cbind(row.names(resSDall), resSDall $log2FoldChange))) 
SED_LFC $V1=sub("-RA","", SED_LFC $V1)
SED=kog.mwu(SED_LFC,gene2kog,Alternative="t")

SS_LFC=as.data.frame((cbind(row.names(resSSall), resSSall$log2FoldChange))) 
SS_LFC $V1=sub("-RA","", SS_LFC $V1)
SS=kog.mwu(SS_LFC,gene2kog,Alternative="t")

ktable=makeDeltaRanksTable(list("C1"=C,"D1"=D,"SS"=SS,"SED"=SED))
#remove nuclear structure because only has 4 genes and values are outliers
ktable <- ktable[-c(4,13), ] #22

heat.colors = colorRampPalette(wes_palette("Zissou1", 30, type = "continuous"), bias=1.0)(100)

pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=heat.colors) 

pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
corrPlot(x="C1",y="D1",ktable)
corrPlot(x="C1",y="SS",ktable)
corrPlot(x="C1",y="SED",ktable)
corrPlot(x="D1",y="SED",ktable)
corrPlot(x="D1",y="SS",ktable)
corrPlot(x="SS",y="SED",ktable)

ktests=rbind(C,D,SS,SED)
ktests$zoox=(c(rep("C",24),rep("D",24),rep("SS",24),rep("SED",24)))
write.csv(ktests,file="kogTestResults_ambo_aten_2022.csv",quote=F, row.names=F)
#Supplementary Data 9

####### look at effects of symbiosis on response to heat in juveniles 
ll=load("AGF_juvenile_2018_WaldGroup2Contrasts_zooxTemp_aten.Rdata") 
gene2kog=read.table("Aten_gene2kogClass.tab",header=T,sep="\t",quote=NULL,fill=T)

D_HA=as.data.frame((cbind(row.names(resD1.HA), resD1.HA$log2FoldChange))) 
D_HA$V1=sub("-RA","", D_HA$V1)
D_HA =kog.mwu(D_HA,gene2kog,Alternative="t")

SS_HA=as.data.frame((cbind(row.names(resSS.HA), resSS.HA$log2FoldChange))) 
SS_HA$V1=sub("-RA","", SS_HA$V1)
SS_HA=kog.mwu(SS_HA,gene2kog,Alternative="t")

SED_HA=as.data.frame((cbind(row.names(resSED.HA), resSED.HA$log2FoldChange))) 
SED_HA$V1=sub("-RA","", SED_HA$V1)
SED_HA=kog.mwu(SED_HA,gene2kog,Alternative="t")

#### load in larval response to heat stress (overall)

ll=load("AGF_larval_2018_WaldCrossTreat.Rdata")

Larva_HA=as.data.frame((cbind(row.names(resHT), resHT$log2FoldChange))) 
Larva_HA $V1=sub("-RA","", Larva_HA$V1)
Larva_HA =kog.mwu(Larva_HA,gene2kog,Alternative="t")

#### load in larval thermal tolerance to heat stress (overall)

ll=load("AGF_larval_2018_PreRes_aten.Rdata")

Larva_Tolerance=as.data.frame((cbind(row.names(resP), resP$log2FoldChange))) 
Larva_Tolerance $V1=sub("-RA","", Larva_Tolerance $V1)
Larva_Tolerance =kog.mwu(Larva_Tolerance,gene2kog,Alternative="t")


ktable=makeDeltaRanksTable(list("D1"=D_HA,"SS"=SS_HA,"SED"=SED_HA,"Larval Response"=Larva_HA, "Larval Tolerance"=Larva_Tolerance))
#remove nuclear structure because only has 4 genes and values are outliers
ktable <- ktable[-c(3,20), ]

heat.colors = colorRampPalette(wes_palette("Zissou1", 30, type = "continuous"), bias=1.8)(100)

pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=heat.colors)

pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
corrPlot(x="D1",y="SS",ktable) 
corrPlot(x="D1",y="SED",ktable) 
corrPlot(x="SS",y="SED",ktable) 
corrPlot(x="SS",y="Larval Response",ktable)
corrPlot(x="SS",y="Larval Tolerance",ktable) 
corrPlot(x="D1",y="Larval Response",ktable) 
corrPlot(x="D1",y="Larval Tolerance",ktable) 
corrPlot(x="SED",y="Larval Response",ktable)
corrPlot(x="SED",y="Larval Tolerance",ktable)
corrPlot(x="Larval Response",y="Larval Tolerance",ktable) 

ktests=rbind(D_HA,SS_HA,SED_HA,Larva_Tolerance,Larva_HA)
ktests$zoox=(c(rep("D",24),rep("SS",24),rep("SED",24), rep("LT",24), rep("LR",24)))

write.csv(ktests,file="kogTestResults_heat_aten.csv",quote=F, row.names=F)
#Supplementary Data 9

### principal coordinate calculation on all samples

All=load("AGF_juvenile_2018_WaldCrossTreatZoox_aten2.Rdata")
colDataH=colData[colData$Treatment.1=="Hot",]
#remove effect of cross
#vstL <- limma::removeBatchEffect(vstCZT.df, colData$Cross)

conditions=colData

#grp=rep("#EBCC2A",ncol(vstCZT.df))
#grp[grep("SED",conditions$Zoox)]="#3B9AB2"
#grp[grep("SS",conditions$Zoox)]="#F21A00"
#grp[grep("D1",conditions$Zoox)]="purple"

grp=rep("#3B9AB2",ncol(vstCZT.df))
grp[grep("Hot",conditions$Treatment.1)]="#F21A00"

fitpc=capscale(dist(t(vstCZT.df),method="manhattan")~1)
#fitpc=capscale(dist(t(vstCZT.df))~1)
summary(eigenvals(fitpc))

plot(fitpc$CA$u,pch=16, col="white",main="AGF Juveniles 2018 All Samples", xlab="PCoA1 10.92%", ylab="PCoA2 6.95%")
abline(h=0,v=0, col="grey")
points(fitpc$CA$u,pch=16, col=grp)
#ordispider(fitpc$CA$u,conditions$Cross,draw="polygon",label=F)
ordiellipse(fitpc$CA$u,conditions$group2,draw="polygon",label=T, kind="sd")

ad=adonis(t(vstCZT.df)~Cross+Treatment.1+Zoox,data=conditions,method="manhattan")
#Call:
#adonis(formula = t(vstCZT.df) ~ Cross + Treatment.1 + Zoox, data = conditions,      method = "manhattan") 
#
#Permutation: free
#Number of permutations: 999
#
#Terms added sequentially (first to last)
#
#             Df  SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#Cross         9  152583453 16953717  2.9492 0.14507  0.001 ***
#Treatment.1   1   56645906 56645906  9.8538 0.05386  0.001 ***
#Zoox          3   83705530 27901843  4.8536 0.07959  0.001 ***
#Residuals   132  758821575  5748648         0.72148           
#Total       145 1051756464                  1.00000           
      


##### PCA for ambo
Amb=load("AGF_juvenile_2018_WaldZooxAmbient_aten.Rdata")
colDataA=colData[colData$Treatment.1=="Ambient",]

#vstL <- limma::removeBatchEffect(vstZ.df, colDataA$Cross)

conditions=colDataA

grp=rep("#EBCC2A",ncol(vstZ.df))
grp[grep("SED",conditions$Zoox)]="#3B9AB2"
grp[grep("SS",conditions$Zoox)]="#F21A00"
grp[grep("D1",conditions$Zoox)]="purple"

fitpc=capscale(dist(t(vstZ.df),method="manhattan")~1)
summary(eigenvals(fitpc))

plot(fitpc$CA$u,pch=16, col="white",main="AGF Juveniles 2018 Ambient Samples", xlab="MDS1 10.63%", ylab="MDS2 6.56%")
abline(h=0,v=0, col="grey")
points(fitpc$CA$u,pch=16, col=grp)
ordiellipse(fitpc$CA$u,conditions$Zoox,draw="polygon",label=T)

ad=adonis(t(vstZ.df)~Cross+Zoox,data=conditions,method="manhattan")

#Call:
#adonis(formula = t(vstZ.df) ~ Cross + Zoox, data = conditions,      method = "manhattan") 
#
#Permutation: free
#Number of permutations: 999
#
#Terms added sequentially (first to last)
#
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#Cross       9 124338372 13815375  2.5366 0.15858  0.001 ***
#Zoox        3  82401261 27467087  5.0432 0.10510  0.001 ***
#Residuals 106 577310883  5446329         0.73632           
#Total     118 784050515                  1.00000           




################### DAPC
library('DESeq2')
library('adegenet')
library('dplyr')
library(ggplot2)
library(wesanderson)
library('stringr')
library(ggridges)
library(cowplot)

ll=load("AGF_juvenile_2018_WaldCrossTreatZoox_aten.Rdata")

colData$group1<-gsub("SED","SD", colData$group1)
colData$group3=factor(paste0(colData$group1,colData$Treatment.1))


##### just symbiosis at ambient, remove effect of cross
ll=load("AGF_juvenile_2018_WaldZooxAmbient_aten.Rdata")
colDataA=colData[colData$Treatment.1=="Ambient",]

colDataA$dam=str_sub(colDataA$Cross, 1,2)
colDataA$sire = str_sub(colDataA$Cross,4,5)
colDataA$group1<-gsub("SED","SD", colDataA$group1)
colDataA$groupall=as.factor(paste0(colDataA$group1, colDataA$Treatment.1))

dim(vstZ.df) #9340  119

vsd2 <- limma::removeBatchEffect(vstZ.df, colDataA$Cross)

dapc0 <- dapc(t(vsd2), colDataA$groupall, n.da=4, n.pca=32)
temp <- optim.a.score(dapc0, n.sim = 5)
#for the vst, they suggest retaining 10 PCs, which explain 0.423 percent of the conserved variance
#for the vst removed effect of cross, retain 13 PCs, which explain 0.461 percent of the conserved variance

dapc <- dapc(t(vsd2), colDataA$groupall, n.da=2, n.pca=13)
#scatter(dapc,bg="white",scree.da=TRUE,scree.pca=TRUE,legend=TRUE,solid=.4)
varexpl <- round((dapc$eig/sum(dapc$eig))[1:2] * 100, 1)

dapc1 <- tibble(sample = rownames(dapc$ind.coord),
               grp = dapc$grp,
               LD1 = dapc$ind.coord[,1],
               LD2 = dapc$ind.coord[,2])
dapc2 <- dapc1 %>%
  group_by(grp) %>%
  summarize(c1 = mean(LD1),
            c2 = mean(LD2)) %>%
  full_join(dapc1)

###symbiosis plus cross
dapc.fig <-   ggplot(dapc2, aes(shape = factor(str_sub(grp, 1,5)), 
		fill = factor(str_sub(grp, 6,7)), color = factor(str_sub(grp,6,7)))) +
		geom_segment(mapping = aes(x = LD1, y = LD2, xend = c1, yend = c2), lwd = 0.25, col = "grey") +
		geom_point(aes(x = c1, y = c2), size = 3) +
		geom_point(aes(x = LD1, y = LD2), size = 2, show.legend = FALSE) +
		scale_color_manual(name="Symbiosis", labels = c("C1", "D1", "Sediments","SS"),values = c("#EBCC2A","purple","#3B9AB2","#F21A00")) +
		scale_fill_manual(name="Symbiosis", labels = c("C1", "D1", "Sediments","SS"),values = c("#EBCC2A","purple","#3B9AB2","#F21A00")) +
		scale_shape_manual(name = "Cross", values=c(21,22,23,24,25,0,1,2,3,8)) +
		guides(fill = guide_legend(override.aes = list(shape = 21, size = 3))) +
		guides(shape = guide_legend(override.aes = list(fill = "black", size = 3))) +
		labs(x = paste0("LD1 [", varexpl[1],"%]"), y = paste0("LD2 [", varexpl[2],"%]")) +
		theme_bw()

### symbiosis
dapc.fig <-   ggplot(dapc2, aes(shape = factor(str_sub(grp, 6,7)), 
		fill = factor(str_sub(grp, 6,7)), color = factor(str_sub(grp,6,7)))) +
		geom_segment(mapping = aes(x = LD1, y = LD2, xend = c1, yend = c2), lwd = 0.25, col = "grey") +
		geom_point(aes(x = c1, y = c2), size = 3) +
		geom_point(aes(x = LD1, y = LD2), size = 2, show.legend = FALSE) +
		scale_color_manual(name="Symbiosis", labels = c("C1", "D1", "Sediments","SS"),values = c("#EBCC2A","purple","#3B9AB2","#F21A00")) +
		scale_fill_manual(name="Symbiosis", labels = c("C1", "D1", "Sediments","SS"),values = c("#EBCC2A","purple","#3B9AB2","#F21A00")) +
		scale_shape_manual(name = "Symbiosis", values=c(21,22,23,24)) +
		guides(fill = guide_legend(override.aes = list(shape = 21, size = 3))) +
		guides(shape = guide_legend(override.aes = list(fill = "black", size = 3))) +
		labs(x = paste0("LD1 [", varexpl[1],"%]"), y = paste0("LD2 [", varexpl[2],"%]")) +
		theme_bw()

############ symbiosis and heat
ll=load("AGF_juvenile_2018_WaldCrossTreatZoox_aten.Rdata")

colData$dam=str_sub(colData$Cross, 1,2)
colData$sire = str_sub(colData$Cross,4,5)
colData$group1<-gsub("SED","SD", colData$group1)
colData$groupall=as.factor(paste0(colData$group1, colData$Treatment.1))

dim(vstCZT.df) #9340  146
colData$group2<-gsub("SED","SD", colData$group2)

vsd2 <- limma::removeBatchEffect(vstCZT.df, colData$Cross)

dapc2 <- dapc(t(vsd2), colData$groupall, n.da=4, n.pca=32)
temp <- optim.a.score(dapc2, n.sim = 5) 
#for the vst, they suggest retaining 11 PCs, which explain 0.436 percent of the conserved variance
#for the vst removed effect of cross, retain 9 PCs, which explain 0.409 percent of the conserved variance

dapc <- dapc(t(vsd2), colData$groupall, n.da=2, n.pca=9)
#scatter(dapc,bg="white",scree.da=TRUE,scree.pca=TRUE,legend=TRUE,solid=.4)
varexpl <- round((dapc$eig/sum(dapc$eig))[1:2] * 100, 1) #31.9 22.8

#scatter(dapc, bg = "white", legend = TRUE, scree.da = FALSE)

dapc1 <- tibble(sample = rownames(dapc$ind.coord),
               grp = dapc$grp,
               LD1 = dapc$ind.coord[,1],
               LD2 = dapc$ind.coord[,2])
dapc2 <- dapc1 %>%
  group_by(grp) %>%
  summarize(c1 = mean(LD1),
            c2 = mean(LD2)) %>%
  full_join(dapc1)

dapc.fig <-   ggplot(dapc2, aes(shape = factor(str_sub(grp, 6,7)), 
		fill = factor(str_sub(grp, 8,10)))) +
		geom_segment(mapping = aes(x = LD1, y = LD2, xend = c1, yend = c2), lwd = 0.25, col = "grey") +
		geom_point(aes(x = c1, y = c2), size = 3) +
		geom_point(aes(x = LD1, y = LD2), size = 2, show.legend = FALSE) +
		scale_shape_manual(name="Symbiosis", labels = c("C1", "D1", "Sediments","SS"), values = c(21,22,23,24)) +
		scale_fill_manual(name="Treatment", values =c("#3B9AB2","#F21A00"))+
		guides(fill = guide_legend(override.aes = list(shape = 21, size = 3))) +
		guides(shape = guide_legend(override.aes = list(fill = "black", size = 3))) +
		labs(x = paste0("LD1 [", varexpl[1],"%]"), y = paste0("LD2 [", varexpl[2],"%]")) +
		theme_bw()

ridges= ggplot(dapc2, aes(x = LD1, y = factor(str_sub(grp, 6,7)), fill = factor(str_sub(grp, 8,10)), height = ..density..)) +
		geom_density_ridges(scale = 1, stat = "density") +
		scale_y_discrete(expand = c(0.01, 0)) +
		scale_x_continuous(expand = c(0.01, 0)) +
		scale_fill_manual(name="Treatment", values =c("#3B9AB2","#F21A00")) +
		theme_ridges() 

fig <- plot_grid(ridges,dapc.fig, ncol = 1, align = "v", axis="b",rel_heights = c(1, 1.5))
