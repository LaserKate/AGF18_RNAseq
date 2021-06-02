
setwd("/Users/mes0192/Dropbox/KQ_GBR_spawnseq/analysis/DESeq2_juves")
setwd("/Users/mariestrader/Dropbox/StuffFromLaptop/Quigley")

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
library(kog.mwu)

gcountsALL=read.table("AGF_juvenile_2018_geneCounts.txt", header=T) 

# clean up file names
colnames(gcountsALL)<-sub("X","", colnames(gcountsALL))
colnames(gcountsALL) <- sub("_[^_]+$", "", colnames(gcountsALL))
colnames(gcountsALL)<-gsub("\\.","_", colnames(gcountsALL))

length(gcountsALL[,1]) #25349
dim(gcountsALL) #25349   146


summary(colSums(gcountsALL))
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 158481  298566  369970  383801  439854  962906 



### REMOVE GENES WITH LOW MEAN COUNTS ###

mns = apply(gcountsALL, 1, mean)

gcounts=gcountsALL[mns>10,] #get rid of genes that show little or no expression


table(mns > 10)
#FALSE  TRUE 
#18712  6637  
dim(gcounts) #6637  146


### BUILD A DATAFRAME ASSOCIATING SAMPLE NAMESWITH TREATMENT CONDITIONS ###
colData=read.csv("J19234meta.csv", header=T)
colData=colData[!duplicated(colData$sample), ]
rownames(colData)<-colData$sample
dim(colData) #146  10
colData$group1=factor(paste0(colData$Cross,colData$Zoox))
colData$group2=factor(paste0(colData$Zoox,colData$Treatment.1))
colData$group3=factor(paste0(colData$Cross,colData$Treatment.1))

purebred<-(colData$Cross=="SB"| colData$dam=="CU" | colData$dam=="LS")
Ndam[Ndam ==T]<- 1
Ndam[Ndam ==F]<- 0
geneL$Ndam = as.factor(Ndam)

all(rownames(colData) %in% colnames(gcounts)) #TRUE
gcounts <- gcounts[, rownames(colData)]
all(rownames(colData) == colnames(gcounts))



### OUTLIERS ### - not going to remove any 

#dds<-DESeqDataSetFromMatrix(countData=gcounts, colData=colData, design= ~ Treatment.1+Cross+Zoox)
#vsd=varianceStabilizingTransformation(dds, blind=T)
#e=ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
#arrayQualityMetrics(e, intgroup=c("Zoox"), force=T, outdir= "report_for_genes_treat")
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
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 1654, 25% 
#LFC < 0 (down)   : 1861, 28% 
#outliers [1]     : 1, 0.015% 
#low counts [2]   : 0, 0% 
table(resHT$padj < 0.05)
#FALSE  TRUE 
# 3576  3060 
res=data.frame(cbind("gene"=row.names(resHT),"stat"= resHT$stat))
head(res)
write.csv(res,file="GO_symb/resJuves.heat_stat.csv",quote=F, row.names=F)


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
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 292, 4.4% 
#LFC < 0 (down)   : 311, 4.7% 
#outliers [1]     : 1, 0.015% 
#low counts [2]   : 0, 0%  
table(resBKBK.HA $padj < 0.05)
#FALSE  TRUE 
# 6271   365  
res=data.frame(cbind("gene"=row.names(resBKBK.HA),"stat"= resBKBK.HA$stat))
head(res)
write.csv(res,file="GO_symb/resBKBK.HA_stat.csv",quote=F, row.names=F)

# DR - Davies = Central
resDRDR.HA<-results(dds, contrast=c('group3', 'DRxDRHot', 'DRxDRAmbient')) #here is where the two contrasting conditions get defined
mcols(resDRDR.HA,use.names=TRUE)
summary(resDRDR.HA)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 522, 7.9% 
#LFC < 0 (down)   : 513, 7.7% 
#outliers [1]     : 1, 0.015% 
table(resDRDR.HA $padj < 0.05)
#FALSE  TRUE 
# 6010   626 
res=data.frame(cbind("gene"=row.names(resDRDR.HA),"stat"= resDRDR.HA$stat))
head(res)
write.csv(res,file="GO_symb/resDRDR.HA_stat.csv",quote=F, row.names=F)

# SB - Sand Bank 7 - North 
resSBSB.HA<-results(dds, contrast=c('group3', 'SBxSBHot', 'SBxSBAmbient')) #here is where the two contrasting conditions get defined
mcols(resSBSB.HA,use.names=TRUE)
summary(resSBSB.HA)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 312, 4.7% 
#LFC < 0 (down)   : 386, 5.8% 
#outliers [1]     : 1, 0.015% 
table(resSBSB.HA $padj < 0.05)
#FALSE  TRUE 
# 6187   449 
res=data.frame(cbind("gene"=row.names(resSBSB.HA),"stat"= resSBSB.HA $stat))
head(res)
write.csv(res,file="GO_symb/resSBSB.HA_stat.csv",quote=F, row.names=F)

# CU - Curd Reef - North 
# no hot juves

# LS - Long Sandy - North 
resLSLS.HA<-results(dds, contrast=c('group3', 'LSxLSHot', 'LSxLSAmbient')) #here is where the two contrasting conditions get defined
mcols(resLSLS.HA,use.names=TRUE)
summary(resLSLS.HA)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 251, 3.8% 
#LFC < 0 (down)   : 276, 4.2% 
#outliers [1]     : 1, 0.015% 
#low counts [2]   : 385, 5.8% 
table(resLSLS.HA $padj < 0.05)
#FALSE  TRUE 
# 5998   253  
res=data.frame(cbind("gene"=row.names(resLSLS.HA),"stat"= resLSLS.HA $stat))
head(res)
write.csv(res,file="GO_symb/resLSLS.HA_stat.csv",quote=F, row.names=F)

## North by North population crosses ##
## CUXSB
resCUSB.HA<-results(dds, contrast=c('group3', 'CUxSBHot', 'CUxSBAmbient')) #here is where the two contrasting conditions get defined
mcols(resCUSB.HA,use.names=TRUE)
summary(resCUSB.HA)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 5, 0.075% 
#LFC < 0 (down)   : 22, 0.33% 
#outliers [1]     : 1, 0.015% 
#low counts [2]   : 385, 5.8%  
table(resCUSB.HA $padj < 0.05)
#FALSE  TRUE 
#6242     9  
res=data.frame(cbind("gene"=row.names(resCUSB.HA),"stat"= resCUSB.HA $stat))
head(res)
write.csv(res,file="GO_symb/resCUSB.HA_stat.csv",quote=F, row.names=F)

## LSXCU
resLSCU.HA<-results(dds, contrast=c('group3', 'LSxCUHot', 'LSxCUAmbient')) #here is where the two contrasting conditions get defined
mcols(resLSCU.HA,use.names=TRUE)
summary(resLSCU.HA)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 62, 0.93% 
#LFC < 0 (down)   : 73, 1.1% 
#outliers [1]     : 1, 0.015% 
#low counts [2]   : 514, 7.7% 
table(resLSCU.HA $padj < 0.05)
#FALSE  TRUE 
#6049    73   
res=data.frame(cbind("gene"=row.names(resLSCU.HA),"stat"= resLSCU.HA $stat))
head(res)
write.csv(res,file="GO_symb/resLSCU.HA_stat.csv",quote=F, row.names=F)


## Central by Central population crosses ## - There are none! 

## Flip Flop Crosses- North Mom X Central Dad ##
##SBXBK
# no juves for this cross

##CUXDR
resCUDR.HA<-results(dds, contrast=c('group3', 'CUxDRHot', 'CUxDRAmbient')) #here is where the two contrasting conditions get defined
mcols(resCUDR.HA,use.names=TRUE)
summary(resCUDR.HA)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 9, 0.14% 
#LFC < 0 (down)   : 28, 0.42% 
#outliers [1]     : 1, 0.015% 
table(resCUDR.HA$padj < 0.05)
#FALSE  TRUE 
# 6630     6  
res=data.frame(cbind("gene"=row.names(resCUDR.HA),"stat"= resCUDR.HA$stat))
head(res)
write.csv(res,file="GO_symb/resCUDR.HA_stat.csv",quote=F, row.names=F)


## Flip Flop Crosses- Central Mom X North Dad ##
##DRXSB
resDRSB.HA<-results(dds, contrast=c('group3', 'DRxSBHot', 'DRxSBAmbient')) #here is where the two contrasting conditions get defined
mcols(resDRSB.HA,use.names=TRUE)
summary(resDRSB.HA)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 290, 4.4% 
#LFC < 0 (down)   : 368, 5.5% 
#outliers [1]     : 1, 0.015% 
#low counts [2]   : 0, 0% 
table(resDRSB.HA $padj < 0.05)
#FALSE  TRUE 
# 6277   359   
res=data.frame(cbind("gene"=row.names(resDRSB.HA),"stat"= resDRSB.HA $stat))
head(res)
write.csv(res,file="GO_symb/resDRSB.HA_stat.csv",quote=F, row.names=F)

save(vstG3,vstG3.df,colData,resBKBK.HA,resCUDR.HA,resCUSB.HA,resDRDR.HA,resDRSB.HA,resLSCU.HA,resLSLS.HA,resSBSB.HA,file="AGF_juvenile_2018_WaldGroup3Contrasts_crossTemp_macbook.Rdata")

#################################### WALD TEST - FULL MODEL - INTERACTIONS ###

dds<-DESeqDataSetFromMatrix(gcounts,
	colData = colData, 
	design = ~group2)

vstG2=vst(dds)
vstG2.df = assay(vstG3)
dds<-DESeq(dds, minReplicatesForReplace=Inf) 

#####SED
resSED.HA<-results(dds, contrast=c('group2', 'SEDHot', 'SEDAmbient')) #here is where the two contrasting conditions get defined
mcols(resSED.HA,use.names=TRUE)
summary(resSED.HA)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 1324, 20% 
#LFC < 0 (down)   : 1417, 21% 
#outliers [1]     : 3, 0.045% 
table(resSED.HA$padj < 0.05)
#FALSE  TRUE 
# 4422  2212    
res=data.frame(cbind("gene"=row.names(resSED.HA),"stat"= resSED.HA$stat))
head(res)
write.csv(res,file="GO_symb/resSED.HA_stat.csv",quote=F, row.names=F)

#####SS
resSS.HA<-results(dds, contrast=c('group2', 'SSHot', 'SSAmbient')) #here is where the two contrasting conditions get defined
mcols(resSS.HA,use.names=TRUE)
summary(resSS.HA)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 702, 11% 
#LFC < 0 (down)   : 833, 13% 
#outliers [1]     : 3, 0.045% 
#low counts [2]   : 0, 0% 
table(resSS.HA$padj < 0.05)
#FALSE  TRUE 
# 5519  1115   
res=data.frame(cbind("gene"=row.names(resSS.HA),"stat"= resSS.HA$stat))
head(res)
write.csv(res,file="GO_symb/resSS.HA_stat.csv",quote=F, row.names=F)

#####D1
resD1.HA<-results(dds, contrast=c('group2', 'D1Hot', 'D1Ambient')) #here is where the two contrasting conditions get defined
mcols(resD1.HA,use.names=TRUE)
summary(resD1.HA)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 727, 11% 
#LFC < 0 (down)   : 838, 13% 
#outliers [1]     : 3, 0.045% 
table(resD1.HA$padj < 0.05)
#FALSE  TRUE 
# 5511  1123    
res=data.frame(cbind("gene"=row.names(resD1.HA),"stat"= resD1.HA$stat))
head(res)
write.csv(res,file="GO_symb/resD1.HA_stat.csv",quote=F, row.names=F)

save(vstG2,vstG2.df,colData,resSED.HA,resSS.HA,resD1.HA,file="AGF_juvenile_2018_WaldGroup2Contrasts_zooxTemp_macbook.Rdata")

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
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 455, 6.9% 
#LFC < 0 (down)   : 407, 6.1% 
#outliers [1]     : 4, 0.06% 
#low counts [2]   : 0, 0% 
table(resCD $padj < 0.05)
#FALSE  TRUE 
# 6094   539    
res=data.frame(cbind("gene"=row.names(resCD),"stat"= resCD$stat))
head(res)
write.csv(res,file="GO_symb/resCD_stat.csv",quote=F, row.names=F)

resCSED<-results(dds, contrast=c('Zoox', 'C1', 'SED')) #here is where the two contrasting conditions get defined
mcols(resCSED,use.names=TRUE)
summary(resCSED)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 1385, 21% 
#LFC < 0 (down)   : 1351, 20% 
#outliers [1]     : 4, 0.06% 
table(resCSED$padj < 0.05)
#FALSE  TRUE 
# 4329  2304 
res=data.frame(cbind("gene"=row.names(resCSED),"stat"= resCSED$stat))
head(res)
write.csv(res,file="GO_symb/resCSED_stat.csv",quote=F, row.names=F)

resCSS<-results(dds, contrast=c('Zoox', 'C1', 'SS')) #here is where the two contrasting conditions get defined
mcols(resCSS,use.names=TRUE)
summary(resCSS)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 529, 8% 
#LFC < 0 (down)   : 470, 7.1% 
#outliers [1]     : 4, 0.06% 
#low counts [2]   : 0, 0% 
table(resCSS$padj < 0.05)
#FALSE  TRUE 
#5918   715
res=data.frame(cbind("gene"=row.names(resCSS),"stat"= resCSS$stat))
head(res)
write.csv(res,file="GO_symb/resCSS_stat.csv",quote=F, row.names=F)

resDSED<-results(dds, contrast=c('Zoox', 'D1', 'SED')) #here is where the two contrasting conditions get defined
mcols(resDSED,use.names=TRUE)
summary(resDSED)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 694, 10% 
#LFC < 0 (down)   : 597, 9% 
#outliers [1]     : 4, 0.06% 
table(resDSED$padj < 0.05)
#FALSE  TRUE 
#5668   965 
res=data.frame(cbind("gene"=row.names(resDSED),"stat"= resDSED$stat))
head(res)
write.csv(res,file="GO_symb/resDSED_stat.csv",quote=F, row.names=F) 

resDSS<-results(dds, contrast=c('Zoox', 'D1', 'SS')) #here is where the two contrasting conditions get defined
mcols(resDSS,use.names=TRUE)
summary(resDSS)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 129, 1.9% 
#LFC < 0 (down)   : 126, 1.9% 
#outliers [1]     : 4, 0.06% 
table(resDSS$padj < 0.05)
#FALSE  TRUE 
# 6440   193
res=data.frame(cbind("gene"=row.names(resDSS),"stat"= resDSS$stat))
head(res)
write.csv(res,file="GO_symb/resDSS_stat.csv",quote=F, row.names=F)

resSSSED<-results(dds, contrast=c('Zoox', 'SS', 'SED')) #here is where the two contrasting conditions get defined
mcols(resSSSED,use.names=TRUE)
summary(resSSSED)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 862, 13% 
#LFC < 0 (down)   : 814, 12% 
#outliers [1]     : 4, 0.06% 
table(resSSSED$padj < 0.05)
#FALSE  TRUE 
# 5326  1307
res=data.frame(cbind("gene"=row.names(resSSSED),"stat"= resSSSED$stat))
head(res)
write.csv(res,file="GO_symb/resSSSED_stat.csv",quote=F, row.names=F)

############# run ambo juve samples for each sym type separately for KOG_MWU
dds<-DESeqDataSetFromMatrix(gcountsA,
	colData = colDataA, 
	design = ~C1)

dds<-DESeq(dds, minReplicatesForReplace=Inf) 

resC1all<-results(dds, contrast=c('C1', '0', '1')) #here is where the two contrasting conditions get defined
mcols(resC1all,use.names=TRUE)
summary(resC1all)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 932, 14% 
#LFC < 0 (down)   : 967, 15% 
#outliers [1]     : 6, 0.09% 

### D1
dds<-DESeqDataSetFromMatrix(gcountsA,
	colData = colDataA, 
	design = ~D1)

dds<-DESeq(dds, minReplicatesForReplace=Inf) 

resD1all<-results(dds, contrast=c('D1', '0', '1')) #here is where the two contrasting conditions get defined
mcols(resD1all,use.names=TRUE)
summary(resD1all)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 88, 1.3% 
#LFC < 0 (down)   : 73, 1.1% 

### SS
dds<-DESeqDataSetFromMatrix(gcountsA,
	colData = colDataA, 
	design = ~SS)

dds<-DESeq(dds, minReplicatesForReplace=Inf) 

resSSall<-results(dds, contrast=c('SS', '0', '1')) #here is where the two contrasting conditions get defined
mcols(resSSall,use.names=TRUE)
summary(resSSall)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 163, 2.5% 
#LFC < 0 (down)   : 210, 3.2% 
#outliers [1]     : 6, 0.09% 


### SD
dds<-DESeqDataSetFromMatrix(gcountsA,
	colData = colDataA, 
	design = ~SD)

dds<-DESeq(dds, minReplicatesForReplace=Inf) 

resSDall<-results(dds, contrast=c('SD', '0', '1')) #here is where the two contrasting conditions get defined
mcols(resSDall,use.names=TRUE)
summary(resSDall)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 1227, 18% 
#LFC < 0 (down)   : 1124, 17% 
#outliers [1]     : 6, 0.09% 


save(vstZ.df, resCD, resCSED, resCSS, resDSED, resDSS, resSSSED, resSSall, resSDall, resD1all, resC1all, colData,file="AGF_juvenile_2018_WaldZooxAmbient.Rdata")


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
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 22, 0.33% 
#LFC < 0 (down)   : 79, 1.2% 
#outliers [1]     : 5, 0.075% 
table(resSSDhot$padj < 0.05)
#FALSE  TRUE 
# 6555    77   
res=data.frame(cbind("gene"=row.names(resSSDhot),"stat"= resSSDhot$stat))
head(res)
write.csv(res,file="GO_symb/resSSDhot_stat.csv",quote=F, row.names=F)

resSSSEDhot<-results(dds, contrast=c('Zoox', 'SS', 'SED')) #here is where the two contrasting conditions get defined
mcols(resSSSEDhot,use.names=TRUE)
summary(resSSSEDhot)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 20, 0.3% 
#LFC < 0 (down)   : 28, 0.42% 
#outliers [1]     : 5, 0.075% 
table(resSSSEDhot$padj < 0.05)
#FALSE  TRUE 
# 6605    27   
res=data.frame(cbind("gene"=row.names(resSSSEDhot),"stat"= resSSSEDhot$stat))
head(res)
write.csv(res,file="GO_symb/resSSSEDhot_stat.csv",quote=F, row.names=F)

resD1SEDhot<-results(dds, contrast=c('Zoox', 'D1', 'SED')) #here is where the two contrasting conditions get defined
mcols(resD1SEDhot,use.names=TRUE)
summary(resD1SEDhot)
#out of 6637 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 105, 1.6% 
#LFC < 0 (down)   : 69, 1% 
#outliers [1]     : 5, 0.075% 
table(resD1SEDhot $padj < 0.05)
#FALSE  TRUE 
# 6518   114   
res=data.frame(cbind("gene"=row.names(resD1SEDhot),"stat"= resD1SEDhot $stat))
head(res)
write.csv(res,file="GO_symb/resD1SEDhot_stat.csv",quote=F, row.names=F)


################ KOG_MWU #######################
library(KOGMWU)

###### look at effects of symbiosis on juveniles in ambient conditions
ll=load("AGF_juvenile_2018_WaldZooxAmbient.Rdata") 
gene2kog=read.table("Amil_gene2kogClass.tab",header=T,sep="\t",quote=NULL,fill=T)

C_LFC=as.data.frame((cbind(row.names(resC1all), resC1all $log2FoldChange))) #collect pvalues for each gene
C_LFC$V1=sub("-RA","", C_LFC$V1)
C=kog.mwu(C_LFC,gene2kog,Alternative="t")

D_LFC=as.data.frame((cbind(row.names(resD1all), resD1all $log2FoldChange))) #collect pvalues for each gene
D_LFC$V1=sub("-RA","", D_LFC$V1)
D=kog.mwu(D_LFC,gene2kog,Alternative="t")

SED_LFC=as.data.frame((cbind(row.names(resSDall), resSDall $log2FoldChange))) #collect pvalues for each gene
SED_LFC $V1=sub("-RA","", SED_LFC $V1)
SED=kog.mwu(SED_LFC,gene2kog,Alternative="t")

SS_LFC=as.data.frame((cbind(row.names(resSSall), resSSall$log2FoldChange))) #collect pvalues for each gene
SS_LFC $V1=sub("-RA","", SS_LFC $V1)
SS=kog.mwu(SS_LFC,gene2kog,Alternative="t")

ktable=makeDeltaRanksTable(list("C1"=C,"D1"=D,"SS"=SS,"SED"=SED))
#remove nuclear structure because only has 4 genes and values are outliers
ktable <- ktable[-c(8), ]

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
ktests$zoox=(c(rep("C",23),rep("D",23),rep("SS",23),rep("SED",23)))
write.csv(ktests,file="kogTestResults_ambo.csv",quote=F, row.names=F)

####### look at effects of symbiosis on response to heat in juveniles 
ll=load("AGF_juvenile_2018_WaldGroup2Contrasts_zooxTemp_macbook.Rdata") 
gene2kog=read.table("Amil_gene2kogClass.tab",header=T,sep="\t",quote=NULL,fill=T)

D_HA=as.data.frame((cbind(row.names(resD1.HA), resD1.HA$log2FoldChange))) #collect pvalues for each gene
D_HA$V1=sub("-RA","", D_HA$V1)
D_HA =kog.mwu(D_HA,gene2kog,Alternative="t")

SS_HA=as.data.frame((cbind(row.names(resSS.HA), resSS.HA$log2FoldChange))) #collect pvalues for each gene
SS_HA$V1=sub("-RA","", SS_HA$V1)
SS_HA=kog.mwu(SS_HA,gene2kog,Alternative="t")

SED_HA=as.data.frame((cbind(row.names(resSED.HA), resSED.HA$log2FoldChange))) #collect pvalues for each gene
SED_HA$V1=sub("-RA","", SED_HA$V1)
SED_HA=kog.mwu(SED_HA,gene2kog,Alternative="t")

#### load in larval response to heat stress (overall)
setwd("/Users/mes0192/Dropbox/KQ_GBR_spawnseq/analysis/DESeq2_larvae")
ll=load("AGF_larval_2018_WaldCrossTreat.Rdata")

Larva_HA=as.data.frame((cbind(row.names(resHT), resHT$log2FoldChange))) #collect pvalues for each gene
Larva_HA $V1=sub("-RA","", Larva_HA$V1)
Larva_HA =kog.mwu(Larva_HA,gene2kog,Alternative="t")

#### load in larval thermal tolerance to heat stress (overall)
setwd("/Users/mes0192/Dropbox/KQ_GBR_spawnseq/analysis/DESeq2_larvae")
ll=load("AGF_larval_2018_PreRes.Rdata")

Larva_Tolerance=as.data.frame((cbind(row.names(resP), resP$log2FoldChange))) #collect pvalues for each gene
Larva_Tolerance $V1=sub("-RA","", Larva_Tolerance $V1)
Larva_Tolerance =kog.mwu(Larva_Tolerance,gene2kog,Alternative="t")


ktable=makeDeltaRanksTable(list("D1"=D_HA,"SS"=SS_HA,"SED"=SED_HA,"Larval Response"=Larva_HA, "Larval Tolerance"=Larva_Tolerance))
#remove nuclear structure because only has 4 genes and values are outliers
ktable <- ktable[-c(17), ]

heat.colors = colorRampPalette(wes_palette("Zissou1", 30, type = "continuous"), bias=0.8)(100)

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

ktests=rbind(C,D,SS,SED,Larva_Tolerance,Larva_HA)
ktests$zoox=(c(rep("C",23),rep("D",23),rep("SS",23),rep("SED",23), rep("LT",23), rep("LR",23)))
write.csv(ktests,file="kogTestResults_heat.csv",quote=F, row.names=F)



### principal coordinate calculation on all samples
All=load("AGF_juvenile_2018_WaldCrossTreatZoox.Rdata")

#remove effect of cross
vstL <- limma::removeBatchEffect(vstCZT.df, colData$Cross)

conditions=colData

grp=rep("#EBCC2A",ncol(vstCZT.df))
grp[grep("SED",conditions$Zoox)]="#3B9AB2"
grp[grep("SS",conditions$Zoox)]="#F21A00"
grp[grep("D1",conditions$Zoox)]="purple"

grp=rep("#3B9AB2",ncol(vstCZT.df))
grp[grep("Hot",conditions$Treatment.1)]="#F21A00"

fitpc=capscale(dist(t(vstL),method="manhattan")~1)
summary(eigenvals(fitpc))

plot(fitpc$CA$u,pch=16, col="white",main="AGF Juveniles 2018 All Samples", xlab="MDS1 11.15%", ylab="MDS2 7.43%")
abline(h=0,v=0, col="grey")
points(fitpc$CA$u,pch=16, col=grp)
ordiellipse(fitpc$CA$u,conditions$group2,draw="polygon",label=T)
#ordispider(fitpc$CA$u,conditions$group2,draw="polygon",label=T)

ad=adonis(t(vstCZT.df)~Cross+Treatment.1+Zoox,data=conditions,method="manhattan")
#adonis(formula = t(vstCZT.df) ~ Cross + Treatment.1 + Zoox, data = conditions,      method = "manhattan") 
#
#Permutation: free
#Number of permutations: 999
#
#Terms added sequentially (first to last)
#
#             Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#Cross         9  64803500  7200389  2.4929 0.12597  0.001 ***
#Treatment.1   1  27720286 27720286  9.5973 0.05388  0.001 ***
#Zoox          3  40667366 13555789  4.6933 0.07905  0.001 ***
#Residuals   132 381261186  2888342         0.74110           
#Total       145 514452338                  1.00000           
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
        
labs=c("Cross","Treatment","Zoox","Residuals")
cols=c("#EBCC2A","#3B9AB2","#F21A00","grey80")
#cols = append(gg_color_hue(3), 'grey')
labs2 = paste(labs, round(ad$aov.tab$R2[1:4]*100, digits=1))
pie(ad$aov.tab$R2[1:4],labels=labs2,col=cols,main="AGF Juveniles 2018")

##### PCA for ambo
Amb=load("AGF_juvenile_2018_WaldZooxAmbient.Rdata")
colDataA=colData[colData$Treatment.1=="Ambient",]

vstL <- limma::removeBatchEffect(vstZ.df, colDataA$Cross)

conditions=colDataA

grp=rep("#EBCC2A",ncol(vstL))
grp[grep("SED",conditions$Zoox)]="#3B9AB2"
grp[grep("SS",conditions$Zoox)]="#F21A00"
grp[grep("D1",conditions$Zoox)]="purple"

fitpc=capscale(dist(t(vstL),method="manhattan")~1)
summary(eigenvals(fitpc))

plot(fitpc$CA$u,pch=16, col="white",main="AGF Juveniles 2018 Ambient Samples", xlab="MDS1 11.5%", ylab="MDS2 6.20%")
abline(h=0,v=0, col="grey")
points(fitpc$CA$u,pch=16, col=grp)
ordiellipse(fitpc$CA$u,conditions$Zoox,draw="polygon",label=T)

ad=adonis(t(vstZ.df)~Cross+Zoox,data=conditions,method="manhattan")

#adonis(formula = t(vstZ.df) ~ Cross + Zoox, data = conditions,      method = "manhattan") 
#
#Permutation: free
#Number of permutations: 999
#
#Terms added sequentially (first to last)
#
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#Cross       9  51759932  5751104  2.1915 0.14039  0.001 ***
#Zoox        3  38744537 12914846  4.9214 0.10509  0.001 ***
#Residuals 106 278170097  2624246         0.75451           
#Total     118 368674566                  1.00000           
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##### PCA for heat
All=load("AGF_juvenile_2018_WaldCrossTreatZoox.Rdata")
colDataH=colData[colData$Treatment.1=="Hot",]

all(rownames(colDataH) %in% colnames(vstCZT.df)) #TRUE
vstCZT.dfH <- vstCZT.df[, rownames(colDataH)]
all(rownames(colDataH) == colnames(vstCZT.dfH))

vstL <- limma::removeBatchEffect(vstCZT.dfH, colDataH$Cross)

conditions=colDataH

grp=rep("#EBCC2A",ncol(vstL))
grp[grep("SED",conditions$Zoox)]="#3B9AB2"
grp[grep("SS",conditions$Zoox)]="#F21A00"
grp[grep("D1",conditions$Zoox)]="purple"

fitpc=capscale(dist(t(vstCZT.dfH),method="manhattan")~1)
summary(eigenvals(fitpc))

plot(fitpc$CA$u,pch=16, col="white",main="AGF Juveniles 2018 Heat Samples", xlab="MDS1 13.41%", ylab="MDS2 9.57%")
abline(h=0,v=0, col="grey")
points(fitpc$CA$u,pch=16, col=grp)
ordiellipse(fitpc$CA$u,conditions$Zoox,draw="polygon",label=T)

ad=adonis(t(vstCZT.dfH)~Cross+Zoox,data=conditions,method="manhattan")
#adonis(formula = t(vstCZT.dfH) ~ Cross + Zoox, data = conditions,      method = "manhattan") 
#
#Permutation: free
#Number of permutations: 999
#
#Terms added sequentially (first to last)
#
#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Cross      7  25472491 3638927  1.6237 0.36077  0.001 ***
#Zoox       2   7035026 3517513  1.5695 0.09964  0.003 ** 
#Residuals 17  38099058 2241121         0.53960           
#Total     26  70606575                 1.00000           
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



