
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
library(gridExtra)
library(wesanderson)

gcountsALL=read.table("AGF_larvae_2018_geneCounts_aten.txt", header=T) 

# clean up file names
colnames(gcountsALL)<-sub("X","", colnames(gcountsALL))
colnames(gcountsALL) <- sub("_[^_]+$", "", colnames(gcountsALL))
colnames(gcountsALL)<-gsub("\\.","_", colnames(gcountsALL))

length(gcountsALL[,1]) #24924
dim(gcountsALL) #24924    96


summary(colSums(gcountsALL))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  135519 1345494 1463876 1444555 1552249 1859418  


### REMOVE GENES WITH LOW MEAN COUNTS ###

mns = apply(gcountsALL, 1, mean)

gcounts=gcountsALL[mns>10,] #get rid of genes that show little or no expression
#gcounts=gcounts[mns>10,] #only highly expressed genes for WGCNA

table(mns > 10)
#FALSE  TRUE 
#13732 11192 
dim(gcounts) #11192    96


### BUILD A DATAFRAME ASSOCIATING SAMPLE NAMESWITH TREATMENT CONDITIONS ###
colData=read.csv("J19188meta.csv", header=T)
colData=colData[!duplicated(colData$sample), ]
rownames(colData)<-colData$sample
dim(colData) #96  3
colData$group=factor(paste0(colData$cross,colData$treatment))
#colData <- colData[order(colData$treatment),]
#colData$ind.n=factor(rep(1:3, each=1))
colData$dam=str_sub(colData$cross, 1,2)
colData$sire = str_sub(colData$cross,4,5)
#colData$ind=factor(rep(1:32, each=3))

all(rownames(colData) %in% colnames(gcounts)) #TRUE
gcounts <- gcounts[, rownames(colData)]
all(rownames(colData) == colnames(gcounts))

colData$treatment=as.factor(colData$treatment)
colData$cross=as.factor(colData$cross)
### OUTLIERS ### - not going to remove any 

#dds<-DESeqDataSetFromMatrix(countData=gcounts, colData=colData, design= ~ treatment+cross)
#vsd=varianceStabilizingTransformation(dds, blind=T)
#e=ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
#arrayQualityMetrics(e, intgroup=c("treatment"), force=T, outdir= "report_for_genes_treat")
#arrayQualityMetrics(e, intgroup=c("cross"), force=T, outdir= "report_for_genes_cross")

### WALD TEST - FULL MODEL ###

ddsCT<-DESeqDataSetFromMatrix(gcounts,

	colData = colData, 

	design = ~cross+treatment)


vstCT=vst(ddsCT)
vstCT.df = assay(vstCT)

### Wald test for treatment in larvae ###

dds<-DESeq(ddsCT, minReplicatesForReplace=Inf) 
resHT<-results(dds, contrast=c('treatment', 'Hot', 'Ambient')) #here is where the two contrasting conditions get defined
mcols(resHT,use.names=TRUE)
summary(resHT)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 3798, 34%
#LFC < 0 (down)     : 4103, 37%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%

table(resHT$padj < 0.05)
#FALSE  TRUE 
#  3811  7381 
res=data.frame(cbind("gene"=row.names(resHT),"stat"= resHT$stat))
head(res)
write.csv(res,file="resHT_stat.csv",quote=F, row.names=F)


resHP<-results(dds, contrast=c('treatment', 'Hot', 'Pre')) #here is where the two contrasting conditions get defined
mcols(resHP,use.names=TRUE)
summary(resHP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 3726, 33%
#LFC < 0 (down)     : 3994, 36%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%

table(resHP$padj < 0.05)
#FALSE  TRUE 
# 4038  7154

resAP<-results(dds, contrast=c('treatment', 'Ambient', 'Pre')) #here is where the two contrasting conditions get defined
mcols(resAP,use.names=TRUE)
summary(resAP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2570, 23%
#LFC < 0 (down)     : 2568, 23%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%

table(resAP$padj < 0.05)
#FALSE  TRUE 
# 6769  4423 

save(vstCT, vstCT.df,resHT,resHP,resAP, colData,file="AGF_larval_2018_WaldCrossTreat.Rdata")


### LRT TEST - FULL MODEL ###
dds<-DESeqDataSetFromMatrix(gcounts,

	colData = colData, 

	design = formula(~ treatment+dam+sire))

vstLRT=vst(dds)
vstLRT.df = assay(vstLRT)

## LRT for treatment
dds <- DESeq(dds, test="LRT", reduced=~dam+sire)
resT <- results(dds)
summary(resT)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 4456, 40%
#LFC < 0 (down)     : 4630, 41%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
table(resT$padj < 0.05)
#FALSE  TRUE 
# 2580  8612

## LRT for dam
dds <- DESeq(dds, test="LRT", reduced=~treatment+sire)
resD <- results(dds)
summary(resD)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 3968, 35%
#LFC < 0 (down)     : 4167, 37%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%

table(resD$padj < 0.05)
#FALSE  TRUE 
#  3877  7315  

## LRT for sire
dds <- DESeq(dds, test="LRT", reduced=~treatment+dam)
resS <- results(dds)
summary(resS)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 3559, 32%
#LFC < 0 (down)     : 3970, 35%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%

table(resS$padj < 0.05)
#FALSE  TRUE 
#4495  6697 

############################################ WALD TEST - FULL MODEL - INTERACTIONS ###

dds<-DESeqDataSetFromMatrix(gcounts,

	colData = colData, 

	design = ~group)


vstG=vst(dds)
vstG.df = assay(vstG)
dds<-DESeq(dds, minReplicatesForReplace=Inf) 


### cross contrasts ###
# population by population 
# BKBK - Backnumbers = Central
resBKBK.AP<-results(dds, contrast=c('group', 'BKxBKAmbient', 'BKxBKPre')) #here is where the two contrasting conditions get defined
mcols(resBKBK.AP,use.names=TRUE)
summary(resBKBK.AP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2339, 21%
#LFC < 0 (down)     : 2090, 19%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%

table(resBKBK.AP$padj < 0.05)
#FALSE  TRUE 
# 7494  3698 
res=data.frame(cbind("gene"=row.names(resBKBK.AP),"stat"=resBKBK.AP$stat))
head(res)
write.csv(res,file="resBKBK.AP_stat.csv",quote=F, row.names=F)


resBKBK.HP<-results(dds, contrast=c('group', 'BKxBKHot', 'BKxBKPre')) #here is where the two contrasting conditions get defined
mcols(resBKBK.HP,use.names=TRUE)
summary(resBKBK.HP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2932, 26%
#LFC < 0 (down)     : 3026, 27%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%

table(resBKBK.HP $padj < 0.05)
#FALSE  TRUE 
# 6011  5181 
res=data.frame(cbind("gene"=row.names(resBKBK.HP),"stat"= resBKBK.HP $stat))
head(res)
write.csv(res,file="resBKBK.HP_stat.csv",quote=F, row.names=F)

resBKBK.HA<-results(dds, contrast=c('group', 'BKxBKHot', 'BKxBKAmbient')) #here is where the two contrasting conditions get defined
mcols(resBKBK.HA,use.names=TRUE)
summary(resBKBK.HA)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2210, 20%
#LFC < 0 (down)     : 2290, 20%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%

table(resBKBK.HA $padj < 0.05)
#FALSE  TRUE 
# 7329  3863  
res=data.frame(cbind("gene"=row.names(resBKBK.HA),"stat"= resBKBK.HA $stat))
head(res)
write.csv(res,file="resBKBK.HA_stat.csv",quote=F, row.names=F)
resBKBK.HA$lp=-log(resBKBK.HA$pvalue,10)
resBKBK.HA$lp[resBKBK.HA$log2FoldChange<0]=-resBKBK.HA$lp[resBKBK.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resBKBK.HA),"lp"=resBKBK.HA$lp))
write.csv(dev_lpv,file="resBKBK.HA_lp.csv",quote=F, row.names=F)

# DR - Davies = Central
resDRDR.AP<-results(dds, contrast=c('group', 'DRxDRAmbient', 'DRxDRPre')) #here is where the two contrasting conditions get defined
mcols(resDRDR.AP,use.names=TRUE)
summary(resDRDR.AP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2175, 19%
#LFC < 0 (down)     : 2419, 22%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
table(resDRDR.AP $padj < 0.05)
#FALSE  TRUE 
# 7198  3994 
res=data.frame(cbind("gene"=row.names(resDRDR.AP),"stat"= resDRDR.AP $stat))
head(res)
write.csv(res,file="resDRDR.AP_stat.csv",quote=F, row.names=F)


resDRDR.HP<-results(dds, contrast=c('group', 'DRxDRHot', 'DRxDRPre')) #here is where the two contrasting conditions get defined
mcols(resDRDR.HP,use.names=TRUE)
summary(resDRDR.HP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2378, 21%
#LFC < 0 (down)     : 2457, 22%

table(resDRDR.HP $padj < 0.05)
#FALSE  TRUE 
# 7095  4097 
res=data.frame(cbind("gene"=row.names(resDRDR.HP),"stat"= resDRDR.HP $stat))
head(res)
write.csv(res,file="resDRDR.HP_stat.csv",quote=F, row.names=F)


resDRDR.HA<-results(dds, contrast=c('group', 'DRxDRHot', 'DRxDRAmbient')) #here is where the two contrasting conditions get defined
mcols(resDRDR.HA,use.names=TRUE)
summary(resDRDR.HA)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2291, 20%
#LFC < 0 (down)     : 2246, 20%

table(resDRDR.HA $padj < 0.05)
#FALSE  TRUE 
# 7275  3917  
res=data.frame(cbind("gene"=row.names(resDRDR.HA),"stat"= resDRDR.HA $stat))
head(res)
write.csv(res,file="resDRDR.HA_stat.csv",quote=F, row.names=F)
resDRDR.HA$lp=-log(resDRDR.HA$pvalue,10)
resDRDR.HA$lp[resDRDR.HA$log2FoldChange<0]=-resDRDR.HA$lp[resDRDR.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resDRDR.HA),"lp"= resDRDR.HA$lp))
write.csv(dev_lpv,file="resDRDR.HA_lp.csv",quote=F, row.names=F)


# SB - Sand Bank 7 - North 
resSBSB.AP<-results(dds, contrast=c('group', 'SBxSBAmbient', 'SBxSBPre')) #here is where the two contrasting conditions get defined
mcols(resSBSB.AP,use.names=TRUE)
summary(resSBSB.AP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 859, 7.7%
#LFC < 0 (down)     : 952, 8.5%

table(resSBSB.AP $padj < 0.05)
#FALSE  TRUE 
# 9742  1450 
res=data.frame(cbind("gene"=row.names(resSBSB.AP),"stat"= resSBSB.AP $stat))
head(res)
write.csv(res,file="resSBSB.AP_stat.csv",quote=F, row.names=F)


resSBSB.HP<-results(dds, contrast=c('group', 'SBxSBHot', 'SBxSBPre')) #here is where the two contrasting conditions get defined
mcols(resSBSB.HP,use.names=TRUE)
summary(resSBSB.HP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2651, 24%
#LFC < 0 (down)     : 2645, 24%

table(resSBSB.HP $padj < 0.05)
#FALSE  TRUE 
# 6586  4606 
res=data.frame(cbind("gene"=row.names(resSBSB.HP),"stat"= resSBSB.HP $stat))
head(res)
write.csv(res,file="resSBSB.HP_stat.csv",quote=F, row.names=F)


resSBSB.HA<-results(dds, contrast=c('group', 'SBxSBHot', 'SBxSBAmbient')) #here is where the two contrasting conditions get defined
mcols(resSBSB.HA,use.names=TRUE)
summary(resSBSB.HA)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2626, 23%
#LFC < 0 (down)     : 2658, 24%

table(resSBSB.HA $padj < 0.05)
#FALSE  TRUE 
# 6612  4580 
res=data.frame(cbind("gene"=row.names(resSBSB.HA),"stat"= resSBSB.HA $stat))
head(res)
write.csv(res,file="resSBSB.HA_stat.csv",quote=F, row.names=F)
resSBSB.HA$lp=-log(resSBSB.HA$pvalue,10)
resSBSB.HA$lp[resSBSB.HA$log2FoldChange<0]=-resSBSB.HA$lp[resSBSB.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resSBSB.HA),"lp"= resSBSB.HA$lp))
write.csv(dev_lpv,file="resSBSB.HA_lp.csv",quote=F, row.names=F)

# CU - Curd Reef - North 
resCUCU.AP<-results(dds, contrast=c('group', 'CUxCUAmbient', 'CUxCUPre')) #here is where the two contrasting conditions get defined
mcols(resCUCU.AP,use.names=TRUE)
summary(resCUCU.AP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 472, 4.2%
#LFC < 0 (down)     : 846, 7.6%
table(resCUCU.AP $padj < 0.05)
#FALSE  TRUE 
#10185  1007 
res=data.frame(cbind("gene"=row.names(resCUCU.AP),"stat"= resCUCU.AP $stat))
head(res)
write.csv(res,file="resCUCU.AP_stat.csv",quote=F, row.names=F)



resCUCU.HP<-results(dds, contrast=c('group', 'CUxCUHot', 'CUxCUPre')) #here is where the two contrasting conditions get defined
mcols(resCUCU.HP,use.names=TRUE)
summary(resCUCU.HP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2334, 21%
#LFC < 0 (down)     : 2417, 22%

table(resCUCU.HP $padj < 0.05)
#FALSE  TRUE 
# 7114  4078  
res=data.frame(cbind("gene"=row.names(resCUCU.HP),"stat"= resCUCU.HP $stat))
head(res)
write.csv(res,file="resCUCU.HP_stat.csv",quote=F, row.names=F)


resCUCU.HA<-results(dds, contrast=c('group', 'CUxCUHot', 'CUxCUAmbient')) #here is where the two contrasting conditions get defined
mcols(resCUCU.HA,use.names=TRUE)
summary(resCUCU.HA)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1893, 17%
#LFC < 0 (down)     : 1983, 18%

table(resCUCU.HA $padj < 0.05)
#FALSE  TRUE 
# 7995  3197   
res=data.frame(cbind("gene"=row.names(resCUCU.HA),"stat"= resCUCU.HA $stat))
head(res)
write.csv(res,file="resCUCU.HA_stat.csv",quote=F, row.names=F)
resCUCU.HA$lp=-log(resCUCU.HA$pvalue,10)
resCUCU.HA$lp[resCUCU.HA$log2FoldChange<0]=-resCUCU.HA$lp[resCUCU.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resCUCU.HA),"lp"= resCUCU.HA$lp))
write.csv(dev_lpv,file="resCUCU.HA_lp.csv",quote=F, row.names=F)


# LS - Long Sandy - North 
resLSLS.AP<-results(dds, contrast=c('group', 'LSxLSAmbient', 'LSxLSPre')) #here is where the two contrasting conditions get defined
mcols(resLSLS.AP,use.names=TRUE)
summary(resLSLS.AP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 371, 3.3%
#LFC < 0 (down)     : 538, 4.8%

table(resLSLS.AP $padj < 0.05)
#FALSE  TRUE 
#10560   632  
res=data.frame(cbind("gene"=row.names(resLSLS.AP),"stat"= resLSLS.AP $stat))
head(res)
write.csv(res,file="resLSLS.AP_stat.csv",quote=F, row.names=F)


resLSLS.HP<-results(dds, contrast=c('group', 'LSxLSHot', 'LSxLSPre')) #here is where the two contrasting conditions get defined
mcols(resLSLS.HP,use.names=TRUE)
summary(resLSLS.HP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1925, 17%
#LFC < 0 (down)     : 2056, 18%

table(resLSLS.HP $padj < 0.05)
#FALSE  TRUE 
#7866  3326  
res=data.frame(cbind("gene"=row.names(resLSLS.HP),"stat"= resLSLS.HP $stat))
head(res)
write.csv(res,file="resLSLS.HP_stat.csv",quote=F, row.names=F)


resLSLS.HA<-results(dds, contrast=c('group', 'LSxLSHot', 'LSxLSAmbient')) #here is where the two contrasting conditions get defined
mcols(resLSLS.HA,use.names=TRUE)
summary(resLSLS.HA)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1941, 17%
#LFC < 0 (down)     : 2086, 19%

table(resLSLS.HA $padj < 0.05)
#FALSE  TRUE 
# 7858  3334  
res=data.frame(cbind("gene"=row.names(resLSLS.HA),"stat"= resLSLS.HA $stat))
head(res)
write.csv(res,file="resLSLS.HA_stat.csv",quote=F, row.names=F)
resLSLS.HA$lp=-log(resLSLS.HA$pvalue,10)
resLSLS.HA$lp[resLSLS.HA$log2FoldChange<0]=-resLSLS.HA$lp[resLSLS.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resLSLS.HA),"lp"= resLSLS.HA$lp))
write.csv(dev_lpv,file="resLSLS.HA_lp.csv",quote=F, row.names=F)


## North by North population crosses ##
## CUXSB
resCUSB.AP<-results(dds, contrast=c('group', 'CUxSBAmbient', 'CUxSBPre')) #here is where the two contrasting conditions get defined
mcols(resCUSB.AP,use.names=TRUE)
summary(resCUSB.AP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 410, 3.7%
#LFC < 0 (down)     : 618, 5.5%

table(resCUSB.AP $padj < 0.05)
#FALSE  TRUE 
#10419   773  
res=data.frame(cbind("gene"=row.names(resCUSB.AP),"stat"= resCUSB.AP $stat))
head(res)
write.csv(res,file="resCUSB.AP_stat.csv",quote=F, row.names=F)


resCUSB.HP<-results(dds, contrast=c('group', 'CUxSBHot', 'CUxSBPre')) #here is where the two contrasting conditions get defined
mcols(resCUSB.HP,use.names=TRUE)
summary(resCUSB.HP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2325, 21%
#LFC < 0 (down)     : 2402, 21%

table(resCUSB.HP $padj < 0.05)
#FALSE  TRUE 
# 7073  4119   
res=data.frame(cbind("gene"=row.names(resCUSB.HP),"stat"= resCUSB.HP $stat))
head(res)
write.csv(res,file="resCUSB.HP_stat.csv",quote=F, row.names=F)


resCUSB.HA<-results(dds, contrast=c('group', 'CUxSBHot', 'CUxSBAmbient')) #here is where the two contrasting conditions get defined
mcols(resCUSB.HA,use.names=TRUE)
summary(resCUSB.HA)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2250, 20%
#LFC < 0 (down)     : 2346, 21%

table(resCUSB.HA $padj < 0.05)
#FALSE  TRUE 
# 7255  3937  
res=data.frame(cbind("gene"=row.names(resCUSB.HA),"stat"= resCUSB.HA $stat))
head(res)
write.csv(res,file="resCUSB.HA_stat.csv",quote=F, row.names=F)
resCUSB.HA$lp=-log(resCUSB.HA$pvalue,10)
resCUSB.HA$lp[resCUSB.HA$log2FoldChange<0]=-resCUSB.HA$lp[resCUSB.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resCUSB.HA),"lp"= resCUSB.HA$lp))
write.csv(dev_lpv,file="resCUSB.HA_lp.csv",quote=F, row.names=F)


## LSXCU
resLSCU.AP<-results(dds, contrast=c('group', 'LSxCUAmbient', 'LSxCUPre')) #here is where the two contrasting conditions get defined
mcols(resLSCU.AP,use.names=TRUE)
summary(resLSCU.AP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 473, 4.2%
#LFC < 0 (down)     : 621, 5.5%

table(resLSCU.AP $padj < 0.05)
#FALSE  TRUE 
# 10377   815  
res=data.frame(cbind("gene"=row.names(resLSCU.AP),"stat"= resLSCU.AP $stat))
head(res)
write.csv(res,file="resLSCU.AP_stat.csv",quote=F, row.names=F)


resLSCU.HP<-results(dds, contrast=c('group', 'LSxCUHot', 'LSxCUPre')) #here is where the two contrasting conditions get defined
mcols(resLSCU.HP,use.names=TRUE)
summary(resLSCU.HP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1981, 18%
#LFC < 0 (down)     : 1991, 18%
table(resLSCU.HP $padj < 0.05)
#FALSE  TRUE 
# 7858  3334   
res=data.frame(cbind("gene"=row.names(resLSCU.HP),"stat"= resLSCU.HP $stat))
head(res)
write.csv(res,file="resLSCU.HP_stat.csv",quote=F, row.names=F)


resLSCU.HA<-results(dds, contrast=c('group', 'LSxCUHot', 'LSxCUAmbient')) #here is where the two contrasting conditions get defined
mcols(resLSCU.HA,use.names=TRUE)
summary(resLSCU.HA)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1639, 15%
#LFC < 0 (down)     : 1734, 15%

table(resLSCU.HA $padj < 0.05)
#FALSE  TRUE 
# 8423  2769  
res=data.frame(cbind("gene"=row.names(resLSCU.HA),"stat"= resLSCU.HA $stat))
head(res)
write.csv(res,file="resLSCU.HA_stat.csv",quote=F, row.names=F)
resLSCU.HA$lp=-log(resLSCU.HA$pvalue,10)
resLSCU.HA$lp[resLSCU.HA$log2FoldChange<0]=-resLSCU.HA$lp[resLSCU.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resLSCU.HA),"lp"= resLSCU.HA$lp))
write.csv(dev_lpv,file="resLSCU.HA_lp.csv",quote=F, row.names=F)


## LSXSB
resLSSB.AP<-results(dds, contrast=c('group', 'LSxSBAmbient', 'LSxSBPre')) #here is where the two contrasting conditions get defined
mcols(resLSSB.AP,use.names=TRUE)
summary(resLSSB.AP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2211, 20%
#LFC < 0 (down)     : 2404, 21%

table(resLSSB.AP $padj < 0.05)
#FALSE  TRUE 
# 7333  3859   
res=data.frame(cbind("gene"=row.names(resLSSB.AP),"stat"= resLSSB.AP $stat))
head(res)
write.csv(res,file="resLSSB.AP_stat.csv",quote=F, row.names=F)


#No hot larvae for this cross

## Central by Central population crosses ## - There are none! 

## Flip Flop Crosses- North Mom X Central Dad ##
##SBXBK
resSBBK.AP<-results(dds, contrast=c('group', 'SBxBKAmbient', 'SBxBKPre')) #here is where the two contrasting conditions get defined
mcols(resSBBK.AP,use.names=TRUE)
summary(resSBBK.AP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1094, 9.8%
#LFC < 0 (down)     : 1341, 12%
table(resSBBK.AP $padj < 0.05)
#FALSE  TRUE 
# 9200  1992  
res=data.frame(cbind("gene"=row.names(resSBBK.AP),"stat"= resSBBK.AP $stat))
head(res)
write.csv(res,file="resSBBK.AP_stat.csv",quote=F, row.names=F)

resSBBK.HP<-results(dds, contrast=c('group', 'SBxBKHot', 'SBxBKPre')) #here is where the two contrasting conditions get defined
mcols(resSBBK.HP,use.names=TRUE)
summary(resSBBK.HP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1892, 17%
#LFC < 0 (down)     : 2039, 18%
table(resSBBK.HP $padj < 0.05)
#FALSE  TRUE 
# 7852  3340  
res=data.frame(cbind("gene"=row.names(resSBBK.HP),"stat"= resSBBK.HP $stat))
head(res)
write.csv(res,file="resSBBK.HP_stat.csv",quote=F, row.names=F)

resSBBK.HA<-results(dds, contrast=c('group', 'SBxBKHot', 'SBxBKAmbient')) #here is where the two contrasting conditions get defined
mcols(resSBBK.HA,use.names=TRUE)
summary(resSBBK.HA)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2028, 18%
#LFC < 0 (down)     : 2075, 19%

table(resSBBK.HA $padj < 0.05)
#FALSE  TRUE 
# 7747  3445  
res=data.frame(cbind("gene"=row.names(resSBBK.HA),"stat"= resSBBK.HA $stat))
head(res)
write.csv(res,file="resSBBK.HA_stat.csv",quote=F, row.names=F)
resSBBK.HA$lp=-log(resSBBK.HA$pvalue,10)
resSBBK.HA$lp[resSBBK.HA$log2FoldChange<0]=-resSBBK.HA$lp[resSBBK.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resSBBK.HA),"lp"= resSBBK.HA$lp))
write.csv(dev_lpv,file="resSBBK.HA_lp.csv",quote=F, row.names=F)


##CUXBK
resCUBK.AP<-results(dds, contrast=c('group', 'CUxBKAmbient', 'CUxBKPre')) #here is where the two contrasting conditions get defined
mcols(resCUBK.AP,use.names=TRUE)
summary(resCUBK.AP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 628, 5.6%
#LFC < 0 (down)     : 657, 5.9%

table(resCUBK.AP $padj < 0.05)
#FALSE  TRUE 
#10245   947  
res=data.frame(cbind("gene"=row.names(resCUBK.AP),"stat"= resCUBK.AP $stat))
head(res)
write.csv(res,file="resCUBK.AP_stat.csv",quote=F, row.names=F)


resCUBK.HP<-results(dds, contrast=c('group', 'CUxBKHot', 'CUxBKPre')) #here is where the two contrasting conditions get defined
mcols(resCUBK.HP,use.names=TRUE)
summary(resCUBK.HP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2539, 23%
#LFC < 0 (down)     : 2339, 21%

table(resCUBK.HP $padj < 0.05)
#FALSE  TRUE 
# 7022  4170  
res=data.frame(cbind("gene"=row.names(resCUBK.HP),"stat"= resCUBK.HP $stat))
head(res)
write.csv(res,file="resCUBK.HP_stat.csv",quote=F, row.names=F)


resCUBK.HA<-results(dds, contrast=c('group', 'CUxBKHot', 'CUxBKAmbient')) #here is where the two contrasting conditions get defined
mcols(resCUBK.HA,use.names=TRUE)
summary(resCUBK.HA)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1961, 18%
#LFC < 0 (down)     : 1896, 17%

table(resCUBK.HA $padj < 0.05)
#FALSE  TRUE 
# 8025  3167   
res=data.frame(cbind("gene"=row.names(resCUBK.HA),"stat"= resCUBK.HA $stat))
head(res)
write.csv(res,file="resCUBK.HA_stat.csv",quote=F, row.names=F)
resCUBK.HA$lp=-log(resCUBK.HA$pvalue,10)
resCUBK.HA$lp[resCUBK.HA$log2FoldChange<0]=-resCUBK.HA$lp[resCUBK.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resCUBK.HA),"lp"= resCUBK.HA$lp))
write.csv(dev_lpv,file="resCUBK.HA_lp.csv",quote=F, row.names=F)


## Flip Flop Crosses- Central Mom X North Dad ##
##DRXSB
resDRSB.AP<-results(dds, contrast=c('group', 'DRxSBAmbient', 'DRxSBPre')) #here is where the two contrasting conditions get defined
mcols(resDRSB.AP,use.names=TRUE)
summary(resDRSB.AP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1254, 11%
#LFC < 0 (down)     : 1453, 13%

table(resDRSB.AP $padj < 0.05)
#FALSE  TRUE 
# 9154  2038   
res=data.frame(cbind("gene"=row.names(resDRSB.AP),"stat"= resDRSB.AP $stat))
head(res)
write.csv(res,file="resDRSB.AP_stat.csv",quote=F, row.names=F)


resDRSB.HP<-results(dds, contrast=c('group', 'DRxSBHot', 'DRxSBPre')) #here is where the two contrasting conditions get defined
mcols(resDRSB.HP,use.names=TRUE)
summary(resDRSB.HP)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1784, 16%
#LFC < 0 (down)     : 1908, 17%
table(resDRSB.HP $padj < 0.05)
#FALSE  TRUE 
# 8107  3085   
res=data.frame(cbind("gene"=row.names(resDRSB.HP),"stat"= resDRSB.HP $stat))
head(res)
write.csv(res,file="resDRSB.HP_stat.csv",quote=F, row.names=F)


resDRSB.HA<-results(dds, contrast=c('group', 'DRxSBHot', 'DRxSBAmbient')) #here is where the two contrasting conditions get defined
mcols(resDRSB.HA,use.names=TRUE)
summary(resDRSB.HA)
#out of 11192 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1964, 18%
#LFC < 0 (down)     : 1932, 17%

table(resDRSB.HA $padj < 0.05)
#FALSE  TRUE 
# 8024  3168   
res=data.frame(cbind("gene"=row.names(resDRSB.HA),"stat"= resDRSB.HA $stat))
head(res)
write.csv(res,file="resDRSB.HA_stat.csv",quote=F, row.names=F)
resDRSB.HA$lp=-log(resDRSB.HA$pvalue,10)
resDRSB.HA$lp[resDRSB.HA$log2FoldChange<0]=-resDRSB.HA$lp[resDRSB.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resDRSB.HA),"lp"= resDRSB.HA$lp))
write.csv(dev_lpv,file="resDRSB.HA_lp.csv",quote=F, row.names=F)


save(vstG,vstG.df,colData,resBKBK.AP,resBKBK.HP,resBKBK.HA,resDRDR.AP,resDRDR.HP,resDRDR.HA,resSBSB.AP,resSBSB.HP,resSBSB.HA,resCUCU.AP,resCUCU.HP,resCUCU.HA,resLSLS.AP,resLSLS.HP,resLSLS.HA,resCUSB.AP,resCUSB.HP,resCUSB.HA,resLSCU.AP,resLSCU.HP,resLSCU.HA,resLSSB.AP,resSBBK.AP,resSBBK.HP,resSBBK.HA,resCUBK.AP,resCUBK.HP,resCUBK.HA,resDRSB.AP,resDRSB.HP,resDRSB.HA,file="AGF_larval_2018_WaldGroupContrasts.Rdata")


##### pick out genes specific to CU dam heat response ############
ll=load("AGF_larval_2018_WaldGroupContrasts.Rdata")
table(resCUBK.HA$padj<= 0.05) #3167
table(resCUCU.HA$padj<= 0.05) #3197
table(resCUSB.HA$padj<= 0.05) #3937

CUgenes=rldPs[rldPs$resCUCU.HA_padj<=0.05 & rldPs$resCUSB.HA_padj<=0.05 & rldPs$resCUBK.HA_padj<=0.05, ] #1179 genes

CUgenes=rldPs[rldPs$resCUCU.HA_padj<=0.05 & rldPs$resCUSB.HA_padj<=0.05 & rldPs$resCUBK.HA_padj<=0.05 & rldPs$resBKBK.HA_padj>0.05 & rldPs$resDRDR.HA_padj>0.05 & rldPs$resSBSB.HA_padj>0.05 & rldPs$resLSLS.HA_padj>0.05 & rldPs$resLSCU.HA_padj>0.05 & rldPs$resSBBK.HA_padj>0.05 & rldPs$resDRSB.HA_padj>0.05, ] #12 genes

CUgenes$gene=row.names(CUgenes)
genes=as.data.frame(CUgenes[-c(97:159) ])
names(genes)=paste0(colData$group, colData$sample)
genes$gene=row.names(genes)

geneL=genes %>%
	gather("sample", "expr", -gene)
geneL$cross = str_sub(geneL$sample, 1,5)
geneL$treat = str_sub(geneL$sample, 6,8)
geneL$dam = str_sub(geneL$sample, 1,2)
CUdam<-(geneL$cross=="CUxSB"| geneL$cross=="CUxCU" | geneL$cross=="CUxBK")
CUdam[CUdam ==T]<- "T"
CUdam[CUdam ==F]<- "F"
geneL$CUdam = as.factor(CUdam)

annot=read.table("aten_annotation_table.tsv",header=T,sep="\t",quote=NULL,fill=T)
annot=as.data.frame(annot[c(1,10) ],)
colnames(annot)=c("transcript","gene") 
genes = geneL$gene
genes2annot = match(genes,annot$transcript)
geneL_annot=data.frame(GeneID=geneL$gene, annot[genes2annot,],geneL)

geneL_annot$treat_f = factor(geneL$treat, levels=c("Pre","Amb","Hot"))
#geneL$cross_f = factor(geneL$cross, levels=c('LSSB','LSLS','DRDR','SBSB','LSCU','DRSB','BKBK','CUCU','CUSB','CUBK','SBBK'))

summ=summarySE(data= geneL_annot,measurevar="expr",groupvars=c("gene","treat_f","CUdam"))
#summ=summarySE(data= geneL,measurevar="expr",groupvars=c("gene","treat_f","cross"))

pd <- position_dodge(0.15)
ggplot(summ,aes(x=treat_f,y=expr, colour= CUdam))+
	geom_point(aes(group= CUdam),position=pd,size=2)+
	geom_line(aes(group= CUdam,linetype= CUdam),position=pd)+
	geom_errorbar(aes(ymin=expr-se,ymax=expr+se),lwd=0.4,width=0.3,position=pd)+
	scale_color_manual(values=c("F"="#3B9AB2","T"="#F21A00"))+
	ggtitle("CUdam specific heat response genes")+
	theme_bw()+
	facet_wrap(~gene,scales="free_y")

#summ$gene
#"Amiloride-sensitive amine oxidase [copper-containing] (DAO) (Diamine oxidase) (EC 1.4.3.22) (Amiloride-binding protein 1) (Amine oxidase copper domain-containing protein 1) (Histaminase)"
#"Beta-hexosaminidase (EC 3.2.1.52) (Beta-N-acetylhexosaminidase) (Chitobiase) (N-acetyl-beta-glucosaminidase)"
#"Coactosin-like protein"
#"Cubilin (460 kDa receptor) (Intestinal intrinsic factor receptor) (Intrinsic factor-cobalamin receptor) (Intrinsic factor-vitamin B12 receptor)"
#"Inhibitor of nuclear factor kappa-B kinase subunit alpha (I-kappa-B kinase alpha) (IKK-A) (IKK-alpha) (IkBKA) (IkappaB kinase) (EC 2.7.11.10) (Conserved helix-loop-helix ubiquitous kinase) (I-kappa-B kinase 1) (IKK1) (Nuclear factor NF-kappa-B inhibitor kinase alpha) (NFKBIKA) (Transcription factor 16) (TCF-16)"
#"Interleukin-1 receptor-associated kinase 1-binding protein 1 (IRAK1-binding protein 1)" 
#"Partner of Y14 and mago (PYM homolog 1 exon junction complex-associated factor) (Protein wibg homolog)"
#"RRP15-like protein (Ribosomal RNA-processing protein 15)"
#"Tetratricopeptide repeat protein 28 (TPR repeat protein 28) (TPR repeat-containing big gene cloned at Keio)"
#"Transposon TX1 uncharacterized 149 kDa protein (ORF 2)"

######################################################## run DESeq2 on only the PRE samples
gcountsALL=read.table("AGF_larvae_2018_geneCounts_aten.txt", header=T) 

# clean up file names
colnames(gcountsALL)<-sub("X","", colnames(gcountsALL))
colnames(gcountsALL) <- sub("_[^_]+$", "", colnames(gcountsALL))
colnames(gcountsALL)<-gsub("\\.","_", colnames(gcountsALL))

length(gcountsALL[,1]) #25142
dim(gcountsALL) #25142    96


summary(colSums(gcountsALL))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  135519 1345494 1463876 1444555 1552249 1859418 


### BUILD A DATAFRAME ASSOCIATING SAMPLE NAMESWITH TREATMENT CONDITIONS ###
colData=read.csv("J19188metaPRE.csv", header=T)
colData=colData[!duplicated(colData$sample), ]
rownames(colData)<-colData$sample
dim(colData) #30  3

all(rownames(colData) %in% colnames(gcountsALL)) #TRUE
gcounts <- gcountsALL[, rownames(colData)]
all(rownames(colData) == colnames(gcounts))
dim(colData) #30  4

### REMOVE GENES WITH LOW MEAN COUNTS ###

mns = apply(gcounts, 1, mean)

gcountsPRE=gcounts[mns>10,] #get rid of genes that show little or no expression
#gcounts=gcounts[mns>10,] #only highly expressed genes for WGCNA

table(mns > 10)
#FALSE  TRUE 
#13835 11089
dim(gcountsPRE) #11089    30

### LRT TEST - FULL MODEL ###
dds<-DESeqDataSetFromMatrix(gcountsPRE,

	colData = colData, 

	design = formula(~ survivalheatPRE))

vstPRE=vst(dds)
vstPRE.df = assay(vstPRE)

## LRT for treatment
dds <- DESeq(dds, test="LRT", reduced=~1)
resP <- results(dds)
summary(resP)
#out of 11089 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 366, 3.3%
#LFC < 0 (down)     : 330, 3%
#outliers [1]       : 5, 0.045%
#low counts [2]     : 0, 0%


table(resP$padj < 0.05)
#FALSE  TRUE 
# 10646   438    
res=data.frame(cbind("gene"=row.names(resP),"stat"= resP$stat))
head(res)
write.csv(res,file="resP_stat.csv",quote=F, row.names=F)

res=data.frame(resP)
res$pvalue_go[res$pvalue<=0.01]<-1
res$pvalue_go[res$pvalue>0.01]<-0 
res$pvalue_go=as.character(res$pvalue_go) 
summary(res$pvalue_go)
#   0    1 NA's 
#8286  534    2
res=res[complete.cases(res),]
res$names=row.names(res)
survival_fishers=data.frame(cbind("gene"= res$names,"pvalue"= res$pvalue_go))
#survival_fishers$pvalue=as.character(survival_fishers$pvalue)
write.csv(survival_fishers,file="survival_LRT_fishers_pval0.01.csv",quote=F, row.names=F)

save(resP,file="AGF_larval_2018_PreRes_aten.Rdata")

### heatmap

vals=(cbind(resP$stat, resP$padj)) #collect pvalues for each gene
colnames(vals)=c("stat_resP","padj_resP") 
rldpvals=as.data.frame(cbind(vstPRE.df,vals)) #combine RLD with pvals

annot=read.table("aten_annotation_table.tsv",header=T,sep="\t",quote=NULL,fill=T)
annot=as.data.frame(annot[c(1,10) ],)
colnames(annot)=c("transcript","gene") 
genes = row.names(rldpvals)
genes2annot = match(genes,annot$transcript)
rldpvals_annot=data.frame(GeneID=row.names(rldpvals), annot[genes2annot,], rldpvals)

rldpvals_annot = rldpvals_annot[rldpvals_annot $padj_resP<= 0.01,]
rldpvals_annot = rldpvals_annot[complete.cases(rldpvals_annot),]
dim(rldpvals_annot) #134  35
write.csv(rldpvals_annot,file="PRE_sig_rldpvals_annot.csv",quote=F, row.names=F)

row.names(rldpvals_annot)=make.names(rldpvals_annot$gene,unique=T)
rldpvals_annot = rldpvals_annot[, 4:33]

means=apply(rldpvals_annot,1,mean) # means of rows
explc= rldpvals_annot-means


heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=1)(100)
heat.colors = colorRampPalette(wes_palette("Zissou1", 30, type = "continuous"), bias=1.1)(100)
pdf("Heatmap_genes_survivalPre.pdf",height=10,width=25)
pheatmap(explc,color=heat.colors,cluster_cols=F,border_color=NA,clustering_distance_rows="correlation")
dev.off()
#pheatmap(explc,color=heat.colors,cluster_cols=T,border_color=NA)



### principal coordinate calculation ###
ll=load("AGF_larval_2018_WaldGroupContrasts.Rdata")

conditions=colData

#conditions$treat_c=scale_colour_manual(values=c("Pre"="#EBCC2A","Ambient"="#3B9AB2", "Hot"="#F21A00"))
grp=rep("#EBCC2A",ncol(vstG.df))
grp[grep("Ambient",conditions$treatment)]="#3B9AB2"
grp[grep("Hot",conditions$treatment)]="#F21A00"

fitpc=capscale(dist(t(vstG.df),method="manhattan")~1)
#fitpc=capscale(dist(t(vstG.df),method="euclidean")~1)
summary(eigenvals(fitpc))

plot(fitpc$CA$u,pch=16, col="white",main="AGF Larvae 2018", xlab="PC1 27.17%", ylab="PC2 20.58%")
abline(h=0,v=0, col="grey")
points(fitpc$CA$u,pch=16, col=grp)
ordiellipse(fitpc$CA$u,conditions$treat,draw="polygon",label=T)

plot(fitpc$CA$u,pch=16, col="white",main="AGF Larvae 2018", xlab="PC1 25.6%", ylab="PC2 14.9%")
abline(h=0,v=0, col="grey")
points(fitpc$CA$u,pch=16, col=grp)
ordispider(fitpc$CA$u,conditions$cross,col="black", label=T)
#ordiellipse(fitpc$CA$u,conditions$cross,draw="polygon",label=T)

ad=adonis(t(vstG.df)~treatment+cross,data=conditions,method="manhattan")
#Call:
#adonis(formula = t(vstG.df) ~ treatment + cross, data = conditions,      method = "manhattan") 
#
#Permutation: free
#Number of permutations: 999
#
#Terms added sequentially (first to last)
#
#          Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)    
#treatment  2 257372697 128686348 29.7889 0.29251  0.001 ***
#cross     10 263941810  26394181  6.1098 0.29998  0.001 ***
#Residuals 83 358555233   4319943         0.40751           
#Total     95 879869740                   1.00000           
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


