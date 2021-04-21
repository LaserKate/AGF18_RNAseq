

setwd("/Users/mes0192/Dropbox/KQ_GBR_spawnseq/analysis/DESeq2_larvae")
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
library(gridExtra)


#Questions we can answer with this dataest
# is "Pre" gene expresison predictive of survival? Correlate pre-gene expression with survival, Dixon et al. 2015, Rachel's papers
# do crosses from northern parents show a different functional response to heat stress than crosses from northern pops? Does this ring true for moms vs. sires from north or central? look for DEG overlap and GO overlap. 
# How does transcriptional plasticity predict survival? 
# Does transcriptional plasticity differ based on where the dam or sire was from? 


gcountsALL=read.table("AGF_larvae_2018_geneCounts.txt", header=T) 

# clean up file names
colnames(gcountsALL)<-sub("X","", colnames(gcountsALL))
colnames(gcountsALL) <- sub("_[^_]+$", "", colnames(gcountsALL))
colnames(gcountsALL)<-gsub("\\.","_", colnames(gcountsALL))

length(gcountsALL[,1]) #25142
dim(gcountsALL) #25142    96


summary(colSums(gcountsALL))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  79710  795300  871500  849100  917300 1102000 


### REMOVE GENES WITH LOW MEAN COUNTS ###

mns = apply(gcountsALL, 1, mean)

gcounts=gcountsALL[mns>10,] #get rid of genes that show little or no expression
#gcounts=gcounts[mns>10,] #only highly expressed genes for WGCNA

table(mns > 10)
#FALSE  TRUE 
#16195  8947 
dim(gcounts) #8947   96


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


rldCT=rlog(ddsCT)
rldCT.df = assay(rldCT)

### Wald test for treatment in adults ###

dds<-DESeq(ddsCT, minReplicatesForReplace=Inf) 
resHT<-results(dds, contrast=c('treatment', 'Hot', 'Ambient')) #here is where the two contrasting conditions get defined
mcols(resHT,use.names=TRUE)
summary(resHT)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 3065, 34%
#LFC < 0 (down)     : 3121, 35%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
table(resHT$padj < 0.05)
#FALSE  TRUE 
# 3212  5735 

resHP<-results(dds, contrast=c('treatment', 'Hot', 'Pre')) #here is where the two contrasting conditions get defined
mcols(resHP,use.names=TRUE)
summary(resHP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2915, 33%
#LFC < 0 (down)     : 2999, 34%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
table(resHP$padj < 0.05)
#FALSE  TRUE 
# 3478  5469

resAP<-results(dds, contrast=c('treatment', 'Ambient', 'Pre')) #here is where the two contrasting conditions get defined
mcols(resAP,use.names=TRUE)
summary(resAP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1857, 21%
#LFC < 0 (down)     : 1891, 21%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
table(resAP$padj < 0.05)
#FALSE  TRUE 
# 5797  3150

save(rldCT,rldCT.df,resHT,resHP,resAP, colData,file="AGF_larval_2018_WaldCrossTreat.Rdata")



### LRT TEST - FULL MODEL ###
dds<-DESeqDataSetFromMatrix(gcounts,

	colData = colData, 

	design = formula(~ treatment+dam+sire))

rldLRT=rlog(dds)
rldLRT.df = assay(rldLRT)

## LRT for treatment
dds <- DESeq(dds, test="LRT", reduced=~dam+sire)
resT <- results(dds)
summary(resT)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 3382, 38%
#LFC < 0 (down)     : 3624, 41%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
table(resT$padj < 0.05)
#FALSE  TRUE 
# 2348  6599

## LRT for dam
dds <- DESeq(dds, test="LRT", reduced=~treatment+sire)
resD <- results(dds)
summary(resD)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2835, 32%
#LFC < 0 (down)     : 3024, 34%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
table(resD$padj < 0.05)
#FALSE  TRUE 
# 3754  5193 

## LRT for sire
dds <- DESeq(dds, test="LRT", reduced=~treatment+dam)
resS <- results(dds)
summary(resS)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2472, 28%
#LFC < 0 (down)     : 2773, 31%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
table(resS$padj < 0.05)
#FALSE  TRUE 
# 4360  4587

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
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1877, 21%
#LFC < 0 (down)     : 1664, 19%
table(resBKBK.AP$padj < 0.05)
#FALSE  TRUE 
# 5987  2960 
#res=data.frame(cbind("gene"=row.names(resBKBK.AP),"stat"=resBKBK.AP$stat))
#head(res)
#write.csv(res,file="GO/resBKBK.AP_stat.csv",quote=F, row.names=F)


resBKBK.HP<-results(dds, contrast=c('group', 'BKxBKHot', 'BKxBKPre')) #here is where the two contrasting conditions get defined
mcols(resBKBK.HP,use.names=TRUE)
summary(resBKBK.HP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2153, 24%
#LFC < 0 (down)     : 2186, 24%
table(resBKBK.HP $padj < 0.05)
#FALSE  TRUE 
# 5131  3816 
#res=data.frame(cbind("gene"=row.names(resBKBK.HP),"stat"= resBKBK.HP $stat))
#head(res)
#write.csv(res,file="GO/resBKBK.HP_stat.csv",quote=F, row.names=F)

resBKBK.HA<-results(dds, contrast=c('group', 'BKxBKHot', 'BKxBKAmbient')) #here is where the two contrasting conditions get defined
mcols(resBKBK.HA,use.names=TRUE)
summary(resBKBK.HA)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1584, 18%
#LFC < 0 (down)     : 1679, 19%
table(resBKBK.HA $padj < 0.05)
#FALSE  TRUE 
# 6124  2823 
#res=data.frame(cbind("gene"=row.names(resBKBK.HA),"stat"= resBKBK.HA $stat))
#head(res)
#write.csv(res,file="GO/resBKBK.HA_stat.csv",quote=F, row.names=F)
resBKBK.HA$lp=-log(resBKBK.HA$pvalue,10)
resBKBK.HA$lp[resBKBK.HA$log2FoldChange<0]=-resBKBK.HA$lp[resBKBK.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resBKBK.HA),"lp"=resBKBK.HA$lp))
write.csv(dev_lpv,file="GO_lpvals/resBKBK.HA_lp.csv",quote=F, row.names=F)

# DR - Davies = Central
resDRDR.AP<-results(dds, contrast=c('group', 'DRxDRAmbient', 'DRxDRPre')) #here is where the two contrasting conditions get defined
mcols(resDRDR.AP,use.names=TRUE)
summary(resDRDR.AP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1497, 17%
#LFC < 0 (down)     : 1726, 19%
table(resDRDR.AP $padj < 0.05)
#FALSE  TRUE 
# 6193  2754 
#res=data.frame(cbind("gene"=row.names(resDRDR.AP),"stat"= resDRDR.AP $stat))
#head(res)
#write.csv(res,file="GO/resDRDR.AP_stat.csv",quote=F, row.names=F)


resDRDR.HP<-results(dds, contrast=c('group', 'DRxDRHot', 'DRxDRPre')) #here is where the two contrasting conditions get defined
mcols(resDRDR.HP,use.names=TRUE)
summary(resDRDR.HP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1670, 19%
#LFC < 0 (down)     : 1796, 20%
table(resDRDR.HP $padj < 0.05)
#FALSE  TRUE 
# 5994  2953 
#res=data.frame(cbind("gene"=row.names(resDRDR.HP),"stat"= resDRDR.HP $stat))
#head(res)
#write.csv(res,file="GO/resDRDR.HP_stat.csv",quote=F, row.names=F)


resDRDR.HA<-results(dds, contrast=c('group', 'DRxDRHot', 'DRxDRAmbient')) #here is where the two contrasting conditions get defined
mcols(resDRDR.HA,use.names=TRUE)
summary(resDRDR.HA)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1596, 18%
#LFC < 0 (down)     : 1579, 18%
table(resDRDR.HA $padj < 0.05)
#FALSE  TRUE 
# 6259  2688 
#res=data.frame(cbind("gene"=row.names(resDRDR.HA),"stat"= resDRDR.HA $stat))
#head(res)
#write.csv(res,file="GO/resDRDR.HA_stat.csv",quote=F, row.names=F)
resDRDR.HA$lp=-log(resDRDR.HA$pvalue,10)
resDRDR.HA$lp[resDRDR.HA$log2FoldChange<0]=-resDRDR.HA$lp[resDRDR.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resDRDR.HA),"lp"= resDRDR.HA$lp))
write.csv(dev_lpv,file="GO_lpvals/resDRDR.HA_lp.csv",quote=F, row.names=F)


# SB - Sand Bank 7 - North 
resSBSB.AP<-results(dds, contrast=c('group', 'SBxSBAmbient', 'SBxSBPre')) #here is where the two contrasting conditions get defined
mcols(resSBSB.AP,use.names=TRUE)
summary(resSBSB.AP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 639, 7.1%
#LFC < 0 (down)     : 701, 7.8%
table(resSBSB.AP $padj < 0.05)
#FALSE  TRUE 
# 7914  1033 
#res=data.frame(cbind("gene"=row.names(resSBSB.AP),"stat"= resSBSB.AP $stat))
#head(res)
#write.csv(res,file="GO/resSBSB.AP_stat.csv",quote=F, row.names=F)


resSBSB.HP<-results(dds, contrast=c('group', 'SBxSBHot', 'SBxSBPre')) #here is where the two contrasting conditions get defined
mcols(resSBSB.HP,use.names=TRUE)
summary(resSBSB.HP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1916, 21%
#LFC < 0 (down)     : 1855, 21%
table(resSBSB.HP $padj < 0.05)
#FALSE  TRUE 
# 5717  3230 
#res=data.frame(cbind("gene"=row.names(resSBSB.HP),"stat"= resSBSB.HP $stat))
#head(res)
#write.csv(res,file="GO/resSBSB.HP_stat.csv",quote=F, row.names=F)


resSBSB.HA<-results(dds, contrast=c('group', 'SBxSBHot', 'SBxSBAmbient')) #here is where the two contrasting conditions get defined
mcols(resSBSB.HA,use.names=TRUE)
summary(resSBSB.HA)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1901, 21%
#LFC < 0 (down)     : 1912, 21%
table(resSBSB.HA $padj < 0.05)
#FALSE  TRUE 
# 5715  3232 
#res=data.frame(cbind("gene"=row.names(resSBSB.HA),"stat"= resSBSB.HA $stat))
#head(res)
#write.csv(res,file="GO/resSBSB.HA_stat.csv",quote=F, row.names=F)
resSBSB.HA$lp=-log(resSBSB.HA$pvalue,10)
resSBSB.HA$lp[resSBSB.HA$log2FoldChange<0]=-resSBSB.HA$lp[resSBSB.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resSBSB.HA),"lp"= resSBSB.HA$lp))
write.csv(dev_lpv,file="GO_lpvals/resSBSB.HA_lp.csv",quote=F, row.names=F)

# CU - Curd Reef - North 
resCUCU.AP<-results(dds, contrast=c('group', 'CUxCUAmbient', 'CUxCUPre')) #here is where the two contrasting conditions get defined
mcols(resCUCU.AP,use.names=TRUE)
summary(resCUCU.AP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 291, 3.3%
#LFC < 0 (down)     : 593, 6.6%
table(resCUCU.AP $padj < 0.05)
#FALSE  TRUE 
# 8278   669 
#res=data.frame(cbind("gene"=row.names(resCUCU.AP),"stat"= resCUCU.AP $stat))
#head(res)
#write.csv(res,file="GO/resCUCU.AP_stat.csv",quote=F, row.names=F)



resCUCU.HP<-results(dds, contrast=c('group', 'CUxCUHot', 'CUxCUPre')) #here is where the two contrasting conditions get defined
mcols(resCUCU.HP,use.names=TRUE)
summary(resCUCU.HP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1686, 19%
#LFC < 0 (down)     : 1814, 20%
table(resCUCU.HP $padj < 0.05)
#FALSE  TRUE 
# 5923  3024  
#res=data.frame(cbind("gene"=row.names(resCUCU.HP),"stat"= resCUCU.HP $stat))
#head(res)
#write.csv(res,file="GO/resCUCU.HP_stat.csv",quote=F, row.names=F)


resCUCU.HA<-results(dds, contrast=c('group', 'CUxCUHot', 'CUxCUAmbient')) #here is where the two contrasting conditions get defined
mcols(resCUCU.HA,use.names=TRUE)
summary(resCUCU.HA)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1322, 15%
#LFC < 0 (down)     : 1424, 16%
table(resCUCU.HA $padj < 0.05)
#FALSE  TRUE 
# 6694  2253  
#res=data.frame(cbind("gene"=row.names(resCUCU.HA),"stat"= resCUCU.HA $stat))
#head(res)
#write.csv(res,file="GO/resCUCU.HA_stat.csv",quote=F, row.names=F)
resCUCU.HA$lp=-log(resCUCU.HA$pvalue,10)
resCUCU.HA$lp[resCUCU.HA$log2FoldChange<0]=-resCUCU.HA$lp[resCUCU.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resCUCU.HA),"lp"= resCUCU.HA$lp))
write.csv(dev_lpv,file="GO_lpvals/resCUCU.HA_lp.csv",quote=F, row.names=F)


# LS - Long Sandy - North 
resLSLS.AP<-results(dds, contrast=c('group', 'LSxLSAmbient', 'LSxLSPre')) #here is where the two contrasting conditions get defined
mcols(resLSLS.AP,use.names=TRUE)
summary(resLSLS.AP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 247, 2.8%
#LFC < 0 (down)     : 374, 4.2%
table(resLSLS.AP $padj < 0.05)
#FALSE  TRUE 
# 7984   442  
#res=data.frame(cbind("gene"=row.names(resLSLS.AP),"stat"= resLSLS.AP $stat))
#head(res)
#write.csv(res,file="GO/resLSLS.AP_stat.csv",quote=F, row.names=F)


resLSLS.HP<-results(dds, contrast=c('group', 'LSxLSHot', 'LSxLSPre')) #here is where the two contrasting conditions get defined
mcols(resLSLS.HP,use.names=TRUE)
summary(resLSLS.HP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1305, 15%
#LFC < 0 (down)     : 1433, 16%
table(resLSLS.HP $padj < 0.05)
#FALSE  TRUE 
# 6704  2243  
#res=data.frame(cbind("gene"=row.names(resLSLS.HP),"stat"= resLSLS.HP $stat))
#head(res)
#write.csv(res,file="GO/resLSLS.HP_stat.csv",quote=F, row.names=F)


resLSLS.HA<-results(dds, contrast=c('group', 'LSxLSHot', 'LSxLSAmbient')) #here is where the two contrasting conditions get defined
mcols(resLSLS.HA,use.names=TRUE)
summary(resLSLS.HA)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1361, 15%
#LFC < 0 (down)     : 1456, 16%
table(resLSLS.HA $padj < 0.05)
#FALSE  TRUE 
# 6709  2238  
#res=data.frame(cbind("gene"=row.names(resLSLS.HA),"stat"= resLSLS.HA $stat))
#head(res)
#write.csv(res,file="GO/resLSLS.HA_stat.csv",quote=F, row.names=F)
resLSLS.HA$lp=-log(resLSLS.HA$pvalue,10)
resLSLS.HA$lp[resLSLS.HA$log2FoldChange<0]=-resLSLS.HA$lp[resLSLS.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resLSLS.HA),"lp"= resLSLS.HA$lp))
write.csv(dev_lpv,file="GO_lpvals/resLSLS.HA_lp.csv",quote=F, row.names=F)


## North by North population crosses ##
## CUXSB
resCUSB.AP<-results(dds, contrast=c('group', 'CUxSBAmbient', 'CUxSBPre')) #here is where the two contrasting conditions get defined
mcols(resCUSB.AP,use.names=TRUE)
summary(resCUSB.AP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 262, 2.9%
#LFC < 0 (down)     : 433, 4.8%
table(resCUSB.AP $padj < 0.05)
#FALSE  TRUE 
# 8425   522  
#res=data.frame(cbind("gene"=row.names(resCUSB.AP),"stat"= resCUSB.AP $stat))
#head(res)
#write.csv(res,file="GO/resCUSB.AP_stat.csv",quote=F, row.names=F)


resCUSB.HP<-results(dds, contrast=c('group', 'CUxSBHot', 'CUxSBPre')) #here is where the two contrasting conditions get defined
mcols(resCUSB.HP,use.names=TRUE)
summary(resCUSB.HP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1666, 19%
#LFC < 0 (down)     : 1779, 20%
table(resCUSB.HP $padj < 0.05)
#FALSE  TRUE 
# 5970  2977  
#res=data.frame(cbind("gene"=row.names(resCUSB.HP),"stat"= resCUSB.HP $stat))
#head(res)
#write.csv(res,file="GO/resCUSB.HP_stat.csv",quote=F, row.names=F)


resCUSB.HA<-results(dds, contrast=c('group', 'CUxSBHot', 'CUxSBAmbient')) #here is where the two contrasting conditions get defined
mcols(resCUSB.HA,use.names=TRUE)
summary(resCUSB.HA)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1586, 18%
#LFC < 0 (down)     : 1746, 20%
table(resCUSB.HA $padj < 0.05)
#FALSE  TRUE 
# 6115  2832  
#res=data.frame(cbind("gene"=row.names(resCUSB.HA),"stat"= resCUSB.HA $stat))
#head(res)
#write.csv(res,file="GO/resCUSB.HA_stat.csv",quote=F, row.names=F)
resCUSB.HA$lp=-log(resCUSB.HA$pvalue,10)
resCUSB.HA$lp[resCUSB.HA$log2FoldChange<0]=-resCUSB.HA$lp[resCUSB.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resCUSB.HA),"lp"= resCUSB.HA$lp))
write.csv(dev_lpv,file="GO_lpvals/resCUSB.HA_lp.csv",quote=F, row.names=F)


## LSXCU
resLSCU.AP<-results(dds, contrast=c('group', 'LSxCUAmbient', 'LSxCUPre')) #here is where the two contrasting conditions get defined
mcols(resLSCU.AP,use.names=TRUE)
summary(resLSCU.AP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 256, 2.9%
#LFC < 0 (down)     : 410, 4.6%
table(resLSCU.AP $padj < 0.05)
#FALSE  TRUE 
# 8455   492  
#res=data.frame(cbind("gene"=row.names(resLSCU.AP),"stat"= resLSCU.AP $stat))
#head(res)
#write.csv(res,file="GO/resLSCU.AP_stat.csv",quote=F, row.names=F)


resLSCU.HP<-results(dds, contrast=c('group', 'LSxCUHot', 'LSxCUPre')) #here is where the two contrasting conditions get defined
mcols(resLSCU.HP,use.names=TRUE)
summary(resLSCU.HP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1378, 15%
#LFC < 0 (down)     : 1450, 16%
table(resLSCU.HP $padj < 0.05)
#FALSE  TRUE 
# 6611  2336  
#res=data.frame(cbind("gene"=row.names(resLSCU.HP),"stat"= resLSCU.HP $stat))
#head(res)
#write.csv(res,file="GO/resLSCU.HP_stat.csv",quote=F, row.names=F)


resLSCU.HA<-results(dds, contrast=c('group', 'LSxCUHot', 'LSxCUAmbient')) #here is where the two contrasting conditions get defined
mcols(resLSCU.HA,use.names=TRUE)
summary(resLSCU.HA)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1152, 13%
#LFC < 0 (down)     : 1253, 14%
table(resLSCU.HA $padj < 0.05)
#FALSE  TRUE 
# 7028  1919  
#res=data.frame(cbind("gene"=row.names(resLSCU.HA),"stat"= resLSCU.HA $stat))
#head(res)
#write.csv(res,file="GO/resLSCU.HA_stat.csv",quote=F, row.names=F)
resLSCU.HA$lp=-log(resLSCU.HA$pvalue,10)
resLSCU.HA$lp[resLSCU.HA$log2FoldChange<0]=-resLSCU.HA$lp[resLSCU.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resLSCU.HA),"lp"= resLSCU.HA$lp))
write.csv(dev_lpv,file="GO_lpvals/resLSCU.HA_lp.csv",quote=F, row.names=F)


## LSXSB
resLSSB.AP<-results(dds, contrast=c('group', 'LSxSBAmbient', 'LSxSBPre')) #here is where the two contrasting conditions get defined
mcols(resLSSB.AP,use.names=TRUE)
summary(resLSSB.AP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1598, 18%
#LFC < 0 (down)     : 1698, 19%
table(resLSSB.AP $padj < 0.05)
#FALSE  TRUE 
# 6150  2797  
#res=data.frame(cbind("gene"=row.names(resLSSB.AP),"stat"= resLSSB.AP $stat))
#head(res)
#write.csv(res,file="GO/resLSSB.AP_stat.csv",quote=F, row.names=F)


#No hot larvae for this cross

## Central by Central population crosses ## - There are none! 

## Flip Flop Crosses- North Mom X Central Dad ##
##SBXBK
resSBBK.AP<-results(dds, contrast=c('group', 'SBxBKAmbient', 'SBxBKPre')) #here is where the two contrasting conditions get defined
mcols(resSBBK.AP,use.names=TRUE)
summary(resSBBK.AP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 802, 9%
#LFC < 0 (down)     : 885, 9.9%
table(resSBBK.AP $padj < 0.05)
#FALSE  TRUE 
# 7563  1384  
#res=data.frame(cbind("gene"=row.names(resSBBK.AP),"stat"= resSBBK.AP $stat))
#head(res)
#write.csv(res,file="GO/resSBBK.AP_stat.csv",quote=F, row.names=F)

resSBBK.HP<-results(dds, contrast=c('group', 'SBxBKHot', 'SBxBKPre')) #here is where the two contrasting conditions get defined
mcols(resSBBK.HP,use.names=TRUE)
summary(resSBBK.HP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1345, 15%
#LFC < 0 (down)     : 1407, 16%
table(resSBBK.HP $padj < 0.05)
#FALSE  TRUE 
# 6672  2275  
#res=data.frame(cbind("gene"=row.names(resSBBK.HP),"stat"= resSBBK.HP $stat))
#head(res)
#write.csv(res,file="GO/resSBBK.HP_stat.csv",quote=F, row.names=F)

resSBBK.HA<-results(dds, contrast=c('group', 'SBxBKHot', 'SBxBKAmbient')) #here is where the two contrasting conditions get defined
mcols(resSBBK.HA,use.names=TRUE)
summary(resSBBK.HA)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1425, 16%
#LFC < 0 (down)     : 1498, 17%
table(resSBBK.HA $padj < 0.05)
#FALSE  TRUE 
# 6526  2421  
#res=data.frame(cbind("gene"=row.names(resSBBK.HA),"stat"= resSBBK.HA $stat))
#head(res)
#write.csv(res,file="GO/resSBBK.HA_stat.csv",quote=F, row.names=F)
resSBBK.HA$lp=-log(resSBBK.HA$pvalue,10)
resSBBK.HA$lp[resSBBK.HA$log2FoldChange<0]=-resSBBK.HA$lp[resSBBK.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resSBBK.HA),"lp"= resSBBK.HA$lp))
write.csv(dev_lpv,file="GO_lpvals/resSBBK.HA_lp.csv",quote=F, row.names=F)


##CUXBK
resCUBK.AP<-results(dds, contrast=c('group', 'CUxBKAmbient', 'CUxBKPre')) #here is where the two contrasting conditions get defined
mcols(resCUBK.AP,use.names=TRUE)
summary(resCUBK.AP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 409, 4.6%
#LFC < 0 (down)     : 419, 4.7%
table(resCUBK.AP $padj < 0.05)
#FALSE  TRUE 
# 8313   634  
#res=data.frame(cbind("gene"=row.names(resCUBK.AP),"stat"= resCUBK.AP $stat))
#head(res)
#write.csv(res,file="GO/resCUBK.AP_stat.csv",quote=F, row.names=F)


resCUBK.HP<-results(dds, contrast=c('group', 'CUxBKHot', 'CUxBKPre')) #here is where the two contrasting conditions get defined
mcols(resCUBK.HP,use.names=TRUE)
summary(resCUBK.HP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1800, 20%
#LFC < 0 (down)     : 1671, 19%
table(resCUBK.HP $padj < 0.05)
#FALSE  TRUE 
# 6071  2876  
#res=data.frame(cbind("gene"=row.names(resCUBK.HP),"stat"= resCUBK.HP $stat))
#head(res)
#write.csv(res,file="GO/resCUBK.HP_stat.csv",quote=F, row.names=F)


resCUBK.HA<-results(dds, contrast=c('group', 'CUxBKHot', 'CUxBKAmbient')) #here is where the two contrasting conditions get defined
mcols(resCUBK.HA,use.names=TRUE)
summary(resCUBK.HA)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1378, 15%
#LFC < 0 (down)     : 1398, 16%
table(resCUBK.HA $padj < 0.05)
#FALSE  TRUE 
# 6735  2212   
#res=data.frame(cbind("gene"=row.names(resCUBK.HA),"stat"= resCUBK.HA $stat))
#head(res)
#write.csv(res,file="GO/resCUBK.HA_stat.csv",quote=F, row.names=F)
resCUBK.HA$lp=-log(resCUBK.HA$pvalue,10)
resCUBK.HA$lp[resCUBK.HA$log2FoldChange<0]=-resCUBK.HA$lp[resCUBK.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resCUBK.HA),"lp"= resCUBK.HA$lp))
write.csv(dev_lpv,file="GO_lpvals/resCUBK.HA_lp.csv",quote=F, row.names=F)


## Flip Flop Crosses- Central Mom X North Dad ##
##DRXSB
resDRSB.AP<-results(dds, contrast=c('group', 'DRxSBAmbient', 'DRxSBPre')) #here is where the two contrasting conditions get defined
mcols(resDRSB.AP,use.names=TRUE)
summary(resDRSB.AP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 838, 9.4%
#LFC < 0 (down)     : 977, 11%
table(resDRSB.AP $padj < 0.05)
#FALSE  TRUE 
# 7559  1388   
#res=data.frame(cbind("gene"=row.names(resDRSB.AP),"stat"= resDRSB.AP $stat))
#head(res)
#write.csv(res,file="GO/resDRSB.AP_stat.csv",quote=F, row.names=F)


resDRSB.HP<-results(dds, contrast=c('group', 'DRxSBHot', 'DRxSBPre')) #here is where the two contrasting conditions get defined
mcols(resDRSB.HP,use.names=TRUE)
summary(resDRSB.HP)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1232, 14%
#LFC < 0 (down)     : 1357, 15%
table(resDRSB.HP $padj < 0.05)
#FALSE  TRUE 
# 6809  2138   
#res=data.frame(cbind("gene"=row.names(resDRSB.HP),"stat"= resDRSB.HP $stat))
#head(res)
#write.csv(res,file="GO/resDRSB.HP_stat.csv",quote=F, row.names=F)


resDRSB.HA<-results(dds, contrast=c('group', 'DRxSBHot', 'DRxSBAmbient')) #here is where the two contrasting conditions get defined
mcols(resDRSB.HA,use.names=TRUE)
summary(resDRSB.HA)
#out of 8947 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1355, 15%
#LFC < 0 (down)     : 1376, 15%
table(resDRSB.HA $padj < 0.05)
#FALSE  TRUE 
# 6754  2193   
#res=data.frame(cbind("gene"=row.names(resDRSB.HA),"stat"= resDRSB.HA $stat))
#head(res)
#write.csv(res,file="GO/resDRSB.HA_stat.csv",quote=F, row.names=F)
resDRSB.HA$lp=-log(resDRSB.HA$pvalue,10)
resDRSB.HA$lp[resDRSB.HA$log2FoldChange<0]=-resDRSB.HA$lp[resDRSB.HA$log2FoldChange<0]
dev_lpv=data.frame(cbind("gene"=row.names(resDRSB.HA),"lp"= resDRSB.HA$lp))
write.csv(dev_lpv,file="GO_lpvals/resDRSB.HA_lp.csv",quote=F, row.names=F)


save(rldG,rldG.df,colData,resBKBK.AP,resBKBK.HP,resBKBK.HA,resDRDR.AP,resDRDR.HP,resDRDR.HA,resSBSB.AP,resSBSB.HP,resSBSB.HA,resCUCU.AP,resCUCU.HP,resCUCU.HA,resLSLS.AP,resLSLS.HP,resLSLS.HA,resCUSB.AP,resCUSB.HP,resCUSB.HA,resLSCU.AP,resLSCU.HP,resLSCU.HA,resLSSB.AP,resSBBK.AP,resSBBK.HP,resSBBK.HA,resCUBK.AP,resCUBK.HP,resCUBK.HA,resDRSB.AP,resDRSB.HP,resDRSB.HA,file="AGF_larval_2018_WaldGroupContrasts_macbook.Rdata")

################ Examine differences in overall L2FC as means of transcriptional plasticity #######################

ll=load("AGF_larval_2018_WaldGroupContrasts_macbook.Rdata")
#"rldG"       "rldG.df"    "colData"    "resBKBK.AP" "resBKBK.HP" "resBKBK.HA" "resDRDR.AP" "resDRDR.HP" "resDRDR.HA" "resSBSB.AP" "resSBSB.HP" "resSBSB.HA" "resCUCU.AP" "resCUCU.HP" "resCUCU.HA" "resLSLS.AP" "resLSLS.HP" "resLSLS.HA" "resCUSB.AP" "resCUSB.HP" "resCUSB.HA" "resLSCU.AP" "resLSCU.HP" "resLSCU.HA" "resLSSB.AP" "resSBBK.AP" "resSBBK.HP" "resSBBK.HA" "resCUBK.AP" "resCUBK.HP" "resCUBK.HA" "resDRSB.AP" "resDRSB.HP" "resDRSB.HA"

FCs=(cbind(resBKBK.AP$stat, resBKBK.AP$log2FoldChange, resBKBK.HP$stat, resBKBK.HP$log2FoldChange, resBKBK.HA$stat, resBKBK.HA$log2FoldChange, resDRDR.AP$stat, resDRDR.AP$log2FoldChange, resDRDR.HP$stat, resDRDR.HP$log2FoldChange, resDRDR.HA$stat, resDRDR.HA$log2FoldChange, resSBSB.AP$stat, resSBSB.AP$log2FoldChange, resSBSB.HP$stat, resSBSB.HP$log2FoldChange, resSBSB.HA$stat, resSBSB.HA$log2FoldChange, resCUCU.AP$stat, resCUCU.AP$log2FoldChange, resCUCU.HP$stat, resCUCU.HP$log2FoldChange, resCUCU.HA$stat, resCUCU.HA$log2FoldChange, resLSLS.AP$stat, resLSLS.AP$log2FoldChange, resLSLS.HP$stat, resLSLS.HP$log2FoldChange, resLSLS.HA$stat, resLSLS.HA$log2FoldChange, resCUSB.AP$stat, resCUSB.AP$log2FoldChange, resCUSB.HP$stat, resCUSB.HP$log2FoldChange, resCUSB.HA$stat, resCUSB.HA$log2FoldChange, resLSCU.AP$stat, resLSCU.AP$log2FoldChange, resLSCU.HP$stat, resLSCU.HP$log2FoldChange, resLSCU.HA$stat, resLSCU.HA$log2FoldChange, resLSSB.AP$stat, resLSSB.AP$log2FoldChange, resSBBK.AP$stat, resSBBK.AP$log2FoldChange, resSBBK.HP$stat, resSBBK.HP$log2FoldChange, resSBBK.HA$stat, resSBBK.HA$log2FoldChange, resCUBK.AP$stat, resCUBK.AP$log2FoldChange, resCUBK.HP$stat, resCUBK.HP$log2FoldChange, resCUBK.HA$stat, resCUBK.HA$log2FoldChange, resDRSB.AP$stat, resDRSB.AP$log2FoldChange, resDRSB.HP$stat, resDRSB.HP$log2FoldChange, resDRSB.HA$stat, resDRSB.HA$log2FoldChange)) #collect pvalues for each gene

colnames(FCs)=c("resBKBK.AP_stat", "resBKBK.AP_log2FoldChange", "resBKBK.HP_stat", "resBKBK.HP_log2FoldChange", "resBKBK.HA_stat", "resBKBK.HA_log2FoldChange", "resDRDR.AP_stat", "resDRDR.AP_log2FoldChange", "resDRDR.HP_stat", "resDRDR.HP_log2FoldChange", "resDRDR.HA_stat", "resDRDR.HA_log2FoldChange", "resSBSB.AP_stat", "resSBSB.AP_log2FoldChange", "resSBSB.HP_stat", "resSBSB.HP_log2FoldChange", "resSBSB.HA_stat", "resSBSB.HA_log2FoldChange", "resCUCU.AP_stat", "resCUCU.AP_log2FoldChange", "resCUCU.HP_stat", "resCUCU.HP_log2FoldChange", "resCUCU.HA_stat", "resCUCU.HA_log2FoldChange", "resLSLS.AP_stat", "resLSLS.AP_log2FoldChange", "resLSLS.HP_stat", "resLSLS.HP_log2FoldChange", "resLSLS.HA_stat", "resLSLS.HA_log2FoldChange", "resCUSB.AP_stat", "resCUSB.AP_log2FoldChange", "resCUSB.HP_stat", "resCUSB.HP_log2FoldChange", "resCUSB.HA_stat", "resCUSB.HA_log2FoldChange", "resLSCU.AP_stat", "resLSCU.AP_log2FoldChange", "resLSCU.HP_stat", "resLSCU.HP_log2FoldChange", "resLSCU.HA_stat", "resLSCU.HA_log2FoldChange", "resLSSB.AP_stat", "resLSSB.AP_log2FoldChange", "resSBBK.AP_stat", "resSBBK.AP_log2FoldChange", "resSBBK.HP_stat", "resSBBK.HP_log2FoldChange", "resSBBK.HA_stat", "resSBBK.HA_log2FoldChange", "resCUBK.AP_stat", "resCUBK.AP_log2FoldChange", "resCUBK.HP_stat", "resCUBK.HP_log2FoldChange", "resCUBK.HA_stat", "resCUBK.HA_log2FoldChange", "resDRSB.AP_stat", "resDRSB.AP_log2FoldChange", "resDRSB.HP_stat", "resDRSB.HP_log2FoldChange", "resDRSB.HA_stat", "resDRSB.HA_log2FoldChange") 
rldpFCs=as.data.frame(cbind(rldG.df,FCs)) #combine RLD with FCs
write.csv(rldpFCs,file="rldpFCs_larvae.csv",quote=F, row.names=F)

Ps=(cbind(resBKBK.AP$pvalue, resBKBK.AP$padj, resBKBK.HP$pvalue, resBKBK.HP$padj, resBKBK.HA$pvalue, resBKBK.HA$padj, resDRDR.AP$pvalue, resDRDR.AP$padj, resDRDR.HP$pvalue, resDRDR.HP$padj, resDRDR.HA$pvalue, resDRDR.HA$padj, resSBSB.AP$pvalue, resSBSB.AP$padj, resSBSB.HP$pvalue, resSBSB.HP$padj, resSBSB.HA$pvalue, resSBSB.HA$padj, resCUCU.AP$pvalue, resCUCU.AP$padj, resCUCU.HP$pvalue, resCUCU.HP$padj, resCUCU.HA$pvalue, resCUCU.HA$padj, resLSLS.AP$pvalue, resLSLS.AP$padj, resLSLS.HP$pvalue, resLSLS.HP$padj, resLSLS.HA$pvalue, resLSLS.HA$padj, resCUSB.AP$pvalue, resCUSB.AP$padj, resCUSB.HP$pvalue, resCUSB.HP$padj, resCUSB.HA$pvalue, resCUSB.HA$padj, resLSCU.AP$pvalue, resLSCU.AP$padj, resLSCU.HP$pvalue, resLSCU.HP$padj, resLSCU.HA$pvalue, resLSCU.HA$padj, resLSSB.AP$pvalue, resLSSB.AP$padj, resSBBK.AP$pvalue, resSBBK.AP$padj, resSBBK.HP$pvalue, resSBBK.HP$padj, resSBBK.HA$pvalue, resSBBK.HA$padj, resCUBK.AP$pvalue, resCUBK.AP$padj, resCUBK.HP$pvalue, resCUBK.HP$padj, resCUBK.HA$pvalue, resCUBK.HA$padj, resDRSB.AP$pvalue, resDRSB.AP$padj, resDRSB.HP$pvalue, resDRSB.HP$padj, resDRSB.HA$pvalue, resDRSB.HA$padj)) #collect pvalues for each gene

colnames(Ps)=c("resBKBK.AP_pvalue", "resBKBK.AP_padj", "resBKBK.HP_pvalue", "resBKBK.HP_padj", "resBKBK.HA_pvalue", "resBKBK.HA_padj", "resDRDR.AP_pvalue", "resDRDR.AP_padj", "resDRDR.HP_pvalue", "resDRDR.HP_padj", "resDRDR.HA_pvalue", "resDRDR.HA_padj", "resSBSB.AP_pvalue", "resSBSB.AP_padj", "resSBSB.HP_pvalue", "resSBSB.HP_padj", "resSBSB.HA_pvalue", "resSBSB.HA_padj", "resCUCU.AP_pvalue", "resCUCU.AP_padj", "resCUCU.HP_pvalue", "resCUCU.HP_padj", "resCUCU.HA_pvalue", "resCUCU.HA_padj", "resLSLS.AP_pvalue", "resLSLS.AP_padj", "resLSLS.HP_pvalue", "resLSLS.HP_padj", "resLSLS.HA_pvalue", "resLSLS.HA_padj", "resCUSB.AP_pvalue", "resCUSB.AP_padj", "resCUSB.HP_pvalue", "resCUSB.HP_padj", "resCUSB.HA_pvalue", "resCUSB.HA_padj", "resLSCU.AP_pvalue", "resLSCU.AP_padj", "resLSCU.HP_pvalue", "resLSCU.HP_padj", "resLSCU.HA_pvalue", "resLSCU.HA_padj", "resLSSB.AP_pvalue", "resLSSB.AP_padj", "resSBBK.AP_pvalue", "resSBBK.AP_padj", "resSBBK.HP_pvalue", "resSBBK.HP_padj", "resSBBK.HA_pvalue", "resSBBK.HA_padj", "resCUBK.AP_pvalue", "resCUBK.AP_padj", "resCUBK.HP_pvalue", "resCUBK.HP_padj", "resCUBK.HA_pvalue", "resCUBK.HA_padj", "resDRSB.AP_pvalue", "resDRSB.AP_padj", "resDRSB.HP_pvalue", "resDRSB.HP_padj", "resDRSB.HA_pvalue", "resDRSB.HA_padj") 
rldPs=as.data.frame(cbind(rldG.df,Ps))

write.csv(rldPs,file="rldpPs_larvae.csv",quote=F, row.names=T)

##### pick out genes specific to CU dam heat response ############
table(resCUBK.HA$padj<= 0.05) #2213
table(resCUCU.HA$padj<= 0.05) #2258
table(resCUSB.HA$padj<= 0.05) #2840

CUgenes=rldPs[rldPs$resCUCU.HA_padj<=0.05 & rldPs$resCUSB.HA_padj<=0.05 & rldPs$resCUBK.HA_padj<=0.05, ] #797 genes

CUgenes=rldPs[rldPs$resCUCU.HA_padj<=0.05 & rldPs$resCUSB.HA_padj<=0.05 & rldPs$resCUBK.HA_padj<=0.05 & rldPs$resBKBK.HA_padj>0.05 & rldPs$resDRDR.HA_padj>0.05 & rldPs$resSBSB.HA_padj>0.05 & rldPs$resLSLS.HA_padj>0.05 & rldPs$resLSCU.HA_padj>0.05 & rldPs$resSBBK.HA_padj>0.05 & rldPs$resDRSB.HA_padj>0.05, ] #6 genes

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
CUdam[CUdam ==T]<- 1
CUdam[CUdam ==F]<- 0
geneL$CUdam = as.factor(CUdam)


geneL$treat_f = factor(geneL$treat, levels=c("Pre","Amb","Hot"))
#geneL$cross_f = factor(geneL$cross, levels=c('LSSB','LSLS','DRDR','SBSB','LSCU','DRSB','BKBK','CUCU','CUSB','CUBK','SBBK'))

summ=summarySE(data= geneL,measurevar="expr",groupvars=c("gene","treat_f","CUdam"))
summ=summarySE(data= geneL,measurevar="expr",groupvars=c("gene","treat_f","cross"))

pd <- position_dodge(0.15)
ggplot(summ,aes(x=treat_f,y=expr, colour= cross))+
	geom_point(aes(group= cross),position=pd,size=2)+
	geom_line(aes(group= cross,linetype= cross),position=pd)+
	geom_errorbar(aes(ymin=expr-se,ymax=expr+se),lwd=0.4,width=0.3,position=pd)+
	ggtitle("CUdam specific heat response genes")+
	theme_bw()+
	facet_wrap(~gene,scales="free_y")



################ plot genes with high Fst SNPs for CU pop
rldPs=read.csv("rldpPs_larvae.csv", header=T)
row.names(rldPs)=rldPs$X

Cands1=rldPs[grepl("Amillepora12232",rownames(rldPs)), ]
Cands2=rldPs[grepl("Amillepora25324",rownames(rldPs)), ]
Cands3=rldPs[grepl("Amillepora33421",rownames(rldPs)), ]

Cands=rbind(Cands1,Cands2,Cands3)

genes=as.data.frame(Cands[-c(1,98:159) ])
names(genes)=paste0(colData$group, colData$sample)
genes$gene=row.names(genes)

geneL=genes %>%
	gather("sample", "expr", -gene)
geneL$cross = str_sub(geneL$sample, 1,5)
geneL$treat = str_sub(geneL$sample, 6,8)
geneL$dam = str_sub(geneL$sample, 1,2)

geneL$CU<- with(geneL, ifelse(geneL$dam == "CU", geneL$cross, "noCU"))


geneL$treat_f = factor(geneL$treat, levels=c("Pre","Amb","Hot"))
#geneL$cross_f = factor(geneL$cross, levels=c('LSSB','LSLS','DRDR','SBSB','LSCU','DRSB','BKBK','CUCU','CUSB','CUBK','SBBK'))

summ=summarySE(data=geneL,measurevar="expr",groupvars=c("gene","treat_f","CU"))
summ=summarySE(data=geneL,measurevar="expr",groupvars=c("gene","treat_f","cross"))

#CUxBK, filled in square, 15, CUCU, filled in diamond 18, CUSB, filled in triangle 17
pd <- position_dodge(0.1)
ggplot(summ,aes(x=treat_f,y=expr))+
	geom_point(aes(shape= CU),size=3)+
	geom_line(aes(group= CU),position=pd)+
	geom_errorbar(aes(ymin=expr-se,ymax=expr+se),lwd=0.4,width=0.3,position=pd)+
	scale_shape_manual(values=c(15,18,17,1))+
	ggtitle("Genes with high Fst SNPs differentiating CU pop")+
	theme_bw()+
	facet_wrap(~gene,scales="free_y")


######################################################## run DESeq2 on only the PRE samples
gcountsALL=read.table("AGF_larvae_2018_geneCounts.txt", header=T) 

# clean up file names
colnames(gcountsALL)<-sub("X","", colnames(gcountsALL))
colnames(gcountsALL) <- sub("_[^_]+$", "", colnames(gcountsALL))
colnames(gcountsALL)<-gsub("\\.","_", colnames(gcountsALL))

length(gcountsALL[,1]) #25142
dim(gcountsALL) #25142    96


summary(colSums(gcountsALL))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  79710  795300  871500  849100  917300 1102000 


### BUILD A DATAFRAME ASSOCIATING SAMPLE NAMESWITH TREATMENT CONDITIONS ###
colData=read.csv("J19188metaPRE.csv", header=T)
colData=colData[!duplicated(colData$sample), ]
rownames(colData)<-colData$sample
dim(colData) #30  3

all(rownames(colData) %in% colnames(gcountsALL)) #TRUE
gcounts <- gcountsALL[, rownames(colData)]
all(rownames(colData) == colnames(gcounts))


### REMOVE GENES WITH LOW MEAN COUNTS ###

mns = apply(gcounts, 1, mean)

gcountsPRE=gcounts[mns>10,] #get rid of genes that show little or no expression
#gcounts=gcounts[mns>10,] #only highly expressed genes for WGCNA

table(mns > 10)
#FALSE  TRUE 
#16320  8822
dim(gcountsPRE) #8822   30

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
out of 8822 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)     : 196, 2.2% 
LFC < 0 (down)   : 154, 1.7% 
outliers [1]     : 2, 0.023% 
low counts [2]   : 0, 0% 

table(resP$padj < 0.05)
#FALSE  TRUE 
# 8613   207   
res=data.frame(cbind("gene"=row.names(resP),"stat"= resP$stat))
head(res)
write.csv(res,file="GO/resP_stat.csv",quote=F, row.names=F)

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
write.csv(survival_fishers,file="GO/survival_LRT_fishers_pval0.01.csv",quote=F, row.names=F)

save(resP,file="AGF_larval_2018_PreRes.Rdata")

### heatmap

vals=(cbind(resP$stat, resP$padj)) #collect pvalues for each gene
colnames(vals)=c("stat_resP","padj_resP") 
rldpvals=as.data.frame(cbind(rldLRT.df,vals)) #combine RLD with pvals

annot=read.table("Amillepora_trinotate_annotation_report_edited.txt",header=F,sep="\t",quote=NULL,fill=T)
annot=as.data.frame(annot[c(2:3) ],)
colnames(annot)=c("transcript","gene") 
genes = row.names(rldpvals)
genes2annot = match(genes,annot$transcript)
rldpvals_annot=data.frame(GeneID=row.names(rldpvals), annot[genes2annot,], rldpvals)

rldpvals_annot = rldpvals_annot[rldpvals_annot $padj_resP<= 0.01,]
rldpvals_annot = rldpvals_annot[complete.cases(rldpvals_annot),]
dim(rldpvals_annot) #207  32
row.names(rldpvals_annot)=paste0(rldpvals_annot$Gene_ID,rldpvals_annot$gene)
rldpvals_annot = rldpvals_annot[, 4:33]

means=apply(rldpvals_annot,1,mean) # means of rows
explc= rldpvals_annot-means


heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=1)(100)
heat.colors = colorRampPalette(wes_palette("Zissou1", 30, type = "continuous"), bias=1.1)(100)
pdf("Heatmap_genes_survivalPre.pdf",height=10,width=12)
pheatmap(explc,color=heat.colors,cluster_cols=F,border_color=NA,clustering_distance_rows="correlation")
dev.off()
pheatmap(explc,color=heat.colors,cluster_cols=T,border_color=NA)



### principal coordinate calculation ###

conditions=colData
col3 <- wes_palette("Zissou1", 3, type = "continuous")
col11 <- wes_palette("FantasticFox1", 11, type = "continuous")

conditions$treat_c=scale_colour_manual(values=c("Pre"="#EBCC2A","Ambient"="#3B9AB2", "Hot"="#F21A00"))
grp=rep("#EBCC2A",ncol(rldG.df))
grp[grep("Ambient",conditions$treatment)]="#3B9AB2"
grp[grep("Hot",conditions$treatment)]="#F21A00"



fitpc=capscale(dist(t(rldG.df),method="manhattan")~1)
summary(eigenvals(fitpc))

plot(fitpc$CA$u,pch=16, col="white",main="AGF Larvae 2018", xlab="MDS1 26.8%", ylab="MDS2 15.9%")
abline(h=0,v=0, col="grey")
points(fitpc$CA$u,pch=16, col=grp)
ordispider(fitpc$CA$u,conditions$treat,col="black", label=T)
#ordihull(fitpc$CA$u,conditions$treat,draw="polygon",label=T, col=grp)

plot(fitpc$CA$u,pch=16, col="white",main="AGF Larvae 2018", xlab="MDS1 26.8%", ylab="MDS2 15.9%")
abline(h=0,v=0, col="grey")
points(fitpc$CA$u,pch=16, col=grp)
ordispider(fitpc$CA$u,conditions$cross,col="black", label=T)
#ordiellipse(fitpc$CA$u,conditions$cross,draw="polygon",label=T)

ad=adonis(t(rldG.df)~treatment+cross,data=conditions,method="manhattan")
#Call:
#adonis(formula = t(rldA.df) ~ treatment + cross, data = conditions,      method = "manhattan") 
#
#Permutation: free
#Number of permutations: 999
#
#Terms added sequentially (first to last)
#
#          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#treatment  2  96503203 48251602 27.2906 0.28829  0.001 ***
#cross     10  91492082  9149208  5.1747 0.27332  0.001 ***
#Residuals 83 146749300  1768064         0.43839           
#Total     95 334744586                  1.00000           



labs=c("Treatment","Genotype","Residuals")
cols=c("coral2","turquoise4","grey80")
#cols = append(gg_color_hue(3), 'grey')
labs2 = paste(labs, round(ad$aov.tab$R2[1:3]*100, digits=1))
pie(ad$aov.tab$R2[1:3],labels=labs2,col=cols,main="AGF larvae 2018")







