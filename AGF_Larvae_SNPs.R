library(vcfR)
library(adegenet)
library(pcadapt)
library(ggplot2)
library(pegas)
library (ape)
library(hierfstat)
library(dartR)
library(wesanderson)
library(pheatmap)
#library('dplyr')
library(stringr)
library(vegan)


#loading all samples together, and all species separate
setwd("/Users/mes0192/Dropbox/KQ_GBR_spawnseq/analysis/DESeq2_larvae")

### analysis with all samples
colData.l=read.csv("J19188meta.csv", header=T)
#colData.l=read.csv("J19188meta_kqedited.csv", header=T) 
colData.l=colData.l[!duplicated(colData.l$sample), ]
dim(colData.l) #96 correct
rownames(colData.l)<-colData.l$sample
colData.l$group=as.factor(paste0(colData.l$cross,colData.l$treatment))
rownames(colData.l)<-sub("_","-", rownames(colData.l))
colData.l$sample<-sub("_","-", colData.l$sample)
dim(colData.l) #96  3


##Adegenet DAPC
aten.l<-read.vcfR('./SNPs/Aten_larval_SNPs_filtered_all_1.recode.vcf')
aten_gl.l=vcfR2genlight(aten.l) #converts the matrix to Genlight format for adegenet
aten_gl.l@ind.names<- sub("_[^_]+$", "", aten_gl.l@ind.names)
aten_gl.l@ind.names<- sub("_[^_]+$", "", aten_gl.l@ind.names)


all(rownames(colData.l) %in% aten_gl.l@ind.names) #TRUE
all(rownames(colData.l) == aten_gl.l@ind.names) #FALSE
colData.l <- colData.l[aten_gl.l@ind.names,]
all(rownames(colData.l) == aten_gl.l@ind.names) #TRUE

hotlarvae=colData.l[colData.l$treatment=="Hot",]

#subset the .gl file for just the hot larvae
hotlarvae.gl <- gl.keep.ind(aten_gl.l, ind.list=row.names(hotlarvae), recalc=T, mono.rm=T)

#call pop for all samples
pop(aten_gl.l) <- as.factor(colData.l$cross)

#call pop for hot samples
pop(hotlarvae.gl) <- as.factor(hotlarvae $cross)

mat.l=as.matrix(aten_gl.l)

#call PC for all larvae
pc.l <- gl.pcoa(aten_gl.l)

#call PC for hot larvae
pc.l <- gl.pcoa(hotlarvae.gl)


scatter(pc.l, posi="bottomright")

pcoa.l=gl.pcoa.plot(pc.l, hotlarvae.gl, labels="pop", xaxis=1, yaxis=2)
pcoa.l+theme_bw()

all(rownames(hotlarvae) %in% hotlarvae.gl@ind.names) #TRUE
all(rownames(hotlarvae) == hotlarvae.gl@ind.names) #FALSE
hotlarvae <- hotlarvae[hotlarvae.gl@ind.names,]
all(rownames(hotlarvae) == hotlarvae.gl@ind.names) #TRUE


dapc1.l<-dapc(hotlarvae.gl, hotlarvae$group, perc.pca=80, n.da=2) # retained 80PCs , Chose 2DF 

varexpl.l <- round((dapc1.l$eig/sum(dapc1.l$eig))[1:2] * 100, 1) #62.4 37.6

scatter(dapc1.l, bg = "white", legend = TRUE, scree.da = FALSE)

larval_load_SNPs=as.data.frame(dapc1.l$ind.coord)
larval_load_SNPs$sample<-row.names(larval_load_SNPs)
colnames(larval_load_SNPs)=c("SNP_dapc1","SNP_dapc2","sample")
save(larval_load_SNPs,file="DAPC_larval_loadingsHOT_SNPs.Rdata")

dapc2.l <- tibble(sample = rownames(dapc1.l$ind.coord),
               grp = dapc1.l$grp,
               LD1 = dapc1.l$ind.coord[,1],
               LD2 = dapc1.l$ind.coord[,2])
dapc3.l <- dapc2.l %>%
  group_by(grp) %>%
  summarize(c1 = mean(LD1),
            c2 = mean(LD2)) %>%
  			full_join(dapc2.l, by="sample")
        #full_join(dapc2.l)

dapc3.l<-dapc3.l%>%
  full_join(dapc2.l)

dapc.fig.l <-   ggplot(dapc3.l, aes(shape = factor(str_sub(grp, 1,5)), color = factor(str_sub(grp, 6,8)),
		fill = factor(str_sub(grp, 6,8)))) +
		geom_segment(mapping = aes(x = LD1, y = LD2, xend = c1, yend = c2), lwd = 0.25, col = "grey") +
		geom_point(aes(x = c1, y = c2), size = 3) +
		#geom_point(aes(x = LD1, y = LD2), size = 2, show.legend = FALSE) +
		scale_shape_manual(name = "Cross", values=c(21,22,23,24,25,0,1,2,3,8,5)) +
		scale_fill_manual(name="Treatment", values = c("#3B9AB2", "#F21A00","#EBCC2A"))+
		scale_color_manual(name="Treatment", values = c("#3B9AB2", "#F21A00","#EBCC2A"))+
		guides(fill = guide_legend(override.aes = list(shape = 21, size = 3))) +
		guides(shape = guide_legend(override.aes = list(fill = "black", size = 3))) +
		labs(x = paste0("LD1 [", varexpl[1],"%]"), y = paste0("LD2 [", varexpl[2],"%]")) +
		theme_bw()

##### Load in larval GE results

ll.larv=load("DAPC_larval_loadings_hotGE.Rdata")
larval_load_GE$sample<-sub("_","-", larval_load_GE$sample)
ll.larv.snps=load("DAPC_larval_loadingsHOT_SNPs.Rdata")


dapc_loadings.l=merge(larval_load_GE,larval_load_SNPs,by="sample")

traits.l=merge(dapc_loadings.l,colData.l,by="sample")

Hot.l=traits.l #30 samples correct
rownames(Hot.l)=make.names(paste0(Hot.l$cross,Hot.l$treatment), unique=T)

Hot.l$perc=(Hot.l$survivalheat)/100
Hot.l <- Hot.l[order(Hot.l$perc),] #order by survival in heat (low survival to high survival)
dim(Hot.l)

gm=glm(survivlheatperc ~ GE_dapc1+ SNP_dapc1,data=Hot,family=binomial(link = "logit"))
summary(gm)
#Warning message: In eval(family$initialize) : non-integer #successes in a binomial glm!
#Call:
#glm(formula = survivlheatperc ~ GE_dapc1 + SNP_dapc1, family = binomial(link = "logit"), 
#    data = Hot)

#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-0.69959  -0.21122   0.01853   0.18521   0.60968  
#
#Coefficients:
#              Estimate Std. Error z value Pr(>|z|)  
#(Intercept) -4.452e+00  2.332e+00  -1.909   0.0563 .
#GE_dapc1     3.540e-01  1.739e-01   2.036   0.0418 *
#SNP_dapc1   -3.411e-17  4.552e-17  -0.749   0.4536 
#survival.pca<- prcomp(Hot[,c(2:5,8)],center=T, scale.=T)
#summary(survival.pca)
#biplot(survival.pca, scale=0, pch=16)

#missing larval samples (not including pre) #should be 189, but only 96 in dataset
larv_ambo_Hot<-Hot.l%>% filter(treatment == c("Hot", "Ambient"))
glm_larv=glmer(perc ~ GE_dapc1+ SNP_dapc1+(1|treatment), data=larv_ambo_Hot,family=binomial(link = "logit")) #scale issues w GE and SNP
#glm_larv2=glmer(perc ~ GE_dapc1+ SNP_dapc1+(1|cross), data=larv_ambo_Hot,family=binomial(link = "logit"))
#glm_larv3=glm(perc ~ GE_dapc1+ SNP_dapc1, data=larv_ambo_Hot,family=binomial(link = "logit"))
summary(glm_larv)
plot(glm_larv)

#rescale
#simple rescale
library(dplyr)
library(lme4)
#########################
#Loadings from PCA made from ambient and hot samples ***
#glm_larv.rescale=glmer(perc ~ log(GE_dapc1+1)+ SNP_dapc1+(1|treatment), data=larv_ambo_Hot,family=binomial(link = "logit"))
larv_ambo_Hot_rescale<-larv_ambo_Hot %>% 
  mutate(zscore.GE = (GE_dapc1 - mean(GE_dapc1))/sd(GE_dapc1))%>% 
  mutate(zscore.SNP = (SNP_dapc1 - mean(SNP_dapc1))/sd(SNP_dapc1))

glm_larv.rescale_final=glmer(perc ~ zscore.GE+ zscore.SNP+(1|treatment), data=larv_ambo_Hot_rescale,family=binomial(link = "logit"))
plot(glm_larv.rescale_final) #much better resids
#########################
#Loadings from PCA made from only hot samples***
larv_onlyHotDAPC_rescale<-Hot.l %>% 
  mutate(zscore.GE = (GE_dapc1 - mean(GE_dapc1))/sd(GE_dapc1))%>% 
  mutate(zscore.SNP = (SNP_dapc1 - mean(SNP_dapc1))/sd(SNP_dapc1))

glm_larv.rescale_onlyhotDAPC=glmer(perc ~ zscore.GE+ zscore.SNP+(1|cross), data=larv_onlyHotDAPC_rescale,family=binomial(link = "logit"))
plot(glm_larv.rescale_onlyhotDAPC) #this looks good, now evidence of heterogen of variance
############################

X.var.larv.resc<- var(as.vector(lme4::fixef(glm_larv.rescale_final) %*% t(glm_larv.rescale_final@pp$X))) #2.237577
#lme4::fixef(glm_larv.rescale_final)
#(Intercept)   zscore.GE  zscore.SNP 
#2.0150666  -1.5174796   0.5371215 

#GE_dapc1 (-(1.5174796/2.237577) *100 =67.81798% of X.var which is mostly fixed effects as R2.marginal shows
#SNP_dapc1(0.5371215 /2.237577)*100 =24.0%

#####
X.var.larv.resc.hot<- var(as.vector(lme4::fixef(glm_larv.rescale_onlyhotDAPC) %*% t(glm_larv.rescale_onlyhotDAPC@pp$X))) #1.636799
#these are the estimates, not stdevs or the variance components
#lme4::fixef(glm_larv.rescale_onlyhotDAPC)
#(Intercept)   zscore.GE  zscore.SNP 
#0.63297291  1.20771478 -0.08952766 
#GE_dapc1 ((1.20771478/1.636799) *100) =73.78516% of X.var which is mostly fixed effects as R2.marginal shows
#SNP_dapc1(-(0.08952766/1.636799)*100 =5.46968%


# Variance-Covariance Matrix of fixed effects: 
vc_fixed.l <- as.matrix(vcov(glm_larv.rescale_onlyhotDAPC))
# Variance of fixed effects: 
var_fixed.l <- diag(vc_fixed.l); var_fixed.l
#(Intercept)   zscore.GE  zscore.SNP 
#0.3649382   0.8542799   0.9470256 
#can also run without intercept at -1
#(0.3649382+0.8542799+0.9470256) =2.166244
#zscore.GE: (0.8542799/2.166244)*100 =39.436%
#zscore.SNP: (0.9470256/2.166244)*100 =43.7174%

## Extract the variance components for the random effects (not including the residuals)
Z.var.larv <- sum(
  sapply(
    VarCorr(glm_larv.rescale_final)[!sapply(unique(unlist(strsplit(names(ranef(glm_larv.rescale_final)),":|/"))), function(l)
      length(unique(glm_larv.rescale_final@frame[,l])) == nrow(glm_larv.rescale_final@frame))],
    function(Sigma) {
      X <- model.matrix(glm_larv.rescale_final)
      Z <- X[,rownames(Sigma)]
      sum(diag(Z %*% Sigma %*% t(Z)))/nrow(X) } ) )
#1.135679e-09
Z.var.larv.heat <- sum(
  sapply(
    VarCorr(glm_larv.rescale_onlyhotDAPC)[!sapply(unique(unlist(strsplit(names(ranef(glm_larv.rescale_onlyhotDAPC)),":|/"))), function(l)
      length(unique(glm_larv.rescale_onlyhotDAPC@frame[,l])) == nrow(glm_larv.rescale_onlyhotDAPC@frame))],
    function(Sigma) {
      X <- model.matrix(glm_larv.rescale_onlyhotDAPC)
      Z <- X[,rownames(Sigma)]
      sum(diag(Z %*% Sigma %*% t(Z)))/nrow(X) } ) )
#Z.var.larv.heat #1.004311

## Extract the variance componts for the residuals
R.var.larv <- attr(lme4::VarCorr(glm_larv.rescale_final), "sc")^2
#1
R.var.larv.heat <- attr(lme4::VarCorr(glm_larv.rescale_onlyhotDAPC), "sc")^2
#1

## The marginal R2 (proportion of variance due to fixed effects)
(R2.marginal.larv <- X.var.larv.resc/(X.var.larv.resc+Z.var.larv+R.var.larv))
# 0.691127 x 100 % of variance due to fixed effects
(R2.marginal.larv.heat <- X.var.larv.resc.hot/(X.var.larv.resc.hot+Z.var.larv.heat+R.var.larv.heat))
#0.4495329

## The proportion of variance due to random effects
(R2.random.larv <- Z.var.larv/(X.var.larv.resc+Z.var.larv+R.var.larv))
# 3.887425e-09 x 100% % (treatment and and obs)
(R2.random.larv.heat <- Z.var.larv.heat/(X.var.larv.resc.hot+Z.var.larv.heat+R.var.larv.heat))
#0.2758256

########
## The proportion of variance due to residuals
VarCorr(glm_larv.rescale_final) 
#only treatment temp
#0.00011219 stdev ==100% of resid bc only one factor here

#So treatmentis 100%% of total 69%% 

VarCorr(glm_larv.rescale_onlyhotDAPC)
#cross is 1.0022 stdev, this is 100% of residual bc only one factor, which is cross
########

(R2.resid.larv <- R.var.larv/(X.var.larv.resc+Z.var.larv+R.var.larv))
# 0.308873 x 100% just due to residuals
(R2.resid.larv.hot <- R.var.larv.heat/(X.var.larv.resc.hot+Z.var.larv.heat+R.var.larv.heat))
#0.2746415 x 100% just due to residuals


## The conditional R2 (proportion of variance due to fixed and random effects)
(R2.conditional.larv <- (X.var.larv.resc+Z.var.larv)/(X.var.larv.resc+Z.var.larv+R.var.larv)) #0.691127% = 69% this is total variabiliy expalined by all
#0.691127% is total from R2.conditional

## The conditional R2 (proportion of variance due to fixed and random effects)
(R2.conditional.larv.hot <- (X.var.larv.resc.hot+Z.var.larv.heat)/(X.var.larv.resc.hot+Z.var.larv.heat+R.var.larv.heat)) #0.691127% = 69% this is total variabiliy expalined by all
#0.7253585% is total from R2.conditional

#F + R =Conditional
#F + R= 73%
#Fixed= 44.9%
#73% (conditional) -44.9% (marginal, only fixed) = 28.1% only random
#73-44.9
#27 residuals 

##########################################################
surv=Hot.l[,c(2:5,8)]
survivalheat=surv$survivalheat

fitpc=rda(surv ~ survivalheat,scale=T)
axes2plot=c(2,3)
plot(fitpc,choices=axes2plot)
str(fitpc)
pc1expl=100*round(summary(eigenvals(fitpc))[2, axes2plot[1]],2)
pc2expl=100*round(summary(eigenvals(fitpc))[2, axes2plot[2]],2)
heat.colors = wes_palette("Zissou1", 11, type = "continuous")


plot(fitpc$CA$u,pch=16, cex=1,col= heat.colors,xlim=c(-0.6,0.6),ylim=c(-0.5,0.7),xlab=paste("PC1 ( ",pc1expl,"% )",sep=""),ylab=paste("PC2 ( ",pc2expl,"% )",sep=""),mgp=c(2.3,1,0),main="Larval Survival in Heat")
arrows(rep(0,4),rep(0,4), fitpc$CA$v[,1]/2, fitpc $CA$v[,2]/1.5,length=0.05)
text(fitpc$CA$v[,1]/1.3, fitpc$CA$v[,2]/1.3+c(0,0,0,0),labels=row.names(fitpc$CA$v) ,cex=0.8,col=c(rep(1,4),rep(4,3)))

summary(lm(surv$GE_dapc1~fitpc$CA$u[,1])) #p=0.000581
summary(lm(surv$SNP_dapc1~fitpc$CA$u[,1])) #p=4.18e-07



######################################### extract SNPs differentiating CU from other pops 
##SNP file for only the purebred larval crosses
aten<-read.vcfR('./SNPs/Aten_larval_SNPs_filtered_1.recode.vcf')

aten_gl=vcfR2genlight(aten) #converts the matrix to Genlight format for adegenet
aten_gl@ind.names<- sub("_[^_]+$", "", aten_gl@ind.names)
aten_gl@ind.names<- sub("_[^_]+$", "", aten_gl@ind.names)
snps=aten_gl@loc.names
write.csv(snps,file="loc.names.csv",quote=F, row.names=F)
aten.gen<-vcfR2genind(aten)
pop(aten.gen) <-c(rep("DR",3),rep("SB",3),rep("BK",3),rep("LS",3),rep("BK",6),rep("DR",3),rep("CU",3),rep("SB",3),rep("LS",3),rep("DR",3),rep("CU",6),rep("LS",3),rep("SB",3))

pc <- gl.pcoa(aten_gl, nfactors=5)

heat.colors = wes_palette("Zissou1", 5, type = "continuous")
pcoa=gl.pcoa.plot(pc, aten_gl, labels="pop", xaxis=1, yaxis=2)
pcoa+theme_bw() + 
	scale_color_manual(values=c("BK"="#3B9AB2", "DR"="#78B7C5", "LS"="#EBCC2A","SB"="#E1AF00","CU"="#F21A00"))

gl.pcoa.scree(pc)

grp <- find.clusters(aten.gen, max.n.clust=44)


dapc1<-dapc(aten.gen, grp$grp,n.pca=40, n.da=10) # retained 40PCs , Chose 2DF 
#$n.pca: 40 first PCs of PCA used
#$n.da: 4 discriminant functions saved
#$var (proportion of conserved variance): 0.939

set.seed(4) 
quantile(dapc1$var.contr,0.99) #0.0004698106
contrib <- loadingplot(dapc1$var.contr, axis=1,thres= 0.0004718359, lab.jitter=1)
sig_snps=contrib$var.values
sig_snps=sort(sig_snps, decreasing=T)
length(sig_snps) #240 for 99 percentile, 5332 for 75 percentile

Amil20691 <- tab(genind2genpop(aten.gen[loc=c("Amillepora20691-RA_311")]),freq=TRUE)
Amil08427 <- tab(genind2genpop(aten.gen[loc=c("Amillepora08427-RA_3148")]),freq=TRUE)
Amil16135 <- tab(genind2genpop(aten.gen[loc=c("Amillepora16135-RA_1788")]),freq=TRUE)
Amil14669 <- tab(genind2genpop(aten.gen[loc=c("Amillepora14669-RA_21")]),freq=TRUE)
Amil04950 <- tab(genind2genpop(aten.gen[loc=c("Amillepora04950-RA_2219")]),freq=TRUE)
Amil35676 <- tab(genind2genpop(aten.gen[loc=c("Amillepora35676-RA_1827")]),freq=TRUE)

par(mfrow=c(2,3), mar=c(5.1,4.1,4.1,.1),las=3) 
matplot(Amil20691, pch=c("1","0"), type="b", xlab="pop",ylab="allele frequency", xaxt="n", cex=1.5, main="Amillepora20691-RA_311") 
axis(side=1, at=1:5, lab=c("DR","SB","BK","LS","CU"))
matplot(Amil08427, pch=c("1","0"), type="b", xlab="pop", ylab="allele frequency", xaxt="n", cex=1.5, main="Amillepora08427-RA_3148") 
axis(side=1, at=1:5, lab=c("DR","SB","BK","LS","CU"))
matplot(Amil16135, pch=c("1","0"), type="b", xlab="pop", ylab="allele frequency", xaxt="n", cex=1.5, main="Amillepora16135-RA_1788") 
axis(side=1, at=1:5, lab=c("DR","SB","BK","LS","CU"))
matplot(Amil14669, pch=c("1","0"), type="b", xlab="pop",ylab="allele frequency", xaxt="n", cex=1.5, main="Amillepora14669-RA_21") 
axis(side=1, at=1:5, lab=c("DR","SB","BK","LS","CU"))
matplot(Amil04950, pch=c("1","0"), type="b", xlab="pop", ylab="allele frequency", xaxt="n", cex=1.5, main="Amillepora04950-RA_2219") 
axis(side=1, at=1:5, lab=c("DR","SB","BK","LS","CU"))
matplot(Amil35676, pch=c("1","0"), type="b", xlab="pop", ylab="allele frequency", xaxt="n", cex=1.5, main="Amillepora35676-RA_1827") 
axis(side=1, at=1:5, lab=c("DR","SB","BK","LS","CU"))

compoplot(dapc1, col=myPal,lab="", txt.leg=paste("group", 1:5), ncol=2)
compoplot(dapc1, col=myPal, txt.leg=paste("group", 1:5), ncol=2, posi="bottomright",xlab="individuals") #test of composition plot similar to structure 
scatter(dapc1, posi.da="bottomright", bg="white", pch=17:22)

tre <- nj(dist(as.matrix(aten_gl)))
plot(tre, typ="fan", cex=0.7, main="NJ distance tree Aten purebred larvae")

x.dist<-dist(aten_gl)
heat.colors = wes_palette("Zissou1", 100, type = "continuous")
pheatmap(as.matrix(x.dist),cex=1.2,color=heat.colors,border_color="NA",cluster_cols=T,cluster_rows=T)



