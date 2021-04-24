library(vcfR)
library(adegenet)
library(pcadapt)
library(ggplot2)
#library(pegas)
#library (ape)
#library(hierfstat)
library(dartR)
library(wesanderson)
#library(pheatmap)
#library('dplyr')
library(stringr)
library(vegan)
library(tidyverse)

colData=read.csv("J19234meta.csv", header=T) 
colData=colData[!duplicated(colData$sample), ]
rownames(colData)<-colData$sample
dim(colData) #146  3
colData$group1=factor(paste0(colData$Cross,colData$Zoox))
colData$group2=factor(paste0(colData$Zoox,colData$Treatment.1))
colData$group3=factor(paste0(colData$group1,colData$Treatment.1))


# PCAadapt
aten<-read.pcadapt('./SNPs/Aten_juves_SNPs_filtered_all_1.recode.vcf', type='vcf')
x <- pcadapt(aten, K = 20) 
plot(x, option='screeplot', K=20) #suggests best k=2
pop.aten=colData$Cross
print(pop.aten)
summary(singular.values(x))
plot(x, option = "scores", pop=colData$Cross,pch=16,col="white",main="A.tenuis larvae", xlab="PC1  21.26%", ylab="PC2 10.56%") #PC1=21.02% PC2=19.98
abline(h=0,v=0, col="grey")
text(x, option="scores", cex=0.9, font=1)
head(x$scores) 

##Adegenet DAPC
aten<-read.vcfR('./SNPs/Aten_juves_SNPs_filtered_all_1.recode.vcf')
aten_gl=vcfR2genlight(aten) #converts the matrix to Genlight format for adegenet
aten_gl@ind.names<- sub("_[^_]+$", "", aten_gl@ind.names)
aten_gl@ind.names<- sub("_[^_]+$", "", aten_gl@ind.names)


all(rownames(colData) %in% aten_gl@ind.names) #TRUE
all(rownames(colData) == aten_gl@ind.names) #FALSE
colData <- colData[aten_gl@ind.names,]
all(rownames(colData) == aten_gl@ind.names) #TRUE

hotjuves=colData[colData$Treatment.1=="Hot",]

#subset the .gl file for just the hot larvae
hotjuves.gl <- gl.keep.ind(aten_gl, ind.list=row.names(hotjuves),recalc=T, mono.rm=T)


pop(hotjuves.gl) <- as.factor(hotjuves$Cross)
mat=as.matrix(hotjuves.gl)

pc <- gl.pcoa(hotjuves.gl)

scatter(pc, posi="bottomright")


pcoa=gl.pcoa.plot(pc, hotjuves.gl, labels="pop", xaxis=1, yaxis=2)
pcoa+theme_bw()


dapc1<-dapc(hotjuves.gl, hotjuves$group3, perc.pca = 80, n.da=2) # retained 80PCs , Chose 2DF 

varexpl <- round((dapc1$eig/sum(dapc1$eig))[1:2] * 100, 1) #43.2 33.5

scatter(dapc1, bg = "white", legend = TRUE, scree.da = FALSE)

juves_load_SNPs=as.data.frame(dapc1$ind.coord)
juves_load_SNPs$sample<-row.names(juves_load_SNPs)
colnames(juves_load_SNPs)=c("SNP_dapc1","SNP_dapc2","sample")
save(juves_load_SNPs,file="DAPC_juvenile_SNP_loadings.Rdata")

dapc2 <- tibble(sample = rownames(dapc1$ind.coord),
               grp = dapc1$grp,
               LD1 = dapc1$ind.coord[,1],
               LD2 = dapc1$ind.coord[,2])
dapc3 <- dapc2 %>%
  group_by(grp) %>%
  summarize(c1 = mean(LD1),
            c2 = mean(LD2)) %>%
  			full_join(dapc2)

dapc.fig <-   ggplot(dapc3, aes(shape = factor(str_sub(grp, 1,5)), color = factor(str_sub(grp, 6,8)),
		fill = factor(str_sub(grp, 6,8)))) +
		geom_segment(mapping = aes(x = LD1, y = LD2, xend = c1, yend = c2), lwd = 0.25, col = "grey") +
		geom_point(aes(x = c1, y = c2), size = 3) +
		#geom_point(aes(x = LD1, y = LD2), size = 2, show.legend = FALSE) +
		scale_shape_manual(name = "Cross", values=c(21,22,23,24,25,0,1,2,3,8,5)) +
		scale_fill_manual(name="Symbiosis", labels = c("C1", "D1", "Sediments","SS"),values = wes_palette("Zissou1", 4, type = "continuous"))+
		scale_color_manual(name="Symbiosis", labels = c("C1", "D1", "Sediments","SS"),values = wes_palette("Zissou1", 4, type = "continuous"))+
		guides(fill = guide_legend(override.aes = list(shape = 21, size = 3))) +
		guides(shape = guide_legend(override.aes = list(fill = "black", size = 3))) +
		labs(x = paste0("LD1 [", varexpl[1],"%]"), y = paste0("LD2 [", varexpl[2],"%]")) +
		theme_bw()

##### Load in juvenile GE results

ll=load("DAPC_juves_loadings_hotGE.Rdata")
ll2=load("DAPC_juvenile_SNP_loadings.Rdata")
#juvenile_load_GE$sample<-sub("_","-", larval_load_GE$sample)

dapc_loadings=merge(juves_load_GE,juves_load_SNPs,by="sample")

traits=merge(dapc_loadings,colData,by="sample")

Hot=traits[complete.cases(traits),]
#rownames(Hot)=make.names(paste0(Hot$groupall), unique=T)
rownames(Hot)=make.names(paste0(Hot$group1), unique=T) #changed l to 1 
Hot <- Hot[order(Hot$X58daysHeatSurv),] #order by survival in heat (low survival to high survival)
Hot$perc =(Hot$X58daysHeatSurv)/100

##################3
surv=Hot[,c(2:5,14)]
survival=surv$X58daysHeatSurv

fitpc=rda(surv ~ survival,scale=T)
axes2plot=c(2,3)
plot(fitpc,choices=axes2plot)
str(fitpc)
pc1expl=100*round(summary(eigenvals(fitpc))[2, axes2plot[1]],2)
pc2expl=100*round(summary(eigenvals(fitpc))[2, axes2plot[2]],2)
heat.colors = wes_palette("Zissou1", 10, type = "continuous")


plot(fitpc$CA$u,pch=16, cex=1,col= heat.colors,xlim=c(-0.7,0.6),ylim=c(-0.6,0.5),main="Juvenile Survival at heat at 58 days")
arrows(rep(0,4),rep(0,4), fitpc$CA$v[,1]/2, fitpc $CA$v[,2]/1.5,length=0.05)
text(fitpc$CA$v[,1]/1.3, fitpc$CA$v[,2]/1.3+c(0,0,0,0),labels=row.names(fitpc$CA$v) ,cex=0.8,col=c(rep(1,4),rep(4,3)))

summary(lm(surv$GE_dapc1~fitpc$CA$u[,1])) #p=0.000581
summary(lm(surv$SNP_dapc1~fitpc$CA$u[,1])) #p=4.18e-07


#######################
Hot_Hot<-Hot%>% filter(Treatment.1 == "Hot") #dim 81 x 18 good
gm_Hot=glm(perc ~ GE_dapc1+ SNP_dapc1,data=Hot_Hot,family=binomial(link = "logit"))
summary(gm_Hot)
plot(gm_Hot)
gm_Hot_v2=glmer(perc ~ GE_dapc1+ SNP_dapc1+(1|Treatment.1), data=Hot,family=binomial(link = "logit"))
summary(gm_Hot_v2)

gm_Hot_v2_onlyHot=glmer(perc ~ GE_dapc1+ SNP_dapc1+(1|Cross), data=Hot_Hot,family=binomial(link = "logit"))
X.var<- var(as.vector(lme4::fixef(gm_Hot_v2) %*% t(gm_Hot_v2@pp$X))) #0.1499041
#lme4::fixef(gm_Hot_v2) 
#(Intercept)    GE_dapc1   SNP_dapc1 
#2.17810577  0.05818529  0.10739495 

#GE_dapc1 (0.05818529/0.1499041) *100 =38.81501% of X.var which is mostly fixed effects as R2.marginal shows
#SNP_dapc1(0.10739495/0.1499041)*100 =71.64244%

## Extract the variance components for the random effects (not including the residuals)
Z.var <- sum(
  sapply(
    VarCorr(gm_Hot_v2)[!sapply(unique(unlist(strsplit(names(ranef(gm_Hot_v2)),":|/"))), function(l)
      length(unique(gm_Hot_v2@frame[,l])) == nrow(gm_Hot_v2@frame))],
    function(Sigma) {
      X <- model.matrix(gm_Hot_v2)
      Z <- X[,rownames(Sigma)]
      sum(diag(Z %*% Sigma %*% t(Z)))/nrow(X) } ) )
## Extract the variance componts for the residuals
R.var <- attr(lme4::VarCorr(gm_Hot_v2), "sc")^2
#1

## The marginal R2 (proportion of variance due to fixed effects)
(R2.marginal <- X.var/(X.var+Z.var+R.var))
# 0.1303623 x 100 % of variance due to fixed effects

## The proportion of variance due to random effects
(R2.random <- Z.var/(X.var+Z.var+R.var))
# 3.212154e-18 x 100% % (treatment and and obs)

## The proportion of variance due to residuals
VarCorr(gm_Hot_v2) 
#only treatment temp
#1.9219e-09 stdev ==100% of resid bc only one factor here

#So treatmentis 100%% of total 13%% 

(R2.resid <- R.var/(X.var+Z.var+R.var))
# 0.869 x 100 just due to residuals

## The conditional R2 (proportion of variance due to fixed and random effects)
(R2.conditional <- (X.var+Z.var)/(X.var+Z.var+R.var)) #0.130% = 13% this is total variabiliy expalined by all

#13% is total from R2.conditional

########################
########################
#now only PCA run on hot samples
Juveniles_onlyHotDAPC_rescale<-Hot %>% 
  mutate(zscore.GE = (GE_dapc1 - mean(GE_dapc1))/sd(GE_dapc1))%>% 
  mutate(zscore.SNP = (SNP_dapc1 - mean(SNP_dapc1))/sd(SNP_dapc1))

gm_Hot_v2_onlyHotPCA <- glmer(perc ~ zscore.GE+ zscore.SNP+(1|Cross), data=Juveniles_onlyHotDAPC_rescale,family=binomial(link = "logit"))
plot(gm_Hot_v2_onlyHotPCA)

X.var.onlyHotPCA<- var(as.vector(lme4::fixef(gm_Hot_v2_onlyHotPCA) %*% t(gm_Hot_v2_onlyHotPCA@pp$X))) #0.0411134 total variance from fixed ok
#lme4::fixef(gm_Hot_v2_onlyHotPCA)  #these are the mean estimates, not Stdev or variance
#Fixed Effects:
#(Intercept)    zscore.GE   zscore.SNP  
#3.27505      0.01823      0.20341 
#GE_dapc1 (0.01823/0.0411134) *100 =44.34077% of X.var which is mostly fixed effects as R2.marginal shows
#SNP_dapc1(0.20341/0.0411134)*100 =3.09018e-14%

# Variance-Covariance Matrix of fixed effects: 
vc_fixed.j <- as.matrix(vcov(gm_Hot_v2_onlyHotPCA))
# Variance of fixed effects: 
var_fixed.j <- diag(vc_fixed.j); var_fixed.j
#(Intercept)   zscore.GE  zscore.SNP 
#1.075270    1.123996    1.381905
#can also run without intercept at -1
#(1.075270+1.123996+1.381905) =3.581171
#zscore.GE: (1.123996/3.581171)*100 =31.38627%
#zscore.SNP: (1.381905/3.581171)*100 =38.58808%


## Extract the variance components for the random effects (not including the residuals)
Z.var.onlyHotPCA <- sum(
  sapply(
    VarCorr(gm_Hot_v2_onlyHotPCA)[!sapply(unique(unlist(strsplit(names(ranef(gm_Hot_v2_onlyHotPCA)),":|/"))), function(l)
      length(unique(gm_Hot_v2_onlyHotPCA@frame[,l])) == nrow(gm_Hot_v2_onlyHotPCA@frame))],
    function(Sigma) {
      X <- model.matrix(gm_Hot_v2_onlyHotPCA)
      Z <- X[,rownames(Sigma)]
      sum(diag(Z %*% Sigma %*% t(Z)))/nrow(X) } ) )
#2.603461e-09

## Extract the variance componts for the residuals
R.var.onlyHotPCA <- attr(lme4::VarCorr(gm_Hot_v2_onlyHotPCA), "sc")^2
#1

## The marginal R2 (proportion of variance due to fixed effects)
(R2.marginal.onlyHotPCA <- X.var.onlyHotPCA/(X.var.onlyHotPCA+Z.var.onlyHotPCA+R.var.onlyHotPCA))
# 0.03948984 x 100 % of variance due to fixed effects
#0.03948984 * 100 =3.94

## The proportion of variance due to random effects
(R2.random.onlyHotPCA <- Z.var.onlyHotPCA/(X.var.onlyHotPCA+Z.var.onlyHotPCA+R.var.onlyHotPCA))
# 2.500651e-09 * 100% % (treatment and and obs)

## The proportion of variance due to residuals
VarCorr(gm_Hot_v2.onlyHotPCA) #this is only random effects bit!
#only treatment temp
#1.9219e-09 stdev ==100% of resid bc only one factor here

#So treatmentis 100%% of total 13%% 

(R2.resid.onlyHotPCA <- R.var.onlyHotPCA/(X.var.onlyHotPCA+Z.var.onlyHotPCA+R.var.onlyHotPCA))
# 0.9605102 * 100 just due to residuals

## The conditional R2 (proportion of variance due to fixed and random effects)
(R2.conditional.onlyHotPCA <- (X.var.onlyHotPCA+Z.var.onlyHotPCA)/(X.var.onlyHotPCA+Z.var.onlyHotPCA+R.var.onlyHotPCA)) #0.130% = 13% this is total variabiliy expalined by all
#0.03948984 * 100 = 3.948984%

#F + R =Conditional
#F + R= 3.948984%
#Fixed= %
#3.94% (conditional) -3.94% (marginal, only fixed) = 0% only random
#3.94-3.94
#0 % due to random effect


######################################
#######################################

#create df
juveniles.perc<-c(38.81501, 71.64, 13.0)
juv.factors<-c("GE_dapc1", "SNP_dapc1", "R2cond")
juv.fact2<-c("Juv", "Juv", "Juv")
juv.relativ.cont <- data.frame(juveniles.perc, juv.factors, juv.fact2, stringsAsFactors=FALSE)

ggplot(juv.relativ.cont, aes(x=juv.factors, y=juveniles.perc))+geom_bar(stat="identity", aes(fill=juv.fact2), width = 0.5)

#stacked barplot, two bars with total r2. conditional per life stage written above each bar
#73% TOTAL% explained in larvae (fixed is 44.9%), 
#3.9% TOTAL% explained in juveniles (fixed is 3.9%)

#larvae: #zscore.GE: (0.8542799/2.166244)*100 =39.436%. Then 39.4% of 44.9% = 17.6 %
#larvae: #zscore.SNP:(0.9470256/2.166244)*100 =43.7174% Then 43.7% of 44.9% = 19.6 %

#juvs: zscore.GE: (1.123996/3.581171)*100 =31.38627%
#juvs: zscore.SNP: (1.381905/3.581171)*100 =38.58808%


F.R.num<-c(44.9, 28.1, 27, 3.9, 0, 96)
F.R.factors.larv.juv<-c("Fixed", "Random", "Residuals", "Fixed", "Random", "Residuals")
fact3<-c("Larvae", "Larvae","Larvae", "Juveniles", "Juveniles","Juveniles")
relativ.cont_totalFix <- data.frame(F.R.num, F.R.factors.larv.juv, fact3, stringsAsFactors=FALSE)
relativ.cont_totalFix$fact3 <- factor(relativ.cont_totalFix$fact3,levels = c("Larvae", "Juveniles"))

#Perc<-c(67.81, 24.0, 38.81501, 71.64)
Perc.onlyhot<-c(39.4, 43.7, 31.4, 38.5) #only fixed

factors.larv.juv<-c("GE_dapc1", "SNP_dapc1","GE_dapc1", "SNP_dapc1")
fact2<-c("Larvae", "Larvae", "Juveniles", "Juveniles")

relativ.cont <- data.frame(Perc.onlyhot, factors.larv.juv, fact2, stringsAsFactors=FALSE)
relativ.cont$fact2 <- factor(relativ.cont$fact2,levels = c("Larvae", "Juveniles"))


Only.fixed<-ggplot(relativ.cont, aes(x=fact2, y=Perc.onlyhot))+geom_bar(stat="identity", aes(fill=factors.larv.juv), width = 0.5)+
  scale_fill_manual(name="Treatment", values = c("#3B9AB2", "#F21A00","#EBCC2A"))+
  scale_color_manual(name="Treatment", values = c("#3B9AB2", "#F21A00","#EBCC2A"))+
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3))) +
  guides(shape = guide_legend(override.aes = list(fill = "black", size = 3))) +
  labs(x = paste0("Life-stage"), y = paste0("Percent contribution to fixed effects [","%]")) +
  theme_bw()

All.effects<-ggplot(relativ.cont_totalFix, aes(x=fact3, y=F.R.num))+geom_bar(stat="identity", aes(fill=F.R.factors.larv.juv), width = 0.5)+
  scale_fill_manual(name="Treatment", values = c("#3B9AB2", "#F21A00","#EBCC2A"))+
  scale_color_manual(name="Treatment", values = c("#3B9AB2", "#F21A00","#EBCC2A"))+
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3))) +
  guides(shape = guide_legend(override.aes = list(fill = "black", size = 3))) +
  labs(x = paste0("Life-stage"), y = paste0("Percent contribution to surival [","%]")) +
  theme_bw()

library(gridExtra)
grid.arrange(All.effects,Only.fixed, nrow=1)



####### extract SNPs differentiating CU from other pops 
aten<-read.vcfR('./SNPs/Aten_larval_SNPs_filtered_1.recode.vcf')

###need to insert all the commands for this from old script

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



