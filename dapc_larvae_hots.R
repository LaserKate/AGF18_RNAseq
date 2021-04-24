
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
#library(tidyr)
library(readr)
#library(Rmisc)
library(ggplot2)
library(ggridges)
library(cowplot)
library(wesanderson)

ll=load("AGF_larval_2018_WaldGroupContrasts.Rdata")

dim(rldG.df)
#order by survival
colData$cross_f = factor(colData$cross, levels=c('LSxSB','LSxLS','DRxDR','SBxSB','LSxCU','DRxSB','BKxBK','CUxCU','CUxSB','CUxBK','SBxBK'))

vsd2 <- limma::removeBatchEffect(rldG.df, colData$cross)

##Calling PCs
#pcp=prcomp(t(rldG.df[degs,]), retx=TRUE, center=TRUE, scale=TRUE)
#scores=pcp$x
#screeplot(pcp,bstick=T)

#How many PCs should be kept?
dapc2 <- dapc(t(rldG.df), colData$group, n.da=4, n.pca=35)
temp <- optim.a.score(dapc2, n.sim = 5) 


#subset for only hot larvae
hotlarvae=colData[colData$treatment=="Hot",]

all(rownames(hotlarvae) %in% colnames(rldG.df)) #TRUE
rldhot <- rldG.df[, rownames(hotlarvae)]
all(rownames(hotlarvae) == colnames(rldhot))


dapc <- dapc(t(rldhot), hotlarvae$cross, n.da=2, n.pca=7)
scatter(dapc,bg="white",scree.da=TRUE,scree.pca=TRUE,legend=TRUE,solid=.4)
varexpl <- round((dapc$eig/sum(dapc$eig))[1:2] * 100, 1) #38.4 26.4

#export loadings for genotype comparison 
larval_load_GE=as.data.frame(dapc$ind.coord)
larval_load_GE$sample<-row.names(larval_load_GE)
colnames(larval_load_GE)=c("GE_dapc1","GE_dapc2","sample")
save(larval_load_GE,file="DAPC_larval_loadings_hotGE.Rdata")


scatter(dapc, bg = "white", legend = TRUE, scree.da = FALSE)

dapc1 <- tibble(sample = rownames(dapc$ind.coord),
               grp = dapc$grp,
               LD1 = dapc$ind.coord[,1],
               LD2 = dapc$ind.coord[,2])
dapc2 <- dapc1 %>%
  group_by(grp) %>%
  summarize(c1 = mean(LD1),
            c2 = mean(LD2)) %>%
  			full_join(dapc1)

dapc.fig <-   ggplot(dapc2, aes(shape = factor(str_sub(grp, 1,5)), color = factor(str_sub(grp, 6,8)),
		fill = factor(str_sub(grp, 6,8)))) +
		geom_segment(mapping = aes(x = LD1, y = LD2, xend = c1, yend = c2), lwd = 0.25, col = "grey") +
		geom_point(aes(x = c1, y = c2), size = 3) +
		geom_point(aes(x = LD1, y = LD2), size = 2, show.legend = FALSE) +
		scale_shape_manual(name = "Cross", values=c(21,22,23,24,25,0,1,2,3,8,5)) +
		scale_fill_manual(name="Treatment", values = c("#3B9AB2", "#F21A00","#EBCC2A"))+
		scale_color_manual(name="Treatment", values = c("#3B9AB2", "#F21A00","#EBCC2A"))+
		guides(fill = guide_legend(override.aes = list(shape = 21, size = 3))) +
		guides(shape = guide_legend(override.aes = list(fill = "black", size = 3))) +
		labs(x = paste0("LD1 [", varexpl[1],"%]"), y = paste0("LD2 [", varexpl[2],"%]")) +
		theme_bw()

ridges= ggplot(dapc2, aes(x = LD1, y = factor(str_sub(grp, 1,5)), fill = factor(str_sub(grp, 6,8)), height = ..density..)) +
		geom_density_ridges(scale = 4, stat = "density") +
		scale_y_discrete(expand = c(0.01, 0)) +
		scale_x_continuous(expand = c(0.01, 0)) +
		scale_fill_manual(name="treatment", values = c("#3B9AB2", "#F21A00","#EBCC2A")) +
		theme_ridges() 

fig <- plot_grid(ridges,dapc.fig, ncol = 1, align = "v", axis="b",rel_heights = c(1, 1))



