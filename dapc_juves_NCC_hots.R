setwd("/Users/mes0192/Dropbox/KQ_GBR_spawnseq/analysis/DESeq2_juves")
setwd("/Users/mariestrader/Dropbox/KQ_GBR_spawnseq/analysis/DESeq2_juves")
setwd("/Users/mariestrader/Dropbox/StuffFromLaptop/Quigley")

library('DESeq2')
library('adegenet')
library('dplyr')
library(ggplot2)
library(wesanderson)
library('stringr')

#library('arrayQualityMetrics')
#library('vegan')
#library('rgl')
#library('ape')
#library('pheatmap')
#library('VennDiagram')
#
#library(tidyr)
#library(readr)
#library(Rmisc)
library(ggridges)
library(cowplot)
#library(wesanderson)
#library('variancePartition')

ll=load("AGF_juvenile_2018_WaldCrossTreatZoox.Rdata")

colData$group1<-gsub("SED","SD", colData$group1)
colData$group3=factor(paste0(colData$group1,colData$Treatment.1))

hotjuves=colData[colData$Treatment.1=="Hot",]

all(rownames(hotjuves) %in% colnames(vstCZT.df)) #TRUE
rldhot <- vstCZT.df[, rownames(hotjuves)]
all(rownames(hotjuves) == colnames(rldhot))

dapc <- dapc(t(rldhot), hotjuves$group3, n.da=2, n.pca=7)
#scatter(dapc,bg="white",scree.da=TRUE,scree.pca=TRUE,legend=TRUE,solid=.4)
varexpl <- round((dapc$eig/sum(dapc$eig))[1:2] * 100, 1) #57.9 24.9


juves_load_GE=as.data.frame(dapc$ind.coord)
juves_load_GE$sample<-row.names(juves_load_GE)
colnames(juves_load_GE)=c("GE_dapc1","GE_dapc2","sample")
save(juves_load_GE,file="DAPC_juves_loadings_hotGE.Rdata")

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

############## cross and symbiosis
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

ridges= ggplot(dapc2, aes(x = LD1, y = factor(str_sub(grp, 1,5)), fill = factor(str_sub(grp, 6,7)), height = ..density..)) +
		geom_density_ridges(scale = 1, stat = "density") +
		scale_y_discrete(expand = c(0.01, 0)) +
		scale_x_continuous(expand = c(0.01, 0)) +
		scale_fill_manual(name="Symbiosis", values = c("#EBCC2A","purple","#3B9AB2","#F21A00")) +
		theme_ridges() 

fig <- plot_grid(ridges,dapc.fig, ncol = 1, align = "v", axis="b",rel_heights = c(1, 1.5))


##### just symbiosis at ambient, remove effect of cross
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
ll=load("AGF_juvenile_2018_WaldCrossTreatZoox.Rdata")

colData$dam=str_sub(colData$Cross, 1,2)
colData$sire = str_sub(colData$Cross,4,5)
colData$group1<-gsub("SED","SD", colData$group1)
colData$groupall=as.factor(paste0(colData$group1, colData$Treatment.1))

dim(vstCZT.df)
colData$group2<-gsub("SED","SD", colData$group2)

vsd2 <- limma::removeBatchEffect(vstCZT.df, colData$Cross)

dapc2 <- dapc(t(vsd2), colData$groupall, n.da=4, n.pca=32)
temp <- optim.a.score(dapc2, n.sim = 5) 
#for the vst, they suggest retaining 12 PCs which explains 0.439 percent of the conserved variance
#for the rld removed batch effects, they suggest retaining 7 PCs which explains 0.439 percent of the conserved variance

dapc <- dapc(t(vsd2), colData$groupall, n.da=2, n.pca=28)
#scatter(dapc,bg="white",scree.da=TRUE,scree.pca=TRUE,legend=TRUE,solid=.4)
varexpl <- round((dapc$eig/sum(dapc$eig))[1:2] * 100, 1) 

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

juves_load_GE=as.data.frame(dapc$ind.coord)
juves_load_GE$sample<-row.names(juves_load_GE)
colnames(juves_load_GE)=c("GE_dapc1","GE_dapc2","sample")
save(juves_load_GE,file="DAPC_juvenile_loadings.Rdata")

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



