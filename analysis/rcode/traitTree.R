# aim of this code is to visualize the phylogeny of spp in bc trait study on a phylogeny
# trying to visualize both height and cue responses on one tree

# Code to add dots to the end of the phylogeny was from chatgpt

rm(list=ls())
options(stringsAsFactors = FALSE)

library(rstan)
library(shinystan)
#library(reshape2)
#library(bayesplot)
library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)
library(phytools)
#library(ggpubr)
library(lattice)
require(cowplot)
require(caper)

 
setwd("~/Documents/github/Treetraits") 

spInfo <- read.csv("analysis/input/species_ring.csv")
spInfo <- spInfo[, c("species.name","species","type","transect","ring.type")]
pheno <- read.csv("analysis/input/phenoDataWChill.csv")
trtPheno <- read.csv("analysis/input/trtPhenoDummy.csv")

htDat <- trtPheno[complete.cases(trtPheno$ht),]
height <- aggregate(htDat["ht"], htDat[c("species")], FUN = mean)

load("analysis/output/htContLatHundowLat_treeShrubDummyInt.Rdata")
sumHt <- summary(mdlhtTS)$summary
postHt <- rstan::extract(mdlhtTS)

b_photoSpHt <- sumHt[grep("alphaPhotoSp\\[", rownames(sumHt)), 1]
b_forceSpHt <- sumHt[grep("alphaForceSp\\[", rownames(sumHt)), 1]
b_chillSpHt <- sumHt[grep("alphaChillSp\\[", rownames(sumHt)), 1]

## Need a nicer phylogy - colour tips based on intercept values:
tree <- read.tree("..//pheno_bc/input/SBphylo_phenobc.tre")

head(tree$tip.label)
length(tree$tip.label) #47
tree$tip.label[tree$tip.label=="Cornus_asperifolia"] <- "Cornus_stolonifera"
tree$tip.label[tree$tip.label=="Alnus_alnobetula"] <- "Alnus_viridis"
tree$tip.label[tree$tip.label== "Fagus_grandifolia_var._caroliniana"] <- "Fagus_grandifolia"
tree$tip.label[tree$tip.label== "Spiraea_alba_var._latifolia"] <- "Spiraea_alba"
tree$tip.label[tree$tip.label== "Rhamnus_arguta"] <- "Rhamnus_frangula"

spFact <- sort(spInfo$species.name)

phylo.dat <- spInfo[order(spInfo$species.name),]

dat.int <- cbind(phylo.dat, height,b_chillSpHt, b_forceSpHt, b_photoSpHt)

# spInfo$species[spInfo$species.name== "Acer_glabrum"] <- "acegla*"
# spInfo$species[spInfo$species.name== "Acer_pensylvanicum"] <- "acepen*"
# spInfo$species[spInfo$species.name== "Acer_rubrum"] <- "acerub*"
# spInfo$species[spInfo$species.name== "Acer_saccharum"] <- "acesac*"
# spInfo$species[spInfo$species.name== "Betula_alleghaniensis"] <- "betall*"
# spInfo$species[spInfo$species.name== "Betula_lenta"] <- "betlen*"
# spInfo$species[spInfo$species.name== "Betula_papyrifera"] <- "betpap*"
# spInfo$species[spInfo$species.name== "Fagus_grandifolia"] <- "faggra*"
# spInfo$species[spInfo$species.name== "Fraxinus_nigra"] <- "franig*"
# spInfo$species[spInfo$species.name== "Nyssa_sylvatica"] <- "nyssyl*"
# spInfo$species[spInfo$species.name== "Populus_balsamifera"] <- "popbal*"
# spInfo$species[spInfo$species.name== "Populus_grandidentata"] <- "popgra*"
# spInfo$species[spInfo$species.name== "Populus_tremuloides"] <- "poptre*"
# spInfo$species[spInfo$species.name== "Prunus_pensylvanica"] <- "prupen*"
# spInfo$species[spInfo$species.name== "Quercus_alba"] <- "quealb*"
# spInfo$species[spInfo$species.name== "Quercus_rubra"] <- "querub*"
# spInfo$species[spInfo$species.name== "Quercus_velutina"] <- "quevel*"


#rownames(dat.int) <- spFact

# tree needs to be rooted:
namesphy <- tree$tip.label
tree$root.edge <- 0
root(tree, outgroup = "Acer_glabrum")

is.rooted(tree)
tree$node.label<-NULL

dataPhy = comparative.data(tree, dat.int, names.col = "species.name", na.omit = F,
  vcv = T, warn.dropped = T)

phyloplot = dataPhy$phy
x = dataPhy$data$b_chillSpHt
y = dataPhy$data$ht
z = dataPhy$data$b_forceSpHt
w = dataPhy$data$b_photoSpHt
names(x)=dataPhy$phy$tip.label
names(y)=dataPhy$phy$tip.label
names(z)=dataPhy$phy$tip.label
names(w)=dataPhy$phy$tip.label

slope <- contMap(tree, y, plot = T)
slopeCol <- setMap(slope, colors=c("cyan","purple4","red"))
h<-max(nodeHeights(slopeCol$tree))

pdf("analysis/figures/phyloIntColor_height.pdf", height = 9, width = 7)
plot(slopeCol,legend = F, lwd=3, ylim=c(1-0.09*(Ntip(slopeCol$tree)),Ntip(slopeCol$tree)), offset = 1.5)
add.color.bar(60, slopeCol$cols, title = "Height (m)", subtitle="", digits = 1,lims = c(0,20),  prompt = F,x=0.2*h, y = -2)

pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)

n <- Ntip(tree)

# Add coloured dots offset from the tree tips
points(
  x = pp$xx[1:n] + 5,
  y = pp$yy[1:n],
  pch = 16,
  col = tip_cols,
  cex = 1.5
)

legend(
  x = 150,
  y = 0.05,
 # legend = names(group_cols),
  legend = c("Tree","Shrub"),
  col = group_cols,
  pch = 16,
  pt.cex = 1.5,
  bty = "n",
  horiz = TRUE
)
dev.off()


slope <- contMap(tree, x, plot = T)
slopeCol <- setMap(slope, colors=c("cyan","purple4","red"))
h<-max(nodeHeights(slopeCol$tree))

pdf("analysis/figures/phyloIntColor_Chilling.pdf", height = 9, width = 7)
plot(slopeCol,legend = F, lwd=3, ylim=c(1-0.09*(Ntip(slopeCol$tree)),Ntip(slopeCol$tree)))
add.color.bar(60, slopeCol$cols, subtitle = "Chill response (z-scored chill portions per day)", title="", digits = 1,lims = c(-30,5),  prompt = F,x=0.2*h, y = -2)
dev.off()

slope <- contMap(tree, z, plot = T)
slopeCol <- setMap(slope, colors=c("cyan","purple4","red"))
h<-max(nodeHeights(slopeCol$tree))

pdf("analysis/figures/phyloIntColor_Forcing.pdf", height = 9, width = 7)
plot(slopeCol,legend = F, lwd=3, ylim=c(1-0.09*(Ntip(slopeCol$tree)),Ntip(slopeCol$tree)))
add.color.bar(60, slopeCol$cols, subtitle = "Forcing response (z-scored forcing per day)", title="", digits = 1,lims = c(-30,5),  prompt = F,x=0.2*h, y = -2)
dev.off()


####################################################
# adding circles to the tips to show when tree vs shrub

group_cols <- c(
  "tree" = "cyan4",
  "shrub" = "maroon"
)

tip_group <- spInfo$type[
  match(tree$tip.label, spInfo$species.name)
]

tip_cols <- group_cols[as.character(tip_group)]

plotTree(tree, 
         fsize = 0.8, 
         ftype = "i",
         offset = 1)

points(
  x = pp$xx[1:n] + 0.5,
  y = pp$yy[1:n],
  pch = 16,
  col = tip_cols,
  cex = 1.5
)

ape::tiplabels(
  pch = 16,
  col = tip_cols,
  cex = 1.5
)

legend(
  "bottomright",
  legend = names(group_cols),
  col = group_cols,
  pch = 16,
  pt.cex = 1.5,
  bty = "n",
  horiz = TRUE
)


# plot(tree, tip.color = dat.int$type, cex = c(1, 1, 1.5))
##### Try to make the tip colour height and tip shape the resp cue
#install.packages("ggimage")
# require(ggtree)
# require(ggimage)
# p <- ggtree(tree) 
# p2 <- p + geom_tiplab(offset = .6, hjust = .5) +
#   geom_tippoint(aes(shape = trophic_habit, color = trophic_habit, 
#     size = mass_in_kg)) + 
#   theme(legend.position = "right") + 
#   scale_size_continuous(range = c(3, 10))
# 
# p2 %<+% df_inode_data + 
#   geom_label(aes(label = vernacularName.y, fill = posterior)) + 
#   scale_fill_gradientn(colors = RColorBrewer::brewer.pal(3, "YlGnBu"))
# 
# 
# ggplot(tree) + geom_tree() + theme_tree()
# 
# require(trait.plot)
# 
# dotTree(tree, dat.int, col = palette("Greens"))
# 
# dat.intC <- dat.int[,1:2]
# phylo.heatmap(tree, dat.intC, standarize =T)
# 
# 
# ##########################################
# obj<-contMap(tree,x,plot=FALSE)
# objy<-contMap(tree,y,plot=FALSE)
# obj<-setMap(obj,invert=TRUE)
# 
# 
# pdf("figures/chillHeightTree.pdf", width =8, height = 6)
# par(mfrow=c(1,2))
# plot(obj,lwd=7,ftype="off",legend=F, ylim=c(1-0.09*(Ntip(slopeCol$tree)),Ntip(slopeCol$tree)))
# add.color.bar(60, obj$cols, title = "Chilling response", subtitle="", digits = 1,lims = c(10,50),  prompt = F,x=0.2*h, y = -2)
# plot(objy,lwd=7,direction="leftwards",ftype="off",legend=F, ylim=c(1-0.09*(Ntip(slopeCol$tree)),Ntip(slopeCol$tree)))
# 
# add.color.bar(60, objy$cols, title = "Height", subtitle="", digits = 1,lims = c(10,50),  prompt = F,x=0.2*h, y = -2)
# dev.off()
# 
# # Forcing
# obj<-contMap(tree,z,plot=FALSE)
# objy<-contMap(tree,y,plot=FALSE)
# obj<-setMap(obj,invert=TRUE)
# 
# 
# pdf("figures/forceHeightTree.pdf", width =8, height = 6)
# par(mfrow=c(1,2))
# plot(obj,lwd=7,ftype="off",legend=F, ylim=c(1-0.09*(Ntip(slopeCol$tree)),Ntip(slopeCol$tree)))
# add.color.bar(60, obj$cols, title = "Forcing response", subtitle="", digits = 1,lims = c(10,50),  prompt = F,x=0.2*h, y = -2)
# plot(objy,lwd=7,direction="leftwards",ftype="off",legend=F, ylim=c(1-0.09*(Ntip(slopeCol$tree)),Ntip(slopeCol$tree)))
# 
# add.color.bar(60, objy$cols, title = "Height", subtitle="", digits = 1,lims = c(10,50),  prompt = F,x=0.2*h, y = -2)
# dev.off()
# 
# # Photoperiod
# obj<-contMap(tree,w,plot=FALSE)
# objy<-contMap(tree,y,plot=FALSE)
# obj<-setMap(obj,invert=TRUE)
# 
# 
# pdf("figures/photoperiodHeightTree.pdf", width =8, height = 6)
# par(mfrow=c(1,2))
# plot(obj,lwd=7,ftype="off",legend=F, ylim=c(1-0.09*(Ntip(slopeCol$tree)),Ntip(slopeCol$tree)))
# add.color.bar(60, obj$cols, title = "Photoperiod response", subtitle="", digits = 1,lims = c(10,50),  prompt = F,x=0.2*h, y = -2)
# plot(objy,lwd=7,direction="leftwards",ftype="off",legend=F, ylim=c(1-0.09*(Ntip(slopeCol$tree)),Ntip(slopeCol$tree)))
# 
# add.color.bar(60, objy$cols, title = "Height", subtitle="", digits = 1,lims = c(10,50),  prompt = F,x=0.2*h, y = -2)
# dev.off()
# 
# nodelabels(text = tree$node.label,
#   frame = "n", cex=0.8, col= "blue")
# 
# # Could use node number as a proxy of evolutionary relatedness, bottom of the tree (populs = small node numbers) top = high
# 
# nodes <- c(73,73, 73,73, 60,60, 68,69,63,63,63,61,56,78,75,79,87,85,85,90,89,92,52,52,52,66,57,58,58,71,89,89,65,75,71,81,68,70,70,70,84,91,91,82,82,82,92)
# 
# dat.int$node <- nodes
# dat.int <- dat.int[,c("species.name","node","ht","b_chillSpHt","b_forceSpHt","b_photoSpHt")]
# 
# htChill <- ggplot(dat.int, aes(x = ht, y = b_chillSpHt, col = node))+
#   geom_point(size =2) +
#   labs( x = "Height", y = "Chilling response", main = NA)+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     panel.background = element_blank(), axis.line = element_line(colour = "black"),
#     legend.key=element_rect(fill="white"), legend.title = element_blank(),legend.position = "none") +
#   scale_color_viridis() 
# 
# htForce <- ggplot(dat.int, aes(x = ht, y = b_forceSpHt, col = node))+
#   geom_point(size =2) +
#   labs( x = "Height", y = "Forcing response", main = NA)+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     panel.background = element_blank(), axis.line = element_line(colour = "black"),
#     legend.key=element_rect(fill="white"), legend.title = element_blank(),legend.position = "none") +
#   scale_color_viridis() 
# 
# 
# htPhoto <- ggplot(dat.int, aes(x = ht, y = b_photoSpHt, col = node))+
#   geom_point(size =2) +
#   labs( x = "Height", y = "Photoperiod response", main = NA)+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     panel.background = element_blank(), axis.line = element_line(colour = "black"),
#     legend.key=element_rect(fill="white"), legend.title = element_blank(),legend.position = "none") +
#   scale_color_viridis() 
# 
# pdf("figures/heightCuePhyloRelatedness.pdf", height =4, width = 12)
# plot_grid( htChill, htForce, htPhoto, ncol = 3, nrow =1,align = "v")
# dev.off()
