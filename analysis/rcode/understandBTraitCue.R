# Started Dec 5, 2023

# why is chilling so variable for SSD and LMA 
# are the estimates being pulled by spp with little data?

rm(list=ls())
options(stringsAsFactors = FALSE)

library(tidyr)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(rstan)
library(bayesplot)# nice posterior check plots 
library(shinystan)
library(stringr)
library(truncnorm)
library(tidybayes)
library(cowplot)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/Treetraits")
}  else{
  setwd("/home/deirdre/Treetraits") # for midge
}

trtPheno <- read.csv("input/trtPhenoDummy.csv")

pheno.t <- read.csv("input/phenoDataWChill.csv")

spInfo <- read.csv("input/species_ring.csv")
sp.name <- sort(spInfo$species)
# #########################################################

# What does the data look like?
ssd <- trtPheno[complete.cases(trtPheno$ssd),]
ssd$count <- 1
ssdNo <- aggregate(ssd["count"], ssd[c("species")], FUN = sum)
ssdSpNo <- paste(ssdNo$species, " n = ",ssdNo$count, sep = "")

lma <- trtPheno[complete.cases(trtPheno$lma),]
lma$count <- 1
lmaNo <- aggregate(lma["count"], lma[c("species")], FUN = sum)
lmaSpNo <- paste(lmaNo$species, " n = ",lmaNo$count, sep = "")

ht <- trtPheno[complete.cases(trtPheno$ht),]
ht$count <- 1
htNo <- aggregate(ht["count"], ht[c("species")], FUN = sum)
htSpNo <- paste(htNo$species, " n = ",htNo$count, sep = "")

## Read in model output
load("output/lmaDummyIntGrandZ25.Rdata")
postLMA<- data.frame(rstan::extract(mdlLMA))

load("output/ssdDummyIntGrandZ25.Rdata")
postSSD<- data.frame(rstan::extract(mdlSSD))

load("output/heightDummyIntGrandZ25.Rdata")
postHt<- data.frame(rstan::extract(mdlHt))

# get posterior for betaTraitChill from each model
bTrtChillLMA <- data.frame(postLMA$betaTraitxChill)
bTrtChillSSD <- data.frame(postSSD$betaTraitxChill)
bTrtChillHt <- data.frame(postHt$betaTraitxChill)

bTChill <- cbind(bTrtChillLMA,bTrtChillSSD,bTrtChillHt)
names(bTChill) <- c("LMA", "SSD","Ht")

bTChillL <- melt(bTChill)

# ggplot() +
#   stat_eye(data = bTChillL, aes(x = variable, y = value, fill = variable, group = variable),  cex = 0.75, position = position_dodge(0.9)) + 
#   #  facet_wrap(~trait, ncol=1, scales = "free") + 
#   theme_classic() + #ylim(-20,20) +
#   # theme(legend.position = "none") +
#   labs( x = "Trait", y = "Trait by cue response", main = NA)+
#   scale_fill_manual(values = c("cyan4", "goldenrod","maroon")) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.title = element_blank(), legend.key=element_rect(fill="white")) # removed grey boxes around legends

## Look at the species level estimates:
aChillSpLMA <- postLMA[,colnames(postLMA) %in% grep( "alphaChillSp", colnames(postLMA), value = TRUE)]
colnames(aChillSpLMA) <- sp.name
aChillSpLMAL <- melt(aChillSpLMA)
colnames(aChillSpLMAL) <- c("species","value")

aChillSpLMAL <- merge(aChillSpLMAL, spInfo, by = "species")
aChillSpLMAL <- merge(aChillSpLMAL, lmaNo, by = "species")
aChillSpLMAL$obs <- paste(aChillSpLMAL$species, " n = ",aChillSpLMAL$count, sep = "")


LMABTC <-ggplot() +
  stat_eye(data = aChillSpLMAL, aes(x = obs, y = value, fill = "cyan4"),  cex = 0.75, position = position_dodge(0.9)) + 
  theme_classic() + #ylim(-20,20) +
  # theme(legend.position = "none") +
  labs( x = "LMA", y = "Trait by Chill response", main = NA)+
  scale_fill_manual(values = c("cyan4")) +
  theme(panel.grid.major = element_blank(), legend.position = "none", panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1),axis.line = element_line(colour = "black"),legend.title = element_blank(), legend.key=element_rect(fill="white")) # removed grey boxes around legends

# SSD output
aChillSpSSD <- postSSD[,colnames(postSSD) %in% grep( "alphaChillSp", colnames(postSSD), value = TRUE)]
colnames(aChillSpSSD) <- sp.name
aChillSpSSDL <- melt(aChillSpSSD)
colnames(aChillSpSSDL) <- c("species","value")

aChillSpSSDL <- merge(aChillSpSSDL, spInfo, by = "species")
aChillSpSSDL <- merge(aChillSpSSDL, ssdNo, by = "species")
aChillSpSSDL$obs <- paste(aChillSpSSDL$species, " n = ",aChillSpSSDL$count, sep = "")


SSDBTC <- ggplot() +
  stat_eye(data = aChillSpSSDL, aes(x = obs, y = value, fill = "maroon"),  cex = 0.75, position = position_dodge(0.9)) + 
  theme_classic() + #ylim(-20,20) +
  # theme(legend.position = "none") +
  labs( x = "SSD", y = "Trait by chill response", main = NA)+
  scale_fill_manual(values = c("maroon")) +
  theme(panel.grid.major = element_blank(), legend.position = "none", panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1),axis.line = element_line(colour = "black"),legend.title = element_blank(), legend.key=element_rect(fill="white")) # removed grey boxes around legends


aChillSpHt <- postHt[,colnames(postHt) %in% grep( "alphaChillSp", colnames(postHt), value = TRUE)]
colnames(aChillSpHt) <- sp.name
aChillSpHtL <- melt(aChillSpHt)
colnames(aChillSpHtL) <- c("species","value")

aChillSpHtL <- merge(aChillSpHtL, spInfo, by = "species")
aChillSpHtL <- merge(aChillSpHtL, htNo, by = "species")
aChillSpHtL$obs <- paste(aChillSpHtL$species, " n = ",aChillSpHtL$count, sep = "")


HtBTC <- ggplot() +
  stat_eye(data = aChillSpHtL, aes(x = obs, y = value, fill = "goldenrod"),  cex = 0.75, position = position_dodge(0.9)) + 
  theme_classic() + #ylim(-20,20) +
  # theme(legend.position = "none") +
  labs( x = "Ht", y = "Trait by chill response", main = NA)+
  scale_fill_manual(values = c("goldenrod")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1),axis.line = element_line(colour = "black"),legend.title = element_blank(), legend.key=element_rect(fill="white")) # removed grey boxes around legends

pdf("figures/alphaChillSp_lma_ssd_ht.pdf", height =12, width = 12)
plot_grid(LMABTC, SSDBTC,HtBTC, nrow = 3, align = "v")#, rel_heights = c(1/4, 1/4, 1.2/3))
dev.off()

############ Compare to Forcing ###############
bTrtForceLMA <- data.frame(postLMA$betaTraitxForce)
bTrtForceSSD <- data.frame(postSSD$betaTraitxForce)
bTrtForceHt <- data.frame(postHt$betaTraitxForce)

bTForce <- cbind(bTrtForceLMA,bTrtForceSSD,bTrtForceHt)
names(bTForce) <- c("LMA", "SSD","Ht")

bTForceL <- melt(bTForce)

# ggplot() +
#   stat_eye(data = bTForceL, aes(x = variable, y = value, fill = variable, group = variable),  cex = 0.75, position = position_dodge(0.9)) + 
#   #  facet_wrap(~trait, ncol=1, scales = "free") + 
#   theme_classic() + #ylim(-20,20) +
#   # theme(legend.position = "none") +
#   labs( x = "Trait", y = "Trait by cue response", main = NA)+
#   scale_fill_manual(values = c("cyan4", "goldenrod","maroon")) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.title = element_blank(), legend.key=element_rect(fill="white")) # removed grey boxes around legends

## Look at the species level estimates:
aForceSpLMA <- postLMA[,colnames(postLMA) %in% grep( "alphaForceSp", colnames(postLMA), value = TRUE)]
colnames(aForceSpLMA) <- sp.name
aForceSpLMAL <- melt(aForceSpLMA)
colnames(aForceSpLMAL) <- c("species","value")

aForceSpLMAL <- merge(aForceSpLMAL, spInfo, by = "species")
aForceSpLMAL <- merge(aForceSpLMAL, lmaNo, by = "species")
aForceSpLMAL$obs <- paste(aForceSpLMAL$species, " n = ",aForceSpLMAL$count, sep = "")


LMABTF <-ggplot() +
  stat_eye(data = aForceSpLMAL, aes(x = obs, y = value, fill = "cyan4"),  cex = 0.75, position = position_dodge(0.9)) + 
  theme_classic() + #ylim(-20,20) +
  # theme(legend.position = "none") +
  labs( x = "LMA", y = "Trait by force response", main = NA)+
  scale_fill_manual(values = c("cyan4")) +
  theme(panel.grid.major = element_blank(), legend.position = "none", panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1),axis.line = element_line(colour = "black"),legend.title = element_blank(), legend.key=element_rect(fill="white")) # removed grey boxes around legends

# SSD output
aForceSpSSD <- postSSD[,colnames(postSSD) %in% grep( "alphaForceSp", colnames(postSSD), value = TRUE)]
colnames(aForceSpSSD) <- sp.name
aForceSpSSDL <- melt(aForceSpSSD)
colnames(aForceSpSSDL) <- c("species","value")

aForceSpSSDL <- merge(aForceSpSSDL, spInfo, by = "species")
aForceSpSSDL <- merge(aForceSpSSDL, ssdNo, by = "species")
aForceSpSSDL$obs <- paste(aForceSpSSDL$species, " n = ",aForceSpSSDL$count, sep = "")


SSDBTF <- ggplot() +
  stat_eye(data = aForceSpSSDL, aes(x = obs, y = value, fill = "maroon"),  cex = 0.75, position = position_dodge(0.9)) + 
  theme_classic() + #ylim(-20,20) +
  # theme(legend.position = "none") +
  labs( x = "ssd", y = "Trait by force response", main = NA)+
  scale_fill_manual(values = c("maroon")) +
  theme(panel.grid.major = element_blank(), legend.position = "none", panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1),axis.line = element_line(colour = "black"),legend.title = element_blank(), legend.key=element_rect(fill="white")) # removed grey boxes around legends


aForceSpHt <- postHt[,colnames(postHt) %in% grep( "alphaForceSp", colnames(postHt), value = TRUE)]
colnames(aForceSpHt) <- sp.name
aForceSpHtL <- melt(aForceSpHt)
colnames(aForceSpHtL) <- c("species","value")

aForceSpHtL <- merge(aForceSpHtL, spInfo, by = "species")
aForceSpHtL <- merge(aForceSpHtL, htNo, by = "species")
aForceSpHtL$obs <- paste(aForceSpHtL$species, " n = ",aForceSpHtL$count, sep = "")


HtBTF <- ggplot() +
  stat_eye(data = aForceSpHtL, aes(x = obs, y = value, fill = "goldenrod"),  cex = 0.75, position = position_dodge(0.9)) + 
  theme_classic() + #ylim(-20,20) +
  # theme(legend.position = "none") +
  labs( x = "ht", y = "Trait by force response", main = NA)+
  scale_fill_manual(values = c("goldenrod")) +
  theme(panel.grid.major = element_blank(), legend.position = "none", panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1),axis.line = element_line(colour = "black"),legend.title = element_blank(), legend.key=element_rect(fill="white")) # removed grey boxes around legends

pdf("figures/alphaForceSp_lma_ssd_ht.pdf", height =12, width = 12)
plot_grid(LMABTC, SSDBTC,HtBTC, nrow = 3, align = "v")#, rel_heights = c(1/4, 1/4, 1.2/3))
dev.off()

## Look at the species trait variation ####
muSpLMA <- postLMA[,colnames(postLMA) %in% grep( "mu_grand_sp", colnames(postLMA), value = TRUE)]
colnames(muSpLMA) <- sp.name
muSpLMAL <- melt(muSpLMA)
colnames(muSpLMAL) <- c("species","value")

muSpLMAL <- merge(muSpLMAL, spInfo, by = "species")
muSpLMAL <- merge(muSpLMAL, lmaNo, by = "species")
muSpLMAL$obs <- paste(muSpLMAL$species, " n = ",muSpLMAL$count, sep = "")


LMAmuSp <-ggplot() +
  stat_eye(data = muSpLMAL, aes(x = obs, y = value, fill = "cyan4"),  cex = 0.75, position = position_dodge(0.9)) + 
  theme_classic() + #ylim(-20,20) +
  # theme(legend.position = "none") +
  labs( x = "LMA", y = "mu_grand_sp", main = NA)+
  scale_fill_manual(values = c("cyan4")) +
  theme(panel.grid.major = element_blank(), legend.position = "none", panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1),axis.line = element_line(colour = "black"),legend.title = element_blank(), legend.key=element_rect(fill="white")) # removed grey boxes around legends

# SSD output
muSpSSD <- postSSD[,colnames(postSSD) %in% grep( "mu_grand_sp", colnames(postSSD), value = TRUE)]
colnames(muSpSSD) <- sp.name
muSpSSDL <- melt(muSpSSD)
colnames(muSpSSDL) <- c("species","value")

muSpSSDL <- merge(muSpSSDL, spInfo, by = "species")
muSpSSDL <- merge(muSpSSDL, ssdNo, by = "species")
muSpSSDL$obs <- paste(muSpSSDL$species, " n = ",muSpSSDL$count, sep = "")


SSDmuSp <- ggplot() +
  stat_eye(data = muSpSSDL, aes(x = obs, y = value, fill = "maroon"),  cex = 0.75, position = position_dodge(0.9)) + 
  theme_classic() + #ylim(-20,20) +
  # theme(legend.position = "none") +
  labs( x = "SSD", y = "mu_grand_sp", main = NA)+
  scale_fill_manual(values = c("maroon")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1), legend.position = "none",axis.line = element_line(colour = "black"),legend.title = element_blank(), legend.key=element_rect(fill="white")) # removed grey boxes around legends


muSpHt <- postHt[,colnames(postHt) %in% grep( "mu_grand_sp", colnames(postHt), value = TRUE)]
colnames(muSpHt) <- sp.name
muSpHtL <- melt(muSpHt)
colnames(muSpHtL) <- c("species","value")

muSpHtL <- merge(muSpHtL, spInfo, by = "species")
muSpHtL <- merge(muSpHtL, htNo, by = "species")
muSpHtL$obs <- paste(muSpHtL$species, " n = ",muSpHtL$count, sep = "")


HtmuSp <- ggplot() +
  stat_eye(data = muSpHtL, aes(x = obs, y = value, fill = "goldenrod"),  cex = 0.75, position = position_dodge(0.9)) + 
  theme_classic() + #ylim(-20,20) +
  # theme(legend.position = "none") +
  labs( x = "Height", y = "mu_grand_sp", main = NA)+
  scale_fill_manual(values = c("goldenrod")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1),axis.line = element_line(colour = "black"),legend.title = element_blank(), legend.key=element_rect(fill="white"), legend.position = "none") # removed grey boxes around legends

pdf("figures/alphmuSp_lma_ssd_ht.pdf", height =12, width = 12)
plot_grid(LMAmuSp, SSDmuSp,HtmuSp, nrow = 3, align = "v")#, rel_heights = c(1/4, 1/4, 1.2/3))
dev.off()


pdf("figures/compareHtLMA.pdf", height =12, width = 20)
plot_grid(LMAmuSp,SSDmuSp, HtmuSp,LMABTF,SSDBTF,HtBTF,LMABTC,SSDBTC,HtBTC, nrow = 3, ncol = 3, align = "v")#, rel_heights = c(1/4, 1/4, 1.2/3))
dev.off()

