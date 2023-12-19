rm(list=ls())
options(stringsAsFactors = FALSE)

library(stringr)
library(plyr)
library(rstan)
library(reshape2)
library(cowplot)
require(dplyr)
require(plotrix)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/Treetraits") 
}  else{
  setwd("/home/deirdre/Treetraits") # for midge
}

spInfo <- read.csv("input/species_ring.csv")
pheno <- read.csv("input/phenoDataWChill.csv")
trtPheno <- read.csv("input/trtPhenoDummy.csv")

load("output/htContLat.Rdata")
sumerht <- summary(mdlHt)$summary
postHt <- rstan::extract(mdlHt)

load("output/lmaContLat.Rdata")
postLMA <- rstan::extract(mdlLMA)
#sumerlma <- summary(mdlLMA)$summary

load("output/dbhContLat.Rdata")
postDBH <- rstan::extract(mdlDBH)

load("output/ssdContLat.Rdata")
postSSD <- rstan::extract(mdlSSD)

load("output/cnContLat.Rdata")
postCN <- rstan::extract(mdlCN)


# load("output/heightDummyIntGrandZ25.Rdata")
# sumerht <- summary(mdlHt)$summary
# postHt <- rstan::extract(mdlHt)

bTrtChillHt <- data.frame(postHt$betaTraitxChill)
bTrtForceHt <- data.frame(postHt$betaTraitxForce)
bTrtPhotoHt <- data.frame(postHt$betaTraitxPhoto)

HtCue <- cbind(bTrtChillHt,bTrtForceHt,bTrtPhotoHt)
names(HtCue) <- c("Chill", "Force","Photo")

HtCueL <- reshape2::melt(HtCue)
HtCueL$trait <- "Height"

# load("output/lmaDummyIntGrandZ25.Rdata")
# postLMA <- rstan::extract(mdlLMA)

bTrtChillLMA <- data.frame(postLMA$betaTraitxChill)
bTrtForceLMA <- data.frame(postLMA$betaTraitxForce)
bTrtPhotoLMA <- data.frame(postLMA$betaTraitxPhoto)

LMACue <- cbind(bTrtChillLMA,bTrtForceLMA,bTrtPhotoLMA)
names(LMACue) <- c("Chill", "Force","Photo")

LMACueL <- reshape2::melt(LMACue)
LMACueL$trait <- "LMA"

# load("output/dbhDummyIntGrandZ25.Rdata")
# mdlDBH <- mdlDBH
# postDBH <- rstan::extract(mdlDBH)

bTrtChillDBH <- data.frame(postDBH$betaTraitxChill)
bTrtForceDBH <- data.frame(postDBH$betaTraitxForce)
bTrtPhotoDBH <- data.frame(postDBH$betaTraitxPhoto)

DBHCue <- cbind(bTrtChillDBH,bTrtForceDBH,bTrtPhotoDBH)
names(DBHCue) <- c("Chill", "Force","Photo")

DBHCueL <- reshape2::melt(DBHCue)
DBHCueL$trait <- "DBH"

# load("output/CNDummyIntGrandZ25.Rdata")
# mdlCN <- mdl
# postCN <- rstan::extract(mdlCN)

bTrtChillCN <- data.frame(postCN$betaTraitxChill)
bTrtForceCN <- data.frame(postCN$betaTraitxForce)
bTrtPhotoCN <- data.frame(postCN$betaTraitxPhoto)

CNCue <- cbind(bTrtChillCN,bTrtForceCN,bTrtPhotoCN)
names(CNCue) <- c("Chill", "Force","Photo")

CNCueL <- reshape2::melt(CNCue)
CNCueL$trait <- "C:N"

# load("output/ssdDummyIntGrandZ25.Rdata")
# postSSD <- rstan::extract(mdlSSD)

bTrtChillSSD <- data.frame(postSSD$betaTraitxChill)
bTrtForceSSD <- data.frame(postSSD$betaTraitxForce)
bTrtPhotoSSD <- data.frame(postSSD$betaTraitxPhoto)

SSDCue <- cbind(bTrtChillSSD,bTrtForceSSD,bTrtPhotoSSD)
names(SSDCue) <- c("Chill", "Force","Photo")

SSDCueL <- reshape2::melt(SSDCue)
SSDCueL$trait <- "SSD"

trtCue <- rbind(HtCueL, LMACueL, DBHCueL, CNCueL, SSDCueL)

pdf("figures/betaTraitxCue.pdf", width = 10, height = 5)
ggplot() +
  stat_eye(data = trtCue, aes(x = trait, y = value, fill = variable, group = variable),  cex = 0.75, position = position_dodge(0.9)) + 
#  facet_wrap(~trait, ncol=1, scales = "free") + 
  theme_classic() + #ylim(-20,20) +
 # theme(legend.position = "none") +
  labs( x = "Trait", y = "Trait by cue response", main = NA)+
  scale_fill_manual(values = c("cyan4", "goldenrod","maroon")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.title = element_blank(), legend.key=element_rect(fill="white")) # removed grey boxes around legends
dev.off()
