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

spInfo <- read.csv("input/species_list.csv")
pheno <- read.csv("input/phenoData.csv")
trtPheno <- read.csv("input/trtData.csv")

load("output/heightDummyIntGrand.Rdata")
sumerht <- summary(mdl)$summary
mdlHt <- mdl
postHt <- rstan::extract(mdlHt)

bTrtChillHt <- data.frame(postHt$betaTraitxChill)
bTrtForceHt <- data.frame(postHt$betaTraitxForce)
bTrtPhotoHt <- data.frame(postHt$betaTraitxPhoto)

HtCue <- cbind(bTrtChillHt,bTrtForceHt,bTrtPhotoHt)
names(HtCue) <- c("Chill", "Force","Photo")

HtCueL <- melt(HtCue)
HtCueL$trait <- "Height"

load("output/lmaDummyInt.Rdata")
postLMA <- rstan::extract(mdlLMA)

bTrtChillLMA <- data.frame(postLMA$betaTraitxChill)
bTrtForceLMA <- data.frame(postLMA$betaTraitxForce)
bTrtPhotoLMA <- data.frame(postLMA$betaTraitxPhoto)

LMACue <- cbind(bTrtChillLMA,bTrtForceLMA,bTrtPhotoLMA)
names(LMACue) <- c("Chill", "Force","Photo")

LMACueL <- melt(LMACue)
LMACueL$trait <- "LMA"

load("output/dbhDummyInt.Rdata")
mdlDBH <- mdl
postDBH <- rstan::extract(mdlDBH)

bTrtChillDBH <- data.frame(postDBH$betaTraitxChill)
bTrtForceDBH <- data.frame(postDBH$betaTraitxForce)
bTrtPhotoDBH <- data.frame(postDBH$betaTraitxPhoto)

DBHCue <- cbind(bTrtChillDBH,bTrtForceDBH,bTrtPhotoDBH)
names(DBHCue) <- c("Chill", "Force","Photo")

DBHCueL <- melt(DBHCue)
DBHCueL$trait <- "DBH"

load("output/CNDummyInt.Rdata")
mdlCN <- mdl
postCN <- rstan::extract(mdlCN)

bTrtChillCN <- data.frame(postCN$betaTraitxChill)
bTrtForceCN <- data.frame(postCN$betaTraitxForce)
bTrtPhotoCN <- data.frame(postCN$betaTraitxPhoto)

CNCue <- cbind(bTrtChillCN,bTrtForceCN,bTrtPhotoCN)
names(CNCue) <- c("Chill", "Force","Photo")

CNCueL <- melt(CNCue)
CNCueL$trait <- "C:N"

load("output/ssdDummyInt.Rdata")
mdlSSD <- mdl
postSSD <- rstan::extract(mdlSSD)

bTrtChillSSD <- data.frame(postSSD$betaTraitxChill)
bTrtForceSSD <- data.frame(postSSD$betaTraitxForce)
bTrtPhotoSSD <- data.frame(postSSD$betaTraitxPhoto)

SSDCue <- cbind(bTrtChillSSD,bTrtForceSSD,bTrtPhotoSSD)
names(SSDCue) <- c("Chill", "Force","Photo")

SSDCueL <- melt(SSDCue)
SSDCueL$trait <- "SSD"

trtCue <- rbind(HtCueL, LMACueL, DBHCueL, CNCueL, SSDCueL)


ggplot() +
  stat_eye(data = trtCue, aes(x = trait, y = value, fill = variable, group = variable), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9)) +
  theme_classic() +
 # theme(legend.position = "none") +
  labs( x = "Trait", y = "Trait by cue response", main = NA)+
  scale_fill_manual(values = c("cyan4", "goldenrod","maroon"))

+
  geom_text(aes(label=species),hjust= 0.5, vjust= 1.5, show.legend = F) +
  geom_errorbar(aes(ymin= bChill25, ymax = bChill75), width= 0) +
  geom_errorbar(aes(xmin= bPhoto25, xmax = bPhoto75), width= 0) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.key=element_rect(fill="white")) # removed grey boxes around legends
