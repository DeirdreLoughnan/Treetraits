# November 10, 1023
# Aim of this code is to create mu plots for the five trait models
rm(list=ls())
options(stringsAsFactors = FALSE)

require(rstan)
library(forcats)
library(ggdist)
library(reshape2)
require(cowplot)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/Treetraits") 
}  else{
  setwd("/home/deirdre/Treetraits") # for midge
}

# want plots with just the cues and plot with just the traitxcue interactions

spInfo <- read.csv("input/species_ring.csv")
pheno <- read.csv("input/phenoDataWChill.csv")
# trtPheno <- read.csv("input/trtData.csv")
trtPheno <- read.csv("input/trtPhenoDummy.csv")

load("output/htContLat6.Rdata")
sumHt <- summary(mdlHt)$summary
postHt <- rstan::extract(mdlHt)

load("output/lmaContLat.Rdata")
postLMA <- rstan::extract(mdlLMA)
sumLMA<- summary(mdlLMA)$summary

load("output/dbhContLat6.Rdata")
postDBH <- rstan::extract(mdlDBH)
sumDBH<- summary(mdlDBH)$summary

load("output/ssdContLat6.Rdata")
postSSD <- rstan::extract(mdlSSD)
sumSSD<- summary(mdlSSD)$summary

load("output/lncContLat.Rdata")
postCN <- rstan::extract(mdlPerN)
sumCN<- summary(mdlPerN)$summary
# load("output/cnContLat.Rdata")
# postCN <- rstan::extract(mdlCN)
# sumCN<- summary(mdlCN)$summary

# load("output/heightDummyIntGrandZ25.Rdata")
# sumHt <- summary(mdlHt)$summary
# postHt <- rstan::extract(mdlHt)
# 
# load("output/lmaDummyIntGrandZ25.Rdata")
# postLMA <- rstan::extract(mdlLMA)
# sumLMA<- summary(mdlLMA)$summary
# 
# load("output/dbhDummyIntGrandZ25.Rdata")
# postDBH <- rstan::extract(mdlDBH)
# sumDBH<- summary(mdlDBH)$summary
# 
# load("output/ssdDummyIntGrandZ25.Rdata")
# postSSD <- rstan::extract(mdlSSD)
# sumSSD<- summary(mdlSSD)$summary
# 
# load("output/cnDummyIntGrandZ25.Rdata")
# postCN <- rstan::extract(mdl)
# sumCN<- summary(mdl)$summary


muForceHt = mean((sumHt[grep("muForceSp", rownames(sumHt)), 1]))
muForceHt5 <- quantile(postHt$muForceSp, c(0.05))
muForceHt95 <- quantile(postHt$muForceSp, c(0.95))
muForceHt25 <- quantile(postHt$muForceSp, c(0.25))
muForceHt75 <- quantile(postHt$muForceSp, c(0.75))
muForceHt <- data.frame(cbind(muForceHt, muForceHt5,muForceHt95, muForceHt25,muForceHt75))
muForceHt$cue <- "Forcing"
names(muForceHt) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

muChillHt <- sumHt[grep("muChillSp", rownames(sumHt)), 1]
muChillHt5 <- quantile(postHt$muChillSp, c(0.05))
muChillHt95 <- quantile(postHt$muChillSp, c(0.95))
muChillHt25 <- quantile(postHt$muChillSp, c(0.25))
muChillHt75 <- quantile(postHt$muChillSp, c(0.75))
muChillHt <- data.frame(cbind(muChillHt, muChillHt5,muChillHt95, muChillHt25,muChillHt75))
muChillHt$cue <- "Chilling"
names(muChillHt) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

muPhotoHt <- sumHt[grep("muPhotoSp", rownames(sumHt)), 1]
muPhotoHt5 <- quantile(postHt$muPhotoSp, c(0.05))
muPhotoHt95 <- quantile(postHt$muPhotoSp, c(0.95))
muPhotoHt25 <- quantile(postHt$muPhotoSp, c(0.25))
muPhotoHt75 <- quantile(postHt$muPhotoSp, c(0.75))
muPhotoHt <- data.frame(cbind(muPhotoHt, muPhotoHt5,muPhotoHt95, muPhotoHt25,muPhotoHt75))
muPhotoHt$cue <- "Photoperiod"
names(muPhotoHt) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")


cueHt <- rbind(muForceHt,muChillHt,muPhotoHt)

cueHeightPlot <- ggplot(cueHt,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "cyan4") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-35,5)+
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, width = 0, color = "cyan4") +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + theme(axis.title = element_text( size=17), axis.text.y=element_text(size = 15)) +
  labs( x = "Estimated change in budburst day", y = "", main = NA) +  annotate("text", x = -30, y = 3.4, label = "a) Height", cex = 10) +
  theme(legend.title = element_blank()) 

########################################################################
muForceLMA = mean((sumLMA[grep("muForceSp", rownames(sumLMA)), 1]))
muForceLMA5 <- quantile(postLMA$muForceSp, c(0.05))
muForceLMA95 <- quantile(postLMA$muForceSp, c(0.95))
muForceLMA25 <- quantile(postLMA$muForceSp, c(0.25))
muForceLMA75 <- quantile(postLMA$muForceSp, c(0.75))
muForceLMA <- data.frame(cbind(muForceLMA, muForceLMA5,muForceLMA95, muForceLMA25,muForceLMA75))
muForceLMA$cue <- "Forcing"
names(muForceLMA) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

muChillLMA <- sumLMA[grep("muChillSp", rownames(sumLMA)), 1]
muChillLMA5 <- quantile(postLMA$muChillSp, c(0.05))
muChillLMA95 <- quantile(postLMA$muChillSp, c(0.95))
muChillLMA25 <- quantile(postLMA$muChillSp, c(0.25))
muChillLMA75 <- quantile(postLMA$muChillSp, c(0.75))
muChillLMA <- data.frame(cbind(muChillLMA, muChillLMA5,muChillLMA95, muChillLMA25,muChillLMA75))
muChillLMA$cue <- "Chilling"
names(muChillLMA) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

muPhotoLMA <- sumLMA[grep("muPhotoSp", rownames(sumLMA)), 1]
muPhotoLMA5 <- quantile(postLMA$muPhotoSp, c(0.05))
muPhotoLMA95 <- quantile(postLMA$muPhotoSp, c(0.95))
muPhotoLMA25 <- quantile(postLMA$muPhotoSp, c(0.25))
muPhotoLMA75 <- quantile(postLMA$muPhotoSp, c(0.75))
muPhotoLMA <- data.frame(cbind(muPhotoLMA, muPhotoLMA5,muPhotoLMA95, muPhotoLMA25,muPhotoLMA75))
muPhotoLMA$cue <- "Photoperiod"
names(muPhotoLMA) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

cueLMA <- rbind(muForceLMA,muChillLMA,muPhotoLMA)

cueLMAPlot <- ggplot(cueLMA,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "darkolivegreen") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-35,5)+
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, width = 0, color = "darkolivegreen") +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  theme(axis.title = element_text( size=17), axis.text.y=element_text(size = 15)) +
  labs( x = "Estimated change in budburst day", y = "", main = NA) +  annotate("text", x = -30, y = 3.4, label = "d) LMA", cex = 10) +
  theme(legend.title = element_blank()) 

########################################################################

muForceDBH = mean((sumDBH[grep("muForceSp", rownames(sumDBH)), 1]))
muForceDBH5 <- quantile(postDBH$muForceSp, c(0.05))
muForceDBH95 <- quantile(postDBH$muForceSp, c(0.95))
muForceDBH25 <- quantile(postDBH$muForceSp, c(0.25))
muForceDBH75 <- quantile(postDBH$muForceSp, c(0.75))
muForceDBH <- data.frame(cbind(muForceDBH, muForceDBH5,muForceDBH95, muForceDBH25,muForceDBH75))
muForceDBH$cue <- "Forcing"
names(muForceDBH) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

muChillDBH <- sumDBH[grep("muChillSp", rownames(sumDBH)), 1]
muChillDBH5 <- quantile(postDBH$muChillSp, c(0.05))
muChillDBH95 <- quantile(postDBH$muChillSp, c(0.95))
muChillDBH25 <- quantile(postDBH$muChillSp, c(0.25))
muChillDBH75 <- quantile(postDBH$muChillSp, c(0.75))
muChillDBH <- data.frame(cbind(muChillDBH, muChillDBH5,muChillDBH95, muChillDBH25,muChillDBH75))
muChillDBH$cue <- "Chilling"
names(muChillDBH) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

muPhotoDBH <- sumDBH[grep("muPhotoSp", rownames(sumDBH)), 1]
muPhotoDBH5 <- quantile(postDBH$muPhotoSp, c(0.05))
muPhotoDBH95 <- quantile(postDBH$muPhotoSp, c(0.95))
muPhotoDBH25 <- quantile(postDBH$muPhotoSp, c(0.25))
muPhotoDBH75 <- quantile(postDBH$muPhotoSp, c(0.75))
muPhotoDBH <- data.frame(cbind(muPhotoDBH, muPhotoDBH5,muPhotoDBH95, muPhotoDBH25,muPhotoDBH75))
muPhotoDBH$cue <- "Photoperiod"
names(muPhotoDBH) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

cueDBH <- rbind(muForceDBH,muChillDBH,muPhotoDBH)

cueDBHPlot <- ggplot(cueDBH,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "goldenrod") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-35,5)+
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, width = 0, color = "goldenrod") +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  theme(axis.title = element_text( size=17), axis.text.y=element_text(size = 15)) +
  labs( x = "Estimated change in budburst day", y = "", main = NA) +  annotate("text", x = -11, y = 3.4, label = "b) DBH", cex = 10) +
  theme(legend.title = element_blank()) 

########################################################################
muForceCN = mean((sumCN[grep("muForceSp", rownames(sumCN)), 1]))
muForceCN5 <- quantile(postCN$muForceSp, c(0.05))
muForceCN95 <- quantile(postCN$muForceSp, c(0.95))
muForceCN25 <- quantile(postCN$muForceSp, c(0.25))
muForceCN75 <- quantile(postCN$muForceSp, c(0.75))
muForceCN <- data.frame(cbind(muForceCN, muForceCN5,muForceCN95, muForceCN25,muForceCN75))
muForceCN$cue <- "Forcing"
names(muForceCN) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

muChillCN <- sumCN[grep("muChillSp", rownames(sumCN)), 1]
muChillCN5 <- quantile(postCN$muChillSp, c(0.05))
muChillCN95 <- quantile(postCN$muChillSp, c(0.95))
muChillCN25 <- quantile(postCN$muChillSp, c(0.25))
muChillCN75 <- quantile(postCN$muChillSp, c(0.75))
muChillCN <- data.frame(cbind(muChillCN, muChillCN5,muChillCN95, muChillCN25,muChillCN75))
muChillCN$cue <- "Chilling"
names(muChillCN) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

muPhotoCN <- sumCN[grep("muPhotoSp", rownames(sumCN)), 1]
muPhotoCN5 <- quantile(postCN$muPhotoSp, c(0.05))
muPhotoCN95 <- quantile(postCN$muPhotoSp, c(0.95))
muPhotoCN25 <- quantile(postCN$muPhotoSp, c(0.25))
muPhotoCN75 <- quantile(postCN$muPhotoSp, c(0.75))
muPhotoCN <- data.frame(cbind(muPhotoCN, muPhotoCN5,muPhotoCN95, muPhotoCN25,muPhotoCN75))
muPhotoCN$cue <- "Photoperiod"
names(muPhotoCN) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

cueCN <- rbind(muForceCN,muChillCN,muPhotoCN)

cueCNPlot <- ggplot(cueCN,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "purple4") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-35,5)+
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, width = 0, color = "purple4") +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  theme(axis.title = element_text( size=17), axis.text.y=element_text(size = 15)) +
  # theme(axis.text.x = element_text( size=17,angle = 78,  hjust=1),
  #       axis.text.y=element_text(size = 15),
  #       axis.title=element_text(size=  17),
  #       legend.position = "none") +
  labs( x = "Estimated change in budburst day", y = "", main = NA) +  annotate("text", x = -11.3, y = 3.4, label = "e) C:N", cex = 10) +
  theme(legend.title = element_blank()) 

########################################################################

muForceSSD = mean((sumSSD[grep("muForceSp", rownames(sumSSD)), 1]))
muForceSSD5 <- quantile(postSSD$muForceSp, c(0.05))
muForceSSD95 <- quantile(postSSD$muForceSp, c(0.95))
muForceSSD25 <- quantile(postSSD$muForceSp, c(0.25))
muForceSSD75 <- quantile(postSSD$muForceSp, c(0.75))
muForceSSD <- data.frame(cbind(muForceSSD, muForceSSD5,muForceSSD95, muForceSSD25,muForceSSD75))
muForceSSD$cue <- "Forcing"
names(muForceSSD) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

muChillSSD <- sumSSD[grep("muChillSp", rownames(sumSSD)), 1]
muChillSSD5 <- quantile(postSSD$muChillSp, c(0.05))
muChillSSD95 <- quantile(postSSD$muChillSp, c(0.95))
muChillSSD25 <- quantile(postSSD$muChillSp, c(0.25))
muChillSSD75 <- quantile(postSSD$muChillSp, c(0.75))
muChillSSD <- data.frame(cbind(muChillSSD, muChillSSD5,muChillSSD95, muChillSSD25,muChillSSD75))
muChillSSD$cue <- "Chilling"
names(muChillSSD) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

muPhotoSSD <- sumSSD[grep("muPhotoSp", rownames(sumSSD)), 1]
muPhotoSSD5 <- quantile(postSSD$muPhotoSp, c(0.05))
muPhotoSSD95 <- quantile(postSSD$muPhotoSp, c(0.95))
muPhotoSSD25 <- quantile(postSSD$muPhotoSp, c(0.25))
muPhotoSSD75 <- quantile(postSSD$muPhotoSp, c(0.75))
muPhotoSSD <- data.frame(cbind(muPhotoSSD, muPhotoSSD5,muPhotoSSD95, muPhotoSSD25,muPhotoSSD75))
muPhotoSSD$cue <- "Photoperiod"
names(muPhotoSSD) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

cueSSD <- rbind(muForceSSD,muChillSSD,muPhotoSSD)

cueSSDPlot <- ggplot(cueSSD,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "maroon") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-35,5)+
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, width = 0, color = "maroon") +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
   theme(axis.title = element_text( size=17), axis.text.y=element_text(size = 15)) +
  labs( x = "Estimated change in budburst day", y = "", main = NA) +  annotate("text", x = -30, y = 3.4, label = "c) SSD", cex = 10) +
  theme(legend.title = element_blank()) 

pdf("figures/muCuePlots1006.pdf", height =5, width = 25)
plot_grid( cueHeightPlot,cueDBHPlot,cueSSDPlot,cueLMAPlot,cueCNPlot , ncol = 5, nrow =1,align = "v")
dev.off()


#########################################
betaTraitxForceHt = mean((sumHt[grep("betaTraitxForce", rownames(sumHt)), 1]))
betaTraitxForceHt5 <- quantile(postHt$betaTraitxForce, c(0.05))
betaTraitxForceHt95 <- quantile(postHt$betaTraitxForce, c(0.95))
betaTraitxForceHt25 <- quantile(postHt$betaTraitxForce, c(0.25))
betaTraitxForceHt75 <- quantile(postHt$betaTraitxForce, c(0.75))
betaTraitxForceHt <- data.frame(cbind(betaTraitxForceHt, betaTraitxForceHt5,betaTraitxForceHt95, betaTraitxForceHt25,betaTraitxForceHt75))
betaTraitxForceHt$cue <- "Trait x Forcing"
names(betaTraitxForceHt) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

betaTraitxChillHt <- sumHt[grep("betaTraitxChill", rownames(sumHt)), 1]
betaTraitxChillHt5 <- quantile(postHt$betaTraitxChill, c(0.05))
betaTraitxChillHt95 <- quantile(postHt$betaTraitxChill, c(0.95))
betaTraitxChillHt25 <- quantile(postHt$betaTraitxChill, c(0.25))
betaTraitxChillHt75 <- quantile(postHt$betaTraitxChill, c(0.75))
betaTraitxChillHt <- data.frame(cbind(betaTraitxChillHt, betaTraitxChillHt5,betaTraitxChillHt95, betaTraitxChillHt25,betaTraitxChillHt75))
betaTraitxChillHt$cue <- "Trait x Chilling"
names(betaTraitxChillHt) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

betaTraitxPhotoHt <- sumHt[grep("betaTraitxPhoto", rownames(sumHt)), 1]
betaTraitxPhotoHt5 <- quantile(postHt$betaTraitxPhoto, c(0.05))
betaTraitxPhotoHt95 <- quantile(postHt$betaTraitxPhoto, c(0.95))
betaTraitxPhotoHt25 <- quantile(postHt$betaTraitxPhoto, c(0.25))
betaTraitxPhotoHt75 <- quantile(postHt$betaTraitxPhoto, c(0.75))
betaTraitxPhotoHt <- data.frame(cbind(betaTraitxPhotoHt, betaTraitxPhotoHt5,betaTraitxPhotoHt95, betaTraitxPhotoHt25,betaTraitxPhotoHt75))
betaTraitxPhotoHt$cue <- "Trait x Photoperiod"
names(betaTraitxPhotoHt) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")


cueTraitHt <- rbind(betaTraitxForceHt,betaTraitxChillHt,betaTraitxPhotoHt)

cueTraitHeightPlot <- ggplot(cueTraitHt,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "cyan4") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-15,40) +
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, width = 0, color = "cyan4") +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + theme(axis.title = element_text( size=17), axis.text.y=element_text(size = 15)) +
  labs( x = "Estimated cue response in budburst days", y = "", main = NA) +  annotate("text", x = -11, y = 3.4, label = "a) Height", cex = 10) +
  theme(legend.title = element_blank()) 

########################################################################
betaTraitxForceLMA = mean((sumLMA[grep("betaTraitxForce", rownames(sumLMA)), 1]))
betaTraitxForceLMA5 <- quantile(postLMA$betaTraitxForce, c(0.05))
betaTraitxForceLMA95 <- quantile(postLMA$betaTraitxForce, c(0.95))
betaTraitxForceLMA25 <- quantile(postLMA$betaTraitxForce, c(0.25))
betaTraitxForceLMA75 <- quantile(postLMA$betaTraitxForce, c(0.75))
betaTraitxForceLMA <- data.frame(cbind(betaTraitxForceLMA, betaTraitxForceLMA5,betaTraitxForceLMA95, betaTraitxForceLMA25,betaTraitxForceLMA75))

betaTraitxForceLMA$cue <- "Trait x Forcing"
names(betaTraitxForceLMA) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")


betaTraitxChillLMA <- sumLMA[grep("betaTraitxChill", rownames(sumLMA)), 1]
betaTraitxChillLMA5 <- quantile(postLMA$betaTraitxChill, c(0.05))
betaTraitxChillLMA95 <- quantile(postLMA$betaTraitxChill, c(0.95))
betaTraitxChillLMA25 <- quantile(postLMA$betaTraitxChill, c(0.25))
betaTraitxChillLMA75 <- quantile(postLMA$betaTraitxChill, c(0.75))
betaTraitxChillLMA <- data.frame(cbind(betaTraitxChillLMA, betaTraitxChillLMA5,betaTraitxChillLMA95, betaTraitxChillLMA25,betaTraitxChillLMA75))

betaTraitxChillLMA$cue <- "Trait x Chilling"
names(betaTraitxChillLMA) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

betaTraitxPhotoLMA <- sumLMA[grep("betaTraitxPhoto", rownames(sumLMA)), 1]
betaTraitxPhotoLMA5 <- quantile(postLMA$betaTraitxPhoto, c(0.05))
betaTraitxPhotoLMA95 <- quantile(postLMA$betaTraitxPhoto, c(0.95))
betaTraitxPhotoLMA25 <- quantile(postLMA$betaTraitxPhoto, c(0.25))
betaTraitxPhotoLMA75 <- quantile(postLMA$betaTraitxPhoto, c(0.75))
betaTraitxPhotoLMA <- data.frame(cbind(betaTraitxPhotoLMA, betaTraitxPhotoLMA5,betaTraitxPhotoLMA95, betaTraitxPhotoLMA25,betaTraitxPhotoLMA75))

betaTraitxPhotoLMA$cue <- "Trait x Photoperiod"
names(betaTraitxPhotoLMA) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

cueTraitLMA <- rbind(betaTraitxForceLMA,betaTraitxChillLMA,betaTraitxPhotoLMA)

cueTraitLMAPlot <- ggplot(cueTraitLMA,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "darkolivegreen") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-1,1) +
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, width = 0, color = "darkolivegreen") +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  theme(axis.title = element_text( size=17), axis.text.y=element_text(size = 15)) +
  labs( x = "Estimated change in budburst day", y = "", main = NA) +  annotate("text", x = -8, y = 3.4, label = "d) LMA", cex = 10) +
  theme(legend.title = element_blank()) 

########################################################################

betaTraitxForceDBH = mean((sumDBH[grep("betaTraitxForce", rownames(sumDBH)), 1]))
betaTraitxForceDBH5 <- quantile(postDBH$betaTraitxForce, c(0.05))
betaTraitxForceDBH95 <- quantile(postDBH$betaTraitxForce, c(0.95))
betaTraitxForceDBH25 <- quantile(postDBH$betaTraitxForce, c(0.25))
betaTraitxForceDBH75 <- quantile(postDBH$betaTraitxForce, c(0.75))
betaTraitxForceDBH <- data.frame(cbind(betaTraitxForceDBH, betaTraitxForceDBH5,betaTraitxForceDBH95, betaTraitxForceDBH25,betaTraitxForceDBH75))
betaTraitxForceDBH$cue <- "Trait x Forcing"
names(betaTraitxForceDBH) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

betaTraitxChillDBH <- sumDBH[grep("betaTraitxChill", rownames(sumDBH)), 1]
betaTraitxChillDBH5 <- quantile(postDBH$betaTraitxChill, c(0.05))
betaTraitxChillDBH95 <- quantile(postDBH$betaTraitxChill, c(0.95))
betaTraitxChillDBH25 <- quantile(postDBH$betaTraitxChill, c(0.25))
betaTraitxChillDBH75 <- quantile(postDBH$betaTraitxChill, c(0.75))
betaTraitxChillDBH <- data.frame(cbind(betaTraitxChillDBH, betaTraitxChillDBH5,betaTraitxChillDBH95, betaTraitxChillDBH25,betaTraitxChillDBH75))
betaTraitxChillDBH$cue <- "Trait x Chilling"
names(betaTraitxChillDBH) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

betaTraitxPhotoDBH <- sumDBH[grep("betaTraitxPhoto", rownames(sumDBH)), 1]
betaTraitxPhotoDBH5 <- quantile(postDBH$betaTraitxPhoto, c(0.05))
betaTraitxPhotoDBH95 <- quantile(postDBH$betaTraitxPhoto, c(0.95))
betaTraitxPhotoDBH25 <- quantile(postDBH$betaTraitxPhoto, c(0.25))
betaTraitxPhotoDBH75 <- quantile(postDBH$betaTraitxPhoto, c(0.75))
betaTraitxPhotoDBH <- data.frame(cbind(betaTraitxPhotoDBH, betaTraitxPhotoDBH5,betaTraitxPhotoDBH95, betaTraitxPhotoDBH25,betaTraitxPhotoDBH75))
betaTraitxPhotoDBH$cue <- "Trait x Photoperiod"
names(betaTraitxPhotoDBH) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

cueTraitDBH <- rbind(betaTraitxForceDBH,betaTraitxChillDBH,betaTraitxPhotoDBH)

cueTraitDBHPlot <- ggplot(cueTraitDBH,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "goldenrod") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-15,40) +
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, width = 0, color = "goldenrod") +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  theme(axis.title = element_text( size=17), axis.text.y=element_text(size = 15)) +
  labs( x = "Estimated change in budburst day", y = "", main = NA) +  annotate("text", x = -11, y = 3.4, label = "b) DBH", cex = 10) +
  theme(legend.title = element_blank()) 

########################################################################
betaTraitxForceCN = mean((sumCN[grep("betaTraitxForce", rownames(sumCN)), 1]))
betaTraitxForceCN5 <- quantile(postCN$betaTraitxForce, c(0.05))
betaTraitxForceCN95 <- quantile(postCN$betaTraitxForce, c(0.95))
betaTraitxForceCN25 <- quantile(postCN$betaTraitxForce, c(0.25))
betaTraitxForceCN75 <- quantile(postCN$betaTraitxForce, c(0.75))
betaTraitxForceCN <- data.frame(cbind(betaTraitxForceCN, betaTraitxForceCN5,betaTraitxForceCN95, betaTraitxForceCN25,betaTraitxForceCN75))
betaTraitxForceCN$cue <- "Trait x Forcing"
names(betaTraitxForceCN) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

betaTraitxChillCN <- sumCN[grep("betaTraitxChill", rownames(sumCN)), 1]
betaTraitxChillCN5 <- quantile(postCN$betaTraitxChill, c(0.05))
betaTraitxChillCN95 <- quantile(postCN$betaTraitxChill, c(0.95))
betaTraitxChillCN25 <- quantile(postCN$betaTraitxChill, c(0.25))
betaTraitxChillCN75 <- quantile(postCN$betaTraitxChill, c(0.75))
betaTraitxChillCN <- data.frame(cbind(betaTraitxChillCN, betaTraitxChillCN5,betaTraitxChillCN95, betaTraitxChillCN25,betaTraitxChillCN75))
betaTraitxChillCN$cue <- "Trait x Chilling"
names(betaTraitxChillCN) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

betaTraitxPhotoCN <- sumCN[grep("betaTraitxPhoto", rownames(sumCN)), 1]
betaTraitxPhotoCN5 <- quantile(postCN$betaTraitxPhoto, c(0.05))
betaTraitxPhotoCN95 <- quantile(postCN$betaTraitxPhoto, c(0.95))
betaTraitxPhotoCN25 <- quantile(postCN$betaTraitxPhoto, c(0.25))
betaTraitxPhotoCN75 <- quantile(postCN$betaTraitxPhoto, c(0.75))
betaTraitxPhotoCN <- data.frame(cbind(betaTraitxPhotoCN, betaTraitxPhotoCN5,betaTraitxPhotoCN95, betaTraitxPhotoCN25,betaTraitxPhotoCN75))
betaTraitxPhotoCN$cue <- "Trait x Photoperiod"
names(betaTraitxPhotoCN) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

cueTraitCN <- rbind(betaTraitxForceCN,betaTraitxChillCN,betaTraitxPhotoCN)

cueTraitCNPlot <- ggplot(cueTraitCN,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "purple4") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-15,40) +
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, width = 0, color = "purple4") +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  theme(axis.title = element_text( size=17), axis.text.y=element_text(size = 15)) +
  # theme(axis.text.x = element_text( size=17,angle = 78,  hjust=1),
  #       axis.text.y=element_text(size = 15),
  #       axis.title=element_text(size=  17),
  #       legend.position = "none") +
  labs( x = "Estimated change in budburst day", y = "", main = NA) +  annotate("text", x = -11, y = 3.4, label = "e) C:N", cex = 10) +
  theme(legend.title = element_blank()) 

########################################################################

betaTraitxForceSSD = mean((sumSSD[grep("betaTraitxForce", rownames(sumSSD)), 1]))
betaTraitxForceSSD5 <- quantile(postSSD$betaTraitxForce, c(0.05))
betaTraitxForceSSD95 <- quantile(postSSD$betaTraitxForce, c(0.95))
betaTraitxForceSSD25 <- quantile(postSSD$betaTraitxForce, c(0.25))
betaTraitxForceSSD75 <- quantile(postSSD$betaTraitxForce, c(0.75))
betaTraitxForceSSD <- data.frame(cbind(betaTraitxForceSSD, betaTraitxForceSSD5,betaTraitxForceSSD95, betaTraitxForceSSD25,betaTraitxForceSSD75))
betaTraitxForceSSD$cue <- "Trait x Forcing"
names(betaTraitxForceSSD) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

betaTraitxChillSSD <- sumSSD[grep("betaTraitxChill", rownames(sumSSD)), 1]
betaTraitxChillSSD5 <- quantile(postSSD$betaTraitxChill, c(0.05))
betaTraitxChillSSD95 <- quantile(postSSD$betaTraitxChill, c(0.95))
betaTraitxChillSSD25 <- quantile(postSSD$betaTraitxChill, c(0.25))
betaTraitxChillSSD75 <- quantile(postSSD$betaTraitxChill, c(0.75))
betaTraitxChillSSD <- data.frame(cbind(betaTraitxChillSSD, betaTraitxChillSSD5,betaTraitxChillSSD95, betaTraitxChillSSD25,betaTraitxChillSSD75))
betaTraitxChillSSD$cue <- "Trait x Chilling"
names(betaTraitxChillSSD) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

betaTraitxPhotoSSD <- sumSSD[grep("betaTraitxPhoto", rownames(sumSSD)), 1]
betaTraitxPhotoSSD5 <- quantile(postSSD$betaTraitxPhoto, c(0.05))
betaTraitxPhotoSSD95 <- quantile(postSSD$betaTraitxPhoto, c(0.95))
betaTraitxPhotoSSD25 <- quantile(postSSD$betaTraitxPhoto, c(0.25))
betaTraitxPhotoSSD75 <- quantile(postSSD$betaTraitxPhoto, c(0.75))
betaTraitxPhotoSSD <- data.frame(cbind(betaTraitxPhotoSSD, betaTraitxPhotoSSD5,betaTraitxPhotoSSD95, betaTraitxPhotoSSD25,betaTraitxPhotoSSD75))
betaTraitxPhotoSSD$cue <- "Trait x Photoperiod"
names(betaTraitxPhotoSSD) <- c("mean", "five", "nintyFive","twentyFive","seventyFive","cue")

cueTraitSSD <- rbind(betaTraitxForceSSD,betaTraitxChillSSD,betaTraitxPhotoSSD)

cueTraitSSDPlot <- ggplot(cueTraitSSD,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "maroon") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-15,40) +
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, color = "maroon", width = 0) +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  theme(axis.title = element_text( size=17), axis.text.y=element_text(size = 15)) +
  labs( x = "Estimated change in budburst day", y = "", main = NA) +  annotate("text", x = -11, y = 3.4, label = "c) SSD", cex = 10) +
  theme(legend.title = element_blank()) 


#####
cueHeightPlot <- ggplot(cueHt,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "cyan4") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-25,5)+
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, width = 0, color = "cyan4") +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + theme(axis.title = element_text( size = 30), axis.text.y=element_text(size = 15), axis.text.x=element_text(size = 15)) +
  labs( x = "", y = "", main = NA) +  annotate("text", x = -16, y = 3.4, label = "a) Height", cex = 10) + scale_y_discrete(limits=rev) +
  theme(legend.title = element_blank()) 

cueLMAPlot <- ggplot(cueLMA,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "darkolivegreen") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-25,5)+
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, width = 0, color = "darkolivegreen") +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  theme(axis.title = element_text( size = 30), axis.text.y=element_text(size = 15), axis.text.x=element_text(size = 15)) +
  labs( x = "", y = "", main = NA) +  annotate("text", x = -19, y = 3.4, label = "d) LMA", cex = 10) + scale_y_discrete(limits=rev) +
  theme(legend.title = element_blank()) 

cueDBHPlot <- ggplot(cueDBH,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "goldenrod") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-25,5)+
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, width = 0, color = "goldenrod") +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  theme(axis.title = element_text( size = 30), axis.text.y=element_text(size = 15), axis.text.x=element_text(size = 15)) +
  labs( x = "", y = "", main = NA) +  annotate("text", x = -19, y = 3.4, label = "b) DBH", cex = 10) + scale_y_discrete(limits=rev) +
  theme(legend.title = element_blank()) 

cueCNPlot <- ggplot(cueCN,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "purple4") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-25,5)+
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, width = 0, color = "purple4") +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  theme(axis.title = element_text( size=30), axis.text.y=element_text(size = 15), axis.text.x=element_text(size = 15)) + scale_y_discrete(limits=rev) +
  labs( x = "", y = "", main = NA) +  annotate("text", x = -19, y = 3.4, label = "e) LNC", cex = 10) +
  theme(legend.title = element_blank()) 


cueSSDPlot <- ggplot(cueSSD,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "maroon") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-25,5)+
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, width = 0, color = "maroon") +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  theme(axis.title = element_text( size=30), axis.text.y=element_text(size = 15), axis.text.x=element_text(size = 15)) +
  labs( x = "Change in budburst day     ", y = "", main = NA) +  annotate("text", x = -19, y = 3.4, label = "c) SSD", cex = 10) + scale_y_discrete(limits=rev) +
  theme(legend.title = element_blank()) 

pdf("figures/muCuePlots1006.pdf", height =5, width = 25)
plot_grid( cueHeightPlot,cueDBHPlot,cueSSDPlot,cueLMAPlot,cueCNPlot , ncol = 5, nrow =1 ,align = "v", rel_widths = c(0.9,0.9,1,0.9,0.9))
dev.off()
#####################################################

cueTraitHeightPlot <- ggplot(cueTraitHt,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "cyan4") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-1,1) +
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, width = 0, color = "cyan4") +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + theme(axis.title = element_text( size=17), axis.text.y=element_text(size = 15), axis.text.x=element_text(size = 15)) +
  labs( x = "", y = "", main = NA) +  annotate("text", x = -0.4, y = 3.4, label = "a) Height", cex = 10) + scale_y_discrete(limits=rev) +
  theme(legend.title = element_blank()) 

cueTraitLMAPlot <- ggplot(cueTraitLMA,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "darkolivegreen") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-5,5) +
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, width = 0, color = "darkolivegreen") +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  theme(axis.title = element_text( size=17), axis.text.y=element_text(size = 15), axis.text.x=element_text(size = 15)) +
  labs( x = "", y = "", main = NA) +  annotate("text", x = -2, y = 3.4, label = "d) LMA", cex = 10) + scale_y_discrete(limits=rev) +
  theme(legend.title = element_blank()) 

##
cueTraitDBHPlot <- ggplot(cueTraitDBH,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "goldenrod") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-1,1) +
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, width = 0, color = "goldenrod") +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  theme(axis.title = element_text( size=17), axis.text.y=element_text(size = 15), axis.text.x=element_text(size = 15)) +
  labs( x = "", y = "", main = NA) +  annotate("text", x = -0.5, y = 3.4, label = "b) DBH", cex = 10) + scale_y_discrete(limits=rev) +
  theme(legend.title = element_blank()) 

cueTraitCNPlot <- ggplot(cueTraitCN,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "purple4") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-5,5) +
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, width = 0, color = "purple4") +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  theme(axis.title = element_text( size=17), axis.text.y=element_text(size = 15), axis.text.x=element_text(size = 15)) +
  labs( x = "", y = "", main = NA) +  annotate("text", x = -2, y = 3.4, label = "e) LNC", cex = 10) + scale_y_discrete(limits=rev) +
  theme(legend.title = element_blank()) 

cueTraitSSDPlot <- ggplot(cueTraitSSD,aes(y= cue, x = mean), size = 7) +
  geom_point(size = 7, color = "maroon") +
  geom_vline(xintercept = 0, linetype='dashed') +
  xlim (-20,20) +
  geom_errorbar(aes(xmin= five, xmax = nintyFive), size = 1, color = "maroon", width = 0) +
  #geom_errorbar(aes(xmin= twentyFive, xmax = seventyFive, ymin= cue, ymax = cue), size =2.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  theme(axis.title = element_text( size=30), axis.text.y=element_text(size = 15), axis.text.x=element_text(size = 15)) +
  labs( x = "Change in cue", y = "", main = NA) +  annotate("text", x = -5, y = 3.4, label = "c) SSD", cex = 10) + scale_y_discrete(limits=rev) +
  theme(legend.title = element_blank()) 

pdf("figures/muCueTraitPlots1006.pdf", height =5, width = 25)
plot_grid( cueTraitHeightPlot,cueTraitDBHPlot,cueTraitSSDPlot,cueTraitLMAPlot,cueTraitCNPlot , ncol = 5, nrow =1,align = "hv")
dev.off()
# , rel_widths = c(1, 0.05, 1)

