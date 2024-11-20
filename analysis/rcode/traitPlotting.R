# started March 18, 2023 by Deirdre

# aim of this code is to make the plots for BC traits projects

#1. Figure comapring the site differences for each trait
#2. Figure comparing the betatraitCues for each figure

rm(list=ls())
options(stringsAsFactors = FALSE)

library(stringr)
library(plyr)
library(rstan)
library(reshape2)
library(cowplot)
require(dplyr)
require(plotrix)
require(ggplot2)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/Treetraits") 
}  else{
  setwd("~/Documents/github/Treetraits") # for midge
}

spInfo <- read.csv("analysis/input/species_ring.csv")
 pheno <- read.csv("analysis/input/phenoDataWChill.csv")
# trtPheno <- read.csv("input/trtData.csv")
trtPheno <- read.csv("analysis/input/trtPhenoDummy.csv")

# Ht and DBH and CN natural, LMA and SSD rescaled by 100
# cues and latitude still z scored by only by 1 sd
load("analysis/output/htContLatHundoLatFinal.Rdata")
sumerht <- summary(mdlHt)$summary
postHt <- rstan::extract(mdlHt)

load("analysis/output/lmaContLatHundoLatFinal.Rdata")
postLMA <- rstan::extract(mdlLMA)
#sumerlma <- summary(mdlLMA)$summary

load("analysis/output/dbhContLatHundoLatFinal.Rdata")
postDBH <- rstan::extract(mdlDBH)

load("analysis/output/ssdContLatHundoLatFinal.Rdata")
postSSD <- rstan::extract(mdlSSD6)

load("analysis/output/lncContLatHundoLatFinal.Rdata")
postCN <- rstan::extract(mdlPerN)

# How do transect effects differ?
sumHt <- summary(mdlHt)$summary

a_spHt = mean((sumHt[grep("mu_grand", rownames(sumHt)), 1]))
a_spHt5 <- quantile(postHt$mu_grand, c(0.05))
a_spHt95 <- quantile(postHt$mu_grand, c(0.95))
a_spHt25 <- quantile(postHt$mu_grand, c(0.25))
a_spHt75 <- quantile(postHt$mu_grand, c(0.75))
a_spHt <- cbind(a_spHt, a_spHt5,a_spHt95, a_spHt25,a_spHt75)
#a_spHt <- a_spHt*100

b_tranHt <- sumHt[grep("b_tranE", rownames(sumHt)), 1]
b_tranHt5 <- quantile(postHt$b_tranE, c(0.05))
b_tranHt95 <- quantile(postHt$b_tranE, c(0.95))
b_tranHt25 <- quantile(postHt$b_tranE, c(0.25))
b_tranHt75 <- quantile(postHt$b_tranE, c(0.75))
b_tranHt <- cbind(b_tranHt, b_tranHt5,b_tranHt95, b_tranHt25,b_tranHt75)
#b_tranHt <- b_tranHt*100

b_tranlatHt <- sumHt[grep("b_tranlat", rownames(sumHt)), 1]
b_tranlatHt5 <- quantile(postHt$b_tranlat, c(0.05))
b_tranlatHt95 <- quantile(postHt$b_tranlat, c(0.95))
b_tranlatHt25 <- quantile(postHt$b_tranlat, c(0.25))
b_tranlatHt75 <- quantile(postHt$b_tranlat, c(0.75))
b_tranlatHt <- cbind(b_tranlatHt, b_tranlatHt5,b_tranlatHt95, b_tranlatHt25,b_tranlatHt75)
#b_tranlatHt <- b_tranlatHt*100
## Simulate interaction with transect and latitude:

eData <- subset(trtPheno, transect == "1" )
wData <- subset(trtPheno, transect == "0" )

# Make the other parameters constant

lati <- seq(40, 60, by = 0.5)
#lati <- (lati/100)
tranW <- 0
tranE <- 1

# plot first for west coast
ht_w = a_spHt[1] + b_tranHt[1] * tranW + b_tranlatHt[1] * (tranW*lati)
ht_e = a_spHt[1] + b_tranHt[1] * tranE + b_tranlatHt[1] * (tranE*lati)

ht_w5 = a_spHt[2] + b_tranHt[2] * tranW + b_tranlatHt[2] * (tranW*lati)
ht_e5 = a_spHt[2] + b_tranHt[2] * tranE + b_tranlatHt[2] * (tranE*lati)

ht_w95 = a_spHt[3] + b_tranHt[3] * tranW + b_tranlatHt[3] * (tranW*lati)
ht_e95 = a_spHt[3] + b_tranHt[3] * tranE + b_tranlatHt[3] * (tranE*lati)

htEW <- data.frame(htw = ht_w, hte = ht_e,  ht_w5 =ht_w5, ht_w95 = ht_w95, ht_e5 = ht_e5, ht_e95 = ht_e95  )

intHt <- ggplot(htEW) +
  geom_line(aes(y = htw, x = lati), col = "cyan4", lty = "dashed") +
  geom_ribbon(data = htEW, aes(ymin = ht_w5, ymax = ht_w95, x= lati), alpha = 0.2, fill = "cyan4") +
  geom_line(aes(y = hte, x = lati), col = "cyan3") + 
  geom_ribbon(data = htEW, aes(ymin = ht_e5, ymax = ht_e95, x= lati), alpha = 0.2, fill = "cyan3") + 
  xlab("Latitude") + ylab("Height (m)") +
  xlim (min(lati), max(lati)) + 
  ylim (-70,70) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 15), axis.title = element_text(size = 20))+
  theme(legend.key=element_blank(),legend.text = element_text(size = 15)) +
  #scale_fill_manual( labels = c("Low force", "High force")) +
  scale_color_manual(values = c("cyan3","cyan4"), labels = c("Eastern", "Western"), name = "") +
  #scale_colour_discrete(labels=c("High forcing","Low forcing"), name = "") +
  theme(legend.title = element_blank()) +  annotate("text", x = 41, y = 70, label = "a)", cex = 10) 
  

################################################
##lma

sumLMA <- summary(mdlLMA)$summary
# postLMA <- rstan::extract(mdlLMAhundo)
# lati <- seq(40, 60, by = 0.5)
# lati <- (lati-mean(lati,na.rm=TRUE))/(sd(lati,na.rm=TRUE))
# tranW <- 0
# tranE <- 1

a_spLMA = mean((sumLMA[grep("mu_grand", rownames(sumLMA)), 1]))
a_spLMA5 <- quantile(postLMA$mu_grand, c(0.05))
a_spLMA95 <- quantile(postLMA$mu_grand, c(0.95))
a_spLMA25 <- quantile(postLMA$mu_grand, c(0.25))
a_spLMA75 <- quantile(postLMA$mu_grand, c(0.75))
a_spLMA <- cbind(a_spLMA, a_spLMA5,a_spLMA95, a_spLMA25,a_spLMA75)
a_spLMA <- a_spLMA/100

b_tranLMA <- sumLMA[grep("b_tranE", rownames(sumLMA)), 1]
b_tranLMA5 <- quantile(postLMA$b_tranE, c(0.05))
b_tranLMA95 <- quantile(postLMA$b_tranE, c(0.95))
b_tranLMA25 <- quantile(postLMA$b_tranE, c(0.25))
b_tranLMA75 <- quantile(postLMA$b_tranE, c(0.75))
b_tranLMA <- cbind(b_tranLMA, b_tranLMA5,b_tranLMA95, b_tranLMA25,b_tranLMA75)
b_tranLMA <- b_tranLMA/100

b_tranlatLMA <- sumLMA[grep("b_tranlat", rownames(sumLMA)), 1]
b_tranlatLMA5 <- quantile(postLMA$b_tranlat, c(0.05))
b_tranlatLMA95 <- quantile(postLMA$b_tranlat, c(0.95))
b_tranlatLMA25 <- quantile(postLMA$b_tranlat, c(0.25))
b_tranlatLMA75 <- quantile(postLMA$b_tranlat, c(0.75))
b_tranlatLMA <- cbind(b_tranlatLMA, b_tranlatLMA5,b_tranlatLMA95, b_tranlatLMA25,b_tranlatLMA75)
b_tranlatLMA <- b_tranlatLMA/100
eData <- subset(trtPheno, transect == "1" )
wData <- subset(trtPheno, transect == "0" )

# Make the other parameters constant

# plot first for west coast
LMA_w = a_spLMA[1] + b_tranLMA[1] * tranW + b_tranlatLMA[1] * (tranW*lati)
LMA_e = a_spLMA[1] + b_tranLMA[1] * tranE + b_tranlatLMA[1] * (tranE*lati)

LMA_w5 = a_spLMA[2] + b_tranLMA[2] * tranW + b_tranlatLMA[2] * (tranW*lati)
LMA_e5 = a_spLMA[2] + b_tranLMA[2] * tranE + b_tranlatLMA[2] * (tranE*lati)

LMA_w95 = a_spLMA[3] + b_tranLMA[3] * tranW + b_tranlatLMA[3] * (tranW*lati)
LMA_e95 = a_spLMA[3] + b_tranLMA[3] * tranE + b_tranlatLMA[3] * (tranE*lati)

LMAEW <- data.frame(LMAw = LMA_w, LMAe = LMA_e,  LMA_w5 =LMA_w5, LMA_w95 = LMA_w95, LMA_e5 = LMA_e5, LMA_e95 = LMA_e95  )


intLMA <- ggplot(LMAEW) +
  geom_line(aes(y = (LMAw), x = lati), col = "darkolivegreen", lty = 2) +
  geom_ribbon(data = LMAEW, aes(ymin = (LMA_w5), ymax = (LMA_w95), x= lati), alpha = 0.2, fill = "darkolivegreen") +
  geom_line(aes(y = (LMAe), x = lati), col = "darkolivegreen3") +
  geom_ribbon(data = LMAEW, aes(ymin = (LMA_e5), ymax = (LMA_e95), x= lati), alpha = 0.2, fill = "darkolivegreen3") +
  scale_color_manual(values = c("darkolivegreen","darkolivegreen3"), labels = c("Western", "Eastern"), name = "") +
  xlab("Latitude") + labs(y = bquote('Leaf mass area '~(g/m^2))) +
  xlim (min(lati), max(lati)) + 
  ylim (-0.4,0.4) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 15), axis.title = element_text(size = 20))+
  theme(legend.key=element_blank(), legend.position=c(.8,.85),legend.text = element_text(size = 15)) +
  #scale_fill_manual( labels = c("Low force", "High force")) +
  #scale_colour_discrete(labels=c("High forcing","Low forcing"), name = "") +
  theme(legend.title = element_blank()) +  annotate("text", x = 41, y = 0.4, label = "d)", cex = 10) 
### DBH


sumDBH <- summary(mdlDBH)$summary

a_spDBH = mean((sumDBH[grep("mu_grand", rownames(sumDBH)), 1]))
a_spDBH5 <- quantile(postDBH$mu_grand, c(0.05))
a_spDBH95 <- quantile(postDBH$mu_grand, c(0.95))
a_spDBH25 <- quantile(postDBH$mu_grand, c(0.25))
a_spDBH75 <- quantile(postDBH$mu_grand, c(0.75))
a_spDBH <- cbind(a_spDBH, a_spDBH5,a_spDBH95, a_spDBH25,a_spDBH75)

b_tranDBH <- sumDBH[grep("b_tranE", rownames(sumDBH)), 1]
b_tranDBH5 <- quantile(postDBH$b_tranE, c(0.05))
b_tranDBH95 <- quantile(postDBH$b_tranE, c(0.95))
b_tranDBH25 <- quantile(postDBH$b_tranE, c(0.25))
b_tranDBH75 <- quantile(postDBH$b_tranE, c(0.75))
b_tranDBH <- cbind(b_tranDBH, b_tranDBH5,b_tranDBH95, b_tranDBH25,b_tranDBH75)

b_tranlatDBH <- sumDBH[grep("b_tranlat", rownames(sumDBH)), 1]
b_tranlatDBH5 <- quantile(postDBH$b_tranlat, c(0.05))
b_tranlatDBH95 <- quantile(postDBH$b_tranlat, c(0.95))
b_tranlatDBH25 <- quantile(postDBH$b_tranlat, c(0.25))
b_tranlatDBH75 <- quantile(postDBH$b_tranlat, c(0.75))
b_tranlatDBH <- cbind(b_tranlatDBH, b_tranlatDBH5,b_tranlatDBH95, b_tranlatDBH25,b_tranlatDBH75)



## Simulate interaction with transect and latitude:

eData <- subset(trtPheno, transect == "1" )
wData <- subset(trtPheno, transect == "0" )

# Make the other parameters constant

# lati <- seq(40, 60, by = 0.5)
# lati <- (lati-mean(lati,na.rm=TRUE))/(sd(lati,na.rm=TRUE)*2)
# tranW <- -0.4406491
# tranE <- 0.5669498


# plot first for west coast
DBH_w = a_spDBH[1] + b_tranDBH[1] * tranW + b_tranlatDBH[1] * (tranW*lati)
DBH_e = a_spDBH[1] + b_tranDBH[1] * tranE + b_tranlatDBH[1] * (tranE*lati)

DBH_w5 = a_spDBH[2] + b_tranDBH[2] * tranW + b_tranlatDBH[2] * (tranW*lati)
DBH_e5 = a_spDBH[2] + b_tranDBH[2] * tranE + b_tranlatDBH[2] * (tranE*lati)

DBH_w95 = a_spDBH[3] + b_tranDBH[3] * tranW + b_tranlatDBH[3] * (tranW*lati)
DBH_e95 = a_spDBH[3] + b_tranDBH[3] * tranE + b_tranlatDBH[3] * (tranE*lati)

DBHEW <- data.frame(DBHw = DBH_w, DBHe = DBH_e,  DBH_w5 =DBH_w5, DBH_w95 = DBH_w95, DBH_e5 = DBH_e5, DBH_e95 = DBH_e95  )


intDBH <- ggplot(DBHEW) +
  geom_line(aes(y = DBHw, x = lati), col = "goldenrod4", lty =2) +
  geom_ribbon(data = DBHEW, aes(ymin = DBH_w5, ymax = DBH_w95, x= lati), alpha = 0.2, fill = "goldenrod4") +
  geom_line(aes(y = DBHe, x = lati), col = "goldenrod") +
  geom_ribbon(data = DBHEW, aes(ymin = DBH_e5, ymax = DBH_e95, x= lati), alpha = 0.2, fill = "goldenrod") +
  #scale_color_manual(values = c("goldenrod4","goldenrod"), labels = c("Western", "Eastern"), name = "") +
  scale_linetype_manual(values = c(1,2), labels = c("Western", "Eastern"), name = "") +
  xlab("Latitude") + ylab("Diameter at breast height (m)") +
  xlim (min(lati), max(lati)) + 
  ylim (-70,70) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 15), axis.title = element_text(size = 20))+
  theme(legend.key=element_blank(), legend.position=c(.8,.85),legend.text = element_text(size = 15)) +
  #scale_fill_manual( labels = c("Low force", "High force")) +
  #scale_colour_discrete(labels=c("High forcing","Low forcing"), name = "") +
  theme(legend.title = element_blank()) +  annotate("text", x = 41, y = 70, label = "b)", cex = 10) 
## C:N

sumCN <- summary(mdlPerN)$summary

a_spCN = mean((sumCN[grep("mu_grand", rownames(sumCN)), 1]))
a_spCN5 <- quantile(postCN$mu_grand, c(0.05))
a_spCN95 <- quantile(postCN$mu_grand, c(0.95))
a_spCN25 <- quantile(postCN$mu_grand, c(0.25))
a_spCN75 <- quantile(postCN$mu_grand, c(0.75))
a_spCN <- cbind(a_spCN, a_spCN5,a_spCN95, a_spCN25,a_spCN75)

b_tranCN <- sumCN[grep("b_tranE", rownames(sumCN)), 1]
b_tranCN5 <- quantile(postCN$b_tranE, c(0.05))
b_tranCN95 <- quantile(postCN$b_tranE, c(0.95))
b_tranCN25 <- quantile(postCN$b_tranE, c(0.25))
b_tranCN75 <- quantile(postCN$b_tranE, c(0.75))
b_tranCN <- cbind(b_tranCN, b_tranCN5,b_tranCN95, b_tranCN25,b_tranCN75)

b_tranlatCN <- sumCN[grep("b_tranlat", rownames(sumCN)), 1]
b_tranlatCN5 <- quantile(postCN$b_tranlat, c(0.05))
b_tranlatCN95 <- quantile(postCN$b_tranlat, c(0.95))
b_tranlatCN25 <- quantile(postCN$b_tranlat, c(0.25))
b_tranlatCN75 <- quantile(postCN$b_tranlat, c(0.75))
b_tranlatCN <- cbind(b_tranlatCN, b_tranlatCN5,b_tranlatCN95, b_tranlatCN25,b_tranlatCN75)

# plot first for west coast
CN_w = a_spCN[1] + b_tranCN[1] * tranW + b_tranlatCN[1] * (tranW*lati)
CN_e = a_spCN[1] + b_tranCN[1] * tranE + b_tranlatCN[1] * (tranE*lati)

CN_w5 = a_spCN[2] + b_tranCN[2] * tranW + b_tranlatCN[2] * (tranW*lati)
CN_e5 = a_spCN[2] + b_tranCN[2] * tranE + b_tranlatCN[2] * (tranE*lati)

CN_w95 = a_spCN[3] + b_tranCN[3] * tranW + b_tranlatCN[3] * (tranW*lati)
CN_e95 = a_spCN[3] + b_tranCN[3] * tranE + b_tranlatCN[3] * (tranE*lati)

CNEW <- data.frame(CNw = CN_w, CN_w5 =CN_w5, CN_w95 = CN_w95, CNe = CN_e, CN_e5 = CN_e5, CN_e95 = CN_e95 )

intNit <- ggplot(CNEW) +
  geom_line(aes(y = (CNw), x = lati), col = "purple4", lty = 2) +
  geom_ribbon(data = CNEW, aes(ymin = (CN_w5), ymax = (CN_w95), x= lati), alpha = 0.2, fill = "purple4") +
  geom_line(aes(y = (CNe), x = lati), col = "purple2") +
  geom_ribbon(data = CNEW, aes(ymin = (CN_e5), ymax = (CN_e95), x= lati), alpha = 0.2, fill = "purple2") +
  xlab("Latitude") + ylab("Leaf nitrogen content (%)") +
  xlim (min(lati), max(lati)) + 
  ylim (-70,70) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 15), axis.title = element_text(size = 20)) +
  theme(legend.position=c(45, 20),legend.text = element_text(size = 15))  +
  scale_fill_manual( labels = c("Low force", "High force")) +
  scale_colour_discrete(labels=c("High forcing","Low forcing"), name = "") +
  theme(legend.title = element_blank()) +  annotate("text", x = 41, y = 70, label = "e)", cex = 10) 
  intNit  

sumSSD <- summary(mdlSSD6)$summary

a_spSSD = mean((sumSSD[grep("mu_grand_sp", rownames(sumSSD)), 1]))
a_spSSD5 <- quantile(postSSD$mu_grand_sp, c(0.05))
a_spSSD95 <- quantile(postSSD$mu_grand_sp, c(0.95))
a_spSSD25 <- quantile(postSSD$mu_grand_sp, c(0.25))
a_spSSD75 <- quantile(postSSD$mu_grand_sp, c(0.75))
a_spSSD <- cbind(a_spSSD, a_spSSD5,a_spSSD95, a_spSSD25,a_spSSD75)
a_spSSD <- a_spSSD

b_tranSSD <- sumSSD[grep("b_tranE", rownames(sumSSD)), 1]
b_tranSSD5 <- quantile(postSSD$b_tranE, c(0.05))
b_tranSSD95 <- quantile(postSSD$b_tranE, c(0.95))
b_tranSSD25 <- quantile(postSSD$b_tranE, c(0.25))
b_tranSSD75 <- quantile(postSSD$b_tranE, c(0.75))
b_tranSSD <- cbind(b_tranSSD, b_tranSSD5,b_tranSSD95, b_tranSSD25,b_tranSSD75)

b_tranlatSSD <- sumSSD[grep("b_tranlat", rownames(sumSSD)), 1]
b_tranlatSSD5 <- quantile(postSSD$b_tranlat, c(0.05))
b_tranlatSSD95 <- quantile(postSSD$b_tranlat, c(0.95))
b_tranlatSSD25 <- quantile(postSSD$b_tranlat, c(0.25))
b_tranlatSSD75 <- quantile(postSSD$b_tranlat, c(0.75))
b_tranlatSSD <- cbind(b_tranlatSSD, b_tranlatSSD5,b_tranlatSSD95, b_tranlatSSD25,b_tranlatSSD75)



## Simulate interaction with transect and latitude:

eData <- subset(trtPheno, transect == "1" )
wData <- subset(trtPheno, transect == "0" )

# Make the other parameters constant



# plot first for west coast
SSD_w = a_spSSD[1] + b_tranSSD[1] * tranW + b_tranlatSSD[1] * (tranW*lati)
SSD_e = a_spSSD[1] + b_tranSSD[1] * tranE + b_tranlatSSD[1] * (tranE*lati)

SSD_w5 = a_spSSD[2] + b_tranSSD[2] * tranW + b_tranlatSSD[2] * (tranW*lati)
SSD_e5 = a_spSSD[2] + b_tranSSD[2] * tranE + b_tranlatSSD[2] * (tranE*lati)

SSD_w95 = a_spSSD[3] + b_tranSSD[3] * tranW + b_tranlatSSD[3] * (tranW*lati)
SSD_e95 = a_spSSD[3] + b_tranSSD[3] * tranE + b_tranlatSSD[3] * (tranE*lati)

SSDEW <- data.frame(SSDw = SSD_w, SSDe = SSD_e,  SSD_w5 =SSD_w5, SSD_w95 = SSD_w95, SSD_e5 = SSD_e5, SSD_e95 = SSD_e95  )


intSSD <- ggplot(SSDEW) +
  geom_line(aes(y = (SSDw), x = lati), color = "maroon4", linetype = "dashed") +
  geom_ribbon(data = (SSDEW), aes(ymin = (SSD_w5), ymax = (SSD_w95), x= lati), alpha = 0.2, fill = "maroon4") +
  geom_line(aes(y = (SSDe), x = lati), col = "maroon") +
  geom_ribbon(data = (SSDEW), aes(ymin = (SSD_e5), ymax = (SSD_e95), x= lati), alpha = 0.2, fill = "maroon") +
  xlab("Latitude") +  labs(y = bquote('Wood specific density'~(g/cm^2)))  +
  xlim (min(lati), max(lati)) + 
  ylim (-7,7) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 15), axis.title = element_text(size = 20))+
  theme(legend.key=element_blank(), legend.position=c(.0,.0),legend.text = element_text(size = 15)) +
  #scale_fill_manual( labels = c("Low force", "High force")) +
  #scale_colour_discrete(labels=c("High forcing","Low forcing"), name = "") +
  theme(legend.title = element_blank()) +  annotate("text", x = 42, y = 7, label = "c)", cex = 10) 

# pdf("figures/intrxnPlotsHundoa.pdf", height =5, width = 5)
# intHt
# dev.off()
# 
# pdf("figures/intrxnPlotsHundob.pdf", height =5, width = 5)
# intDBH
# dev.off()
# 
# pdf("figures/intrxnPlotsHundoc.pdf", height =5, width = 5)
# intSSD
# dev.off()
# 
# pdf("figures/intrxnPlotsHundod.pdf", height =5, width = 5)
# intLMA
# dev.off()
# 
# pdf("figures/intrxnPlotsHundoedfd.pdf", height =5, width = 5)
# intNit
# dev.off()

pdf("analysis/figures/intrxnPlotsHundo.pdf", height =5, width = 25)
plot_grid( intHt, intDBH, intSSD, intLMA, intNit , ncol = 5, nrow =1,align = "v")
dev.off()


# #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# I think what we want is a loop that goes through each iteration of the posteriors and calculates the bb, but using 20 for forcing, 12 for photoperiod, 75 (75/10 when rescaled), and smithers to start
# 

#If we are using the old model, we will use the z-scored values for the parameters
# photo <- -0.5033863 #8 h photo
# siteSM <- 0
# force <- -0.3568628 #5/15 C trt
# chill <- -0.3546922 # low chill
# 
# m <- matrix(nrow = 1000, ncol = 47)
# 
# for(sp in 1:47){
#   for (it in 1:nrow(m)){
#     m[it,sp] <- post$a_sp[it,sp]+ post$b_site2[it] * siteSM + post$b_site3[it] * siteSM + post$b_site4[it] * siteSM + 
#       post$b_warm[it,sp] * force + post$b_photo[it, sp] * photo + post$b_chill[it,sp] * chill +
#       post$b_inter_wp[it,sp] * (force*photo) + post$b_inter_wc1[it,sp] * (force*chill) + post$b_inter_pc1[it,sp] * (photo*chill) +
#       post$b_inter_s2c1[it,sp] * (chill*siteSM) + post$b_inter_ws2[it,sp] * (force*siteSM) + post$b_inter_ps2[it,sp] * (photo*siteSM) +
#       post$b_inter_s3c1[it,sp] * (chill*siteSM) + post$b_inter_ws3[it,sp] * (force*siteSM) + post$b_inter_ps3[it,sp] * (photo*siteSM) +
#       post$b_inter_s4c1[it,sp] * (chill*siteSM) + post$b_inter_ws4[it,sp] * (force*siteSM) + post$b_inter_ps4[it,sp] * (photo*siteSM)
#   }
# }

############ SHRUB VS TREE ##############################
spInfo <- read.csv("analysis/input/species_ring.csv")
colnames(spInfo)[colnames(spInfo) == "X"] <- "ringType"

spInfo <- spInfo[, 1:5]
head(spInfo)

## Start with height:
fit <- rstan::extract(mdlHt)

chillB <- data.frame(fit$betaChillSp)

colnames(chillB) <- sort(spInfo$species.name)

longChill <- reshape2::melt(chillB)
colnames(longChill) <- c("species.name","betaCueSp")
longChill$cue <- "Chilling"

head(longChill)

###################################################################
# Photoperiod:
photoB <- data.frame(fit$betaPhotoSp)

colnames(photoB) <- sort(spInfo$species.name)

longPhoto <- reshape2::melt(photoB)
colnames(longPhoto) <- c("species.name","betaCueSp")
longPhoto$cue <- "Photoperiod"
head(longPhoto)

###################################################################
# Forcing:
forceB <- data.frame(fit$betaForceSp)

colnames(forceB) <- sort(spInfo$species.name)

longForce <- reshape2::melt(forceB)
colnames(longForce) <- c("species.name","betaCueSp")
longForce$cue <- "Forcing"

head(longForce)

longCues <- rbind(longForce, longChill, longPhoto)
longCues <- merge(longCues, spInfo, by = "species.name")

ring <- subset(longCues, ringType == "Ring ")
diffuse <- subset(longCues, ringType == "Diffuse")
difRing <- subset(longCues, ringType == "Diffuse/semi-ring")
semi <- subset(longCues, ringType == "Semi-ring" )

longCuesRing <- subset(longCues, ringType != "")

# ggplot() + 
#   geom_violin(data = longCuesRing, aes(x = as.factor(cue), y = betaCueSp, col = factor(cue))) +
#   facet_grid(col = vars(ringType), scales = "free_y") + theme(strip.background = element_blank(), strip.text.y = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.title=element_blank()) +
#    ylab ("Cue response") + xlab ("Cue")
# 
# 
# ggplot() + 
#   geom_violin(data = longCues, aes(x = as.factor(cue), y = betaCueSp, col = factor(cue))) +
#   facet_grid(col = vars(type), scales = "free_y") + theme(strip.background = element_blank(), strip.text.y = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.title=element_blank()) +
#   ylab ("Cue response") + xlab ("Cue")
# 
#   facet_grid(col = vars(type))
#               
#               +
#   stat_pointinterval(data = longest, aes(x = as.factor(cue), y = value, fill = factor(site, level = siteOrder)), .width = c(.5, .95) ,position = position_dodge(0.9)) +
#   theme_classic() +   
#   theme(legend.position = "right", 
#         legend.title = element_blank(),
#         axis.text.x = element_text( size= 16),
#         axis.text.y = element_text( size= 12),
#         axis.title=element_text(size = 14)) +
#   labs( x = "Treatment cue", y = "Cue response (days/standardized unit)", main = NA) +
#   scale_color_manual(values = c("Smithers" = "deepskyblue3",
#                                 "Manning Park" = "palegreen4", 
#                                 "St. Hippolyte"="darkorchid3", 
#                                 "Harvard Forest" = "tomato3"))+
#   scale_fill_manual(values = c("Smithers" = "deepskyblue3",
#                                "Manning Park" = "palegreen4", 
#                                "St. Hippolyte"="darkorchid3", 
#                                "Harvard Forest" = "tomato3"))
#   
#   
#   # compare points:
#   # chilling

longChill <- merge(longChill, spInfo, by = "species.name")
  
meanChill <- aggregate(longChill[c("betaCueSp")], longChill[c("ringType")], FUN = mean)
names(meanChill) <- c( "ringType", "betaCueSp")

meanError75 <- aggregate(longChill[c("betaCueSp")], longChill[c( "ringType")], FUN = function(i) quantile(i, probs = 0.75, na.rm = T))
names(meanError75) <- c( "ringType", "error75")

meanError25 <- aggregate(longChill[c("betaCueSp")], longChill[c( "ringType")], FUN = function(i) quantile(i, probs = 0.25, na.rm = T))
names(meanError25) <- c( "ringType", "error25")

meanError95 <- aggregate(longChill[c("betaCueSp")], longChill[c("ringType")], FUN = function(i) quantile(i, probs = 0.95, na.rm = T))
names(meanError95) <- c( "ringType", "error95")

meanError05 <- aggregate(longChill[c("betaCueSp")], longChill[c( "ringType")], FUN = function(i) quantile(i, probs = 0.05, na.rm = T))
names(meanError05) <- c( "ringType", "error05")


meanChill2 <- merge(meanChill, meanError05, by = c( "ringType"))
meanChill2 <- merge(meanChill2, meanError25, by = c("ringType"))
meanChill2 <- merge(meanChill2, meanError75, by = c("ringType"))
meanChill2 <- merge(meanChill2, meanError95, by = c("ringType"))

# require(dplyr)
# longCues %>% group_by(species.name) %>%
#   summarise(p90 = quantile(betaCueSp, probs=0.9, na.rm=TRUE))
meanChill2 <- subset(meanChill2, ringType !="")

ringChill <- ggplot(meanChill2,aes(y= betaCueSp, x = ringType)) +
  geom_point(size = 7,  color = "cyan4") +
  ylim (-30,3) +
  geom_errorbar(aes(ymin= error05, ymax = error95,xmin= ringType, xmax = ringType), 
                width= 0, linewidth = 0.5, color = "cyan4") +
  geom_errorbar(aes(ymin= error25, ymax = error75,xmin= ringType, xmax = ringType), 
                width= 0, linewidth = 1.5, color = "cyan4") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none")  + 
  labs( x = "Pore distribution type", y = "Chilling response (days/m)", main = NA) +
  theme(legend.title = element_blank(), axis.text.x = element_text(size = 20, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 18), axis.title = element_text(size = 20)) +  
  annotate("text", x = 0.75, y = 3, label = "a)", cex = 10) 
ringChill
# theme(axis.text.x = element_text( size=17,angle = 78,  hjust=1),
#       axis.text.y=element_text(size = 15),
#       axis.title=element_text(size=  17),
#       legend.position = "none") +
#   labs( x = "Ring Type", y = "Chilling response (days/standardized unit)", main = NA) +
#   theme(legend.title = element_blank()) 
  
longforce <- merge(longForce, spInfo, by = "species.name")

meanforce <- aggregate(longforce[c("betaCueSp")], longforce[c("ringType")], FUN = mean)
names(meanforce) <- c( "ringType", "betaCueSp")

meanError75 <- aggregate(longforce[c("betaCueSp")], longforce[c( "ringType")], FUN = function(i) quantile(i, probs = 0.75, na.rm = T))
names(meanError75) <- c( "ringType", "error75")

meanError25 <- aggregate(longforce[c("betaCueSp")], longforce[c( "ringType")], FUN = function(i) quantile(i, probs = 0.25, na.rm = T))
names(meanError25) <- c( "ringType", "error25")

meanError95 <- aggregate(longforce[c("betaCueSp")], longforce[c("ringType")], FUN = function(i) quantile(i, probs = 0.95, na.rm = T))
names(meanError95) <- c( "ringType", "error95")

meanError05 <- aggregate(longforce[c("betaCueSp")], longforce[c( "ringType")], FUN = function(i) quantile(i, probs = 0.05, na.rm = T))
names(meanError05) <- c( "ringType", "error05")


meanforce2 <- merge(meanforce, meanError05, by = c( "ringType"))
meanforce2 <- merge(meanforce2, meanError25, by = c("ringType"))
meanforce2 <- merge(meanforce2, meanError75, by = c("ringType"))
meanforce2 <- merge(meanforce2, meanError95, by = c("ringType"))

meanforce2 <- subset(meanforce2, ringType !="")

ringForce <- ggplot(meanforce2,aes(y= betaCueSp, x = ringType), size = 7) +
  geom_point(size = 7, color = "goldenrod") +
  ylim (-30,3) +
  geom_errorbar(aes(ymin= error05, ymax = error95,xmin= ringType, xmax = ringType), width= 0, linewidth = 0.5, color = "goldenrod") +
  geom_errorbar(aes(ymin= error25, ymax = error75,xmin= ringType, xmax = ringType), width= 0, linewidth = 1.5, color = "goldenrod") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +  
  annotate("text", x = 0.75, y = 3, label = "b)", cex = 10) +
  labs( x = "Pore distribution type", y = "Forcing response (days/m)", main = NA) +
  theme(legend.title = element_blank(), axis.text.x = element_text(size = 20, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 18), axis.title = element_text(size = 20)) 

longphoto <- merge(longPhoto, spInfo, by = "species.name")

meanphoto <- aggregate(longphoto[c("betaCueSp")], longphoto[c("ringType")], FUN = mean)
names(meanphoto) <- c( "ringType", "betaCueSp")

meanError75 <- aggregate(longphoto[c("betaCueSp")], longphoto[c( "ringType")], FUN = function(i) quantile(i, probs = 0.75, na.rm = T))
names(meanError75) <- c( "ringType", "error75")

meanError25 <- aggregate(longphoto[c("betaCueSp")], longphoto[c( "ringType")], FUN = function(i) quantile(i, probs = 0.25, na.rm = T))
names(meanError25) <- c( "ringType", "error25")

meanError95 <- aggregate(longphoto[c("betaCueSp")], longphoto[c("ringType")], FUN = function(i) quantile(i, probs = 0.95, na.rm = T))
names(meanError95) <- c( "ringType", "error95")

meanError05 <- aggregate(longphoto[c("betaCueSp")], longphoto[c( "ringType")], FUN = function(i) quantile(i, probs = 0.05, na.rm = T))
names(meanError05) <- c( "ringType", "error05")


meanphoto2 <- merge(meanphoto, meanError05, by = c( "ringType"))
meanphoto2 <- merge(meanphoto2, meanError25, by = c("ringType"))
meanphoto2 <- merge(meanphoto2, meanError75, by = c("ringType"))
meanphoto2 <- merge(meanphoto2, meanError95, by = c("ringType"))

# require(dplyr)
# longCues %>% group_by(species.name) %>%
#   summarise(p90 = quantile(betaCueSp, probs=0.9, na.rm=TRUE))
meanphoto2 <- subset(meanphoto2, ringType !="")

ringPhoto <- ggplot(meanphoto2,aes(y= betaCueSp, x = ringType), size = 7) +
  geom_point(size = 7, color = "maroon") +
  ylim (-30,3) +
  geom_errorbar(aes(ymin= error05, ymax = error95,xmin= ringType, xmax = ringType), width= 0, linewidth = 0.5, color = "maroon") +
  geom_errorbar(aes(ymin= error25, ymax = error75,xmin= ringType, xmax = ringType), width= 0, linewidth = 1.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  annotate("text", x = 0.75, y = 3, label = "c)", cex = 10) +
  labs( x = "Pore distribution type", y = "Photoperiod response (days/m)", main = NA) +
  theme(legend.title = element_blank(), axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18), axis.title = element_text(size = 20)) 

# pdf("figures/ringPorosityHeightHundoa.pdf", width = 6, height = 8)
# ringChill
# dev.off()
# 
# pdf("figures/ringPorosityHeightHundob.pdf", width = 6, height = 8)
# ringForce
# dev.off()
# 
# pdf("figures/ringPorosityHeightHundoc.pdf", width = 6, height = 8)
# ringPhoto
# dev.off()


pdf("analysis/figures/ringPorosityHeightHundo.pdf", width = 12, height = 6)
plot_grid(ringChill, ringForce, ringPhoto, nrow = 1, ncol = 3, align = "v")
dev.off()


## Old code
# OLD CODE

sumLMA <- summary(mdlLMA)$summary
muGrand = (sumLMA[grep("mu_grand", rownames(sumLMA)), 1])
b_trtSpLMA = (sumLMA[grep("b_muSp", rownames(sumLMA)), 1])
a_trtSpLMA = mean((sumLMA[grep("mu_grand", rownames(sumLMA)), 1]))
b_tranELMA = sumLMA[grep("b_tranE", rownames(sumLMA)), 1]
b_tranlatLMA = sumLMA[grep("b_tranlat", rownames(sumLMA)), 1]

b_phenoSpLMA = (sumLMA[grep("alphaPhenoSp", rownames(sumLMA)), 1])
a_phenoSpLMA = (sumLMA[grep("muPhenoSp", rownames(sumLMA)), 1])

a_chillSpLMA = sumLMA[grep("alphaChillSp", rownames(sumLMA)), 1]
a_forceSpLMA = sumLMA[grep("alphaForceSp", rownames(sumLMA)), 1]
a_photoSpLMA = sumLMA[grep("alphaPhotoSp", rownames(sumLMA)), 1]

b_photoSpLMA = sumLMA[grep("muPhotoSp", rownames(sumLMA)), 1]
b_forceSpLMA = sumLMA[grep("muForceSp", rownames(sumLMA)), 1]
b_chillSpLMA = sumLMA[grep("muChillSp", rownames(sumLMA)), 1]

bTrtChillLMA = sumLMA[grep("betaTraitxChill", rownames(sumLMA)), 1]
bTrtForceLMA = sumLMA[grep("betaTraitxForce", rownames(sumLMA)), 1]
bTrtPhotoLMA = sumLMA[grep("betaTraitxPhoto", rownames(sumLMA)), 1]

# a_trtsp5LMA <- vector()
# for(i in 1:ncol(postLMA$muSp)){  
#   quantU <- quantile(postLMA$mu_grand[,i], c(0.05, 0.95, 0.25, 0.75))
#   a_trtsp5LMA <- rbind(a_trtsp5LMA, quantU)
# }
# colnames(a_trtsp5LMA) <- c("Int5","Int95","Int25","Int75")

a_trtsp5LMA <- quantile(postLMA$mu_grand, c(0.05, 0.95, 0.25, 0.75))
b_tran5LMA <- quantile(postLMA$b_tranE, c(0.05, 0.95, 0.25, 0.75))
b_tranlat5LMA <- quantile(postLMA$b_tranlat, c(0.05, 0.95, 0.25, 0.75))

# b_chill5 <- vector()
# for(i in 1:ncol(postLMA$betaChillSp)){
#   quantU <- round(quantile(postLMA$betaChillSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
#   b_chill5 <- rbind(b_chill5, quantU)
# }
# colnames(b_chill5) <- c("chill5","chill95","chill25","chill75")
# 
# b_force5 <- vector()
# for(i in 1:ncol(postLMA$betaForceSp)){
#   quantU <- round(quantile(postLMA$betaForceSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
#   b_force5 <- rbind(b_force5, quantU)
# }
# colnames(b_force5) <- c("force5","force95","force25","force75")
# 
# b_photo5 <- vector()
# for(i in 1:ncol(postLMA$betaPhotoSp)){
#   quantU <- round(quantile(postLMA$betaPhotoSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
#   b_photo5 <- rbind(b_photo5, quantU)
# }
# colnames(b_photo5) <- c("photo5","photo95","photo25","photo75")

## Simulate interaction with transect and latitude:

eData <- subset(trtPheno, transect == "1" )
wData <- subset(trtPheno, transect == "0" )

# Make the other parameters constant

# plot first for west coast
LMA_w = a_trtSpLMA + b_tranELMA * tranW + b_tranlatLMA * (tranW*lati)
LMA_e = a_trtSpLMA + b_tranELMA * tranE + b_tranlatLMA * (tranE*lati)

LMA_w5 = mean(a_trtsp5LMA[1]) + b_tran5LMA[1] * tranW + b_tranlat5LMA[1] * (tranW*lati)
LMA_e5 = mean(a_trtsp5LMA[1]) + b_tran5LMA[1] * tranE + b_tranlat5LMA[1] * (tranE*lati)

LMA_w95 = mean(a_trtsp5LMA[2]) + b_tran5LMA[2] * tranW + b_tranlat5LMA[2] * (tranW*lati)
LMA_e95 = mean(a_trtsp5LMA[2]) + b_tran5LMA[2] * tranE + b_tranlat5LMA[2] * (tranE*lati)

intLMA <- data.frame(LMAw = c(LMA_w ),LMAw5 = c(LMA_w5),LMAw95 = c(LMA_w95), LMAe = c(LMA_e), LMAe5 = c(LMA_e5), LMAe95 = c(LMA_e95))

ggplot(intLMA) +
  geom_line(aes(y = LMAw, x = lati, col = "cyan4")) +
  geom_ribbon( aes(ymin = LMAw5, ymax = LMAw95, x= lati), alpha = 0.2, fill = "cyan4") +
  geom_line(aes(y = LMAe, x = lati, col = "cyan3")) + 
  geom_ribbon(data = LAEW, aes(ymin = ht_e5, ymax = ht_e95, x= lati), alpha = 0.2, fill = "cyan3") + 
  xlab("Standardized latitude") + ylab("Standardized height") +
  xlim (-0.8,0.8) + 
  ylim (-0.5,0.5) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 15), axis.title = element_text(size = 20))+
  theme(legend.key=element_blank(), legend.position=c(.8,.85),legend.text = element_text(size = 15)) +
  #scale_fill_manual( labels = c("Low force", "High force")) +
  scale_color_manual(values = c("cyan3","cyan4"), labels = c("Eastern", "Western"), name = "") +
  #scale_colour_discrete(labels=c("High forcing","Low forcing"), name = "") +
  theme(legend.title = element_blank())# +  annotate("text", x = -5.4, y = 125, label = "a)", cex = 10) 


# par(mfrow = c(1,1))
# plot(0, type = "n", xlim = c(25,60), ylim = c(-1,1),
#      xlab = "Latitude",
#      ylab = "Trait")
# abline(lm(LMA_w ~ lati), col = "darkslategray", lwd = 3, lty = 2)
# abline(lm(LMA_e ~lati), col = "darkslategray", lwd = 3, lty =1)
#############################################
##dbh

sumDBH <- summary(mdlDBH)$summary
muGrand = (sumDBH[grep("mu_grand", rownames(sumDBH)), 1])
b_trtSpDBH = (sumDBH[grep("b_muSp", rownames(sumDBH)), 1])
a_trtSpDBH = mean((sumDBH[grep("mu_grand_sp", rownames(sumDBH)), 1]))
b_tranEDBH = sumDBH[grep("b_tranE", rownames(sumDBH)), 1]
b_tranlatDBH = sumDBH[grep("b_tranlat", rownames(sumDBH)), 1]

b_phenoSpDBH = (sumDBH[grep("alphaPhenoSp", rownames(sumDBH)), 1])
a_phenoSpDBH = (sumDBH[grep("muPhenoSp", rownames(sumDBH)), 1])

a_chillSpDBH = sumDBH[grep("alphaChillSp", rownames(sumDBH)), 1]
a_forceSpDBH = sumDBH[grep("alphaForceSp", rownames(sumDBH)), 1]
a_photoSpDBH = sumDBH[grep("alphaPhotoSp", rownames(sumDBH)), 1]

b_photoSpDBH = sumDBH[grep("muPhotoSp", rownames(sumDBH)), 1]
b_forceSpDBH = sumDBH[grep("muForceSp", rownames(sumDBH)), 1]
b_chillSpDBH = sumDBH[grep("muChillSp", rownames(sumDBH)), 1]

bTrtChillDBH = sumDBH[grep("betaTraitxChill", rownames(sumDBH)), 1]
bTrtForceDBH = sumDBH[grep("betaTraitxForce", rownames(sumDBH)), 1]
bTrtPhotoDBH = sumDBH[grep("betaTraitxPhoto", rownames(sumDBH)), 1]

a_trtsp5DBH <- vector()
for(i in 1:ncol(postDBH$muSp)){  
  quantU <- round(quantile(postDBH$muSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  a_trtsp5DBH <- rbind(a_trtsp5DBH, quantU)
}
colnames(a_trtsp5DBH) <- c("Int5","Int95","Int25","Int75")

b_tran5DBH <- round(quantile(postDBH$b_tranE, c(0.05, 0.95, 0.25, 0.75)),1)
b_tranlat5DBH <- round(quantile(postDBH$b_tranlat, c(0.05, 0.95, 0.25, 0.75)),1)

b_chill5 <- vector()
for(i in 1:ncol(postDBH$betaChillSp)){
  quantU <- round(quantile(postDBH$betaChillSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_chill5 <- rbind(b_chill5, quantU)
}
colnames(b_chill5) <- c("chill5","chill95","chill25","chill75")

b_force5 <- vector()
for(i in 1:ncol(postDBH$betaForceSp)){
  quantU <- round(quantile(postDBH$betaForceSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_force5 <- rbind(b_force5, quantU)
}
colnames(b_force5) <- c("force5","force95","force25","force75")

b_photo5 <- vector()
for(i in 1:ncol(postDBH$betaPhotoSp)){
  quantU <- round(quantile(postDBH$betaPhotoSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_photo5 <- rbind(b_photo5, quantU)
}
colnames(b_photo5) <- c("photo5","photo95","photo25","photo75")

## Simulate interaction with transect and latitude:

eData <- subset(trtPheno, transect == "1" )
wData <- subset(trtPheno, transect == "0" )


# plot first for west coast
DBH_w = a_trtSpDBH + b_tranEDBH * tranW + b_tranlatDBH * (tranW*lati)
DBH_e = a_trtSpDBH + b_tranEDBH * tranE + b_tranlatDBH * (tranE*lati)

par(mfrow = c(1,1))
plot(0, type = "n", xlim = c(25,60), ylim = c(-10,10),
  xlab = "Latitude",
  ylab = "Trait")
abline(lm(DBH_w ~ lati), col = "darkslategray", lwd = 3, lty = 2)
abline(lm(DBH_e ~lati), col = "darkslategray", lwd = 3, lty =1)


#############################################
sumSSD <- summary(mdlSSD)$summary
muGrand = (sumSSD[grep("mu_grand", rownames(sumSSD)), 1])
b_trtSpSSD = (sumSSD[grep("b_muSp", rownames(sumSSD)), 1])
a_trtSpSSD = mean((sumSSD[grep("mu_grand_sp", rownames(sumSSD)), 1]))
b_tranESSD = sumSSD[grep("b_tranE", rownames(sumSSD)), 1]
b_tranlatSSD = sumSSD[grep("b_tranlat", rownames(sumSSD)), 1]

b_phenoSpSSD = (sumSSD[grep("alphaPhenoSp", rownames(sumSSD)), 1])
a_phenoSpSSD = (sumSSD[grep("muPhenoSp", rownames(sumSSD)), 1])

a_chillSpSSD = sumSSD[grep("alphaChillSp", rownames(sumSSD)), 1]
a_forceSpSSD = sumSSD[grep("alphaForceSp", rownames(sumSSD)), 1]
a_photoSpSSD = sumSSD[grep("alphaPhotoSp", rownames(sumSSD)), 1]

b_photoSpSSD = sumSSD[grep("muPhotoSp", rownames(sumSSD)), 1]
b_forceSpSSD = sumSSD[grep("muForceSp", rownames(sumSSD)), 1]
b_chillSpSSD = sumSSD[grep("muChillSp", rownames(sumSSD)), 1]

bTrtChillSSD = sumSSD[grep("betaTraitxChill", rownames(sumSSD)), 1]
bTrtForceSSD = sumSSD[grep("betaTraitxForce", rownames(sumSSD)), 1]
bTrtPhotoSSD = sumSSD[grep("betaTraitxPhoto", rownames(sumSSD)), 1]

a_trtsp5SSD <- vector()
for(i in 1:ncol(postSSD$muSp)){  
  quantU <- round(quantile(postSSD$muSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  a_trtsp5SSD <- rbind(a_trtsp5SSD, quantU)
}
colnames(a_trtsp5SSD) <- c("Int5","Int95","Int25","Int75")

b_tran5SSD <- round(quantile(postSSD$b_tranE, c(0.05, 0.95, 0.25, 0.75)),1)
b_tranlat5SSD <- round(quantile(postSSD$b_tranlat, c(0.05, 0.95, 0.25, 0.75)),1)

b_chill5 <- vector()
for(i in 1:ncol(postSSD$betaChillSp)){
  quantU <- round(quantile(postSSD$betaChillSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_chill5 <- rbind(b_chill5, quantU)
}
colnames(b_chill5) <- c("chill5","chill95","chill25","chill75")

b_force5 <- vector()
for(i in 1:ncol(postSSD$betaForceSp)){
  quantU <- round(quantile(postSSD$betaForceSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_force5 <- rbind(b_force5, quantU)
}
colnames(b_force5) <- c("force5","force95","force25","force75")

b_photo5 <- vector()
for(i in 1:ncol(postSSD$betaPhotoSp)){
  quantU <- round(quantile(postSSD$betaPhotoSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_photo5 <- rbind(b_photo5, quantU)
}
colnames(b_photo5) <- c("photo5","photo95","photo25","photo75")

## Simulate interaction with transect and latitude:

eData <- subset(trtPheno, transect == "1" )
wData <- subset(trtPheno, transect == "0" )


# plot first for west coast
SSD_w = a_trtSpSSD + b_tranESSD * tranW + b_tranlatSSD * (tranW*lati)
SSD_e = a_trtSpSSD + b_tranESSD * tranE + b_tranlatSSD * (tranE*lati)

par(mfrow = c(1,1))
plot(0, type = "n", xlim = c(25,60), ylim = c(-100,100),
  xlab = "Latitude",
  ylab = "Trait")
abline(lm(SSD_w ~ lati), col = "darkslategray", lwd = 3, lty = 2)
abline(lm(SSD_e ~lati), col = "darkslategray", lwd = 3, lty =1)

#############################################
sumCN <- summary(mdlPerN)$summary
muGrand = (sumCN[grep("mu_grand", rownames(sumCN)), 1])
b_trtSpCN = (sumCN[grep("b_muSp", rownames(sumCN)), 1])
a_trtSpCN = mean((sumCN[grep("mu_grand_sp", rownames(sumCN)), 1]))
b_tranECN = sumCN[grep("b_tranE", rownames(sumCN)), 1]
b_tranlatCN = sumCN[grep("b_tranlat", rownames(sumCN)), 1]

b_phenoSpCN = (sumCN[grep("alphaPhenoSp", rownames(sumCN)), 1])
a_phenoSpCN = (sumCN[grep("muPhenoSp", rownames(sumCN)), 1])

a_chillSpCN = sumCN[grep("alphaChillSp", rownames(sumCN)), 1]
a_forceSpCN = sumCN[grep("alphaForceSp", rownames(sumCN)), 1]
a_photoSpCN = sumCN[grep("alphaPhotoSp", rownames(sumCN)), 1]

b_photoSpCN = sumCN[grep("muPhotoSp", rownames(sumCN)), 1]
b_forceSpCN = sumCN[grep("muForceSp", rownames(sumCN)), 1]
b_chillSpCN = sumCN[grep("muChillSp", rownames(sumCN)), 1]

bTrtChillCN = sumCN[grep("betaTraitxChill", rownames(sumCN)), 1]
bTrtForceCN = sumCN[grep("betaTraitxForce", rownames(sumCN)), 1]
bTrtPhotoCN = sumCN[grep("betaTraitxPhoto", rownames(sumCN)), 1]

a_trtsp5CN <- vector()
for(i in 1:ncol(postCN$muSp)){  
  quantU <- round(quantile(postCN$muSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  a_trtsp5CN <- rbind(a_trtsp5CN, quantU)
}
colnames(a_trtsp5CN) <- c("Int5","Int95","Int25","Int75")

b_tran5CN <- round(quantile(postCN$b_tranE, c(0.05, 0.95, 0.25, 0.75)),1)
b_tranlat5CN <- round(quantile(postCN$b_tranlat, c(0.05, 0.95, 0.25, 0.75)),1)

b_chill5 <- vector()
for(i in 1:ncol(postCN$betaChillSp)){
  quantU <- round(quantile(postCN$betaChillSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_chill5 <- rbind(b_chill5, quantU)
}
colnames(b_chill5) <- c("chill5","chill95","chill25","chill75")

b_force5 <- vector()
for(i in 1:ncol(postCN$betaForceSp)){
  quantU <- round(quantile(postCN$betaForceSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_force5 <- rbind(b_force5, quantU)
}
colnames(b_force5) <- c("force5","force95","force25","force75")

b_photo5 <- vector()
for(i in 1:ncol(postCN$betaPhotoSp)){
  quantU <- round(quantile(postCN$betaPhotoSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_photo5 <- rbind(b_photo5, quantU)
}
colnames(b_photo5) <- c("photo5","photo95","photo25","photo75")

## Simulate interaction with transect and latitude:

eData <- subset(trtPheno, transect == "1" )
wData <- subset(trtPheno, transect == "0" )


# plot first for west coast
CN_w = a_trtSpCN + b_tranECN * tranW + b_tranlatCN * (tranW*lati)
CN_e = a_trtSpCN + b_tranECN * tranE + b_tranlatCN * (tranE*lati)

par(mfrow = c(1,1))
plot(0, type = "n", xlim = c(25,60), ylim = c(-10,10),
  xlab = "Latitude",
  ylab = "Trait")
abline(lm(CN_w ~ lati), col = "darkslategray", lwd = 3, lty = 2)
abline(lm(CN_e ~lati), col = "darkslategray", lwd = 3, lty =1)
#################################
pdf("figures/transectIntrxnZ25.pdf", width =13, height =3)
par(mfrow = c(1,5))
plot(0, type = "n", xlim = c(40,55), ylim = c(-2,2),
  xlab = "Latitude",
  ylab = "Height", cex.lab = 1.3)

#abline(lm(ht_w ~ lati), col = "darkslategray", lwd = 3, lty = 2)
abline(lm(ht_e ~lati), col = "darkslategray4", lwd = 3, lty =1)
abline(lm(ht_w ~lati), col = "darkslategray4", lwd = 3, lty =2)
text(40.8, 2, label = "a)", cex = 1.25)

plot(0, type = "n", xlim = c(40,55), ylim = c(-2,2),
  xlab = "Latitude",
  ylab = "Diameter at breast height", cex.lab = 1.3)
abline(lm(DBH_e ~lati), col = "goldenrod", lwd = 3, lty =1)
abline(lm(DBH_w ~lati), col = "goldenrod", lwd = 3, lty =2)
text(40.8, 2, label = "b)", cex = 1.25)

plot(0, type = "n", xlim = c(40,55), ylim = c(-2,2),
  xlab = "Latitude",
  ylab = "Leaf mass area", cex.lab = 1.3)
abline(lm(LMA_e ~lati), col = "darkolivegreen", lwd = 3, lty =1)
abline(lm(LMA_w ~lati), col = "darkolivegreen", lwd = 3, lty =2)
text(40.8,2, label = "c)", cex = 1.25)


plot(0, type = "n", xlim = c(40,55), ylim = c(-2,2),
  xlab = "Latitude",
  ylab = "Stem specific density", cex.lab = 1.3)
abline(lm(SSD_e ~lati), col = "maroon", lwd = 3, lty =1)
abline(lm(SSD_w ~lati), col = "maroon", lwd = 3, lty =2)
text(40.8,2, label = "d)", cex = 1.25)

plot(0, type = "n", xlim = c(40,55), ylim = c(-2,2),
  xlab = "Latitude",
  ylab = "Carbon:Nitrogen", cex.lab = 1.3)
abline(lm(CN_e ~lati), col = "purple4", lwd = 3, lty =1)
abline(lm(CN_w ~lati), col = "purple4", lwd = 3, lty =2)
text(40.8,2, label = "e)", cex = 1.25)

legend("topright",legend = c(expression("Western"),
  expression("Eastern")
),
  col = c("black", "black"),
  lty = c(2,1), lwd = 1, cex= 1.25, bty = "n")

dev.off()