# aim of this code is to generate the values referenced in the sweave file:

library(rstan)
library(tidybayes)
library(reshape2)


########## General values ###############################
spInfo <- read.csv("..//analysis/input/species_ring.csv")

perTreeRing <- round((nrow(subset(spInfo, !is.na(ring.type)))/47)*100,1)
nRing <- nrow(subset(spInfo, !is.na(ring.type)))
spInfoE <- subset(spInfo, transect != "west")
spInfoW <- subset(spInfo, transect != "east")

spInfoEShrub <- subset(spInfoE, type == "shrub")
spInfoWShrub <- subset(spInfoW, type == "shrub")

spInfoETree <- subset(spInfoE, type == "tree")
spInfoWTree <- subset(spInfoW, type == "tree")

nSpp <- length(unique(spInfo$species.name))
nSppE <- length(unique(spInfoE$species.name))
nSppW <- length(unique(spInfoW$species.name))
nSppES <- length(unique(spInfoEShrub$species.name))
nSppWS <- length(unique(spInfoWShrub$species.name))
nSppET <- length(unique(spInfoETree$species.name))
nSppWT <- length(unique(spInfoWTree$species.name))


pheno <- read.csv("..//analysis/input/phenoDataWChill.csv")
# trtPheno <- read.csv("input/trtData.csv")
trtPheno <- read.csv("..//analysis/input/trtPhenoDummy.csv")

height <- trtPheno[complete.cases(trtPheno$ht),]
nHt <- nrow(height)

leafMA <- trtPheno[complete.cases(trtPheno$lma),]
nLMA <- nrow(leafMA)

diam <- trtPheno[complete.cases(trtPheno$dbh),] 
nDBH <- nrow(diam)

stem <- trtPheno[complete.cases(trtPheno$ssd),]
nSSD <- nrow(stem)

carbNit <- trtPheno[complete.cases(trtPheno$C.N),]
nCNit <- nrow(carbNit)

nit <- trtPheno[complete.cases(trtPheno$per.N),]
nLNC <- nrow(nit)

temp <- trtPheno[,c("sample","species","site")]
temp$count <- 1
noObsTrait <- aggregate(temp["count"], trtPheno[,c("species","site")], FUN = sum)

noObsTrtMin <- min(noObsTrait$count)
noObsTrtMin <- min(noObsTrait$count)

eTrt <- subset(trtPheno, transect == "1")
wTrt <- subset(trtPheno, transect == "0")
noIndivE <- length(unique(eTrt$sample))
noIndivW<- length(unique(wTrt$sample))
noIndiv <- length(unique(trtPheno$sample))

pheno$count <- 1
noObsPheno <- aggregate(pheno["count"], pheno[,c("species","population")], FUN = sum)
noIndivPheno <- nrow(pheno)


# col4fig <- c("mean","sd","25%","50%","75%","Rhat")
# col4table <- c("mean","sd","2.5%","50%","97.5%","Rhat")
# 
# mu_params_4 <- c( 
#   #"a_z",
#   # "lam_interceptsa",
#   "mu_b_warm",
#   "mu_b_photo",
#   "mu_b_chill1",
#   "b_site2",
#   "b_site3",
#   "b_site4",
#   "mu_b_inter_wp",
#   "mu_b_inter_wc1",
#   "mu_b_inter_pc1",
#   "mu_b_inter_ws2",
#   "mu_b_inter_ps2",
#   "mu_b_inter_s2c1",
#   "mu_b_inter_ws3",
#   "mu_b_inter_ps3",
#   "mu_b_inter_s3c1",
#   "mu_b_inter_ws4",
#   "mu_b_inter_ps4",
#   "mu_b_inter_s4c1")
# 
# meanz4 <- sum[mu_params_4, col4table]
# 
# rownames(meanz4) = c( 
#   #"Root trait intercept", "Lambda",
#   "Forcing",
#   "Photoperiod",
#   "Chilling",
#   "Manning Park",
#   "Harvard Forest",
#   "St. Hippolyte",
#   "Forcing x photoperiod",
#   "Forcing x chilling",
#   "Photoperiod x chilling",
#   "Forcing x Manning Park",
#   "Photoperiod x Manning Park",
#   "Chilling x Manning Park",
#   "Forcing x Harvard Forest",
#   "Photoperiod x Harvard Forest",
#   "Chilling x Harvard Forest",
#   "Forcing x St. Hippolyte",
#   "Photoperiod x St. Hippolyte",
#   "Chilling x St. Hippolyte"            
# )
# 
# meanz4.table <- sum[mu_params_4, col4table]
# row.names(meanz4.table) <- row.names(meanz4)
# head(meanz4.table)
# #write.table(meanzew.table , "output/term.mdl.esti.dldf.csv", sep = ",", row.names = FALSE)
# 
# # call table values:
# 
# chillCue <- round(meanz4["Chilling",1],1)
# chillCueU <- round(quantile(fit$mu_b_chill1, c(0.95)),1)
# chillCueL <- round(quantile(fit$mu_b_chill1, c(0.05)),1)


load("..//analysis/output/htContLatHundowLat.Rdata")

htModelFit <- rstan::extract(mdlHt)
muSp <- data.frame(htModelFit$mu_grand_sp)
muSpMean <- colMeans(muSp)
htFold <- -1*round(max(muSpMean)/min(muSpMean),0)

htLat<- as.numeric(round(mean(htModelFit$b_lat),1))


htTran <- as.numeric(round(mean(htModelFit$b_tranE),1))
lower_htTran <- format(round(quantile(htModelFit$b_tranE, prob = 0.05),1), nsmall =1)
upper_htTran <- round(quantile(htModelFit$b_tranE, prob = 0.95),1)

htLatTran <- as.numeric(round(mean(htModelFit$b_tranlat),1))
lower_htLatTran <- format(round(quantile(htModelFit$b_tranlat, prob = 0.05),1), nsmall =1)
upper_htLatTran <- format(round(quantile(htModelFit$b_tranlat, prob = 0.95),1), nsmall =1)

htBFSpMean <- as.numeric(round(mean(htModelFit$betaTraitxForce),1))
lower_htBFSpMean <- format(round(quantile(htModelFit$betaTraitxForce, prob = 0.05),1), nsmall =1)
upper_htBFSpMean <- round(quantile(htModelFit$betaTraitxForce, prob = 0.95),1)
htBFSpMean; lower_htBFSpMean; upper_htBFSpMean

htBCSpMean <- as.numeric(round(mean(htModelFit$betaTraitxChill),1))
lower_htBCSpMean <- round(quantile(htModelFit$betaTraitxChill, prob = 0.05),1)
upper_htBCSpMean <- round(quantile(htModelFit$betaTraitxChill, prob = 0.95),1)
htBCSpMean; lower_htBCSpMean; upper_htBCSpMean

htBPSpMean <- as.numeric(round(mean(htModelFit$betaTraitxPhoto),1))
lower_htBPSpMean <- round(quantile(htModelFit$betaTraitxPhoto, prob = 0.05),1)
upper_htBPSpMean <- round(quantile(htModelFit$betaTraitxPhoto, prob = 0.95),1)
htBPSpMean; lower_htBPSpMean; upper_htBPSpMean

# muForceSp
htMuFSpMean <- as.numeric(round(mean(htModelFit$muForceSp),1))
lower_htMuFSpMean <- format(round(quantile(htModelFit$muForceSp, prob = 0.05),1), nsmall =1)
upper_htMuFSpMean <- round(quantile(htModelFit$muForceSp, prob = 0.95),1)
htMuFSpMean; lower_htMuFSpMean; upper_htMuFSpMean

#muChillSP
htMuCSpMean <- as.numeric(round(mean(htModelFit$muChillSp),1))
lower_htMuCSpMean <- format(round(quantile(htModelFit$muChillSp, prob = 0.05),1), nsmall =1)
upper_htMuCSpMean <- round(quantile(htModelFit$muChillSp, prob = 0.95),1)
htMuCSpMean; lower_htMuCSpMean; upper_htMuCSpMean

#muPhotoSp
htMuPSpMean <- as.numeric(round(mean(htModelFit$muPhotoSp),1))
lower_htMuPSpMean <- format(round(quantile(htModelFit$muPhotoSp, prob = 0.05),1), nsmall =1)
upper_htMuPSpMean <- round(quantile(htModelFit$muPhotoSp, prob = 0.95),1)
htMuPSpMean; lower_htMuPSpMean; upper_htMuPSpMean

htchill <- data.frame(htModelFit$alphaChillSp)
htChillSpMean <- colMeans(htchill)
htChillMax <- max(htChillSpMean)

htsigmaSp <- as.numeric(round(mean(htModelFit$sigma_sp),1))
lower_htsigmaSp <- format(round(quantile(htModelFit$sigma_sp, prob = 0.05),1), nsmall =1)
upper_htsigmaSp <- round(quantile(htModelFit$sigma_sp, prob = 0.95),1)
######## Leaf mass area###############################

load("..//analysis/output/lmaContLatHundowLat.Rdata")
lmaModelFit <- rstan::extract(mdlLMA)
muSp <- data.frame(lmaModelFit$mu_grand_sp)
muSpMean <- colMeans(muSp)/100
lmaFold <- round(max(muSpMean)/min(muSpMean),0)

lmaBFSpMean <- as.numeric(round(mean(lmaModelFit$betaTraitxForce),1))
lower_lmaBFSpMean <- format(round(quantile(lmaModelFit$betaTraitxForce, prob = 0.05),1), nsmall =1)
upper_lmaBFSpMean <- round(quantile(lmaModelFit$betaTraitxForce, prob = 0.95),1)
lmaBFSpMean; lower_lmaBFSpMean; upper_lmaBFSpMean

lmaBCSpMean <- as.numeric(round(mean(lmaModelFit$betaTraitxChill),1))
lower_lmaBCSpMean <- round(quantile(lmaModelFit$betaTraitxChill, prob = 0.05),1)
upper_lmaBCSpMean <- round(quantile(lmaModelFit$betaTraitxChill, prob = 0.95),1)
lmaBCSpMean; lower_lmaBCSpMean; upper_lmaBCSpMean

lmaBPSpMean <- as.numeric(round(mean(lmaModelFit$betaTraitxPhoto),1))
lower_lmaBPSpMean <- round(quantile(lmaModelFit$betaTraitxPhoto, prob = 0.05),1)
upper_lmaBPSpMean <- format(round(quantile(lmaModelFit$betaTraitxPhoto, prob = 0.95),1), nsmall =1)
lmaBPSpMean; lower_lmaBPSpMean; upper_lmaBPSpMean

lmachill <- data.frame(lmaModelFit$alphaChillSp)
lmaChillSpMean <- colMeans(lmachill)
lmaChillMax <- max(lmaChillSpMean)

# muForceSp
lmaMuFSpMean <- as.numeric(round(mean(lmaModelFit$muForceSp),1))
lower_lmaMuFSpMean <- format(round(quantile(lmaModelFit$muForceSp, prob = 0.05),1), nsmall =1)
upper_lmaMuFSpMean <- round(quantile(lmaModelFit$muForceSp, prob = 0.95),1)
lmaMuFSpMean; lower_lmaMuFSpMean; upper_lmaMuFSpMean

#muChillSP
lmaMuCSpMean <- as.numeric(round(mean(lmaModelFit$muChillSp),1))
lower_lmaMuCSpMean <- format(round(quantile(lmaModelFit$muChillSp, prob = 0.05),1), nsmall =1)
upper_lmaMuCSpMean <- format(round(quantile(lmaModelFit$muChillSp, prob = 0.95),1), nsmall =1)
lmaMuCSpMean; lower_lmaMuCSpMean; upper_lmaMuCSpMean

#muPhotoSp
lmaMuPSpMean <- format(as.numeric(round(mean(lmaModelFit$muPhotoSp),1)), nsmall = 1)
lower_lmaMuPSpMean <- format(round(quantile(lmaModelFit$muPhotoSp, prob = 0.05),1), nsmall =1)
upper_lmaMuPSpMean <- round(quantile(lmaModelFit$muPhotoSp, prob = 0.95),1)
lmaMuPSpMean; lower_lmaMuPSpMean; upper_lmaMuPSpMean

lmachill <- data.frame(lmaModelFit$alphaChillSp)
lmaChillSpMean <- colMeans(lmachill)
lmaChillMax <- max(lmaChillSpMean)

lmaLat <- as.numeric(round(mean(lmaModelFit$b_lat),1))


lmaLatTran <- as.numeric(round(mean(lmaModelFit$b_tranlat),1))
lower_lmaLatTran <- format(round(quantile(lmaModelFit$b_tranlat, prob = 0.05),1), nsmall =1)
upper_lmaLatTran <- round(quantile(lmaModelFit$b_tranlat, prob = 0.95),1)

lmasigmaSp <- as.numeric(round(mean(lmaModelFit$sigma_sp),1))
lower_lmasigmaSp <- format(round(quantile(lmaModelFit$sigma_sp, prob = 0.05),1), nsmall =1)
upper_lmasigmaSp <- round(quantile(lmaModelFit$sigma_sp, prob = 0.95),1)
##### C:N ###############################

load("..//analysis/output/lncContLatHundowLat.Rdata")
lncModelFit <- rstan::extract(mdlPerN)
muSp <- data.frame(lncModelFit$mu_grand_sp)
muSpMean <- colMeans(muSp)
lncFold <- round(max(muSpMean)/min(muSpMean),0)

lncBFSpMean <- as.numeric(round(mean(lncModelFit$betaTraitxForce),1))
lower_lncBFSpMean <- format(round(quantile(lncModelFit$betaTraitxForce, prob = 0.05),1), nsmall =1)
upper_lncBFSpMean <- round(quantile(lncModelFit$betaTraitxForce, prob = 0.95),1)
lncBFSpMean; lower_lncBFSpMean; upper_lncBFSpMean

lncBCSpMean <- as.numeric(round(mean(lncModelFit$betaTraitxChill),1))
lower_lncBCSpMean <- round(quantile(lncModelFit$betaTraitxChill, prob = 0.05),1)
upper_lncBCSpMean <- round(quantile(lncModelFit$betaTraitxChill, prob = 0.95),1)
lncBCSpMean; lower_lncBCSpMean; upper_lncBCSpMean

lncBPSpMean <- as.numeric(round(mean(lncModelFit$betaTraitxPhoto),1))
lower_lncBPSpMean <- round(quantile(lncModelFit$betaTraitxPhoto, prob = 0.05),1)
upper_lncBPSpMean <- format(round(quantile(lncModelFit$betaTraitxPhoto, prob = 0.95),1), nsmall =1)
lncBPSpMean; lower_lncBPSpMean; upper_lncBPSpMean

lncchill <- data.frame(lncModelFit$alphaChillSp)
lncChillSpMean <- colMeans(lncchill)
lncChillMax <- max(lncChillSpMean)

# muForceSp
lncMuFSpMean <- as.numeric(round(mean(lncModelFit$muForceSp),1))
lower_lncMuFSpMean <- format(round(quantile(lncModelFit$muForceSp, prob = 0.05),1), nsmall =1)
upper_lncMuFSpMean <- round(quantile(lncModelFit$muForceSp, prob = 0.95),1)
lncMuFSpMean; lower_lncMuFSpMean; upper_lncMuFSpMean

#muChillSP
lncMuCSpMean <- as.numeric(round(mean(lncModelFit$muChillSp),1))
lower_lncMuCSpMean <- format(round(quantile(lncModelFit$muChillSp, prob = 0.05),1), nsmall =1)
upper_lncMuCSpMean <- round(quantile(lncModelFit$muChillSp, prob = 0.95),1)
lncMuCSpMean; lower_lncMuCSpMean; upper_lncMuCSpMean

#muPhotoSp
lncMuPSpMean <- as.numeric(round(mean(lncModelFit$muPhotoSp),1))
lower_lncMuPSpMean <- format(round(quantile(lncModelFit$muPhotoSp, prob = 0.05),1), nsmall =1)
upper_lncMuPSpMean <- round(quantile(lncModelFit$muPhotoSp, prob = 0.95),1)
lncMuPSpMean; lower_lncMuPSpMean; upper_lncMuPSpMean

lncchill <- data.frame(lncModelFit$alphaChillSp)
lncChillSpMean <- colMeans(lncchill)
lncChillMax <- max(lncChillSpMean)

lncLatTran <- as.numeric(round(mean(lncModelFit$b_tranlat),1))
lower_lncLatTran <- format(round(quantile(lncModelFit$b_tranlat, prob = 0.05),1), nsmall =1)
upper_lncLatTran <- round(quantile(lncModelFit$b_tranlat, prob = 0.95),1)

lncLat <- as.numeric(round(mean(lncModelFit$b_lat),1))
lower_lncLat <- format(round(quantile(lncModelFit$b_lat, prob = 0.05),1), nsmall =1)
upper_lncLat <- format(-0.0, nsmall =1)#format(round(quantile(lncModelFit$b_lat, prob = 0.95),1), nsmall =1)
#lncmuSp <- apply(posterior_ssd$muSp, MARGIN = 2, FUN = mean)

lncTran <- as.numeric(round(mean(lncModelFit$b_tranE),1))
lower_lncTran <- format(round(quantile(lncModelFit$b_tranE, prob = 0.05),1), nsmall =1)
upper_lncTran <- round(quantile(lncModelFit$b_tranE, prob = 0.95),1)
#lncmuSp <- apply(posterior_ssd$muSp, MARGIN = 2, FUN = mean)

lncsigmaSp <- as.numeric(round(mean(lncModelFit$sigma_sp),1))
lower_lncsigmaSp <- format(round(quantile(lncModelFit$sigma_sp, prob = 0.05),1), nsmall =1)
upper_lncsigmaSp <- round(quantile(lncModelFit$sigma_sp, prob = 0.95),1)
##### SSD  ###############################
load("..//analysis/output/ssdContLatHundowLat10.Rdata")
ssdModelFit <- rstan::extract(mdlSSD)
muSp <- data.frame(ssdModelFit$mu_grand_sp)
muSpMean <- colMeans(muSp)
ssdFold <- round(max(muSpMean)/min(muSpMean),0)

ssdBFSpMean <- as.numeric(round(mean(ssdModelFit$betaTraitxForce),1))
lower_ssdBFSpMean <- format(round(quantile(ssdModelFit$betaTraitxForce, prob = 0.05),1), nsmall =1)
upper_ssdBFSpMean <- round(quantile(ssdModelFit$betaTraitxForce, prob = 0.95),1)
ssdBFSpMean; lower_ssdBFSpMean; upper_ssdBFSpMean

ssdBCSpMean <- as.numeric(round(mean(ssdModelFit$betaTraitxChill),1))
lower_ssdBCSpMean <- round(quantile(ssdModelFit$betaTraitxChill, prob = 0.05),1)
upper_ssdBCSpMean <- round(quantile(ssdModelFit$betaTraitxChill, prob = 0.95),1)
ssdBCSpMean; lower_ssdBCSpMean; upper_ssdBCSpMean

ssdBPSpMean <- as.numeric(round(mean(ssdModelFit$betaTraitxPhoto),1))
lower_ssdBPSpMean <- round(quantile(ssdModelFit$betaTraitxPhoto, prob = 0.05),1)
upper_ssdBPSpMean <- format(round(quantile(ssdModelFit$betaTraitxPhoto, prob = 0.95),1), nsmall =1)
ssdBPSpMean; lower_ssdBPSpMean; upper_ssdBPSpMean

ssdchill <- data.frame(ssdModelFit$alphaChillSp)
ssdChillSpMean <- colMeans(ssdchill)
ssdChillMax <- max(ssdChillSpMean)

# muForceSp
ssdMuFSpMean <- as.numeric(round(mean(ssdModelFit$muForceSp),1))
lower_ssdMuFSpMean <- format(round(quantile(ssdModelFit$muForceSp, prob = 0.05),1), nsmall =1)
upper_ssdMuFSpMean <- round(quantile(ssdModelFit$muForceSp, prob = 0.95),1)
ssdMuFSpMean; lower_ssdMuFSpMean; upper_ssdMuFSpMean

#muChillSP
ssdMuCSpMean <- as.numeric(round(mean(ssdModelFit$muChillSp),1))
lower_ssdMuCSpMean <- format(round(quantile(ssdModelFit$muChillSp, prob = 0.05),1), nsmall =1)
upper_ssdMuCSpMean <- round(quantile(ssdModelFit$muChillSp, prob = 0.95),1)
ssdMuCSpMean; lower_ssdMuCSpMean; upper_ssdMuCSpMean

#muPhotoSp
ssdMuPSpMean <- as.numeric(round(mean(ssdModelFit$muPhotoSp),1))
lower_ssdMuPSpMean <- format(round(quantile(ssdModelFit$muPhotoSp, prob = 0.05),1), nsmall =1)
upper_ssdMuPSpMean <- round(quantile(ssdModelFit$muPhotoSp, prob = 0.95),1)
ssdMuPSpMean; lower_ssdMuPSpMean; upper_ssdMuPSpMean

ssdchill <- data.frame(ssdModelFit$alphaChillSp)
ssdChillSpMean <- colMeans(ssdchill)
ssdChillMax <- max(ssdChillSpMean)

ssdLat <- as.numeric(round(mean(ssdModelFit$b_lat),1))

ssdLatTran <- as.numeric(round(mean(ssdModelFit$b_tranlat),1))
lower_ssdLatTran <- format(round(quantile(ssdModelFit$b_tranlat, prob = 0.05),1), nsmall =1)
upper_ssdLatTran <- round(quantile(ssdModelFit$b_tranlat, prob = 0.95),1)

ssdsigmaSp <- as.numeric(round(mean(ssdModelFit$sigma_sp),1))
lower_ssdsigmaSp <- format(round(quantile(ssdModelFit$sigma_sp, prob = 0.05),1), nsmall =1)
upper_ssdsigmaSp <- round(quantile(ssdModelFit$sigma_sp, prob = 0.95),1)
########## DBH ###############################
load("..//analysis/output/dbhContLatHundowLat.Rdata")
dbhModelFit <- rstan::extract(mdlDBH)
muSp <- data.frame(dbhModelFit$mu_grand_sp)
muSpMean <- colMeans(muSp)
dbhFold <- round(max(muSpMean)/min(muSpMean),0)

dbhBFSpMean <- as.numeric(round(mean(dbhModelFit$betaTraitxForce),1))
lower_dbhBFSpMean <- format(round(quantile(dbhModelFit$betaTraitxForce, prob = 0.05),1), nsmall =1)
upper_dbhBFSpMean <- round(quantile(dbhModelFit$betaTraitxForce, prob = 0.95),1)
dbhBFSpMean; lower_dbhBFSpMean; upper_dbhBFSpMean

dbhBCSpMean <- as.numeric(round(mean(dbhModelFit$betaTraitxChill),1))
lower_dbhBCSpMean <- round(quantile(dbhModelFit$betaTraitxChill, prob = 0.05),1)
upper_dbhBCSpMean <- round(quantile(dbhModelFit$betaTraitxChill, prob = 0.95),1)
dbhBCSpMean; lower_dbhBCSpMean; upper_dbhBCSpMean

dbhBPSpMean <- as.numeric(round(mean(dbhModelFit$betaTraitxPhoto),1))
lower_dbhBPSpMean <- round(quantile(dbhModelFit$betaTraitxPhoto, prob = 0.05),1)
upper_dbhBPSpMean <- format(round(quantile(dbhModelFit$betaTraitxPhoto, prob = 0.95),1), nsmall =1)
dbhBPSpMean; lower_dbhBPSpMean; upper_dbhBPSpMean

dbhchill <- data.frame(dbhModelFit$alphaChillSp)
dbhChillSpMean <- colMeans(dbhchill)
dbhChillMax <- max(dbhChillSpMean)

# muForceSp
dbhMuFSpMean <- as.numeric(round(mean(dbhModelFit$muForceSp),1))
lower_dbhMuFSpMean <- format(round(quantile(dbhModelFit$muForceSp, prob = 0.05),1), nsmall =1)
upper_dbhMuFSpMean <- round(quantile(dbhModelFit$muForceSp, prob = 0.95),1)
dbhMuFSpMean; lower_dbhMuFSpMean; upper_dbhMuFSpMean

#muChillSP
dbhMuCSpMean <- as.numeric(round(mean(dbhModelFit$muChillSp),1))
lower_dbhMuCSpMean <- format(round(quantile(dbhModelFit$muChillSp, prob = 0.05),1), nsmall =1)
upper_dbhMuCSpMean <- round(quantile(dbhModelFit$muChillSp, prob = 0.95),1)
dbhMuCSpMean; lower_dbhMuCSpMean; upper_dbhMuCSpMean

#muPhotoSp
dbhMuPSpMean <- as.numeric(round(mean(dbhModelFit$muPhotoSp),1))
lower_dbhMuPSpMean <- format(round(quantile(dbhModelFit$muPhotoSp, prob = 0.05),1), nsmall =1)
upper_dbhMuPSpMean <- round(quantile(dbhModelFit$muPhotoSp, prob = 0.95),1)
dbhMuPSpMean; lower_dbhMuPSpMean; upper_dbhMuPSpMean

dbhchill <- data.frame(dbhModelFit$alphaChillSp)
dbhChillSpMean <- colMeans(dbhchill)
dbhChillMax <- max(dbhChillSpMean)

dbhLatTran <- as.numeric(round(mean(dbhModelFit$b_tranlat),1))
lower_dbhLatTran <- format(round(quantile(dbhModelFit$b_tranlat, prob = 0.05),1), nsmall =1)
upper_dbhLatTran <- format(round(quantile(dbhModelFit$b_tranlat, prob = 0.95),1), nsmall = 1)

dbhsigmaSp <- as.numeric(round(mean(dbhModelFit$sigma_sp),1))
lower_dbhsigmaSp <- format(round(quantile(dbhModelFit$sigma_sp, prob = 0.05),1), nsmall =1)
upper_dbhsigmaSp <- round(quantile(dbhModelFit$sigma_sp, prob = 0.95),1)

tot <- read.csv("..//analysis/input/phenoMini.csv")

#Only want the total number that reached 7:
below <- subset(tot, bbch.t < 8)
totalObs <- nrow(below)


totalDays <- max(tot$day)

#######################################
sumHt <- summary(mdlHt)$summary
postHt <- rstan::extract(mdlHt)

a_sp = (sumHt[grep("mu_grand_sp", rownames(sumHt)), 1])
b_photo = sumHt[grep("betaPhotoSp\\[", rownames(sumHt)), 1]
b_chill = sumHt[grep("betaChillSp\\[", rownames(sumHt)), 1]
b_force = sumHt[grep("betaForceSp\\[", rownames(sumHt)), 1]

a_sp5 <- vector()
for(i in 1:ncol(postHt$mu_grand_sp)){
  quantU <- round(quantile(postHt$mu_grand_sp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  a_sp5 <- rbind(a_sp5, quantU)
}
colnames(a_sp5) <- c("Int5","Int95","Int25","Int75")
#a_sp5 <- a_sp5/100

b_chill5 <- vector()
for(i in 1:ncol(postHt$betaChillSp)){
  quantU <- round(quantile(postHt$betaChillSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_chill5 <- rbind(b_chill5, quantU)
}
colnames(b_chill5) <- c("chill5","chill95","chill25","chill75")

b_force5 <- vector()
for(i in 1:ncol(postHt$betaForceSp)){
  quantU <- round(quantile(postHt$betaForceSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_force5 <- rbind(b_force5, quantU)
}
colnames(b_force5) <- c("force5","force95","force25","force75")

b_photo5 <- vector()
for(i in 1:ncol(postHt$betaPhotoSp)){
  quantU <- round(quantile(postHt$betaPhotoSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_photo5 <- rbind(b_photo5, quantU)
}
colnames(b_photo5) <- c("photo5","photo95","photo25","photo75")

# #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#

#If we are using the old model, we will use the z-scored values for the parameters
photo <- -0.5033863 #8 h photo
force <- -0.3568628 #5/15 C trt
chill <- -0.3546922 # low chill

m <- matrix(nrow = 1000, ncol = 47)

for(sp in 1:47){
  for (it in 1:nrow(m)){
    m[it,sp] <- postHt$mu_grand_sp[it,sp]+  
      postHt$betaForceSp[it,sp] * force + 
      postHt$betaPhotoSp[it, sp] * photo + 
      postHt$betaChillSp[it,sp] * chill 
  }
}

photoHigh <- 0.4965051 #8 h photo
siteSM <- 0
forceHigh <- 0.5877121 
chillHigh <- 0.3660412 # high chill for Smithers

mHigh <- matrix(nrow = 1000, ncol = 47)

for(sp in 1:47){
  for (it in 1:nrow(mHigh)){
    mHigh[it,sp] <- postHt$mu_grand_sp[it,sp]+  
      postHt$betaForceSp[it,sp] * forceHigh + 
      postHt$betaPhotoSp[it, sp] * photoHigh + 
      postHt$betaChillSp[it,sp] * chillHigh 
  }
}

spInfo <- spInfo[order(spInfo$species),]
spInfo$meanBB <- colMeans(m)
colnames(m) <- spInfo$species.name

spInfo$meanBBHigh <- colMeans(mHigh)
colnames(mHigh) <- spInfo$species.name

spInfo$Int <- a_sp
spInfo <- cbind(spInfo, a_sp5,b_force5, b_chill5,b_photo5)

spInfo$force <- b_force
spInfo$chill <- b_chill
spInfo$photo <- b_photo

quantile595 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.05, 0.95, 0.25,0.75))
  return(returnQuanilte)
}

bb_quan <- apply(m, 2, quantile595)
bb_t <- t(bb_quan)
bb_df <- data.frame(bb_t)
colnames(bb_df)[colnames(bb_df) == "X5."] <- "bb5"
colnames(bb_df)[colnames(bb_df) == "X95."] <- "bb95"
colnames(bb_df)[colnames(bb_df) == "X25."] <- "bb25"
colnames(bb_df)[colnames(bb_df) == "X75."] <- "bb75"


bb_quanHigh <- apply(mHigh, 2, quantile595)
bb_tHigh <- t(bb_quanHigh)
bb_dfHigh <- data.frame(bb_tHigh)
colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X5."] <- "bb5High"
colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X95."] <- "bb95High"
colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X25."] <- "bb25High"
colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X75."] <- "bb75High"

spInfo <- cbind(spInfo, bb_df)
spInfo$value <- spInfo$meanBB

spInfo <- cbind(spInfo, bb_dfHigh)
spInfo$valueHigh <- spInfo$meanBBHigh

m <- data.frame(m)

long <- reshape2::melt(m)
names(long) <- c("species.name", "valueLow")
 
mHigh <- data.frame(mHigh)

longHigh <- reshape2::melt(mHigh)
names(longHigh) <- c("species.name", "valueHigh")

long <- cbind(long, longHigh[,2])

long <- merge(long,spInfo, by = "species.name")

spOrderData <- spInfo[order(spInfo$meanBB),]
spOrder <- as.factor(spOrderData$species.name)
# 
# spOrderHtData <- spInfo[order(spInfo$ht),]
# spOrderHt <- as.factor(spOrderData$species.name)
# 
# long <- long[order(long$species),]

# longPhotoInfo$mean <- rowMeans(longPhotoInfo[,c("Site1","Site2","Site3","Site4")], na.rm=TRUE)

bChill <- data.frame(postHt$betaChillSp[1:1000,])
colnames(bChill) <- (spInfo$species.name)
longChill <- reshape2::melt(bChill)
names(longChill) <- c("species.name", "chill")

long <- cbind(long, longChill$chill)

# Add forcing
bForce <- data.frame(postHt$betaForceSp[1:1000,])
colnames(bForce) <- (spInfo$species.name)
longForce <- reshape2::melt(bForce)
names(longForce) <- c("species.name", "force")

long <- cbind(long, longForce$force)

# photoperiod
bPhoto <- data.frame(postHt$betaPhotoSp[1:1000,])
colnames(bPhoto) <- (spInfo$species.name)
longPhoto <- reshape2::melt(bPhoto)
names(longPhoto) <- c("species.name", "photo")

long <- cbind(long, longPhoto$photo)

# intercept
aSp <- data.frame(postHt$mu_grand_sp[1:1000,])
colnames(aSp) <- (spInfo$species.name)
longInt <- reshape2::melt(aSp)
names(longInt) <- c("species.name", "int")

long <- cbind(long, longInt$int)

data <- long[order(long$meanBB),]

shrub <- subset(data, type == "shrub")
meanShrubBB <- round(mean(shrub$meanBB), 1)

tree <- subset(data, type == "tree")
meanTreeBB <- round(mean(tree$meanBB), 1)

