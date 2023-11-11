# aim of this code is to generate the values referenced in the sweave file:

library(rstan)
require(shinystan)
library(tidybayes)
library(reshape2)

########## Height ###############################

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


load("output/heightDummyIntGrandZ.Rdata")
htModelFit <- rstan::extract(mdlHt)

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

htchill <- data.frame(htModelFit$alphaChillSp)
htChillSpMean <- colMeans(htchill)
htChillMax <- max(htChillSpMean)

######## Leaf mass area###############################

load("output/lmaDummyIntGrandZ.Rdata")
lmaModelFit <- rstan::extract(mdlLMA)

lmaBFSpMean <- format((round(mean(lmaModelFit$betaTraitxForce),1)), nsmall =1)
lower_lmaBFSpMean <- round(quantile(lmaModelFit$betaTraitxForce, prob = 0.05),1)
upper_lmaBFSpMean <- round(quantile(lmaModelFit$betaTraitxForce, prob = 0.95),1)
lmaBFSpMean; lower_lmaBFSpMean; upper_lmaBFSpMean

lmaBCSpMean <- format((round(mean(lmaModelFit$betaTraitxChill),1)), nsmall =1)
lower_lmaBCSpMean <- round(quantile(lmaModelFit$betaTraitxChill, prob = 0.05),1)
upper_lmaBCSpMean <- round(quantile(lmaModelFit$betaTraitxChill, prob = 0.95),1)
lmaBCSpMean; lower_lmaBCSpMean; upper_lmaBCSpMean

lmaBPSpMean <- format((round(mean(lmaModelFit$betaTraitxPhoto),1)), nsmall =1)
lower_lmaBPSpMean <- round(quantile(lmaModelFit$betaTraitxPhoto, prob = 0.05),1)
upper_lmaBPSpMean <- round(quantile(lmaModelFit$betaTraitxPhoto, prob = 0.95),1)
lmaBPSpMean; lower_lmaBPSpMean; upper_lmaBPSpMean

##### C:N ###############################
load("output/cnDummyIntGrandZ.Rdata")
cnModelFit <- rstan::extract(mdlCN)

cnBFSpMean <- as.numeric(round(mean(cnModelFit$betaTraitxForce),1))
lower_cnBFSpMean <- round(quantile(cnModelFit$betaTraitxForce, prob = 0.05),1)
upper_cnBFSpMean <- round(quantile(cnModelFit$betaTraitxForce, prob = 0.95),1)
cnBFSpMean; lower_cnBFSpMean; upper_cnBFSpMean

cnBCSpMean <- as.numeric(round(mean(cnModelFit$betaTraitxChill),1))
lower_cnBCSpMean <- round(quantile(cnModelFit$betaTraitxChill, prob = 0.05),1)
upper_cnBCSpMean <- round(quantile(cnModelFit$betaTraitxChill, prob = 0.95),1)
cnBCSpMean; lower_cnBCSpMean; upper_cnBCSpMean

cnBPSpMean <- as.numeric(round(mean(cnModelFit$betaTraitxPhoto),1))
lower_cnBPSpMean <- round(quantile(cnModelFit$betaTraitxPhoto, prob = 0.05),1)
upper_cnBPSpMean <- round(quantile(cnModelFit$betaTraitxPhoto, prob = 0.95),1)
cnBPSpMean; lower_cnBPSpMean; upper_cnBPSpMean

#cnmuSp <- apply(posterior_ssd$muSp, MARGIN = 2, FUN = mean)

##### SSD  ###############################
load("output/ssdDummyIntGrandZ.Rdata")
ssdModelFit <- rstan::extract(mdlSSD)

ssdBFSpMean <- as.numeric(round(mean(ssdModelFit$betaTraitxForce),1))
lower_ssdBFSpMean <- round(quantile(ssdModelFit$betaTraitxForce, prob = 0.05),1)
upper_ssdBFSpMean <- round(quantile(ssdModelFit$betaTraitxForce, prob = 0.95),1)
ssdBFSpMean; lower_ssdBFSpMean; upper_ssdBFSpMean

ssdBCSpMean <- as.numeric(round(mean(ssdModelFit$betaTraitxChill),1))
lower_ssdBCSpMean <- round(quantile(ssdModelFit$betaTraitxChill, prob = 0.05),1)
upper_ssdBCSpMean <- round(quantile(ssdModelFit$betaTraitxChill, prob = 0.95),1)
ssdBCSpMean; lower_ssdBCSpMean; upper_ssdBCSpMean

ssdBPSpMean <- as.numeric(round(mean(ssdModelFit$betaTraitxPhoto),1))
lower_ssdBPSpMean <- round(quantile(ssdModelFit$betaTraitxPhoto, prob = 0.05),1)
upper_ssdBPSpMean <- round(quantile(ssdModelFit$betaTraitxPhoto, prob = 0.95),1)
ssdBPSpMean; lower_ssdBPSpMean; upper_ssdBPSpMean

ssdchill <- data.frame(ssdModelFit$alphaChillSp)
ssdChillSpMean <- colMeans(ssdchill)
ssdChillMax <- max(ssdChillSpMean)

########## DBH ###############################
load("output/dbhDummyIntGrandZ.Rdata")
dbhModelFit <- rstan::extract(mdlDBH)

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

