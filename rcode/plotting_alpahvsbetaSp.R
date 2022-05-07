#started May 6 2022 by Deirdre
# aim of this code is to generate the plots comparing the alpha and beta parameters in the joing model

rm(list=ls())
options(stringsAsFactors = FALSE)

require(dplyr)
library(stringr)
library(plyr)
library(rstan)


if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  else{
  setwd("/home/deirdre/pheno_bc") # for midge
}

# Start by getting the pheno data
dl <- read.csv("input/dl_allbb.csv")

temp <- str_split_fixed(dl$trt, "_", 3); head(temp)
dl$chill<- temp[,1]
dl$photo <- temp[,2]
dl$force <- temp[,3]

dl.chill <- read.csv("input/chilling_values_Hope_Smithers.csv")

dl.wchill <- merge(dl, dl.chill, by = c("population","chill"))
dl.wchill$lab3 <- dl.wchill$lab2
dl.wchill$lab2 <- paste(dl.wchill$species, dl.wchill$population, dl.wchill$rep, sep = "_")

df <- read.csv("input/df_dxb_prepped_data.csv")
df.chill <- read.csv("input/chilling_values_eastern.csv")
df.wchill <- merge(df, df.chill, by =c("population","chill"))
df.wchill <- df.wchill[, c("population", "chill","force","photo","lab2", "bb","species", "treatment","Chill_portions","Utah_Model")]

# mergeing the my data with DF
pheno <- rbind.fill(dl.wchill, df.wchill)
pheno$force.n <- pheno$force
pheno$force.n[pheno$force.n == "HF"] <- "1"
pheno$force.n[pheno$force.n == "LF"] <- "0"
pheno$force.n <- as.numeric(pheno$force.n)

pheno$photo.n <- pheno$photo
pheno$photo.n[pheno$photo.n == "HP"] <- "1"
pheno$photo.n[pheno$photo.n == "LP"] <- "0"
pheno$photo.n <- as.numeric(pheno$photo.n)

pheno$transect <- pheno$population
pheno$transect[pheno$transect == "sm"] <- "0"
pheno$transect[pheno$transect == "mp"] <- "0"
pheno$transect[pheno$transect == "kl"] <- "0"
pheno$transect[pheno$transect == "af"] <- "0"
pheno$transect[pheno$transect == "GR"] <- "1"
pheno$transect[pheno$transect == "HF"] <- "1"
pheno$transect[pheno$transect == "SH"] <- "1"
pheno$transect[pheno$transect == "WM"] <- "1"

pheno$force.z2 <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE)*2)
pheno$photo.z2 <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE)*2)
pheno$chillport.z2 <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE)*2)

pheno.term <- pheno[,c("bb", "force.z2", "photo.z2", "population", "species", "lab2","Utah_Model","Chill_portions","chillport.z2","transect")] #"site2.z2", "site3.z2","site4.z2")]
pheno.t <- pheno.term[complete.cases(pheno.term$bb), ] # none,great!
length(unique(pheno.t$species))

pheno.t$species <- tolower(pheno.t$species)
pheno.t$species.fact <- as.numeric(as.factor(pheno.t$species))
sort(unique(pheno.t$species.fact)) # 47 bc two species occur in both transects

# Now get the trait data and subset to only include spp we have pheno data for:
setwd("..//Treetraits")
trtData <- read.csv("data/allTrt.csv", stringsAsFactors = FALSE)
head(trtData)

phenoSp <- sort(unique(pheno.t$species))

trtPheno <- trtData[trtData$species %in% phenoSp, ]
length(unique(trtPheno$species))

trtPheno$transect <- trtPheno$site
trtPheno$transect[trtPheno$transect == "sm"] <- "0"
trtPheno$transect[trtPheno$transect == "mp"] <- "0"
trtPheno$transect[trtPheno$transect == "kl"] <- "0"
trtPheno$transect[trtPheno$transect == "af"] <- "0"
trtPheno$transect[trtPheno$transect == "GR"] <- "1"
trtPheno$transect[trtPheno$transect == "HF"] <- "1"
trtPheno$transect[trtPheno$transect == "SH"] <- "1"
trtPheno$transect[trtPheno$transect == "WM"] <- "1"

files <- list.files(path = "output", pattern ="_stanfit.RDS" )
files
#for (i in 1:length(files)){
i <- 5
  Model <- readRDS(paste("output/", files[i], sep = ""))
  
  ModelFit <- rstan::extract(Model)
  
  muGrandSp <- data.frame(ModelFit$mu_grand_sp)
  muGrandSpMean <- colMeans(muGrandSp)
  
  betaForceSp <- data.frame(ModelFit$betaForceSp)
  alphaForceSp <- data.frame(ModelFit$alphaForceSp)
  
  betaForceSpMean <- colMeans(betaForceSp)
  alphaForceSpMean <- colMeans(alphaForceSp)
  diff <- betaForceSp - alphaForceSp
  diffForceMean <- colMeans(diff)
  
  quantile2575 <- function(x){
    returnQuanilte <- quantile(x, prob = c(0.25, 0.75))
    return(returnQuanilte)
  }
  
  diff_quan <- apply(diff, 2, quantile2575)
  alpha_quan <- apply(alphaForceSp, 2, quantile2575)
  beta_quan <- apply(betaForceSp, 2, quantile2575)
  mugrand_quan <- apply(muGrandSp, 2, quantile2575)
  
  df <- rbind(diffForceMean, diff_quan)
  df_t <- t(df)
  df_df <- data.frame(df_t)
  colnames(df_df)[colnames(df_df) == "X25."] <- "force25"
  colnames(df_df)[colnames(df_df) == "X75."] <- "force75"
  
  alpha <- rbind(alphaForceSpMean, alpha_quan)
  a_t <- t(alpha)
  a_df <- data.frame(a_t)
  colnames(a_df)[colnames(a_df) == "X25."] <- "force25"
  colnames(a_df)[colnames(a_df) == "X75."] <- "force75"
  
  beta <- rbind(betaForceSpMean, beta_quan)
  b_t <- t(beta)
  b_df <- data.frame(b_t)
  colnames(b_df)[colnames(b_df) == "X25."] <- "force25"
  colnames(b_df)[colnames(b_df) == "X75."] <- "force75"
  
  mg<- rbind(muGrandSpMean, mugrand_quan)
  mg_t <- t(mg)
  mg_df <- data.frame(mg_t)
  colnames(mg_df)[colnames(mg_df) == "X25."] <- "trait25"
  colnames(mg_df)[colnames(mg_df) == "X75."] <- "trait75"
  
  
  muForceSp <- data.frame(ModelFit$muForceSp)
  muForceSpMean <- colMeans(muForceSp)
  
  betaTraitxForce<- data.frame(ModelFit$betaTraitxForce)
  betaTraitxForceMean <- colMeans(betaTraitxForce)
  
  pdf(paste("figures/force_bdecomp", files[i], ".pdf", sep = ""), height = 5, width =15)
  par(mfrow = c(1,3))
  plot( x= mg_df$muGrandSpMean, y = df_df$diffForceMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(df_df$force25), max(df_df$force75)), xlab = "trait effect", ylab = "alpha-beta") # blank plot with x range 
  # 3 columns, mean, quantile
  # min and max defined by quantiles
  arrows(
    mg_df[,"muGrandSpMean"], # x mean
    df_df[,"force25"], # y 25
    mg_df[,"muGrandSpMean"],
    df_df[,"force75"],
    length = 0
  )
  
  arrows(
    mg_df[,"trait25"], # x mean
    df_df[,"diffForceMean"], # y 25
    mg_df[,"trait75"], # x mean
    df_df[,"diffForceMean"],
    length = 0
  )
  
  for(r in 1:length(betaTraitxForce[,1])){
    abline(a = 0, b = betaTraitxForce[r,], col=alpha("lightpink", 0.015))
  }
  abline(a = 0, b=betaTraitxForceMean, col = "grey")
  # dev.off()
  #________________________________________________________________#
  plot( x= mg_df$muGrandSpMean, y = a_df$alphaForceSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(a_df$force25), max(a_df$force75))) # blank plot with x range 
  
  # 3 columns, mean, quantile
  # min and max defined by quantiles
  arrows(
    mg_df[,"muGrandSpMean"], # x mean
    a_df[,"force25"], # y 25
    mg_df[,"muGrandSpMean"],
    a_df[,"force75"],
    length = 0
  )
  
  arrows(
    mg_df[,"trait25"], # x mean
    a_df[,"alphaForceSpMean"], # y 25
    mg_df[,"trait75"], # x mean
    a_df[,"alphaForceSpMean"],
    length = 0
  )
  
  #________________________________________________________________#
  plot( x= mg_df$muGrandSpMean, y = b_df$betaForceSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(b_df$force25), 5)) # blank plot with x range 
  # 3 columns, mean, quantile
  # min and max defined by quantiles
  arrows(
    mg_df[,"muGrandSpMean"], # x mean
    b_df[,"force25"], # y 25
    mg_df[,"muGrandSpMean"],
    b_df[,"force75"],
    length = 0
  )
  
  arrows(
    mg_df[,"trait25"], # x mean
    b_df[,"betaForceSpMean"], # y 25
    mg_df[,"trait75"], # x mean
    b_df[,"betaForceSpMean"],
    length = 0
  )
  for(r in 1:length(muForceSp[,1])){
    abline(a = muForceSp[r,], b = betaTraitxForceMean, col=alpha("lightpink", 0.015))
  }
  abline(a=muForceSpMean, b=betaTraitxForceMean, col = "grey")
  dev.off()
  
  #------------------------------------------------------------------------------#
  betaChillSp <- data.frame(ModelFit$betaChillSp)
  alphaChillSp <- data.frame(ModelFit$alphaChillSp)
  
  betaChillSpMean <- colMeans(betaChillSp)
  alphaChillSpMean <- colMeans(alphaChillSp)
  diff <- betaChillSp - alphaChillSp
  diffChillMean <- colMeans(diff)
  
  quantile2575 <- function(x){
    returnQuanilte <- quantile(x, prob = c(0.25, 0.75))
    return(returnQuanilte)
  }
  
  diff_quan <- apply(diff, 2, quantile2575)
  alpha_quan <- apply(alphaChillSp, 2, quantile2575)
  beta_quan <- apply(betaChillSp, 2, quantile2575)
  mugrand_quan <- apply(muGrandSp, 2, quantile2575)
  
  df <- rbind(diffChillMean, diff_quan)
  df_t <- t(df)
  df_df <- data.frame(df_t)
  colnames(df_df)[colnames(df_df) == "X25."] <- "chill25"
  colnames(df_df)[colnames(df_df) == "X75."] <- "chill75"
  
  alpha <- rbind(alphaChillSpMean, alpha_quan)
  a_t <- t(alpha)
  a_df <- data.frame(a_t)
  colnames(a_df)[colnames(a_df) == "X25."] <- "chill25"
  colnames(a_df)[colnames(a_df) == "X75."] <- "chill75"
  
  beta <- rbind(betaChillSpMean, beta_quan)
  b_t <- t(beta)
  b_df <- data.frame(b_t)
  colnames(b_df)[colnames(b_df) == "X25."] <- "chill25"
  colnames(b_df)[colnames(b_df) == "X75."] <- "chill75"
  
  mg<- rbind(muGrandSpMean, mugrand_quan)
  mg_t <- t(mg)
  mg_df <- data.frame(mg_t)
  colnames(mg_df)[colnames(mg_df) == "X25."] <- "trait25"
  colnames(mg_df)[colnames(mg_df) == "X75."] <- "trait75"
  
  
  muChillSp <- data.frame(ModelFit$muChillSp)
  muChillSpMean <- colMeans(muChillSp)
  
  betaTraitxChill<- data.frame(ModelFit$betaTraitxChill)
  betaTraitxChillMean <- colMeans(betaTraitxChill)
  
  pdf(paste("figures/chill_bdecomp", files[i], ".pdf", sep = ""), height = 5, width =15)
  par(mfrow = c(1,3))
  plot( x= mg_df$muGrandSpMean, y = df_df$diffChillMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(df_df$chill25), max(df_df$chill75))) # blank plot with x range 
  # 3 columns, mean, quantile
  # min and max defined by quantiles
  arrows(
    mg_df[,"muGrandSpMean"], # x mean
    df_df[,"chill25"], # y 25
    mg_df[,"muGrandSpMean"],
    df_df[,"chill75"],
    length = 0
  )
  
  arrows(
    mg_df[,"trait25"], # x mean
    df_df[,"diffChillMean"], # y 25
    mg_df[,"trait75"], # x mean
    df_df[,"diffChillMean"],
    length = 0
  )
  
  for(r in 1:length(betaTraitxChill[,1])){
    abline(a = 0, b = betaTraitxChill[r,], col=alpha("lightpink", 0.015))
  }
  abline(a = 0, b=betaTraitxChillMean, col = "grey")
  # dev.off()
  #________________________________________________________________#
  plot( x= mg_df$muGrandSpMean, y = a_df$alphaChillSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(a_df$chill25), max(a_df$chill75))) # blank plot with x range 
  
  # 3 columns, mean, quantile
  # min and max defined by quantiles
  arrows(
    mg_df[,"muGrandSpMean"], # x mean
    a_df[,"chill25"], # y 25
    mg_df[,"muGrandSpMean"],
    a_df[,"chill75"],
    length = 0
  )
  
  arrows(
    mg_df[,"trait25"], # x mean
    a_df[,"alphaChillSpMean"], # y 25
    mg_df[,"trait75"], # x mean
    a_df[,"alphaChillSpMean"],
    length = 0
  )
  
  #________________________________________________________________#
  plot( x= mg_df$muGrandSpMean, y = b_df$betaChillSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(b_df$chill25), 5)) # blank plot with x range 
  # 3 columns, mean, quantile
  # min and max defined by quantiles
  arrows(
    mg_df[,"muGrandSpMean"], # x mean
    b_df[,"chill25"], # y 25
    mg_df[,"muGrandSpMean"],
    b_df[,"chill75"],
    length = 0
  )
  
  arrows(
    mg_df[,"trait25"], # x mean
    b_df[,"betaChillSpMean"], # y 25
    mg_df[,"trait75"], # x mean
    b_df[,"betaChillSpMean"],
    length = 0
  )
  for(r in 1:length(muChillSp[,1])){
    abline(a = muChillSp[r,], b = betaTraitxChillMean, col=alpha("lightpink", 0.015))
  }
  abline(a=muChillSpMean, b=betaTraitxChillMean, col = "grey")  
  dev.off()
  
  #------------------------------------------------------------------------------#
  betaPhotoSp <- data.frame(ModelFit$betaPhotoSp)
  alphaPhotoSp <- data.frame(ModelFit$alphaPhotoSp)
  
  betaPhotoSpMean <- colMeans(betaPhotoSp)
  alphaPhotoSpMean <- colMeans(alphaPhotoSp)
  diff <- betaPhotoSp - alphaPhotoSp
  diffPhotoMean <- colMeans(diff)
  
  quantile2575 <- function(x){
    returnQuanilte <- quantile(x, prob = c(0.25, 0.75))
    return(returnQuanilte)
  }
  
  diff_quan <- apply(diff, 2, quantile2575)
  alpha_quan <- apply(alphaPhotoSp, 2, quantile2575)
  beta_quan <- apply(betaPhotoSp, 2, quantile2575)
  mugrand_quan <- apply(muGrandSp, 2, quantile2575)
  
  df <- rbind(diffPhotoMean, diff_quan)
  df_t <- t(df)
  df_df <- data.frame(df_t)
  colnames(df_df)[colnames(df_df) == "X25."] <- "photo25"
  colnames(df_df)[colnames(df_df) == "X75."] <- "photo75"
  
  alpha <- rbind(alphaPhotoSpMean, alpha_quan)
  a_t <- t(alpha)
  a_df <- data.frame(a_t)
  colnames(a_df)[colnames(a_df) == "X25."] <- "photo25"
  colnames(a_df)[colnames(a_df) == "X75."] <- "photo75"
  
  beta <- rbind(betaPhotoSpMean, beta_quan)
  b_t <- t(beta)
  b_df <- data.frame(b_t)
  colnames(b_df)[colnames(b_df) == "X25."] <- "photo25"
  colnames(b_df)[colnames(b_df) == "X75."] <- "photo75"
  
  mg<- rbind(muGrandSpMean, mugrand_quan)
  mg_t <- t(mg)
  mg_df <- data.frame(mg_t)
  colnames(mg_df)[colnames(mg_df) == "X25."] <- "trait25"
  colnames(mg_df)[colnames(mg_df) == "X75."] <- "trait75"
  
  
  muPhotoSp <- data.frame(ModelFit$muPhotoSp)
  muPhotoSpMean <- colMeans(muPhotoSp)
  
  betaTraitxPhoto<- data.frame(ModelFit$betaTraitxPhoto)
  betaTraitxPhotoMean <- colMeans(betaTraitxPhoto)
  
  pdf(paste("figures/photo_bdecomp", files[i], ".pdf", sep = ""), height = 5, width =15)
  par(mfrow = c(1,3))
  plot( x= mg_df$muGrandSpMean, y = df_df$diffPhotoMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(df_df$photo25), max(df_df$photo75))) # blank plot with x range 
  # 3 columns, mean, quantile
  # min and max defined by quantiles
  arrows(
    mg_df[,"muGrandSpMean"], # x mean
    df_df[,"photo25"], # y 25
    mg_df[,"muGrandSpMean"],
    df_df[,"photo75"],
    length = 0
  )
  
  arrows(
    mg_df[,"trait25"], # x mean
    df_df[,"diffPhotoMean"], # y 25
    mg_df[,"trait75"], # x mean
    df_df[,"diffPhotoMean"],
    length = 0
  )
  
  for(r in 1:length(betaTraitxPhoto[,1])){
    abline(a = 0, b = betaTraitxPhoto[r,], col=alpha("lightpink", 0.015))
  }
  abline(a=0, b=betaTraitxPhotoMean, col = "grey")
  # dev.off()
  #________________________________________________________________#
  plot( x= mg_df$muGrandSpMean, y = a_df$alphaPhotoSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(a_df$photo25), max(a_df$photo75))) # blank plot with x range 
  
  # 3 columns, mean, quantile
  # min and max defined by quantiles
  arrows(
    mg_df[,"muGrandSpMean"], # x mean
    a_df[,"photo25"], # y 25
    mg_df[,"muGrandSpMean"],
    a_df[,"photo75"],
    length = 0
  )
  
  arrows(
    mg_df[,"trait25"], # x mean
    a_df[,"alphaPhotoSpMean"], # y 25
    mg_df[,"trait75"], # x mean
    a_df[,"alphaPhotoSpMean"],
    length = 0
  )
  
  #________________________________________________________________#
  plot( x= mg_df$muGrandSpMean, y = b_df$betaPhotoSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(b_df$photo25), max(a_df$photo75))) # blank plot with x range 
  # 3 columns, mean, quantile
  # min and max defined by quantiles
  arrows(
    mg_df[,"muGrandSpMean"], # x mean
    b_df[,"photo25"], # y 25
    mg_df[,"muGrandSpMean"],
    b_df[,"photo75"],
    length = 0
  )
  
  arrows(
    mg_df[,"trait25"], # x mean
    b_df[,"betaPhotoSpMean"], # y 25
    mg_df[,"trait75"], # x mean
    b_df[,"betaPhotoSpMean"],
    length = 0
  )
  for(r in 1:length(muPhotoSp[,1])){
    abline(a = muPhotoSp[r,], b = betaTraitxPhotoMean, col=alpha("lightpink", 0.015))
  }
  abline(a=muPhotoSpMean, b=betaTraitxPhotoMean, col = "grey")
  dev.off()
  # }
  # 
  # ###################################
  # # getting plots of mean traits and cue slopes:
  # traitsMean$sp.fact = as.numeric(as.factor(traitsMean$speciesname))
  # htMean <- subset(traitsMean, traitname == "Plant_height_vegetative")
  # htDf <- cbind(htMean, bfs_df)
  # 
  # pdf("figures/heightvs")
  # plot(betaForceSpMean ~ traitvalue, data = htDf, pch = 21, col = "black", bg = "purple", xlab = "Height")
  