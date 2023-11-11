# One more thesis chapter to go!

# Started Oct 10 2023 by D Loughnan
# aim of this code is to run the model with latitude as a continuous variable on the real trait data for my bc trait chapter

rm(list=ls())
options(stringsAsFactors = FALSE)

## Load libraries
library(rstan)
require(shinystan)
require(stringr)
library(plyr)
library(dplyr)


library(ggplot2)
library(reshape2)
# library(viridis)
#library(bayesplot)
#library(tidybayes)
library(gridExtra) 

options(mc.cores = 4)

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
# pheno$force.n[pheno$force.n == "HF"] <- "1"
# pheno$force.n[pheno$force.n == "LF"] <- "0"
# pheno$force.n <- as.numeric(pheno$force.n)

pheno$force.n[pheno$force.n == "HF" & pheno$population == "mp"] <- "15"
pheno$force.n[pheno$force.n == "HF" & pheno$population == "sm"] <- "15"
pheno$force.n[pheno$force.n == "LF" & pheno$population == "mp"] <- "10"
pheno$force.n[pheno$force.n == "LF" & pheno$population == "sm"] <- "10"

# (8*20+16*10)/24 #13.33
# (8*15+16*5)/24 #8.33

pheno$force.n[pheno$force.n == "HF" & pheno$population == "HF"] <- "13.33"
pheno$force.n[pheno$force.n == "HF" & pheno$population == "SH"] <- "13.33"
pheno$force.n[pheno$force.n == "LF" & pheno$population == "HF"] <- "8.33"
pheno$force.n[pheno$force.n == "LF" & pheno$population == "SH"] <- "8.33"
pheno$force.n <- as.numeric(pheno$force.n)
# pheno$photo.n <- pheno$photo
# pheno$photo.n[pheno$photo.n == "HP"] <- "1"
# pheno$photo.n[pheno$photo.n == "LP"] <- "0"
# pheno$photo.n <- as.numeric(pheno$photo.n)

pheno$photo.n <- pheno$photo
pheno$photo.n[pheno$photo.n == "HP"] <- "12"
pheno$photo.n[pheno$photo.n == "LP"] <- "8"
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

# head(pheno)
# #add dummy/ site level effects:
# pheno <- pheno %>%
#   mutate ( site2 = if_else(site.n == 2, 1, 0),
#            site3 = if_else(site.n == 3, 1, 0),
#            site4 = if_else(site.n == 4, 1, 0))

pheno$force.z2 <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE))
pheno$photo.z2 <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE))
pheno$chillport.z2 <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE))

pheno.term <- pheno[,c("bb", "force.z2", "photo.z2", "population", "species", "lab2","Utah_Model","Chill_portions","chillport.z2","transect")] #"site2.z2", "site3.z2","site4.z2")]
pheno.t <- pheno.term[complete.cases(pheno.term$bb), ] # none,great!
length(unique(pheno.t$species))
#pheno.t <- pheno.term[complete.cases(pheno.term$bb), ] # 1780 rows data 
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

trtPheno$latitude <- trtPheno$site
trtPheno$latitude[trtPheno$latitude == "sm"] <- 54.7824
trtPheno$latitude[trtPheno$latitude == "mp"] <- 49.0646
trtPheno$latitude[trtPheno$latitude == "kl"] <- 52.1417
trtPheno$latitude[trtPheno$latitude == "af"] <- 50.8837
trtPheno$latitude[trtPheno$latitude == "GR"] <- 43.99837498
trtPheno$latitude[trtPheno$latitude == "HF"] <- 42.5315
trtPheno$latitude[trtPheno$latitude == "SH"] <- 45.9310
trtPheno$latitude[trtPheno$latitude == "WM"] <- 44.92466697

trtPheno$latitude <- as.numeric(trtPheno$latitude)
trtPheno$latZ <- (trtPheno$latitude-mean(trtPheno$latitude,na.rm=TRUE))/(sd(trtPheno$latitude,na.rm=TRUE))

specieslist <- sort(unique(trtPheno$species))
sitelist <- sort(unique(trtPheno$transect))

# write.csv(trtPheno ,"input/trtPhenoZScore.csv", row.names = F)
# write.csv(pheno.t,"input/bbPhenoZScore.csv", row.names = F)

######################################################################
# 1. Height 
height <- trtPheno[complete.cases(trtPheno$ht),]

ht.data <- list(yTraiti = height$ht, 
  N = nrow(height),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(height$species)),
  n_tran = length(unique(height$transect)),
  lati = height$latZ,
  tranE = as.numeric(height$transect),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z2,
  chilli = pheno.t$chillport.z2,
  photoi = pheno.t$photo.z2
)

length(ht.data$tranE)
head(ht.data$tranE)


# mdl <- stan("stan/heightDummyInt.stan",
#   data = ht.data,
#   iter = 4000, warmup = 3000, chains=4,
#   include = FALSE, pars = c("y_hat")
# )

mdlHt <- stan("stan/heightDummyIntGrand.stan",
  data = ht.data,
  iter = 4000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
)

save(mdl, file="output/heightDummyInt.Rdata")
save(mdlHt, file="output/heightDummyIntGrandZ.Rdata")
# ssm <- as.shinystan(mdl)
# launch_shinystan(ssm)

sumer <- summary(mdl)$summary

######################################################################
# 2. Leaf mass area
leafMA <- trtPheno[complete.cases(trtPheno$lma),]

lma.data <- list(yTraiti = leafMA$lma, 
  N = nrow(leafMA),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(leafMA$species)),
  n_tran = length(unique(leafMA$transect)),
  lati = leafMA$latitude,
  tranE = as.numeric(leafMA$transect),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z2,
  chilli = pheno.t$chillport.z2,
  photoi = pheno.t$photo.z2
)

lmaZ.data <- list(yTraiti = leafMA$lma, 
  N = nrow(leafMA),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(leafMA$species)),
  n_tran = length(unique(leafMA$transect)),
  #lati = leafMA$latitude,
  lati = leafMA$latZ,
  tranE = as.numeric(leafMA$transect),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z2,
  chilli = pheno.t$chillport.z2,
  photoi = pheno.t$photo.z2
)

mdlLMA <- stan("stan/lmaDummyInt.stan",
  data = lma.data,
  iter = 4000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
)

mdlLMA <- stan("stan/lmaDummyIntGrand.stan",
  data = lma.data,
  iter = 4000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
)

mdlLMA <- stan("stan/lmaDummyIntGrand.stan",
  data = lmaZ.data,
  iter = 4000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
)

sumer <- data.frame(summary(mdlLMA)$summary[c("mu_grand","b_tranE","b_tranlat", "muForceSp", "muChillSp", "muPhotoSp","muPhenoSp","betaTraitxForce", "betaTraitxChill","betaTraitxPhoto","sigma_traity" ,"sigma_sp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigmaPhenoSp","sigmapheno_y"),c("mean","2.5%","25%","50%", "75%","97.5%")])

bMuSp <- summary(mdl)$summary[grep("b_muSp\\["),c("mean","2.5%","25%","50%", "75%","97.5%")]

save(mdlLMA, file="output/lmaDummyIntGrandZ.Rdata")

sum <-summary(mdlLMA)$summary

# what if we just ran a simple linear model?
# require(lme4)
# simpLin <- lmer(lma ~ transect + transect*latitude + (1|species), leafMA)
# summary(simpLin)

# Random effects:
#   Groups   Name        Variance  Std.Dev.
# species  (Intercept) 9.064e-05 0.00952 
# Residual             2.591e-04 0.01610 
# Number of obs: 1345, groups:  species, 47
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)         0.1233246  0.0143426   8.599
# transect1          -0.2310547  0.0282366  -8.183
# latitude           -0.0016195  0.0002748  -5.894
# transect1:latitude  0.0050203  0.0006139   8.177

# mdl <- stan("stan/justDummyIntTrt.stan",
#   data = lma.data,
#   iter = 4000, warmup = 3000, chains=4,
#   include = FALSE, pars = c("y_hat")
# )
# 
# save(mdl, file="output/lmaDummyInt.Rdata")
# 
sumerTrt <- summary(mdl)$summary
# 
# require(rstanarm)
# 
# armMdl <- stan_lmer(lma ~ transect + transect*latitude + (1|species), data = leafMA)
# 
# sumerTrt <- summary(armMdl)
# head(sumerTrt)
# 
# bMuSpArm <- data.frame(sumerTrt)
# bMuSpArm <- bMuSpArm[5:51,]
# 
# pdf("figures/rstanArmBmuSpesti.pdf")
# hist(bMuSpArm$mean)
# dev.off()
# mean         mcse           sd          10%          50%          90% n_eff      Rhat
# (Intercept)                    0.115272446 2.550651e-04 0.0141260473  0.097737604  0.115132643  0.133558985  3067 0.9993842
# transect1                     -0.198132898 5.700797e-04 0.0263763362 -0.231257354 -0.198734331 -0.163476550  2141 1.0010089
# latitude                      -0.001466237 4.896403e-06 0.0002724896 -0.001821922 -0.001465479 -0.001124786  3097 0.9996112
# transect1:latitude             0.004305947 1.224876e-05 0.0005755368  0.003548270  0.004307676  0.005030400  2208 1.0010652

# stan_lmer
# family:       gaussian [identity]
# formula:      lma ~ transect + transect * latitude + (1 | species)
# observations: 1345
# ------
#   Median MAD_SD
# (Intercept)         0.1    0.0  
# transect1          -0.2    0.0  
# latitude            0.0    0.0  
# transect1:latitude  0.0    0.0  
# 
# Auxiliary parameter(s):
#   Median MAD_SD
# sigma 0.0    0.0   
# 
# Error terms:
#   Groups   Name        Std.Dev.
# species  (Intercept) 0.0098  
# Residual             0.0161  
# Num. levels: species 47 

######################################################################
# 3. diameter at breast height
diam <- trtPheno[complete.cases(trtPheno$dbh),] 

dbh.data <- list(yTraiti = diam$dbh, 
  N = nrow(diam),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(diam$species)),
  n_tran = length(unique(diam$transect)),
  lati = diam$latZ,
  tranE = as.numeric(diam$transect),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z2,
  chilli = pheno.t$chillport.z2,
  photoi = pheno.t$photo.z2
)

mdl <- stan("stan/diamDummyInt.stan",
  data = dbh.data,
  iter = 4000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
)

mdlDBH <- stan("stan/diamDummyIntGrand.stan",
  data = dbh.data,
  iter = 4000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
)

save(mdlDBH, file="output/dbhDummyIntGrandZ.Rdata")

######################################################################
# 4. stem specific density
stem <- trtPheno[complete.cases(trtPheno$ssd),]

ssd.data <- list(yTraiti = stem$ssd, 
  N = nrow(stem),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(stem$species)),
  n_tran = length(unique(stem$transect)),
  lati = stem$latitude,
  tranE = as.numeric(stem$transect),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z2,
  chilli = pheno.t$chillport.z2,
  photoi = pheno.t$photo.z2
)

mdl <- stan("stan/ssdDummyInt.stan",
  data = ssd.data,
  iter = 4000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
)

mdlSSD <- stan("stan/ssdDummyIntGrand.stan",
  data = ssd.data,
  iter = 4000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
)
save(mdlSSD, file="output/ssdDummyIntGrand.Rdata")

######################################################################
# 5. carbon to nitrogen ratio
carbNit <- trtPheno[complete.cases(trtPheno$C.N),]
carbNit$cn.z <- (carbNit$C.N-mean(carbNit$C.N,na.rm=TRUE))/(sd(carbNit$C.N,na.rm=TRUE))
carbNit$cn.port <- (carbNit$C.N/100)

cnZ.data <- list(yTraiti = carbNit$cn.z, 
  N = nrow(carbNit),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(carbNit$species)),
  n_tran = length(unique(carbNit$transect)),
  lati = carbNit$latZ,
  tranE = as.numeric(carbNit$transect),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z2,
  chilli = pheno.t$chillport.z2,
  photoi = pheno.t$photo.z2
)

mdl <- stan("stan/cnDummyInt.stan",
  data = cnZ.data,
  iter = 4000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
  #,control = list(adapt_delta = 0.99, max_treedepth =12)
)

mdlCN <- stan("stan/cnDummyIntGrand.stan",
  data = cnZ.data,
  iter = 4000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
  #,control = list(adapt_delta = 0.99, max_treedepth =12)
)

save(mdlCN, file="output/cnDummyIntGrandZboth.Rdata")

load("output/cnDummyInt.Rdata")
 ssm <- as.shinystan(mdl)
 launch_shinystan(ssm)

sumer <- summary(mdl)$summary

ModelFit <- rstan::extract(mdl)

muSp <- data.frame(ModelFit$b_muSp)
muSpMean <- colMeans(muSp)

betaForceSp <- data.frame(ModelFit$betaForceSp)
betaForceSpMean <- colMeans(betaForceSp)

quantile2575 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.25, 0.75))
  return(returnQuanilte)
}

bf_quan <- apply(betaForceSp, 2, quantile2575) 
mu_quan <- apply(muSp, 2, quantile2575)

bfs <- rbind(betaForceSpMean, bf_quan)
bfs_t <- t(bfs)
bfs_df <- data.frame(bfs_t)
colnames(bfs_df)[colnames(bfs_df) == "X25."] <- "force25"
colnames(bfs_df)[colnames(bfs_df) == "X75."] <- "force75"
bfs_df$species <- specieslist

mg<- rbind(muSpMean, mu_quan)
mg_t <- t(mg)
mg_df <- data.frame(mg_t)
colnames(mg_df)[colnames(mg_df) == "X25."] <- "trait25"
colnames(mg_df)[colnames(mg_df) == "X75."] <- "trait75"
mg_df$species <- specieslist

muForceSp <- data.frame(ModelFit$muForceSp)
muForceSpMean <- colMeans(muForceSp)

betaTraitxForce <- data.frame(ModelFit$betaTraitxForce)
betaTraitxForceMean <- colMeans(betaTraitxForce)

# mg_df_east <- mg_df[mg_df$species %in% eastSp, ]
# mg_df_west <- mg_df[mg_df$species %in% westSp, ]
# 
# bfs_df_east <- bfs_df[bfs_df$species %in% eastSp, ]
# bfs_df_west <- bfs_df[bfs_df$species %in% westSp, ]

plot( x= mg_df$muSpMean, y = bfs_df$betaForceSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bfs_df$force25), max(bfs_df$force75)), ylab = "Species level forcing slope", xlab = "Trait value", cex.lab = 1.5) # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles
arrows(
  mg_df[,"muSpMean"], # x mean
  bfs_df[,"force25"], # y 25
  mg_df[,"muSpMean"],
  bfs_df[,"force75"],
  length = 0, col= "#218380", lwd = 2
)

arrows(
  mg_df[,"trait25"], # x mean
  bfs_df[,"betaForceSpMean"], # y 25
  mg_df[,"trait75"], # x mean
  bfs_df[,"betaForceSpMean"],
  length = 0, col = "#218380", lwd = 2
)


mtext(side = 3, text = "C:N, Forcing", adj = 0, cex = 1.25)
for(j in 1:length(muForceSp[,1])){
  abline(a = muForceSp[j,], b = betaTraitxForceMean, col=alpha("#73d2de", 0.085))
}
abline(a=muForceSpMean, b=betaTraitxForceMean, col = "black")

