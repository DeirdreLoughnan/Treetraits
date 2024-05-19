# Started Aug 22 2023: DL
#the aim of this test data is to get the model working for my BC traits chapter

# the model: a joint model with trait effects then modeled with phenological cues  ---similar to traitors model but with some key differences
#1. site as a summy varaible---not study---only one intercept to fit

# first run: 7 div transitions, rhat 1.14, 3210 trans>max treedepth
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

require(rstan)

setwd("~/Documents/github/Treetraits")

# Nrep <- 20# rep per trait
# Npop <- 8
# Ntran <- 2
# Nspp <- 25# number of species with traits (making this 20 just for speed for now)
# 
# # First making a data frame for the test trait data
# Ntrt <- Nspp * Npop * Nrep * Ntran# total number of traits observations
# Ntrt
#  
# #make a dataframe for height
# trt.dat <- data.frame(matrix(NA, Ntrt, 1))
# names(trt.dat) <- c("rep")
# trt.dat$rep <- c(1:Nrep)
# trt.dat$species <- rep(1:Nspp, each = Nrep)
# trt.dat$pop <- rep(1:Npop, each = Nspp*Nrep*Ntran)
# trt.dat$tran <- rep(1:Ntran, each = Npop*Nrep*Nspp)
# 
# lati <- rnorm(8, 50, 5)
# trt.dat$lat <- rep(lati, each = Nrep*Nspp)
# trt.dat$lat <- as.numeric(trt.dat$lat)
# 
# mu.tranE <- 4
# trt.dat$dumE <- as.numeric(ifelse(trt.dat$tran == "1","0","1"))
# trt.dat$mutranE <- mu.tranE*trt.dat$dumE
# 
# mu.tranlat  <- 2
# sigma.tranlat = 1
# alpha.tranlat <- rnorm(Npop,mu.tranlat, sigma.tranlat)
# trt.dat$alpha.tranlat <- rep(alpha.tranlat, each = Nrep*Nspp)
# 
# # trt.dat$species <- rep(1:Nspp, Nstudy)
# 
# # now generating the species trait data, here it is for height
# mu.grand <- 10
# sigma.species <- 5 # we want to keep the variaiton across spp. high
# 
# #the alphaTraitSp in Faiths original code:
# mu.trtsp <- rnorm(Nspp, 0, sigma.species)
# trt.dat$mu.trtsp <- rep(mu.trtsp, each = Nrep) #adding ht data for ea. sp
# 
# # general variance
# trt.var <- 1 #sigma_traity in the stan code
# trt.dat$trt.er <- rnorm(Ntrt, 0, trt.var)
# 
# # generate yhat - heights -  for this first trt model
# for (i in 1:Ntrt){
#   trt.dat$yTraiti[i] <- mu.grand + 
#   trt.dat$mutranE[i]+ trt.dat$mu.trtsp[i] + trt.dat$trt.er[i] + trt.dat$alpha.tranlat[i] * (trt.dat$dumE[i]*trt.dat$lat[i])  
# }
# 
# # all.data <- list(lday = trt.dat$yTraiti,
# #                  N = Ntrt,
# #                  n_sp = Nspp,
# #                  sp = as.numeric(as.factor(trt.dat$species)),
# #                  #n_tran = Ntran,
# #                  warm = trt.dat$lat,
# #                  site2 = as.numeric(trt.dat$dumE),
# #                  trait_transect = as.numeric(as.factor(trt.dat$tran)))
# 
# all.data <- list(yTraiti = trt.dat$yTraiti,
#                  N = Ntrt,
#                  n_spec = Nspp,
#                  trait_species = as.numeric(as.factor(trt.dat$species)),
#                  n_tran = Ntran,
#                  lati = trt.dat$lat,
#                  tranE = as.numeric(trt.dat$dumE),
#                  trait_transect = as.numeric(as.factor(trt.dat$tran)))
# 
# 
# mdl <- stan("stan/justDummyInt.stan",
#             data = all.data,
#             iter = 4000, warmup = 3000, chains=4,
#             include = FALSE, pars = c("y_hat")
# )
# # save(mdl, file="output/testTraitLatitudeSept20NoGrand.Rdata")
# sumer <- summary(mdl)$summary
# 
# bTranE <- sumer[grep("b_tranE", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
# # muTranLat <- sumer[grep("sigma_tranlat", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
# # muTranLat <- sumer[grep("mu_tranlat", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
# # #muGrand <- sumer[grep("mu_grand", rownames(sumer)), c("mean","2.5%","97.5%")]
# # 
# # muForce <- sumer[grep("muForceSp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
# # muChill <- sumer[grep("muChillSp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
# # muPhoto <- sumer[grep("muPhotoSp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
# #sigma_tran <- sumer[grep("sigma_tran", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
# sigma_sp <- sumer[grep("sigma_sp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
# sigma_traity <- sumer[grep("sigma_traity", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
# 
# bTranE
# sigma_sp
# sigma_traity

##################################################################
# What about just the mu_sp and mu grand
# Nrep <- 40# rep per trait
# Nspp <- 40# number of species with traits (making this 20 just for speed for now)
# 
# # First making a data frame for the test trait data
# Ntrt <- Nspp * Nrep # total number of traits observations
# Ntrt
# 
# #make a dataframe for height
# trt.dat <- data.frame(matrix(NA, Ntrt, 1))
# names(trt.dat) <- c("rep")
# trt.dat$rep <- c(1:Nrep)
# trt.dat$species <- rep(1:Nspp, each = Nrep)
# 
# # trt.dat$species <- rep(1:Nspp, Nstudy)
# 
# # now generating the species trait data, here it is for height
# mu.grand <- 10
# sigma.species <- 5 # we want to keep the variaiton across spp. high
# 
# #the alphaTraitSp in Faiths original code:
# mu.trtsp <- rnorm(Nspp, 0, sigma.species)
# trt.dat$mu.trtsp <- rep(mu.trtsp, each = Nrep) #adding ht data for ea. sp
# 
# # general variance
# trt.var <- 1 #sigma_traity in the stan code
# trt.dat$trt.er <- rnorm(Ntrt, 0, trt.var)
# 
# # generate yhat - heights -  for this first trt model
# for (i in 1:Ntrt){
#   trt.dat$yTraiti[i] <- mu.grand + trt.dat$mu.trtsp[i] + trt.dat$trt.er[i]
# }
# 
# all.data <- list(yTraiti = trt.dat$yTraiti,
#                                   N = Ntrt,
#                                   n_spec = Nspp,
#                                   trait_species = as.numeric(as.factor(trt.dat$species)),
#                                   n_tran = Ntran,
#                                   lati = trt.dat$lat,
#                                   tranE = as.numeric(trt.dat$dumE),
#                                   trait_transect = as.numeric(as.factor(trt.dat$tran)))
# 
# 
# mdl <- stan("stan/justDummy.stan",
#             data = all.data,
#             iter = 4000, warmup = 3000, chains=4,
#             include = FALSE, pars = c("y_hat")
# )
# 
# sumer <- summary(mdl)$summary
# 
# muGrand <- sumer[grep("mu_grand", rownames(sumer)), c("mean","2.5%","97.5%")]
# sigma_sp <- sumer[grep("sigma_sp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
# sigma_traity <- sumer[grep("sigma_traity", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
# 
# muGrand
# sigma_sp
# sigma_traity

##################################################################
# musp+ muGrand + Dummy

Nrep <- 60# rep per trait
#Npop <- 8
Ntran <- 2
Nspp <- 60# number of species with traits (making this 20 just for speed for now)

# First making a data frame for the test trait data
Ntrt <- Nspp * Nrep * Ntran #* Npop# total number of traits observations
Ntrt

#make a dataframe for height
trt.dat <- data.frame(matrix(NA, Ntrt, 1))
names(trt.dat) <- c("rep")
trt.dat$rep <- c(1:Nrep)
trt.dat$species <- rep(1:Nspp, each = Nrep)
#trt.dat$pop <- rep(1:Npop, each = Nspp*Nrep*Ntran)
trt.dat$tran <- rep(1:Ntran, each = Nrep*Nspp)

# lati <- rnorm(8, 50, 5)
# trt.dat$lat <- rep(lati, each = Nrep*Nspp)
# trt.dat$lat <- as.numeric(trt.dat$lat)

mu.tranE <- 4
trt.dat$dumE <- as.numeric(ifelse(trt.dat$tran == "1","0","1"))
trt.dat$mutranE <- mu.tranE*trt.dat$dumE

# mu.tranlat  <- 2
# sigma.tranlat = 1
# alpha.tranlat <- rnorm(Npop,mu.tranlat, sigma.tranlat)
# trt.dat$alpha.tranlat <- rep(alpha.tranlat, each = Nrep*Nspp)

# trt.dat$species <- rep(1:Nspp, Nstudy)

# now generating the species trait data, here it is for height
#mu.grand <- 10
sigma.species <- 5 # we want to keep the variaiton across spp. high

#the alphaTraitSp in Faiths original code:

mu.trtsp <- rnorm(Nspp, 10, sigma.species)
trt.dat$mu.trtsp <- rep(mu.trtsp, each = Nrep) #adding ht data for ea. sp

# general variance
trt.var <- 1 #sigma_traity in the stan code
trt.dat$trt.er <- rnorm(Ntrt, 0, trt.var)

# generate yhat - heights -  for this first trt model
for (i in 1:Ntrt){
  trt.dat$yTraiti[i] <-  trt.dat$mu.trtsp[i] + trt.dat$trt.er[i] + trt.dat$mutranE[i]
}

all.data <- list(yTraiti = trt.dat$yTraiti,
                 N = Ntrt,
                 n_spec = Nspp,
                 trait_species = as.numeric(as.factor(trt.dat$species)),
                 n_tran = Ntran,
                 lati = trt.dat$lat,
                 tranE = as.numeric(trt.dat$dumE),
                 trait_transect = as.numeric(as.factor(trt.dat$tran)))


mdl <- stan("stan/justDummy2.stan",
            data = all.data,
            iter = 4000, warmup = 3000, chains=4,
            include = FALSE, pars = c("y_hat")
)

sumer <- summary(mdl)$summary

#muGrand <- sumer[grep("mu_grand", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
bTranE <- sumer[grep("b_tranE", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
sigma_sp <- sumer[grep("sigma_sp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
sigma_traity <- sumer[grep("sigma_traity", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]

mdl.out <- data.frame( "Parameter" = c("bTranE","sigma_sp","sigma_traity"),  
                       "Test.data.values" = c( mu.tranE,  sigma.species, trt.var) ,
                       "Estiamte"= c(bTranE[1], sigma_sp[1], sigma_traity[1]),
                       "2.5"= c(bTranE[2],   sigma_sp[2], sigma_traity[2]),
                       "25"= c( bTranE[3],   sigma_sp[3], sigma_traity[3]),
                       "50"= c( bTranE[4],   sigma_sp[4], sigma_traity[4]),
                       "75"= c( bTranE[5],   sigma_sp[5], sigma_traity[5]),
                       "97.5"= c(bTranE[6],  sigma_sp[6], sigma_traity[6]))

mdl.out

muTraitSp <- sumer[grep("muSp", rownames(sumer))]

plot(muTraitSp ~ mu.trtsp, xlab = "simulated muTraitSp", ylab = "mdl estimated muTraitSp")
abline(0,1)

##################################################
# muSp +  latitude as a slope
Nrep <- 20# rep per trait
Npop <- 8
Ntran <- 2
Nspp <- 25# number of species with traits (making this 20 just for speed for now)

# First making a data frame for the test trait data
Ntrt <- Nspp * Npop * Nrep# total number of traits observations
Ntrt

#make a dataframe for height
trt.dat <- data.frame(matrix(NA, Ntrt, 1))
names(trt.dat) <- c("rep")
trt.dat$rep <- c(1:Nrep)
trt.dat$species <- rep(1:Nspp, each = Nrep)
trt.dat$pop <- rep(1:Npop, each = Nspp*Nrep)
trt.dat$tran <- rep(1:Ntran, each = 4*Nrep*Nspp)


lati <- rnorm(8, 50, 5)
trt.dat$lat <- rep(lati, each = Nrep*Nspp)
trt.dat$lat <- as.numeric(trt.dat$lat)

mu.tranE <- 4
trt.dat$dumE <- as.numeric(ifelse(trt.dat$tran == "1","0","1"))
trt.dat$mutranE <- mu.tranE*trt.dat$dumE

mu.tranlat  <- 2
# sigma.tranlat = 1
#alpha.tranlat <- rnorm(Npop,mu.tranlat, sigma.tranlat)
#trt.dat$alpha.tranlat <- rep(alpha.tranlat, each = Nrep*Nspp)
trt.dat$alpha.tranlat <- mu.tranlat*(trt.dat$dumE*trt.dat$lat)

mu.lat <- 2
# trt.dat$species <- rep(1:Nspp, Nstudy)

# now generating the species trait data, here it is for height
#mu.grand <- 10
sigma.species <- 5 # we want to keep the variaiton across spp. high

#the alphaTraitSp in Faiths original code:
mu.trtsp <- rnorm(Nspp, 10, sigma.species)
trt.dat$mu.trtsp <- rep(mu.trtsp, each = Nrep) #adding ht data for ea. sp

# general variance
trt.var <- 1 #sigma_traity in the stan code
trt.dat$trt.er <- rnorm(Ntrt, 0, trt.var)

# generate yhat - heights -  for this first trt model
for (i in 1:Ntrt){
  trt.dat$yTraiti[i] <- 
    trt.dat$mu.trtsp[i] + trt.dat$trt.er[i] + mu.lat * (trt.dat$lat[i])  
}
hist(trt.dat$yTraiti)

all.data <- list(yTraiti = trt.dat$yTraiti,
                 N = Ntrt,
                 n_spec = Nspp,
                 trait_species = as.numeric(as.factor(trt.dat$species)),
                 lati = trt.dat$lat,
                 tranE = as.numeric(trt.dat$dumE)
)


mdl <- stan("stan/justSlope.stan",
            data = all.data,
            iter = 4000, warmup = 3000, chains=4,
            include = FALSE, pars = c("y_hat")
)

sumer <- summary(mdl)$summary

#bTranE <- sumer[grep("b_tranE", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
bTranLat <- sumer[grep("b_tranlat", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
bmuSp <- sumer[grep("b_muSp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
sigma_sp <- sumer[grep("sigma_sp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
sigma_traity <- sumer[grep("sigma_traity", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]


mdl.out <- data.frame( "Parameter" = c("bTranLat","sigma_sp","sigma_traity"),  
                       "Test.data.values" = c(  mu.tranlat, sigma.species, trt.var) ,
                       "Estiamte"= c( bTranLat[1],  sigma_sp[1], sigma_traity[1]),
                       "2.5"= c( bTranLat[2],  sigma_sp[2], sigma_traity[2]),
                       "25"= c(  bTranLat[3],  sigma_sp[3], sigma_traity[3]),
                       "50"= c(  bTranLat[4],  sigma_sp[4], sigma_traity[4]),
                       "75"= c(  bTranLat[5],  sigma_sp[5], sigma_traity[5]),
                       "97.5"= c( bTranLat[6],  sigma_sp[6], sigma_traity[6]))

mdl.out

muTraitSp <- sumer[grep("b_muSp", rownames(sumer))]

plot(muTraitSp ~ mu.trtsp, xlab = "simulated muTraitSp", ylab = "mdl estimated muTraitSp")
abline(0,1)

post <- rstan::extract(mdl)

plot(post$b_tranE~post$b_tranlat)
plot(post$muSp[,1]~post$sigma_sp)


##################################################
# muSp + dummy + just interaction
Nrep <- 20# rep per trait
Npop <- 8
Ntran <- 2
Nspp <- 25# number of species with traits (making this 20 just for speed for now)

# First making a data frame for the test trait data
Ntrt <- Nspp * Npop * Nrep# total number of traits observations
Ntrt

#make a dataframe for height
trt.dat <- data.frame(matrix(NA, Ntrt, 1))
names(trt.dat) <- c("rep")
trt.dat$rep <- c(1:Nrep)
trt.dat$species <- rep(1:Nspp, each = Nrep)
trt.dat$pop <- rep(1:Npop, each = Nspp*Nrep)
trt.dat$tran <- rep(1:Ntran, each = 4*Nrep*Nspp)


lati <- rnorm(8, 50, 5)
trt.dat$lat <- rep(lati, each = Nrep*Nspp)
trt.dat$lat <- as.numeric(trt.dat$lat)

mu.tranE <- 4
trt.dat$dumE <- as.numeric(ifelse(trt.dat$tran == "1","0","1"))
trt.dat$mutranE <- mu.tranE*trt.dat$dumE

mu.tranlat  <- 2
# sigma.tranlat = 1
#alpha.tranlat <- rnorm(Npop,mu.tranlat, sigma.tranlat)
#trt.dat$alpha.tranlat <- rep(alpha.tranlat, each = Nrep*Nspp)
trt.dat$alpha.tranlat <- mu.tranlat*(trt.dat$dumE*trt.dat$lat)

# trt.dat$species <- rep(1:Nspp, Nstudy)

# now generating the species trait data, here it is for height
#mu.grand <- 10
sigma.species <- 5 # we want to keep the variaiton across spp. high

#the alphaTraitSp in Faiths original code:
mu.trtsp <- rnorm(Nspp, 10, sigma.species)
trt.dat$mu.trtsp <- rep(mu.trtsp, each = Nrep) #adding ht data for ea. sp

# general variance
trt.var <- 1 #sigma_traity in the stan code
trt.dat$trt.er <- rnorm(Ntrt, 0, trt.var)

# generate yhat - heights -  for this first trt model
for (i in 1:Ntrt){
  trt.dat$yTraiti[i] <- 
    trt.dat$mu.trtsp[i] + trt.dat$trt.er[i] + trt.dat$mutranE[i]+ mu.tranlat * (trt.dat$dumE[i]*trt.dat$lat[i])  
}
# for (i in 1:Ntrt){
#   trt.dat$yTraiti[i] <- 
#     trt.dat$mu.trtsp[i] + trt.dat$trt.er[i] + trt.dat$mutranE[i]+ trt.dat$alpha.tranlat[i] * (trt.dat$dumE[i]*trt.dat$lat[i])  
# }

all.data <- list(yTraiti = trt.dat$yTraiti,
                 N = Ntrt,
                 n_spec = Nspp,
                 trait_species = as.numeric(as.factor(trt.dat$species)),
                 lati = trt.dat$lat,
                 tranE = as.numeric(trt.dat$dumE)
                 )


mdl <- stan("stan/justDummyIntTrait.stan",
            data = all.data,
            iter = 3000, warmup = 2000, chains=4,
            include = FALSE, pars = c("y_hat")
)

save(mdl, file="output/testLatitudeTraitOnly.Rdata")

sumer <- summary(mdl)$summary

bTranE <- sumer[grep("b_tranE", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
bTranLat <- sumer[grep("b_tranlat", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]

sigma_sp <- sumer[grep("sigma_sp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
sigma_traity <- sumer[grep("sigma_traity", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]


mdl.out <- data.frame( "Parameter" = c("bTranE","bTranLat","sigma_sp","sigma_traity"),  
                       "Test.data.values" = c( mu.tranE, mu.tranlat, sigma.species, trt.var) ,
                       "Estiamte"= c(bTranE[1], bTranLat[1],  sigma_sp[1], sigma_traity[1]),
                       "2.5"= c(bTranE[2], bTranLat[2],  sigma_sp[2], sigma_traity[2]),
                       "25"= c( bTranE[3], bTranLat[3],  sigma_sp[3], sigma_traity[3]),
                       "50"= c( bTranE[4], bTranLat[4],  sigma_sp[4], sigma_traity[4]),
                       "75"= c( bTranE[5], bTranLat[5],  sigma_sp[5], sigma_traity[5]),
                       "97.5"= c(bTranE[6], bTranLat[6],  sigma_sp[6], sigma_traity[6]))

mdl.out

muTraitSp <- sumer[grep("b_muSp", rownames(sumer))]

plot(muTraitSp ~ mu.trtsp, xlab = "simulated muTraitSp", ylab = "mdl estimated muTraitSp")
abline(0,1)

post <- rstan::extract(mdl)

plot(post$b_tranE~post$b_tranlat)
plot(post$muSp[,1]~post$sigma_sp)

table(trt.dat$pop, trt.dat$tranE)
