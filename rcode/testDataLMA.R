# Started Aug 22 2023: DL
#the aim of this test data is to get the model working for my BC traits chapter

# the model: a joint model with trait effects then modeled with phenological cues  ---similar to traitors model but with some key differences
#1. site as a dummy varaible---not study---only one intercept to fit

# Modified Nov 22 - still having issues with LMA, SSD, CN --- not able to estimate the beta parameters for these models---just fails---tried just logging the trait value, z-scoring everything, but it just didn't work
# Taking a step back and trying to replicate the issue in the test data


# Making this test data more similar to the lma data
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

require(rstan)
#require(truncnorm)

# setwd("~/Documents/github/Treetraits")
# trtPheno <- read.csv("input/trtPhenoDummy.csv")
# hist(trtPheno$lma)
# hist(log10(trtPheno$lma))
#pheno.t <- read.csv("input/phenoDataWChill.csv")

Nrep <- 20# rep per trait
Npop <- 8
Ntran <- 2
Nspp <- 15# number of species with traits (making this 20 just for speed for now)

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
trt.dat$latZ <- (trt.dat$lat-mean(trt.dat$lat,na.rm=TRUE))/(sd(trt.dat$lat,na.rm=TRUE))

mu.tranE <- 0.5
trt.dat$dumE <- as.numeric(ifelse(trt.dat$tran == "1","0","1"))
trt.dat$mutranE <- mu.tranE*trt.dat$dumE

mu.tranlat  <- 0.5
# sigma.tranlat = 1
#alpha.tranlat <- rnorm(Npop,mu.tranlat, sigma.tranlat)
#trt.dat$alpha.tranlat <- rep(alpha.tranlat, each = Nrep*Nspp)
trt.dat$alpha.tranlat <- mu.tranlat*(trt.dat$dumE*trt.dat$latZ)
# now generating the species trait data, here it is for height

#sigma.species <- 0.5 # we want to keep the variation across spp. high
#mu.grand <- 0.01

# we want the test data to be more like the lma data---most data is btwn 0-0.1, but some data as big as 0.35

# Option 1: could try using the beta distribution - get pretty goodd data shape, but can't figure out how to estimate sigma.sp, everything else gets good estimates 

# mu.trtsp <- rbeta(Nspp, 0.2, 8) # trying the beta distribution
# hist(mu.trtsp)

# perfectMatch <- c(9.055796e-03, 9.314260e-05, 6.424912e-02, 1.019174e-02, 7.409561e-03, 1.999574e-04, 2.244918e-03, 1.111483e-01, 5.114561e-02, 9.392062e-02, 1.589031e-06, 3.281615e-02, 3.265666e-03, 2.883886e-04, 5.830078e-02, 1.302096e-02,6.416374e-03, 1.610923e-05, 1.429036e-08, 6.533133e-03, 5.890707e-04, 3.188054e-01, 3.917652e-02, 1.490867e-03, 1.751903e-03) # randomly generated using rbeta(Nspp, 0.2,8), but matches the real data really well
# hist(perfectMatch)

# Option 2: Faith suggested doing a log transformation:
# log10 gives a mean around -1.5, with values that range from --2.5584884 to -0.4757686, all are negative
sigma.species <- 0.5 # we want to keep the variaiton across spp. high

#the alphaTraitSp in Faiths original code:
mu.trtsp <- rnorm(Nspp, -1.5, sigma.species)
trt.dat$mu.trtsp <- rep(mu.trtsp, each = Nrep) #adding ht data for ea. sp
#hist(mu.trtsp)
# general variance
trt.var <- 0.005 #sigma_traity in the stan code
trt.dat$trt.er <- rnorm(Ntrt, 0, trt.var)

# generate yhat - heights -  for this first trt model
for (i in 1:Ntrt){
  trt.dat$yTraiti[i] <- 
    trt.dat$mu.trtsp[i] + trt.dat$trt.er[i] + trt.dat$mutranE[i] + mu.tranlat * (trt.dat$dumE[i]*trt.dat$latZ[i])  
}


#hist(trt.dat$yTraiti)
# hist( trt.dat$mu.trtsp+ trt.dat$trt.er+ trt.dat$mutranE)
# hist(trt.dat$yTraiti)

# for (i in 1:Ntrt){
#   trt.dat$yTraiti[i] <- 
#     trt.dat$mu.trtsp[i] + trt.dat$trt.er[i] + trt.dat$mutranE[i]+ trt.dat$alpha.tranlat[i] * (trt.dat$dumE[i]*trt.dat$lat[i])  
# }

all.data <- list(yTraiti = trt.dat$yTraiti,
  N = Ntrt,
  n_spec = Nspp,
  trait_species = as.numeric(as.factor(trt.dat$species)),
  lati = trt.dat$latZ,
  tranE = as.numeric(trt.dat$dumE)
)


# mdl <- stan("stan/justDummyIntTrait.stan",
#   data = all.data,
#   iter = 4000, warmup = 3000, chains=4,
#   include = FALSE, pars = c("y_hat")
# )
# 
# sumer <- summary(mdl)$summary
# 
# bTranE <- sumer[grep("b_tranE", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
# bTranLat <- sumer[grep("b_tranlat", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
# 
# sigma_sp <- sumer[grep("sigma_sp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
# sigma_traity <- sumer[grep("sigma_traity", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
# 
# 
# mdl.out <- data.frame( "Parameter" = c("bTranE","bTranLat","sigma_sp","sigma_traity"),
#   "Test.data.values" = c( mu.tranE, mu.tranlat, sigma.species, trt.var) ,
#   "Estiamte"= c(bTranE[1], bTranLat[1],  sigma_sp[1], sigma_traity[1]),
#   "2.5"= c(bTranE[2], bTranLat[2],  sigma_sp[2], sigma_traity[2]),
#   "25"= c( bTranE[3], bTranLat[3],  sigma_sp[3], sigma_traity[3]),
#   "50"= c( bTranE[4], bTranLat[4],  sigma_sp[4], sigma_traity[4]),
#   "75"= c( bTranE[5], bTranLat[5],  sigma_sp[5], sigma_traity[5]),
#   "97.5"= c(bTranE[6], bTranLat[6],  sigma_sp[6], sigma_traity[6]))
# 
# mdl.out

# Parameter Test.data.values    Estiamte        X2.5         X25         X50         X75     X97.5
# 1       bTranE            0.400 0.400078204 0.399762494 0.399972208 0.400081657 0.400186009 0.4003788
# 2     bTranLat            0.500 0.500094131 0.499897676 0.500026294 0.500091901 0.500161112 0.5002968
# 3     sigma_sp            0.500 0.450210001 0.335431760 0.398142985 0.436145667 0.484332891 0.5996707
# 4 sigma_traity            0.005 0.004916394 0.004809158 0.004881562 0.004915106 0.004952073 0.0050256
##### Phenology test data ###########################

Nchill <- 2 # sm high low, mp high low
Nphoto <- 2 # high low
Nforce <- 2 # high amd low

Nrep <- 20
Nspp <- 15
Npop <- 2
Nph <- Nrep*Npop*Nspp*Nchill#*Nphoto*Nforce*Nphoto - bc N force and N chill are completely confounded

pheno.dat <- data.frame(matrix(NA, Nph, 2))
names(pheno.dat) <- c("rep","species")
pheno.dat$rep <- c(1:Nrep)
pheno.dat$species <- rep(c(1:Nspp), each = Nrep)
pheno.dat$pop <- rep(1:Npop, each = Nspp*Nrep)

chillTrt <- c("LC", "HC")
chilli <- rnorm(Nchill, 1, 5)
pheno.dat$chill <- rep(rep(chillTrt, each = Npop*Nrep*Nspp))
pheno.dat$chilli <- rep(rep(chilli, each = Npop*Nrep*Nspp))

# pheno.dat$chill <- rep(rep(chillTrt, each = Nphoto*Npop*Nrep*Nspp), times = Nforce)
# pheno.dat$chilli <- rep(rep(chilli, each = Nphoto*Npop*Nrep*Nspp), times = Nforce)

# forceTrt <- c("LF", "HF")
# forcei <- rnorm(Nforce, 1, 2)
# pheno.dat$force <- rep(forceTrt, each = Nchill*Nphoto*Npop*Nrep*Nspp)
# pheno.dat$forcei <- rep(forcei, each = Nchill*Nphoto*Npop*Nrep*Nspp)
# 
# photoTrt <- c("LP", "HP")
# photoi <- rnorm(Nphoto, 1, 2)
# pheno.dat$photo <- rep(rep(photoTrt, each = Npop*Nrep*Nspp), times = Nforce*Nchill)
# pheno.dat$photoi <- rep(rep(photoi, each = Npop*Nrep*Nspp), times = Nforce*Nchill)

# sanity check that all treatments are the same
# pheno.dat$label <- paste(pheno.dat$chill, pheno.dat$force, pheno.dat$photo, sep ="_")
# pheno.dat$label <- paste(pheno.dat$chilli, pheno.dat$forcei, pheno.dat$photoi, sep ="_")

mu.chill = -14
sigma.chill = 1
alpha.chill.sp <- rnorm(Nspp, mu.chill, sigma.chill)
pheno.dat$alphaChillSp <- rep(alpha.chill.sp, each = Nrep)

# mu.force = -10 
# sigma.force = 1
# alpha.force.sp <- rnorm(Nspp, mu.force, sigma.force)
# pheno.dat$alphaForceSp <- rep(alpha.force.sp, each = Nrep)
# 
# mu.photo = -5
# sigma.photo = 1
# alpha.photo.sp <- rnorm(Nspp, mu.photo, sigma.photo)
# pheno.dat$alphaPhotoSp <- rep(alpha.photo.sp, each = Nrep)

mu.pheno.sp = 80
sigma.pheno.sp = 30
alphaPhenoSp <- rnorm(Nspp, mu.pheno.sp, sigma.pheno.sp)
pheno.dat$alphaPhenoSp <- rep(alphaPhenoSp, each = Nrep)

sigma.pheno.y = 3
pheno.dat$ePhen <- rnorm(Nph, 0, sigma.pheno.y)

# betaTraitxForce <- 0.3 
# betaTraitxPhoto <- -0.2
betaTraitxChill <- -0.1

pheno.datTrait <- merge(pheno.dat, unique(trt.dat[,c("species","mu.trtsp")]), by = "species")
#head(pheno.datTrait,50)

for (i in 1:Nph){
  # pheno.datTrait$betaForceSp[i] <-  pheno.datTrait$alphaForceSp[i] + (betaTraitxForce *  pheno.datTrait$mu.trtsp[i])
  # 
  # pheno.datTrait$betaPhotoSp[i]<- pheno.datTrait$alphaPhotoSp[i] + (betaTraitxPhoto*  pheno.datTrait$mu.trtsp[i])
  # 
  pheno.datTrait$betaChillSp[i] <-pheno.datTrait$alphaChillSp[i] + (betaTraitxChill* pheno.datTrait$mu.trtsp[i])
}

for (i in 1:Nph){
  pheno.datTrait$yMu[i] <-  pheno.datTrait$alphaPhenoSp[i] + pheno.datTrait$betaChillSp[i] * pheno.datTrait$chilli[i] 
  #+  pheno.datTrait$betaForceSp[i] * pheno.datTrait$forcei[i] +  pheno.datTrait$betaPhotoSp[i] * pheno.datTrait$photoi[i] 
}

pheno.datTrait$yPhenoi <- pheno.datTrait$yMu + pheno.datTrait$ePhen
dim(pheno.dat)

hist(pheno.datTrait$yMu)
  
all.data <- list(yTraiti = trt.dat$yTraiti,
                 N = Ntrt,
                 n_spec = Nspp,
                 trait_species = as.numeric(as.factor(trt.dat$species)),
                 n_tran = Ntran,
                 lati = trt.dat$latZ,
                 tranE = as.numeric(trt.dat$dumE),
                 trait_transect = as.numeric(as.factor(trt.dat$tran)),
                 Nph = nrow(pheno.datTrait),
                 phenology_species = as.numeric(as.factor(pheno.datTrait$species)),
                 species2 = as.numeric(as.factor(pheno.datTrait$species)),
                 yPhenoi = pheno.datTrait$yPhenoi,
                 #forcei = pheno.datTrait$forcei,
                 chilli = pheno.datTrait$chilli)
                # photoi = pheno.datTrait$photoi)


# mdl <- stan("stan/testTraitPhenoLatitudeTraitOnly.stan",
#                data = all.data,
#                iter = 4000, warmup = 3000, chains=4,
#             include = FALSE, pars = c("y_hat")
# )

mdl <- stan("stan/justDummyIntZ.stan",
            data = all.data,
            iter = 3000, warmup = 2000, chains=4,
            include = FALSE, pars = c("y_hat")
)
save(mdl, file="output/lmaTestDataZ.Rdata")
sumer <- summary(mdl)$summary

muTran <- sumer[grep("b_tran", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
muTranLat <- sumer[grep("b_tranlat", rownames(sumer)), c("mean","2.5%","97.5%")]

muForce <- sumer[grep("muForceSp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
muChill <- sumer[grep("muChillSp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
muPhoto <- sumer[grep("muPhotoSp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]

sigma_sp <- sumer[grep("sigma_sp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
sigma_traity <- sumer[grep("sigma_traity", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]

sigma_force <- sumer[grep("sigmaForceSp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
sigma_chill <- sumer[grep("sigmaChillSp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
sigma_photo <- sumer[grep("sigmaPhotoSp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]

sigma_phenoy <- sumer[grep("sigmapheno_y", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]

betaFT <- sumer[grep("betaTraitxForce", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
betaCT <- sumer[grep("betaTraitxChill", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
betaPT <- sumer[grep("betaTraitxPhoto", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]

# mdl.out <- data.frame( "Parameter" = c("mu.tran","sigma.tran","sigma.species","trt.var"),
#                        "Test.data.values" = c( mu.tranE, sigma.tran, sigma.species, trt.var) ,
#                        "Estiamte"= c(muTran[1], sigma_tran[1], sigma_sp[1], sigma_traity[1]),
#                        "2.5"= c(muTran[2], sigma_tran[2], sigma_sp[2], sigma_traity[1]),
#                        "97.5"= c( muTran[3],sigma_tran[3], sigma_sp[3], sigma_traity[3]) )
# mdl.out

mdl.out <- data.frame( "Parameter" = c("mutran","mu_forcesp","mu_chillsp","mu_photosp","sigma_tran","sigma_forcesp","sigma_chillsp","sigma_photosp", "sigma_phenoy","betaF", "betaC", "betaP"),  
                       "Test.data.values" = c( mu.tranE, mu.force, mu.chill, mu.photo,sigma.tranE, sigma.force, sigma.chill, sigma.photo,sigma.pheno.y,betaTraitxForce,betaTraitxChill,betaTraitxPhoto) ,
                       "Estiamte"= c(muTran[1], muForce[1], muChill[1], muPhoto[1],sigma_tran[1], sigma_force[1], sigma_chill[1], sigma_photo[1], sigma_phenoy[1], betaFT[1], betaCT[1],betaPT[1]),
                       "2.5"= c(muTran[2],muForce[2], muChill[2], muPhoto[2], sigma_tran[2],sigma_force[2], sigma_chill[2], sigma_photo[2], sigma_phenoy[2],betaFT[2], betaCT[2],betaPT[2]),
                       "55"= c( muTran[3],muForce[3], muChill[3], muPhoto[3],sigma_tran[3], sigma_force[3], sigma_chill[3], sigma_photo[3],  sigma_phenoy[3], betaFT[3], betaCT[3],betaPT[3]), "50"= c(muTran[4],muForce[4], muChill[4], muPhoto[4], sigma_tran[4],sigma_force[4], sigma_chill[4], sigma_photo[4], sigma_phenoy[4],betaFT[4], betaCT[4],betaPT[4]),
                       "75"= c( muTran[5],muForce[5], muChill[5], muPhoto[5],sigma_tran[5], sigma_force[5], sigma_chill[5], sigma_photo[5],  sigma_phenoy[5], betaFT[5], betaCT[5],betaPT[5]),
                       "97.5"= c(muTran[6],muForce[6], muChill[6], muPhoto[6], sigma_tran[6],sigma_force[6], sigma_chill[6], sigma_photo[6], sigma_phenoy[6],betaFT[6], betaCT[6],betaPT[6]))

mdl.out

postLMA<- data.frame(rstan::extract(mdl))

pdf("betaTraitChillPostPrior.pdf")
hist(postLMA$betaTraitxChill, main ="betaTraitxChill", col=rgb(0,0,1,1/4),  xlim = c(-50,50))
hist(rnorm(1000, 0,1), col=rgb(1,0,1,1/4), add = T)
abline(v =0, col="red", lwd=3, lty=2)
dev.off()

muTraitSp <- sumer[grep("muSp", rownames(sumer))]
pdf("latitudeMdlEstivsSimNoGrand.pdf")
plot(muTraitSp ~ mu.trtsp, xlab = "simulated muTraitSp", ylab = "mdl estimated muTraitSp")
abline(0,1)
dev.off()


# muPhenoSp <- sumer[grep("alphaPhenoSp", rownames(sumer))]
# muTraitSp <- sumer[grep("muSp", rownames(sumer))]
# muStudyEsti <- sumer[grep("muStudy", rownames(sumer))]
# muGrandSp <- sumer[grep("mu_grand_sp", rownames(sumer))]

# pdf("GeoffsMdl_5_5_comparisons.pdf",width = 6, height = 6)
# par(mfrow = c(2,2))
# plot(muTraitSp ~ alphaTraitSp, xlab = "simulated muTraitSp", ylab = "mdl estimated muTraitSp")
# abline(0,1)
# plot(muStudyEsti ~ muStudy, xlab = "simulated muStudy", ylab = "mdl estimated muStudy")
# abline(0,1)
# plot(muPhenoSp ~ alphaPhenoSp, xlab = "simulated muPhenoSp", ylab = "mdl estimated muPhenoSp")
# abline(0,1)
# plot(muGrandSp ~ unique(trt.dat$mu_grand_sp), xlab = "simulated muGrandSp", ylab = "mdl estimated muGrandSp")
# abline(0,1)
# dev.off()
# 
# # More complex model runs:
# mdl <- stan("stan/phenoBc_mdl_simp.stan",
#             data = pheno.data,
#             iter = 4000, warmup = 3000, chains=4,
#             include = FALSE, pars = c("y_hat")
# )
# mu_grand <- sumer[grep("mu_grand", rownames(sumer)),c("mean","2.5%","97.5%")]
# sigma_sp <- sumer[grep("sigma_sp", rownames(sumer)), c("mean","2.5%","97.5%")]
# popTwo <- sumer[grep("b_pop2", rownames(sumer)), c("mean","2.5%","97.5%")]
# popThree <- sumer[grep("b_pop3", rownames(sumer)), c("mean","2.5%","97.5%")]
# popFour <- sumer[grep("b_pop4", rownames(sumer)), c("mean","2.5%","97.5%")]
# sigma_traity <- sumer[grep("sigma_traity", rownames(sumer)), c("mean","2.5%","97.5%")]

# mu_chillsp <- sumer[grep("muChillSp", rownames(sumer)), c("mean","2.5%","97.5%")]
# sigma_chillsp <- sumer[grep("sigmaChillSp", rownames(sumer)), c("mean","2.5%","97.5%")]
# beta_tc <- sumer[grep("betaTraitxChill", rownames(sumer)), c("mean","2.5%","97.5%")]
# 
# mu_photosp <- sumer[grep("muPhotoSp", rownames(sumer)),c("mean","2.5%","97.5%")]
# sigma_photosp <- sumer[grep("sigmaPhotoSp", rownames(sumer)),c("mean","2.5%","97.5%")]
# beta_tp <- sumer[grep("betaTraitxPhoto", rownames(sumer)), c("mean","2.5%","97.5%")]
# 
# mu_forcesp <- sumer[grep("muForceSp", rownames(sumer)), c("mean","2.5%","97.5%")]
# mu_phenosp <- sumer[grep("muPhenoSp", rownames(sumer)), c("mean","2.5%","97.5%")]
# alpha.forcingsp <- sumer[grep("alphaForcingSp", rownames(sumer)), c("mean","2.5%","97.5%")]
# sigma_forcesp <- sumer[grep("sigmaForceSp", rownames(sumer)), c("mean","2.5%","97.5%")]
# sigma_phenosp <- sumer[grep("sigmaPhenoSp", rownames(sumer)), c("mean","2.5%","97.5%")]
# sigma_phenoy <- sumer[grep("sigmapheno_y", rownames(sumer)), c("mean","2.5%","97.5%")]
# beta_tf <- sumer[grep("betaTraitxForce", rownames(sumer)),c("mean","2.5%","97.5%")]
# beta_tc <- sumer[grep("betaTraitxChill", rownames(sumer)),c("mean","2.5%","97.5%")]
# beta_tp <- sumer[grep("betaTraitxPhoto", rownames(sumer)),c("mean","2.5%","97.5%")]
# 
# mdl.out <- data.frame( "Parameter" = c("mu_grand","sigma_sp","pop2","pop3","pop4","sigma_traity","mu_forcesp","mu_chillsp","mu_photosp","mu_phenosp","sigma_forcesp","sigma_chillsp","sigma_photosp", "sigma_phenosp", "sigma_phenoy", "beta_tf", "beta_tc","beta_tp"),  
#                        "Test.data.values" = c(mu.grand, sigma.species, mu.pop2,mu.pop3,mu.pop4, trt.var, mu.force, mu.chill, mu.photo, mu.pheno.sp, sigma.force, sigma.chill, sigma.photo, sigma.pheno.sp,sigma.pheno.y, betaTraitxForce, betaTraitxChill, betaTraitxPhoto) ,
#                        "Estiamte"= c(mu_grand[1,1], sigma_sp[1], popTwo[1], popThree[1], popFour[1], sigma_traity[1],  mu_forcesp[1], mu_chillsp[1], mu_photosp[1], mu_phenosp[1], sigma_forcesp[1], sigma_chillsp[1], sigma_photosp[1], sigma_phenosp[1], sigma_phenoy[1], beta_tf[1], beta_tc[1],  beta_tp[1]),
#                        "2.5"= c(mu_grand[1,2], sigma_sp[2], popTwo[2], popThree[2], popFour[2], sigma_traity[2],  mu_forcesp[2], mu_chillsp[2], mu_photosp[2], mu_phenosp[2], sigma_forcesp[2], sigma_chillsp[2], sigma_photosp[2], sigma_phenosp[2], sigma_phenoy[2], beta_tf[2], beta_tc[2],  beta_tp[2]),
#                        "97.5"= c(mu_grand[1,3], sigma_sp[3], popTwo[3], popThree[3], popFour[3], sigma_traity[3],  mu_forcesp[3], mu_chillsp[3], mu_photosp[3], mu_phenosp[3], sigma_forcesp[3], sigma_chillsp[3], sigma_photosp[3], sigma_phenosp[3], sigma_phenoy[3], beta_tf[3], beta_tc[3],  beta_tp[3]) )
# 
# mdl.out
# 
# muPhenoSp <- sumer[grep("alphaPhenoSp", rownames(sumer))]
# muTraitSp <- sumer[grep("muSp", rownames(sumer))]
# muStudyEsti <- sumer[grep("muStudy", rownames(sumer))]
# muGrandSp <- sumer[grep("mu_grand_sp", rownames(sumer))]
# 
# pdf("GeoffsMdl_5_5_comparisons.pdf",width = 6, height = 6)
# par(mfrow = c(2,2))
# plot(muTraitSp ~ alphaTraitSp, xlab = "simulated muTraitSp", ylab = "mdl estimated muTraitSp")
# abline(0,1)
# plot(muStudyEsti ~ muStudy, xlab = "simulated muStudy", ylab = "mdl estimated muStudy")
# abline(0,1)
# plot(muPhenoSp ~ alphaPhenoSp, xlab = "simulated muPhenoSp", ylab = "mdl estimated muPhenoSp")
# abline(0,1)
# plot(muGrandSp ~ unique(trt.dat$mu_grand_sp), xlab = "simulated muGrandSp", ylab = "mdl estimated muGrandSp")
# abline(0,1)
# dev.off()


# trt.out  <- data.frame( 
#   "Parameter" = c("mu_grand","sigma_sp","pop2","pop3","pop4","sigma_traity"),  
#   "Test.data.values" = c(mu.grand, sigma.species, mu.pop2,mu.pop3,mu.pop4, trt.var) ,
#   "Estiamte"= c(mu_grand[1,1], sigma_sp[1], popTwo[1], popThree[1], popFour[1], sigma_traity[1]),
#   "2.5"= c(mu_grand[1,2], sigma_sp[2], popTwo[2], popThree[2], popFour[2], sigma_traity[2]),
#   "97.5"= c(mu_grand[1,3], sigma_sp[3], popTwo[3], popThree[3], popFour[3], sigma_traity[3])
#   
# )
# trt.out
# 
# 
# muSpEsti <- sumer[grep("muSp\\[", rownames(sumer))]
# 
# plot(muSpEsti ~ mu.trtsp, xlab = "simulated muSpecies", ylab = "mdl estimated muSpecies")
# abline(0,1)
