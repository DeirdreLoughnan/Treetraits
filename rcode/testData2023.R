# Started Aug 22 2023: DL
#the aim of this test data is to get the model working for my BC traits chapter

# the model: a joint model with trait effects then modeled with phenological cues  ---similar to traitors model but with some key differences
#1. site as a summy varaible---not study---only one intercept to fit

# first run: 7 div transitions, rhat 1.14, 3210 trans>max treedepth
setwd("~/Documents/github/Treetraits") 

library(rstan)
require(rstanarm)
require(shinystan)
require(bayesplot)

rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

Nrep <- 5 # rep per trait
Npop <- 4 
Nspp <- 20 # number of species with traits (making this 20 just for speed for now)

# First making a data frame for the test trait data
Ntrt <- Nspp * Npop * Nrep # total number of traits observations
Ntrt

#make a dataframe for height
trt.dat <- data.frame(matrix(NA, Ntrt, 1))
names(trt.dat) <- c("rep")
trt.dat$rep <- c(1:Nrep)
trt.dat$species <- rep(1:Nspp, each = Nrep)
trt.dat$pop <- rep(1:Npop, each = Nspp*Nrep)
  
trt.dat$pop2 <- trt.dat$pop
trt.dat$pop3 <- trt.dat$pop
trt.dat$pop4 <- trt.dat$pop

trt.dat$pop2 <- ifelse(trt.dat$pop == "2", 1,0)
trt.dat$pop3 <- ifelse(trt.dat$pop == "3", 1,0)
trt.dat$pop4 <- ifelse(trt.dat$pop == "4", 1,0)

mu.pop2 <- 2
mu.pop3 <- 3
mu.pop4 <- 4

trt.dat$mu.pop2 <- mu.pop2*trt.dat$pop2 
trt.dat$mu.pop3 <-  mu.pop3*trt.dat$pop3 
trt.dat$mu.pop4 <-  mu.pop4*trt.dat$pop4 


# trt.dat$species <- rep(1:Nspp, Nstudy)

# now generating the species trait data, here it is for height
mu.grand <- 10 # the grand mean of the height model
sigma.species <- 5 # we want to keep the variaiton across spp. high

#the alphaTraitSp in Faiths original code:
mu.trtsp <- rnorm(Nspp, 0, sigma.species)
trt.dat$mu.trtsp <- rep(mu.trtsp, each = Nrep) #adding ht data for ea. sp

# general variance
trt.var <- 1 #sigma_traity in the stan code
trt.dat$trt.er <- rnorm(Ntrt, 0, trt.var)

# generate yhat - heights -  for this first trt model
for (i in 1:Ntrt){
  trt.dat$yTraiti[i] <-  mu.grand + trt.dat$mu.trtsp[i] + trt.dat$mu.pop2[i] + trt.dat$mu.pop3[i] + trt.dat$mu.pop4[i] + trt.dat$trt.er[i]
}

for (i in 1:Ntrt){
  trt.dat$muGrandSp[i] <-  trt.dat$mu.trtsp[i] +  mu.grand
}

##### Phenology test data ###########################

Nchill <- 4 # sm high low, mp high low
Nphoto <- 2 # high low
Nforce <- 2 # high amd low

Nrep <- 5
Nspp <- 20
Npop <- 2
Nph <- Nchill*Nphoto*Nrep*Npop*Nspp #*Nforce - bc N force and N chill are completely confounded

pheno.dat <- data.frame(matrix(NA, Nph, 2))
names(pheno.dat) <- c("rep","species")
pheno.dat$rep <- c(1:Nrep)
pheno.dat$species <- rep(c(1:Nspp), each = Nrep)
pheno.dat$pop <- rep(1:Npop, each = Nspp*Nrep)

forcei <- rnorm(Nforce, 1, 2)
pheno.dat$forcei <- rep(forcei, each = Npop*Nrep*Nspp)

chilli <- rnorm(Nchill, 1, 5)
pheno.dat$chilli <- rep(chilli, each = Nrep*Nspp)

photoi <- rnorm(Nphoto, 1, 5)
pheno.dat$photoi <- rep(photoi, each = Npop*Nrep*Nspp)

mu.force = -10 
sigma.force = 1
alpha.force.sp <- rnorm(Nspp, mu.force, sigma.force)
pheno.dat$alphaForceSp <- rep(alpha.force.sp, each = Nrep)

mu.photo = -15
sigma.photo = 1
alpha.photo.sp <- rnorm(Nspp, mu.photo, sigma.photo)
pheno.dat$alphaPhotoSp <- rep(alpha.photo.sp, each = Nrep)

mu.chill = -14
sigma.chill = 1
alpha.chill.sp <- rnorm(Nspp, mu.chill, sigma.chill)
pheno.dat$alphaChillSp <- rep(alpha.chill.sp, each = Nrep)

mu.pheno.sp = 80
sigma.pheno.sp = 30
alphaPhenoSp <- rnorm(Nspp, mu.pheno.sp, sigma.pheno.sp)
pheno.dat$alphaPhenoSp <- rep(alphaPhenoSp, each = Nrep)

sigma.pheno.y = 3
pheno.dat$ePhen <- rnorm(Nph, 0, sigma.pheno.y)

betaTraitxForce <- 0.3 
betaTraitxPhoto <- -0.2
betaTraitxChill <- -0.4

pheno.datTrait <- merge(pheno.dat, unique(trt.dat[,c("species","muGrandSp")]), by = "species")
head(pheno.datTrait,50)

for (i in 1:Nph){
  pheno.datTrait$betaForceSp[i] <-  pheno.datTrait$alphaForceSp[i] + (betaTraitxForce *  pheno.datTrait$muGrandSp[i])
  
  pheno.datTrait$betaPhotoSp[i]<- pheno.datTrait$alphaPhotoSp[i] + (betaTraitxPhoto*  pheno.datTrait$muGrandSp[i])
  
  pheno.datTrait$betaChillSp[i] <-pheno.datTrait$alphaChillSp[i] + (betaTraitxChill* pheno.datTrait$muGrandSp[i])
}

for (i in 1:Nph){
  pheno.datTrait$yMu[i] <-  pheno.datTrait$alphaPhenoSp[i] +  pheno.datTrait$betaForceSp[i] * pheno.datTrait$forcei[i] +  pheno.datTrait$betaPhotoSp[i] * pheno.datTrait$photoi[i] + pheno.datTrait$betaChillSp[i] * pheno.datTrait$chilli[i]
}

pheno.datTrait$yPhenoi <- pheno.datTrait$yMu + pheno.datTrait$ePhen


all.data <- list(yTraiti = trt.dat$yTraiti,
                 N = Ntrt,
                 n_spec = Nspp,
                 trait_species = as.numeric(as.factor(trt.dat$species)),
                 species = as.numeric(as.factor(trt.dat$species)),
                 n_pop = Npop,
                 pop2 = trt.dat$pop2,
                 pop3 = trt.dat$pop3,
                 pop4 = trt.dat$pop4,
                 Nph = nrow(pheno.datTrait),
                 phenology_species = as.numeric(as.factor(pheno.datTrait$species)),
                 species2 = as.numeric(as.factor(pheno.datTrait$species)),
                 yPhenoi = pheno.datTrait$yPhenoi,
                 forcei = pheno.datTrait$forcei,
                 chilli = pheno.datTrait$chilli,
                 photoi = pheno.datTrait$photoi)
                 
mdl <- stan("stan/testTraitPhenoDummy.stan",
               data = all.data,
               iter = 4000, warmup = 3000, chains=4,
            include = FALSE, pars = c("y_hat")
)
save(mdl, file="output/testTraitPheno_Nrep5Nspp20Npop4.Rdata")
sumer <- summary(mdl)$summary

mu_grand <- sumer[grep("mu_grand", rownames(sumer)),c("mean","2.5%","97.5%")]
sigma_sp <- sumer[grep("sigma_sp", rownames(sumer)), c("mean","2.5%","97.5%")]
popTwo <- sumer[grep("b_pop2", rownames(sumer)), c("mean","2.5%","97.5%")]
popThree <- sumer[grep("b_pop3", rownames(sumer)), c("mean","2.5%","97.5%")]
popFour <- sumer[grep("b_pop4", rownames(sumer)), c("mean","2.5%","97.5%")]
sigma_traity <- sumer[grep("sigma_traity", rownames(sumer)), c("mean","2.5%","97.5%")]

mu_chillsp <- sumer[grep("muChillSp", rownames(sumer)), c("mean","2.5%","97.5%")]
sigma_chillsp <- sumer[grep("sigmaChillSp", rownames(sumer)), c("mean","2.5%","97.5%")]
beta_tc <- sumer[grep("betaTraitxChill", rownames(sumer)), c("mean","2.5%","97.5%")]

mu_photosp <- sumer[grep("muPhotoSp", rownames(sumer)),c("mean","2.5%","97.5%")]
sigma_photosp <- sumer[grep("sigmaPhotoSp", rownames(sumer)),c("mean","2.5%","97.5%")]
beta_tp <- sumer[grep("betaTraitxPhoto", rownames(sumer)), c("mean","2.5%","97.5%")]

mu_forcesp <- sumer[grep("muForceSp", rownames(sumer)), c("mean","2.5%","97.5%")]
mu_phenosp <- sumer[grep("muPhenoSp", rownames(sumer)), c("mean","2.5%","97.5%")]
alpha.forcingsp <- sumer[grep("alphaForcingSp", rownames(sumer)), c("mean","2.5%","97.5%")]
sigma_forcesp <- sumer[grep("sigmaForceSp", rownames(sumer)), c("mean","2.5%","97.5%")]
sigma_phenosp <- sumer[grep("sigmaPhenoSp", rownames(sumer)), c("mean","2.5%","97.5%")]
sigma_phenoy <- sumer[grep("sigmapheno_y", rownames(sumer)), c("mean","2.5%","97.5%")]
beta_tf <- sumer[grep("betaTraitxForce", rownames(sumer)),c("mean","2.5%","97.5%")]
beta_tc <- sumer[grep("betaTraitxChill", rownames(sumer)),c("mean","2.5%","97.5%")]
beta_tp <- sumer[grep("betaTraitxPhoto", rownames(sumer)),c("mean","2.5%","97.5%")]

mdl.out <- data.frame( "Parameter" = c("mu_grand","sigma_sp","pop2","pop3","pop4","sigma_traity","mu_forcesp","mu_chillsp","mu_photosp","mu_phenosp","sigma_forcesp","sigma_chillsp","sigma_photosp", "sigma_phenosp", "sigma_phenoy", "beta_tf", "beta_tc","beta_tp"),  
                       "Test.data.values" = c(mu.grand, sigma.species, mu.pop2,mu.pop3,mu.pop4, trt.var, mu.force, mu.chill, mu.photo, mu.pheno.sp, sigma.force, sigma.chill, sigma.photo, sigma.pheno.sp,sigma.pheno.y, betaTraitxForce, betaTraitxChill, betaTraitxPhoto) ,
                       "Estiamte"= c(mu_grand[1,1], sigma_sp[1], popTwo[1], popThree[1], popFour[1], sigma_traity[1],  mu_forcesp[1], mu_chillsp[1], mu_photosp[1], mu_phenosp[1], sigma_forcesp[1], sigma_chillsp[1], sigma_photosp[1], sigma_phenosp[1], sigma_phenoy[1], beta_tf[1], beta_tc[1],  beta_tp[1]),
                       "2.5"= c(mu_grand[1,2], sigma_sp[2], popTwo[2], popThree[2], popFour[2], sigma_traity[2],  mu_forcesp[2], mu_chillsp[2], mu_photosp[2], mu_phenosp[2], sigma_forcesp[2], sigma_chillsp[2], sigma_photosp[2], sigma_phenosp[2], sigma_phenoy[2], beta_tf[2], beta_tc[2],  beta_tp[2]),
                       "97.5"= c(mu_grand[1,3], sigma_sp[3], popTwo[3], popThree[3], popFour[3], sigma_traity[3],  mu_forcesp[3], mu_chillsp[3], mu_photosp[3], mu_phenosp[3], sigma_forcesp[3], sigma_chillsp[3], sigma_photosp[3], sigma_phenosp[3], sigma_phenoy[3], beta_tf[3], beta_tc[3],  beta_tp[3]) )

mdl.out

muPhenoSp <- sumer[grep("alphaPhenoSp", rownames(sumer))]
muTraitSp <- sumer[grep("muSp", rownames(sumer))]
muStudyEsti <- sumer[grep("muStudy", rownames(sumer))]
muGrandSp <- sumer[grep("mu_grand_sp", rownames(sumer))]

pdf("GeoffsMdl_5_5_comparisons.pdf",width = 6, height = 6)
par(mfrow = c(2,2))
plot(muTraitSp ~ alphaTraitSp, xlab = "simulated muTraitSp", ylab = "mdl estimated muTraitSp")
abline(0,1)
plot(muStudyEsti ~ muStudy, xlab = "simulated muStudy", ylab = "mdl estimated muStudy")
abline(0,1)
plot(muPhenoSp ~ alphaPhenoSp, xlab = "simulated muPhenoSp", ylab = "mdl estimated muPhenoSp")
abline(0,1)
plot(muGrandSp ~ unique(trt.dat$mu_grand_sp), xlab = "simulated muGrandSp", ylab = "mdl estimated muGrandSp")
abline(0,1)
dev.off()


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
