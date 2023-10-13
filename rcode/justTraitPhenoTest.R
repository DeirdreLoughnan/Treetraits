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

# muSp + dummy + just interaction
Nrep <- 10# rep per trait
Npop <- 8
Ntran <- 2
Nspp <- 10# number of species with traits (making this 20 just for speed for now)

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

#####################################################################
##### Phenology test data ###########################

Nchill <- 2 # sm high low, mp high low
Nphoto <- 2 # high low
Nforce <- 2 # high amd low

Nrep <- 10
Nspp <- 10
Npop <- 2
Nph <- Nchill*Nphoto*Nforce*Nrep*Npop*Nspp #*Nforce - bc N force and N chill are completely confounded

pheno.dat <- data.frame(matrix(NA, Nph, 2))
names(pheno.dat) <- c("rep","species")
pheno.dat$rep <- c(1:Nrep)
pheno.dat$species <- rep(c(1:Nspp), each = Nrep)
pheno.dat$pop <- rep(1:Npop, each = Nspp*Nrep)

forceTrt <- c("LF", "HF")
forcei <- rnorm(Nforce, 1, 2)
pheno.dat$force <- rep(forceTrt, each = Nchill*Nphoto*Npop*Nrep*Nspp)
pheno.dat$forcei <- rep(forcei, each = Nchill*Nphoto*Npop*Nrep*Nspp)

chillTrt <- c("LC", "HC")
chilli <- rnorm(Nchill, 1, 5)
pheno.dat$chill <- rep(rep(chillTrt, each = Nphoto*Npop*Nrep*Nspp), times = Nforce)
pheno.dat$chilli <- rep(rep(chilli, each = Nphoto*Npop*Nrep*Nspp), times = Nforce)

photoTrt <- c("LP", "HP")
photoi <- rnorm(Nphoto, 1, 2)
pheno.dat$photo <- rep(rep(photoTrt, each = Npop*Nrep*Nspp), times = Nforce*Nchill)
pheno.dat$photoi <- rep(rep(photoi, each = Npop*Nrep*Nspp), times = Nforce*Nchill)

# sanity check that all treatments are the same
pheno.dat$label <- paste(pheno.dat$chill, pheno.dat$force, pheno.dat$photo, sep ="_")
pheno.dat$label <- paste(pheno.dat$chilli, pheno.dat$forcei, pheno.dat$photoi, sep ="_")
temp <- subset(pheno.dat, label == "HC_HF_HP"); dim(temp)

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

pheno.datTrait <- merge(pheno.dat, unique(trt.dat[,c("species","mu.trtsp")]), by = "species")
#head(pheno.datTrait,50)

for (i in 1:Nph){
  pheno.datTrait$betaForceSp[i] <-  pheno.datTrait$alphaForceSp[i] + (betaTraitxForce *  pheno.datTrait$mu.trtsp[i])
  
  pheno.datTrait$betaPhotoSp[i]<- pheno.datTrait$alphaPhotoSp[i] + (betaTraitxPhoto*  pheno.datTrait$mu.trtsp[i])
  
  pheno.datTrait$betaChillSp[i] <-pheno.datTrait$alphaChillSp[i] + (betaTraitxChill* pheno.datTrait$mu.trtsp[i])
}

for (i in 1:Nph){
  pheno.datTrait$yMu[i] <-  pheno.datTrait$alphaPhenoSp[i] +  pheno.datTrait$betaForceSp[i] * pheno.datTrait$forcei[i] +  pheno.datTrait$betaPhotoSp[i] * pheno.datTrait$photoi[i] + pheno.datTrait$betaChillSp[i] * pheno.datTrait$chilli[i]
}

pheno.datTrait$yPhenoi <- pheno.datTrait$yMu + pheno.datTrait$ePhen
dim(pheno.dat)

###############################################################
all.data <- list(yTraiti = trt.dat$yTraiti,
                 N = Ntrt,
                 n_spec = Nspp,
                 trait_species = as.numeric(as.factor(trt.dat$species)),
                 lati = trt.dat$lat,
                 tranE = as.numeric(trt.dat$dumE),
                 Nph = nrow(pheno.datTrait),
                 phenology_species = as.numeric(as.factor(pheno.datTrait$species)),
                 species2 = as.numeric(as.factor(pheno.datTrait$species)),
                 yPhenoi = pheno.datTrait$yPhenoi,
                 forcei = pheno.datTrait$forcei,
                 chilli = pheno.datTrait$chilli,
                 photoi = pheno.datTrait$photoi)


mdl <- stan("stan/justDummyInt.stan",
            data = all.data,
            iter = 4000, warmup = 3000, chains=4,
            include = FALSE, pars = c("y_hat")
)

save(mdl, file="output/testLatitude.Rdata")

sumer <- summary(mdl)$summary

muTraitSp <- sumer[grep("b_muSp", rownames(sumer))]

pdf("speciesEstiVsSim.pdf", width = 5, height = 5)
plot(muTraitSp ~ mu.trtsp, xlab = "simulated muTraitSp", ylab = "mdl estimated muTraitSp")
abline(0,1)
dev.off()

bTranE <- sumer[grep("b_tranE", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
bTranLat <- sumer[grep("b_tranlat", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
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


mdl.out <- data.frame( "Parameter" = c("bTranE","bTranLat","mu_forcesp","mu_chillsp","mu_photosp","sigma_sp","sigma_traity","sigma_forcesp","sigma_chillsp","sigma_photosp", "sigma_phenoy","betaF", "betaC", "betaP"),  
                       "Test.data.values" = c( mu.tranE, mu.tranlat, mu.force, mu.chill, mu.photo, sigma.species, trt.var, sigma.force, sigma.chill, sigma.photo,sigma.pheno.y,betaTraitxForce,betaTraitxChill,betaTraitxPhoto) ,
                       "Estiamte"= c(bTranE[1], bTranLat[1], muForce[1], muChill[1], muPhoto[1],  sigma_sp[1], sigma_traity[1], sigma_force[1], sigma_chill[1], sigma_photo[1], sigma_phenoy[1], betaFT[1], betaCT[1],betaPT[1]),
                       "2.5"= c(bTranE[2], bTranLat[2], muForce[2], muChill[2], muPhoto[2],  sigma_sp[2], sigma_traity[2],sigma_force[2], sigma_chill[2], sigma_photo[2], sigma_phenoy[2],betaFT[2], betaCT[2],betaPT[2]),
                       "25"= c( bTranE[3], bTranLat[3], muForce[3], muChill[3], muPhoto[3],  sigma_sp[3], sigma_traity[3], sigma_force[3], sigma_chill[3], sigma_photo[3],  sigma_phenoy[3], betaFT[3], betaCT[3],betaPT[3]),
                       "50"= c( bTranE[4], bTranLat[4], muForce[4], muChill[4], muPhoto[4],  sigma_sp[4], sigma_traity[4], sigma_force[4], sigma_chill[4], sigma_photo[4],  sigma_phenoy[4], betaFT[4], betaCT[4],betaPT[4]),
                       "75"= c( bTranE[5], bTranLat[5], muForce[5], muChill[5], muPhoto[5],  sigma_sp[5], sigma_traity[5], sigma_force[5], sigma_chill[5], sigma_photo[5],  sigma_phenoy[5], betaFT[5], betaCT[5],betaPT[5]),
                       "97.5"= c(bTranE[6], bTranLat[6], muForce[6], muChill[6], muPhoto[6],  sigma_sp[6], sigma_traity[6], sigma_force[6], sigma_chill[6], sigma_photo[6],  sigma_phenoy[6], betaFT[6], betaCT[6],betaPT[6]))

mdl.out



# post <- rstan::extract(mdl)
# 
# plot(post$b_tranE~post$b_tranlat)
# plot(post$muSp[,1]~post$sigma_sp)
# 
# table(trt.dat$pop, trt.dat$tranE)
