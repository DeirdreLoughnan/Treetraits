# Started December 15, 2024: DL
#the aim of this test data is to double check that we are confident in the SSD values

rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

require(rstan)
#require(truncnorm)

setwd("~/Documents/github/Treetraits")

Nrep <- 10# rep per trait
Npop <- 8
Ntran <- 2
Nspp <- 80# number of species with traits (making this 20 just for speed for now)

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

mu.tranE <- 5
trt.dat$dumE <- as.numeric(ifelse(trt.dat$tran == "1","0","1"))
trt.dat$mutranE <- mu.tranE*trt.dat$dumE

mu.tranlat  <- 2
trt.dat$alpha.tranlat <- mu.tranlat*(trt.dat$dumE*trt.dat$latZ)

sigma.species <- 10 # we want to keep the variation across spp. high

mu.trtsp <- rnorm(Nspp, 50, sigma.species)
trt.dat$mu.trtsp <- rep(mu.trtsp, each = Nrep) #adding ht data for ea. sp
#hist(mu.trtsp)
# general variance
trt.var <- 5 #sigma_traity in the stan code
trt.dat$trt.er <- rnorm(Ntrt, 0, trt.var)

# generate yhat - heights -  for this first trt model
for (i in 1:Ntrt){
  trt.dat$yTraiti[i] <- 
    trt.dat$mu.trtsp[i] + trt.dat$trt.er[i] + trt.dat$mutranE[i] + mu.tranlat * (trt.dat$dumE[i]*trt.dat$latZ[i])  
}

all.data <- list(yTraiti = trt.dat$yTraiti,
  N = Ntrt,
  n_spec = Nspp,
  trait_species = as.numeric(as.factor(trt.dat$species)),
  lati = trt.dat$latZ,
  tranE = as.numeric(trt.dat$dumE)
)


# mdl <- stan("stan/modelDevelopment/justDummyIntTrait.stan",
#   data = all.data,
#   iter = 4000, warmup = 3000, chains=4,
#   include = FALSE, pars = c("y_hat")
# )

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

# 20 reps, 40 spp
# Parameter Test.data.values   Estiamte       X2.5        X25        X50        X75     X97.5
# 1       bTranE              5.0 4.99576221 4.97106368 4.98668986 4.99569333 5.00464410 5.0206824
# 2     bTranLat              2.0 1.99877138 1.98263847 1.99282010 1.99871669 2.00444280 2.0157245
# 3     sigma_sp              0.1 0.09051132 0.06862914 0.08174978 0.08946797 0.09844331 0.1175569
# 4 sigma_traity              0.5 0.50371058 0.49502363 0.50045952 0.50361544 0.50684106 0.5126634

##### Phenology test data ###########################

Nchill <- 8 # sm high low, mp high low
Nphoto <- 2 # high low
Nforce <- 2 # high amd low

# Nrep <- 20
# Nspp <- 15
Npop <- 2
Nph <- Nrep*Npop*Nspp*Nphoto*Nforce*Nphoto #bc N force and N chill are completely confounded

pheno.dat <- data.frame(matrix(NA, Nph, 2))
names(pheno.dat) <- c("rep","species")
pheno.dat$rep <- c(1:Nrep)
pheno.dat$species <- rep(c(1:Nspp), each = Nrep)
pheno.dat$pop <- rep(1:Npop, each = Nspp*Nrep)

# Growth chamber study, so the cue values are constants - unique level chilling for each pop
chillTrt <- c("LC", "HC")
chilli <- rnorm(Nchill, 0, 1)

# pheno.dat$chill <- rep(rep(chillTrt, each = Npop*Nrep*Nspp))
# pheno.dat$chilli <- rep(rep(chilli, each = Npop*Nrep*Nspp))

pheno.dat$chill <- rep(rep(chillTrt, each = Nspp*Nrep))
pheno.dat$chilli <- rep(rep(chilli, each = Nspp*Nrep))

# 
forceTrt <- c("LF", "HF") 
forcei <- c(-0.6688487,  1.2154703)
pheno.dat$force <- rep(forceTrt, each = Nphoto*Npop*Nrep*Nspp)
pheno.dat$forcei <- rep(forcei, each = Nphoto*Npop*Nrep*Nspp)

photoTrt <- c("LP", "HP") # 8 or 12 h
photoi <- c(0.9916049, -1.0082267)
pheno.dat$photo <- rep(rep(photoTrt, each = Npop*Nrep*Nspp), times = Nforce)
pheno.dat$photoi <- rep(rep(photoi, each = Npop*Nrep*Nspp), times = Nforce)

# sanity check that all treatments are the same
# pheno.dat$label <- paste(pheno.dat$chill, pheno.dat$force, pheno.dat$photo, sep ="_")
# pheno.dat$label <- paste(pheno.dat$chilli, pheno.dat$forcei, pheno.dat$photoi, sep ="_")

mu.chill = -14
sigma.chill = 1
alpha.chill.sp <- rnorm(Nspp, mu.chill, sigma.chill)
pheno.dat$alphaChillSp <- rep(alpha.chill.sp, each = Nrep)

mu.force = -10
sigma.force = 1
alpha.force.sp <- rnorm(Nspp, mu.force, sigma.force)
pheno.dat$alphaForceSp <- rep(alpha.force.sp, each = Nrep)

mu.photo = -5
sigma.photo = 1
alpha.photo.sp <- rnorm(Nspp, mu.photo, sigma.photo)
pheno.dat$alphaPhotoSp <- rep(alpha.photo.sp, each = Nrep)

mu.pheno.sp = 80
sigma.pheno.sp = 10
alphaPhenoSp <- rnorm(Nspp, mu.pheno.sp, sigma.pheno.sp)
pheno.dat$alphaPhenoSp <- rep(alphaPhenoSp, each = Nrep)

sigma.pheno.y = 3
pheno.dat$ePhen <- rnorm(Nph, 0, sigma.pheno.y)

betaTraitxForce <- 3
betaTraitxPhoto <- -2
betaTraitxChill <- -4

pheno.datTrait <- merge(pheno.dat, unique(trt.dat[,c("species","mu.trtsp")]), by = "species")
#head(pheno.datTrait,50)

for (i in 1:Nph){
  pheno.datTrait$betaForceSp[i] <-  pheno.datTrait$alphaForceSp[i] + (betaTraitxForce *  pheno.datTrait$mu.trtsp[i])

  pheno.datTrait$betaPhotoSp[i]<- pheno.datTrait$alphaPhotoSp[i] + (betaTraitxPhoto*  pheno.datTrait$mu.trtsp[i])

  pheno.datTrait$betaChillSp[i] <-pheno.datTrait$alphaChillSp[i] + (betaTraitxChill* pheno.datTrait$mu.trtsp[i])
}

for (i in 1:Nph){
  pheno.datTrait$yMu[i] <-  pheno.datTrait$alphaPhenoSp[i] + pheno.datTrait$betaChillSp[i] * pheno.datTrait$chilli[i] +  pheno.datTrait$betaForceSp[i] * pheno.datTrait$forcei[i] +  pheno.datTrait$betaPhotoSp[i] * pheno.datTrait$photoi[i] 
}

pheno.datTrait$yPhenoi <- pheno.datTrait$yMu + pheno.datTrait$ePhen
dim(pheno.dat)

hist(pheno.datTrait$yMu) # these values look reasonable

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
                 forcei = pheno.datTrait$forcei,
                 chilli = pheno.datTrait$chilli,
                 photoi = pheno.datTrait$photoi)

mdl <- stan("stan/modelDevelopment/justDummyIntZ.stan",
            data = all.data,
            iter = 4000, warmup = 3000, chains=4,
            include = FALSE, pars = c("y_hat")
)

save(mdl, file="analysis/output/ssdTestDataZ.Rdata")
sumer <- summary(mdl)$summary

bTranE <- sumer[grep("b_tranE", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
bTranLat <- sumer[grep("b_tranlat", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]

muForce <- sumer[grep("muForceSp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
muChill <- sumer[grep("muChillSp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
muPhoto <- sumer[grep("muPhotoSp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
muPheno <- sumer[grep("muPhenoSp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]

sigma_sp <- sumer[grep("sigma_sp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
sigma_traity <- sumer[grep("sigma_traity", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]

sigma_force <- sumer[grep("sigmaForceSp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
sigma_chill <- sumer[grep("sigmaChillSp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
sigma_photo <- sumer[grep("sigmaPhotoSp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]

sigma_phenoy <- sumer[grep("sigmapheno_y", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]

betaFT <- sumer[grep("betaTraitxForce", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
betaCT <- sumer[grep("betaTraitxChill", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
betaPT <- sumer[grep("betaTraitxPhoto", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]

mdl.out <- data.frame( "Parameter" = 
    c("mutran","mutranLat","mu_forcesp","mu_chillsp","mu_photosp","mu_phenosp",
      "sigma_traity", "sigma_sp","sigma_forcesp","sigma_chillsp","sigma_photosp", 
      "sigma_phenoy","betaF", "betaC", "betaP"),  
                       "Test.data.values" = c( mu.tranE, mu.tranlat ,mu.force, mu.chill, mu.photo, mu.pheno.sp, trt.var, sigma.species, sigma.force, sigma.chill, sigma.photo,sigma.pheno.y,betaTraitxForce,betaTraitxChill,betaTraitxPhoto) ,
                       "Estiamte"= c(bTranE[1], bTranLat[1], muForce[1], muChill[1], muPhoto[1], muPheno[1], sigma_traity[1], sigma_sp[1], sigma_force[1], sigma_chill[1], sigma_photo[1], sigma_phenoy[1], betaFT[1], betaCT[1],betaPT[1]),
                       "2.5"= c(bTranE[2],bTranLat[2],muForce[2], muChill[2], muPhoto[2],muPheno[2], sigma_traity[2], sigma_sp[2], sigma_force[2], sigma_chill[2], sigma_photo[2], sigma_phenoy[2],betaFT[2], betaCT[2],betaPT[2]),
                       "55"= c( bTranE[3],bTranLat[3],muForce[3], muChill[3], muPhoto[3],muPheno[3], sigma_traity[3], sigma_sp[3], sigma_force[3], sigma_chill[3], sigma_photo[3],  sigma_phenoy[3], betaFT[3], betaCT[3],betaPT[3]), "50"= c(bTranE[4],bTranLat[4],muForce[4], muChill[4], muPhoto[4],muPheno[4], sigma_traity[4], sigma_sp[4], sigma_force[4], sigma_chill[4], sigma_photo[4], sigma_phenoy[4],betaFT[4], betaCT[4],betaPT[4]),
                       "75"= c( bTranE[5],bTranLat[5],muForce[5], muChill[5], muPhoto[5],muPheno[5], sigma_traity[5], sigma_sp[5], sigma_force[5], sigma_chill[5], sigma_photo[5],  sigma_phenoy[5], betaFT[5], betaCT[5],betaPT[5]),
                       "97.5"= c(bTranE[6],bTranLat[6],muForce[6], muChill[6], muPhoto[6],muPheno[6], sigma_traity[6], sigma_sp[6],sigma_force[6], sigma_chill[6], sigma_photo[6], sigma_phenoy[6],betaFT[6], betaCT[6],betaPT[6]))

mdl.out

# Parameter Test.data.values    Estiamte        X2.5         X55         X50         X75      X97.5
# 1         mutran                5   4.7180702   4.1419345   4.5146800   4.7135841   4.9221097   5.296150
# 2      mutranLat                2   2.2290689   1.5452195   1.9988924   2.2325799   2.4692511   2.892048
# 3     mu_forcesp              -10 -10.0168913 -12.1073613 -10.7438022 -10.0506820  -9.2684616  -7.945515
# 4     mu_chillsp              -14 -14.5300854 -16.9042080 -15.3608747 -14.5675223 -13.6647878 -12.109585
# 5     mu_photosp               -5  -5.4247022  -6.8634576  -5.9516105  -5.4564071  -4.8874566  -3.895768
# 6     mu_phenosp               80  80.3587164  77.9590554  79.5359570  80.3723926  81.1714601  82.732046
# 7   sigma_traity                5   4.9130902   4.8270385   4.8829401   4.9132615   4.9420710   5.002455
# 8       sigma_sp               10  10.7699891   9.2335413  10.1295446  10.7084796  11.3305787  12.674339
# 9  sigma_forcesp                1   1.0349833   0.7823935   0.9384598   1.0308520   1.1249996   1.309711
# 10 sigma_chillsp                1   0.8102007   0.3266061   0.6802449   0.8198023   0.9574059   1.204739
# 11 sigma_photosp                1   0.8482172   0.6871503   0.7842336   0.8428351   0.9076939   1.038829
# 12  sigma_phenoy                3   3.0212577   2.9838130   3.0084186   3.0208992   3.0343502   3.058598
# 13         betaF                3   3.0042213   2.9635365   2.9907854   3.0046815   3.0178851   3.044582
# 14         betaC               -4  -3.9976865  -4.0426522  -4.0135536  -3.9975950  -3.9818976  -3.950773
# 15         betaP               -2  -1.9979220  -2.0258830  -2.0081868  -1.9974009  -1.9876068  -1.971152

# postLMA<- data.frame(rstan::extract(mdl))
# 
# pdf("betaTraitChillPostPrior.pdf")
# hist(postLMA$betaTraitxForce, main ="betaTraitxChill", col=rgb(0,0,1,1/4),  xlim = c(-5,5))
# hist(rnorm(1000, 0,1), col=rgb(1,0,1,1/4), add = T)
# abline(v =0, col="red", lwd=3, lty=2)
# dev.off()
# 
# muTraitSp <- sumer[grep("muSp", rownames(sumer))]
# pdf("latitudeMdlEstivsSimNoGrand.pdf")
# plot(muTraitSp ~ mu.trtsp, xlab = "simulated muTraitSp", ylab = "mdl estimated muTraitSp")
# abline(0,1)
# dev.off()
# 
# 
