# Started Aug 22 2023: DL
#the aim of this test data is to get the model working for my BC traits chapter

# the model: a joint model with trait effects then modeled with phenological cues  ---similar to traitors model but with some key differences
#1. site as a summy varaible---not study---only one intercept to fit
setwd("~/Documents/github/Treetraits") 

library(rstan)
require(rstanarm)
require(shinystan)
require(bayesplot)

rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

Nrep <- 8 # rep per trait
Npop <- 4 
Nspp <- 50 # number of species with traits (making this 20 just for speed for now)

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

#constant population effects

#now generating the effects of study
# sigma.pop2 <- 2
# sigma.pop3 <- 3
# sigma.pop4 <- 4
# 
# mu.pop1 <- rnorm(Npop, 0, sigma.pop1)
# mu.pop2 <- rnorm(Npop, 0, sigma.pop2)
# mu.pop3 <- rnorm(Npop, 0, sigma.pop3)
# mu.pop4 <- rnorm(Npop, 0, sigma.pop4)

#intercept for each study
#trt.dat$mu.pop1 <- rep(mu.pop1, each = Nspp) # generate data for ea study
# trt.dat$mu.pop2 <- rep(mu.pop2, each = Nspp) # generate data for ea study
# trt.dat$mu.pop3 <- rep(mu.pop3, each = Nspp) # generate data for ea study
# trt.dat$mu.pop4 <- rep(mu.pop4, each = Nspp) # generate data for ea study

# general variance
trt.var <- 1 #sigmaTrait_y in the stan code
trt.dat$trt.er <- rnorm(Ntrt, 0, trt.var)

# generate yhat - heights -  for this first trt model
for (i in 1:Ntrt){
  trt.dat$yTraiti[i] <-  mu.grand + trt.dat$mu.trtsp[i] + trt.dat$mu.pop2[i] + trt.dat$mu.pop3[i] + trt.dat$mu.pop4[i] + trt.dat$trt.er[i]
}

all.data <- list(yTraiti = trt.dat$yTraiti,
                 N = Ntrt,
                 n_spec = Nspp,
                 trait_species = as.numeric(as.factor(trt.dat$species)),
                 species = as.numeric(as.factor(trt.dat$species)),
                 n_pop = Npop,
                 pop2 = trt.dat$pop2,
                 pop3 = trt.dat$pop3,
                 pop4 = trt.dat$pop4)
                 
mdl <- stan("stan/testTraitOnlyDummy.stan",
               data = all.data,
               iter = 4000, warmup = 3000, chains=4,
            include = FALSE, pars = c("y_hat")
)
save(mdl, file="output/testNrep20Nspp50Npop4.Rdata")
sumer <- summary(mdl)$summary

mu_grand <- sumer[grep("mu_grand", rownames(sumer)),c("mean","2.5%","97.5%")]
sigma_sp <- sumer[grep("sigma_sp", rownames(sumer)), c("mean","2.5%","97.5%")]
popTwo <- sumer[grep("b_pop2", rownames(sumer)), c("mean","2.5%","97.5%")]
popThree <- sumer[grep("b_pop3", rownames(sumer)), c("mean","2.5%","97.5%")]
popFour <- sumer[grep("b_pop4", rownames(sumer)), c("mean","2.5%","97.5%")]
sigma_traity <- sumer[grep("sigma_traity", rownames(sumer)), c("mean","2.5%","97.5%")]

trt.out  <- data.frame( 
  "Parameter" = c("mu_grand","sigma_sp","pop2","pop3","pop4","sigma_traity"),  
  "Test.data.values" = c(mu.grand, sigma.species, mu.pop2,mu.pop3,mu.pop4, trt.var) ,
  "Estiamte"= c(mu_grand[1,1], sigma_sp[1], popTwo[1], popThree[1], popFour[1], sigma_traity[1]),
  "2.5"= c(mu_grand[1,2], sigma_sp[2], popTwo[2], popThree[2], popFour[2], sigma_traity[2]),
  "97.5"= c(mu_grand[1,3], sigma_sp[3], popTwo[3], popThree[3], popFour[3], sigma_traity[3])
  
)
trt.out


muSpEsti <- sumer[grep("muSp\\[", rownames(sumer))]

plot(muSpEsti ~ mu.trtsp, xlab = "simulated muSpecies", ylab = "mdl estimated muSpecies")
abline(0,1)
