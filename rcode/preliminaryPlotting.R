# Started Feb 4, 2022

# The purpose of this code is to plot the trait values for each species against the the slope estimates from the phenology model:

# Get the model output from pheno_bc repo:
rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(rstan)
library(shinystan)
#library(reshape2)
library(bayesplot)
library(ggplot2)
#library(RColorBrewer)
library(dplyr)
library(plyr)


if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  else{
  setwd("/home/deirdre/pheno_bc") # for midge
}

# get the output from my data
load("output/tbb_ncp_chillportions_zsc_dl.Rda")

sumt <- summary(mdl.t)$summary
bforce_dl <- sumt[grep("b_force", rownames(sumt)), "mean"]; bforce_dl
bchill_dl <- sumt[grep("b_chill", rownames(sumt)), "mean"]; bchill_dl
bphoto_dl <- sumt[grep("b_photo", rownames(sumt)), "mean"]; bphoto_dl <- bphoto_dl[20:38]

# get the output from my data
load("output/termMdlEstiContChillStndDFlynn.Rda")

sumdf <- summary(mdl.std)$summary
bforce_df <- sumdf[grep("b_force", rownames(sumdf)), "mean"]; bforce_df
bchill_df <- sumdf[grep("b_chill", rownames(sumdf)), "mean"]; bchill_df
bphoto_df <- sumdf[grep("b_photo", rownames(sumdf)), "mean"]; bphoto_df <- bphoto_df[29:56]


##################################################
setwd("..//Treetraits")
trtData <- read.csv("data/allTrt.csv")
head(trtData)

unique(trtData$C.N)

spp <- c("acegla", "acepen", "acerub", "acesac", "alninc","alnvir", "amealn", "aromel", "betall", "betlen", "betpap",
         "corcor", "corsto", "faggra","franig", "hamvir", "ilemuc", "kalang", "loncan", "loninv", "lyolig", "menfer", "nyssyl", "popbal","popgra", "poptre", "prupen", "quealb", "querub", "quevel", "rhafra", "rhoalb", "rhopri", "riblac", "rubpar", "samrac", "shecan","sorsco", "spialb", "spibet", "spipyr", "symalb", "vacmem", "vacmyr", "vibcas", "vibedu", "viblan")

dlspp <- c("acegla", "alninc","alnvir", "amealn", "betpap",
           "corsto", "loninv", "popbal", "poptre", 
           "riblac", "rubpar", "samrac", "shecan",
           "sorsco", "spibet", "spipyr", "symalb",
           "vacmem", "vibedu")

dfspp <- c("acepen", "acerub", "acesac", "alninc","aromel", "betall", "betlen", "betpap",
           "corcor","faggra","franig", "hamvir", "ilemuc", "kalang", "loncan", "lyolig",  
           "nyssyl", "popgra",  "prupen", "quealb", "querub", "quevel", "rhafra", "rhopri",
           "spialb", "vacmyr", "vibcas", "viblan")



#trtDataSpp <- trtData[trtData$species %in% dfspp,]
trtDataSpp <- trtData
trtDataSpp$species.fact <- as.numeric(as.factor(trtDataSpp$species))

# now run simple trait model to get trait effects for each species:

lma_data  <- trtDataSpp[complete.cases(trtDataSpp$lma),]
lma_data$species.fact <- as.numeric(as.factor(lma_data$species))

lma_datalist <- list(yTraiti = lma_data$lma, 
                   N = nrow(lma_data), 
                   n_spec = length(unique(lma_data$species)), 
                   species = lma_data$species.fact, 
                   prior_mu_grand_mu = 0.5,
                   prior_mu_grand_sigma = 1,
                   prior_sigma_sp_mu = 1,
                   prior_sigma_sp_sigma = 1,
                   prior_sigma_traity_mu = 1,
                   prior_sigma_traity_sigma = 1
) 


mdl.lma <- stan('stan/bc_trait_only_2.stan',
                  data = lma_datalist,
                  iter = 6000,
                  warmup = 4000,
                  chains = 4,
                  include = FALSE,
                  pars = c("mu_y","y_hat"))
save(mdl.lma, file = "output_lma_traitonly.Rda")

sum.lma <- summary(mdl.lma)$summary

lmaTrt_df <- sum.lma[grep("muSp", rownames(sum.lma)), "mean"]; lmaTrt_df
########################################################
ht_data  <- trtDataSpp[complete.cases(trtDataSpp$ht),]
ht_datalist <- list(yTraiti = ht_data$ht, 
                     N = nrow(ht_data), 
                     n_spec = length(unique(ht_data$species)), 
                     species = as.numeric(as.factor(ht_data$species)), 
                     prior_mu_grand_mu = 15,
                     prior_mu_grand_sigma = 10,
                     prior_sigma_sp_mu = 5,
                     prior_sigma_sp_sigma = 1,
                     prior_sigma_traity_mu = 5,
                     prior_sigma_traity_sigma = 1
) 

mdl.ht <- stan('stan/bc_trait_only_2.stan',
                data = ht_datalist,
                iter = 6000,
                warmup = 4000,
                chains = 4,
                include = FALSE,
                pars = "mu_y")
save(mdl.ht, file = "output_ht_traitonly_df.Rda")

sum.ht <- summary(mdl.ht)$summary

htTrt_df <- sum.ht[grep("muSp", rownames(sum.ht)), "mean"]; htTrt_df
######################################################
dbh_data  <- trtDataSpp[complete.cases(trtDataSpp$dbh),]
dbh_data$species.fact <- as.numeric(as.factor(dbh_data$species))

dbh_datalist <- list(yTraiti = dbh_data$dbh, 
                     N = nrow(dbh_data), 
                     n_spec = length(unique(dbh_data$species)), 
                     species = dbh_data$species.fact, 
                     prior_mu_grand_mu = 0.5,
                     prior_mu_grand_sigma = 1,
                     prior_sigma_sp_mu = 1,
                     prior_sigma_sp_sigma = 1,
                     prior_sigma_traity_mu = 1,
                     prior_sigma_traity_sigma = 1
) 


mdl.dbh <- stan('stan/bc_trait_only_2.stan',
                data = dbh_datalist,
                iter = 6000,
                warmup = 4000,
                chains = 4,
                include = FALSE,
                pars = "mu_y")
save(mdl.dbh, file = "output_dbh_traitonly_df.Rda")

sum.dbh <- summary(mdl.dbh)$summary

dbhTrt_df <- sum.dbh[grep("muSp", rownames(sum.dbh)), "mean"]; dbhTrt_df
########################################################
ssd_data  <- trtDataSpp[complete.cases(trtDataSpp$ssd),]
ssd_data$species.fact <- as.numeric(as.factor(ssd_data$species))

ssd_datalist <- list(yTraiti = ssd_data$ssd, 
                     N = nrow(ssd_data), 
                     n_spec = length(unique(ssd_data$species)), 
                     species = ssd_data$species.fact, 
                     prior_mu_grand_mu = 0.5,
                     prior_mu_grand_sigma = 1,
                     prior_sigma_sp_mu = 1,
                     prior_sigma_sp_sigma = 1,
                     prior_sigma_traity_mu = 1,
                     prior_sigma_traity_sigma = 1
) 


mdl.ssd <- stan('stan/bc_trait_only_2.stan',
                data = ssd_datalist,
                iter = 6000,
                warmup = 4000,
                chains = 4,
                include = FALSE,
                pars = "mu_y")
save(mdl.ssd, file = "output_ssd_traitonly_df.Rda")

sum.ssd <- summary(mdl.ssd)$summary

ssdTrt_df <- sum.ssd[grep("muSp", rownames(sum.ssd)), "mean"]; ssdTrt_df
########################################################
perCN_data  <- trtDataSpp[complete.cases(trtDataSpp$C.N),]
perCN_data$species.fact <- as.numeric(as.factor(perCN_data$species))

perCN_datalist <- list(yTraiti = perCN_data$C.N, 
                     N = nrow(perCN_data), 
                     n_spec = length(unique(perCN_data$species)), 
                     species = perCN_data$species.fact, 
                     prior_mu_grand_mu = 0.5,
                     prior_mu_grand_sigma = 1,
                     prior_sigma_sp_mu = 1,
                     prior_sigma_sp_sigma = 1,
                     prior_sigma_traity_mu = 1,
                     prior_sigma_traity_sigma = 1
) 


mdl.perCN <- stan('stan/bc_trait_only_2.stan',
                data = perCN_datalist,
                iter = 6000,
                warmup = 4000,
                chains = 4,
                include = FALSE,
                pars = "mu_y")
save(mdl.perCN, file = "output_perCN_traitonly_df.Rda")

sum.perCN <- summary(mdl.perCN)$summary
perCNTrt_df <- sum.perCN[grep("muSp", rownames(sum.perCN)), "mean"]; perCNTrt_df
########################################################

# get the phylogeny data:
load("output/dl_df_allbb_4sites.Rda")

sum.pheno <- summary(mdl.t)$summary

spp <- sort(unique(trtData$species))

sppPheno <- tolower(sppPheno)
bforce <- data.frame(sum.pheno[grep("^b_warm", rownames(sum.pheno)), "mean"])
bforce$spp <- sppPheno
colnames(bforce) <- c("muSp","species")

bchill <- data.frame(sum.pheno[grep("^b_chill", rownames(sum.pheno)), "mean"])
bchill$spp <- sppPheno
colnames(bchill) <- c("muSp","species")

bphoto <- data.frame(sum.pheno[grep("^b_photo", rownames(sum.pheno)), "mean"])
bphoto$spp <- sppPheno
colnames(bphoto) <- c("muSp","species")
bphoto<- bphoto[50:98,]

perCNTrt_df <- data.frame(sum.perCN[grep("muSp", rownames(sum.perCN)), "mean"])
perCNTrt_df$spp <- dfspp
colnames(perCNTrt_df) <- c("muSp","species")

perCNTrt_dl <- data.frame(sum.perCN.dl[grep("muSp", rownames(sum.perCN.dl)), "mean"])
perCNTrt_dl$spp <- dlspp
colnames(perCNTrt_dl) <- c("muSp","species")

trtOut <- rbind(perCNTrt_df, perCNTrt_dl)
trtOut <- aggregate(trtOut["muSp"], trtOut[c("species")], FUN = mean)

head(trtOut)
head(bphoto)

bforce <- subset(bforce, species != "menfer")
bforce <- subset(bforce, species != "rhoalb")

plot(bforce$muSp ~ trtOut$muSp, col = "maroon", pch =19)
dfPt <- bforce[bforce$species %in% dfspp, ]
dftrt <- trtOut[trtOut$species %in% dfspp, ]
points(dfPt$muSp ~ dftrt$muSp, col = "darkslategray4", pch =19)

bchill <- subset(bchill, species != "menfer")
bchill <- subset(bchill, species != "rhoalb")

plot(bchill$muSp ~ trtOut$muSp, col = "maroon", pch =19)
dfPt <- bchill[bchill$species %in% dfspp, ]
dftrt <- trtOut[trtOut$species %in% dfspp, ]
points(dfPt$muSp ~ dftrt$muSp, col = "darkslategray4", pch =19)

bphoto <- subset(bphoto, species != "menfer")
bphoto <- subset(bphoto, species != "rhoalb")

plot(bphoto$muSp ~ trtOut$muSp, col = "maroon", pch =19)
dfPt <- bphoto[bphoto$species %in% dfspp, ]
dftrt <- trtOut[trtOut$species %in% dfspp, ]
points(dfPt$muSp ~ dftrt$muSp, col = "darkslategray4", pch =19)

########################################################
# get dl mdl output
load("output/output_lma_traitonly.Rda")
sum.lma <- summary(mdl.lma)$summary
lmaTrt <- sum.lma[grep("muSp", rownames(sum.lma)), "mean"]; 

load("output/output_ht_traitonly.Rda")
sum.ht <- summary(mdl.ht)$summary
htTrt <- sum.ht[grep("muSp", rownames(sum.ht)), "mean"];

load("output/output_ssd_traitonly.Rda")
sum.ssd <- summary(mdl.ssd)$summary
ssdTrt <- sum.ssd[grep("muSp", rownames(sum.ssd)), "mean"];

load("output/output_dbh_traitonly.Rda")
sum.dbh <- summary(mdl.dbh)$summary
dbhTrt <- sum.dbh[grep("muSp", rownames(sum.dbh)), "mean"]

load("output/output_perCN_traitonly_dl.Rda")
sum.perCN.dl<- summary(mdl.perCN)$summary
perCNTrt_dl <- sum.perCN.dl[grep("muSp", rownames(sum.perCN.dl)), "mean"]; 

perCNTrt_df

pdf("figures/muSpvscue.pdf", width = 15, height = 25)
par(mfrow = c(6, 3))
plot(bforce_dl ~ ssdTrt, col = "maroon", pch =19, ylim = c(-15, 5), xlim = c(-0.3,0.3))
points(bforce_df ~ssdTrt_df, col = "darkslategray4", pch =19)

plot(bchill_dl ~ ssdTrt, col = "maroon", pch =19, ylim = c(-10, 30), xlim = c(-0.1, 0.2))
points(bchill_df ~ ssdTrt_df, col = "darkslategray4", pch =19)

plot(bphoto_dl ~ ssdTrt, col = "maroon", pch =19, ylim = c(-15, 0), xlim = c(-0.25, 0.2))
points(bphoto_df ~ ssdTrt_df, col = "darkslategray4", pch =19)

plot(bforce_dl ~ perCNTrt, col = "maroon", pch =19, ylim = c(-15, 4), xlim = c(37,47))
points(bforce_df ~ perCTrt_df, col = "darkslategray4", pch =19)

plot(bchill_dl ~ perCTrt, col = "maroon", pch =19, ylim = c(-15, 30), xlim = c(38, 53))
points(bchill_df ~ perCTrt_df, col = "darkslategray4", pch =19)

plot(bphoto_dl ~ perCTrt, col = "maroon", pch =19, ylim = c(-15, 0), xlim = c(38, 53))
points(bphoto_df ~ perCTrt_df, col = "darkslategray4", pch =19)

plot(bforce_dl ~ perNTrt, col = "maroon", pch =19, ylim = c(-15, 0), xlim = c(4,-1))
points(bforce_df ~ perNTrt_df, col = "darkslategray4", pch =19)

plot(bchill_dl ~ perNTrt, col = "maroon", pch =19, ylim = c(-15, 30), xlim = c(4,-1))
points(bchill_df ~ perNTrt_df, col = "darkslategray4", pch =19)

plot(bphoto_dl ~ perNTrt, col = "maroon", pch =19, ylim = c(-15, 1), xlim = c(4,-1))
points(bphoto_df ~ perNTrt_df, col = "darkslategray4", pch =19)

plot(bforce_dl ~ htTrt, col = "maroon", pch =19, ylim = c(-15, 5), xlim = c(-5,15))
points(bforce_df ~ htTrt_df, col = "darkslategray4", pch =19)

plot(bchill_dl ~ htTrt, col = "maroon", pch =19, ylim = c(-15, 35), xlim = c(-5,15))
points(bchill_df ~ htTrt_df, col = "darkslategray4", pch =19)

plot(bphoto_dl ~ htTrt, col = "maroon", pch =19, ylim = c(-15, 1), xlim = c(-5,15))
points(bphoto_df ~ htTrt_df, col = "darkslategray4", pch =19)

plot(bforce_dl ~ lmaTrt, col = "maroon", pch =19, ylim = c(-15, 5), xlim = c(-0.05,0.05))
points(bforce_df ~ lmaTrt_df, col = "darkslategray4", pch =19)

plot(bchill_dl ~ lmaTrt, col = "maroon", pch =19, ylim = c(-15, 30), xlim = c(-0.05,0.05))
points(bchill_df ~ lmaTrt_df, col = "darkslategray4", pch =19)

plot(bphoto_dl ~ lmaTrt, col = "maroon", pch =19)
plot(bphoto_df ~ lmaTrt_df, col = "darkslategray4", pch =19)


plot(bforce_dl ~ dbhTrt, col = "maroon", pch =19)
plot(bforce_df ~ dbhTrt_df, col = "darkslategray4", pch =19)

plot(bchill_dl ~ dbhTrt, col = "maroon", pch =19)
plot(bchill_df ~ dbhTrt_df, col = "darkslategray4", pch =19)

plot(bphoto_dl ~ dbhTrt, col = "maroon", pch =19)
plot(bphoto_df ~ dbhTrt_df, col = "darkslategray4", pch =19)

dev.off()

