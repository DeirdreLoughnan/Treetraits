#Running a posterior predictive check on the fit traitors model during teh November 2021 retreit
#The plot is to check the relationship bewteen model predicted day of year and the intercept day of year for each species (alphaPhenoSp)

rm(list = ls())
dd
## Load libraries
library(rstan)
library(ggplot2)
library(reshape2)
library(viridis)
library(bayesplot)
library(tidybayes)
library(gridExtra) # for arranging plots 
library(patchwork) # another way of arranging plots 


#set the traits we are assessing
#this code will only work if you have the model outputs saved locally 

if(length(grep("deirdreloughnan", getwd())>0)) {  setwd("~/Documents/github/pheno_bc")
} else if
(length(grep("Lizzie", getwd())>0)) {   setwd("~/Documents/git/projects/treegarden/budreview/ospree/analyses/traits")
}


#Make a dataframe for saving traiit estimates for results section
traits <- c("height","lma", "cn", "dbh","ssd")
traitsDF <- data.frame(matrix(NA, 4,18))
# names(traitsDF) <- c("Trait", "GrandMean", "GrandMean_upper", "GrandMean_lower", 
#   "SpeciesSigma",  "SpeciesSigma_upper", "SpeciesSigma_lower", 
#   "StudySigma",  "StudySigma_upper", "StudySigma_lower", 
#   "MaxValue",  "MaxValue_upper", "MaxValue_lower", "MaxValueSp", 
#   "MinValue", "MinValueSp", "MinValue_upper", "MinValue_lower")
# traitsDF$Trait <- traits

#Load Trait Data
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

pheno$site.n <- pheno$population
pheno$site.n[pheno$site.n == "sm"] <- "1"
pheno$site.n[pheno$site.n == "mp"] <- "2"
pheno$site.n[pheno$site.n == "HF"] <- "3"
pheno$site.n[pheno$site.n == "SH"] <- "4"
pheno$site.n <- as.numeric(pheno$site.n)

head(pheno)
#add dummy/ site level effects:
pheno <- pheno %>%
  mutate ( site2 = if_else(site.n == 2, 1, 0),
    site3 = if_else(site.n == 3, 1, 0),
    site4 = if_else(site.n == 4, 1, 0))

# standardize the 0/1 and standardize sites? 
pheno$force.z2 <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE)*2)
pheno$photo.z2 <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE)*2)
pheno$chillport.z2 <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE)*2)

pheno$site2.z2 <- (pheno$site2-mean(pheno$site2,na.rm=TRUE))/(sd(pheno$site2,na.rm=TRUE)*2)
pheno$site3.z2 <- (pheno$site3-mean(pheno$site3,na.rm=TRUE))/(sd(pheno$site3,na.rm=TRUE)*2)
pheno$site4.z2 <- (pheno$site4-mean(pheno$site4,na.rm=TRUE))/(sd(pheno$site4,na.rm=TRUE)*2)

#going to split it into analysis of terminal bb and lateral bb
# Starting with the terminal buds:
#pheno.term <- pheno[,c("tbb", "chill.n", "force.n", "photo.n", "site.n", "species", "lab2")]
pheno.term <- pheno[,c("bb", "force.z2", "photo.z2", "population", "species", "lab2","Utah_Model","Chill_portions","chillport.z2", "site2.z2", "site3.z2","site4.z2")]
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

head(trtPheno)

trtPheno$transect <- trtPheno$site
trtPheno$transect[trtPheno$transect == "GR"] <- "east"
trtPheno$transect[trtPheno$transect == "WM"] <- "east"
trtPheno$transect[trtPheno$transect == "SH"] <- "east"
trtPheno$transect[trtPheno$transect == "HF"] <- "east"

trtPheno$transect[trtPheno$transect == "af"] <- "west"
trtPheno$transect[trtPheno$transect == "kl"] <- "west"
trtPheno$transect[trtPheno$transect == "sm"] <- "west"
trtPheno$transect[trtPheno$transect == "mp"] <- "west"

filePathData <- "output/"
traitModelNames <- grep(".Rdata", list.files(filePathData), value = TRUE) 


traitPlotList <- list() #Make a list to save trait plots so we can make 1 pannel with all 4 plots at the end 

  
  traiti <- 4
  
  #Load lma model fit
  lmaModel <- load(paste(filePathData,traitModelNames[traiti], sep = "/"))
  traitName <- gsub(".Rdata", "", traitModelNames[traiti])
  lmaModelFit <- rstan::extract(mdlLMA)
  
  #sensible cue values
  # #-------------------------------------
  # forcingValue <- 0.85 # 20 degrees C
  # chillinValue <- 1 #coudl go up to 2 or 3 
  # photoValue <- -0.25 # about 12 Or 0.5(about 16)
  
  #Extracting  postreior values 
  #----------------------------------
  
  #meanInterceptValues
  # alphaPhenoSpdf <- data.frame(lmaModelFit$alphaPhenoSp)
  # alphaPhenoSpMean <- colMeans(alphaPhenoSpdf)
  # 
  # #Forcing slope values 
  # betaForceSpdf <- data.frame(lmaModelFit$betaForceSp)
  # betaForceSpMean <- colMeans(betaForceSpdf)
  # 
  # #Chilling slope values 
  # betaChillSpdf <- data.frame(lmaModelFit$betaChillSp)
  # betaChillSpMean <- colMeans(betaChillSpdf)
  # 
  
  #Photoperiod slope values 
  # betaPhotoSpdf <- data.frame(lmaModelFit$betaPhotoSp)
  # betaPhotoSpMean <- colMeans(betaPhotoSpdf)
  # 
  # #Overall model 
  # sigmapheno_yMean <- mean(lmaModelFit$sigmapheno_y)
  # 
  # #Predict DOY based on model estimated parameters, One DOY value per species
  # yPhenoi <- vector()
  # 
  # for(ip in 1:length(betaPhotoSpMean)){
  #   yPhenoi[ip] <- alphaPhenoSpMean[ip] + betaForceSpMean[ip] * forcingValue + betaPhotoSpMean[ip] * photoValue + betaChillSpMean[ip]* chillinValue
  # }
  
  #Trait data 
  #--------------

    specieslist <- sort(unique(trtPheno$species))
    trtPheno2 <- trtPheno[,2:14]
    trtPhenoLong <- melt(trtPheno2, id = c("species", "site", "transect"))
    trtPhenoLong <- trtPhenoLong[complete.cases(trtPhenoLong$value),]
    trtPhenoLong$value <- as.numeric(trtPhenoLong$value)
   
    meanRealTrait <- aggregate(trtPhenoLong["value"], trtPhenoLong[c("species","variable")], FUN = mean)
    names(meanRealTrait) <- c("species","trait","meanTrait")

    meanRealLMA <- subset(meanRealTrait, trait == "lma")
  
  
  #Trait predicted means
  #--------------------
  
  mu_Df <- data.frame(lmaModelFit$b_muSp)
  colnames(mu_Df) <- specieslist
  longMeans <- melt(mu_Df)
  colnames(longMeans) <- c("speciesname", "traitMean")
  
  mu_mean <- colMeans(mu_Df)
  
  color_scheme_set("viridis")
  lmaData <- subset(trtPhenoLong, variable == "lma")

    
    traitFit <- ggplot(data = lmaData, aes(y = species, x = value, colour = "black"))+
      stat_eye(data = longMeans, aes(y = speciesname, x = traitMean))+
      geom_point( alpha = 0.5, size = 1.2, aes(colour = "red"))+
      theme_classic() +  
      xlim (0,0.2) +
      scale_y_discrete(limits=rev) +
      theme(text = element_text(size=16))+
      geom_point(data = meanRealLMA, aes(x = meanTrait,y = species, colour = "purple"), shape = 8, size = 3)+
      labs(title = "Leaf mass area", y = "Species", x ="Trait Value")+ 
      scale_color_identity(name = "Model fit",
        breaks = c("black", "red", "purple"),
        labels = c("Model Posterior", "Raw Data", "Data Mean"),
        guide = guide_legend(override.aes = list(
          linetype = c(NA, NA, NA),
          shape = c(19, 20, 8)))) + 
      theme(legend.title = element_blank())
  
  
  #######################################################################
    traiti <- 3
    htModel <- load(paste(filePathData,traitModelNames[traiti], sep = "/"))
    traitName <- gsub(".Rdata", "", traitModelNames[traiti])
    htModelFit <- rstan::extract(mdl)
    
    meanRealHt <- subset(meanRealTrait, trait == "ht")
    
        #Trait predicted means
    #--------------------
    
    mu_Df <- data.frame(htModelFit$b_muSp)
    colnames(mu_Df) <- specieslist
    longMeans <- melt(mu_Df)
    colnames(longMeans) <- c("speciesname", "traitMean")
    
    mu_mean <- colMeans(mu_Df)
    
    color_scheme_set("viridis")
    htData <- subset(trtPhenoLong, variable == "ht")

    
   ggplot(data = htData, aes(y = species, x = value, colour = "black"))+
      stat_eye(data = longMeans, aes(y = speciesname, x = traitMean))+
      geom_point( alpha = 0.5, size = 1.2, aes(colour = "red"))+
      theme_classic() +  
      scale_y_discrete(limits=rev) +
      theme(text = element_text(size=16))+
      geom_point(data = meanRealHt, aes(x = meanTrait,y = species, colour = "purple"), shape = 8, size = 3)+
      labs(title = "Height", y = "Species", x ="Trait Value")+ 
      scale_color_identity(name = "Model fit",
        breaks = c("black", "red", "purple"),
        labels = c("Model Posterior", "Raw Data", "Data Mean"),
        guide = guide_legend(override.aes = list(
          linetype = c(NA, NA, NA),
          shape = c(19, 20, 8)))) + 
      theme(legend.title = element_blank())
   
   #######################################################################
   traiti <- 1
   cnModel <- load(paste(filePathData,traitModelNames[traiti], sep = "/"))
   traitName <- gsub(".Rdata", "", traitModelNames[traiti])
   cnModelFit <- rstan::extract(mdl)
   
   meanRealCN <- subset(meanRealTrait, trait == "C.N")
   
   #Trait predicted means
   #--------------------
   
   mu_Df <- data.frame(cnModelFit$b_muSp)
   colnames(mu_Df) <- specieslist
   longMeans <- melt(mu_Df)
   colnames(longMeans) <- c("speciesname", "traitMean")
   
   mu_mean <- colMeans(mu_Df)
   
   color_scheme_set("viridis")
   cnData <- subset(trtPhenoLong, variable == "C.N")
   
   ggplot(data = cnData, aes(y = species, x = (value/100), colour = "black"))+
     stat_eye(data = longMeans, aes(y = speciesname, x = traitMean))+
     geom_point( alpha = 0.5, size = 1.2, aes(colour = "red"))+
     theme_classic() +  
     scale_y_discrete(limits=rev) +
     theme(text = element_text(size=16))+
     geom_point(data = meanRealCN, aes(x = (meanTrait/100),y = species, colour = "purple"), shape = 8, size = 3)+
     labs(title = "Carbon:Nitrogen", y = "Species", x ="Trait Value")+ 
     scale_color_identity(name = "Model fit",
       breaks = c("black", "red", "purple"),
       labels = c("Model Posterior", "Raw Data", "Data Mean"),
       guide = guide_legend(override.aes = list(
         linetype = c(NA, NA, NA),
         shape = c(19, 20, 8)))) + 
     theme(legend.title = element_blank())
  
   #######################################################################
   traiti <- 2
   dbhModel <- load(paste(filePathData,traitModelNames[traiti], sep = "/"))
   traitName <- gsub(".Rdata", "", traitModelNames[traiti])
   dbhModelFit <- rstan::extract(mdl)
   
   meanRealdbh <- subset(meanRealTrait, trait == "dbh")
   
   #Trait predicted means
   #--------------------
   
   mu_Df <- data.frame(dbhModelFit$b_muSp)
   colnames(mu_Df) <- specieslist
   longMeans <- melt(mu_Df)
   colnames(longMeans) <- c("speciesname", "traitMean")
   
   mu_mean <- colMeans(mu_Df)
   
   color_scheme_set("viridis")
   dbhData <- subset(trtPhenoLong, variable == "dbh")
   
   ggplot(data = dbhData, aes(y = species, x = (value), colour = "black"))+
     stat_eye(data = longMeans, aes(y = speciesname, x = traitMean))+
     geom_point( alpha = 0.5, size = 1.2, aes(colour = "red"))+
     theme_classic() +  
     scale_y_discrete(limits=rev) +
     theme(text = element_text(size=16))+
     geom_point(data = meanRealdbh, aes(x = (meanTrait),y = species, colour = "purple"), shape = 8, size = 3)+
     labs(title = "Diameter at breast height", y = "Species", x ="Trait Value")+ 
     scale_color_identity(name = "Model fit",
       breaks = c("black", "red", "purple"),
       labels = c("Model Posterior", "Raw Data", "Data Mean"),
       guide = guide_legend(override.aes = list(
         linetype = c(NA, NA, NA),
         shape = c(19, 20, 8)))) + 
     theme(legend.title = element_blank())
   
   #######################################################################
   traiti <- 5
   ssdModel <- load(paste(filePathData,traitModelNames[traiti], sep = "/"))
   traitName <- gsub(".Rdata", "", traitModelNames[traiti])
   ssdModelFit <- rstan::extract(mdl)
   
   meanRealssd <- subset(meanRealTrait, trait == "ssd")
   
   #Trait predicted means
   #--------------------
   
   mu_Df <- data.frame(ssdModelFit$b_muSp)
   colnames(mu_Df) <- specieslist
   longMeans <- melt(mu_Df)
   colnames(longMeans) <- c("speciesname", "traitMean")
   
   mu_mean <- colMeans(mu_Df)
   
   color_scheme_set("viridis")
   ssdData <- subset(trtPhenoLong, variable == "ssd")
   
   ggplot(data = ssdData, aes(y = species, x = (value), colour = "black"))+
     stat_eye(data = longMeans, aes(y = speciesname, x = traitMean))+
     geom_point( alpha = 0.5, size = 1.2, aes(colour = "red"))+
     theme_classic() +  
     scale_y_discrete(limits=rev) +
     theme(text = element_text(size=16))+
     geom_point(data = meanRealssd, aes(x = (meanTrait),y = species, colour = "purple"), shape = 8, size = 3)+
     labs(title = "Stem specific density", y = "Species", x ="Trait Value")+ 
     scale_color_identity(name = "Model fit",
       breaks = c("black", "red", "purple"),
       labels = c("Model Posterior", "Raw Data", "Data Mean"),
       guide = guide_legend(override.aes = list(
         linetype = c(NA, NA, NA),
         shape = c(19, 20, 8)))) + 
     theme(legend.title = element_blank())
   
  #ggsave(paste(paste0("figures/traitFit_", traitName), "png", sep = "."), traitFit,   width = 10, height = 6,  units = "in")
#   traitPlotList[[traiti]] <-traitFit 
#   
#   #saving model values in a table
#   
#   traitOutput <- data.frame(summary(lmaModel))
#   speciesMeans <-  aggregate(longMeans$traitMean, by = list(longMeans$speciesname), FUN = mean)
#   names(speciesMeans) <- c("speciesname", "traitMean")
#   
#   #Get Trait Model output for the results section. 
#   
#   
#   
#   #Mean trait value
#   traitsDF$GrandMean[traiti] <- mean(lmaModelFit$mu_grand)
#   #Uncertainty around mu grand
#   traitsDF$GrandMean_upper[traiti] <- HPDI( as.vector(lmaModelFit$mu_grand) , prob=0.90 )[2]
#   traitsDF$GrandMean_lower[traiti] <- HPDI( as.vector(lmaModelFit$mu_grand) , prob=0.90 )[1]
#   
#   #speciesSigma
#   traitsDF$SpeciesSigma[traiti] <- mean(lmaModelFit$sigma_sp)
#   #Uncertainty around speciesSigma
#   traitsDF$SpeciesSigma_upper[traiti] <- HPDI( as.vector(lmaModelFit$sigma_sp) , prob=0.90 )[2]
#   traitsDF$SpeciesSigma_lower[traiti] <- HPDI( as.vector(lmaModelFit$sigma_sp) , prob=0.90 )[1]
#   
#   #studySigma
#   traitsDF$StudySigma[traiti] <- mean(lmaModelFit$sigma_study)
#   #Uncertainty around studySigma
#   traitsDF$StudySigma_upper[traiti] <- HPDI( as.vector(lmaModelFit$sigma_study) , prob=0.90 )[2]
#   traitsDF$StudySigma_lower[traiti] <- HPDI( as.vector(lmaModelFit$sigma_study) , prob=0.90 )[1]
#   
#   
#   
#   #Max trait value
#   traitsDF$MaxValue[traiti] <- max(speciesMeans$traitMean)
#   #max trait value species id
#   traitsDF$MaxValueSp[traiti] <- as.character(speciesMeans$speciesname[speciesMeans$traitMean == max(speciesMeans$traitMean)])
#   #Uncertainty around max species
#   traitsDF$MaxValue_upper[traiti] <- HPDI( as.vector(longMeans$traitMean[longMeans$speciesname == traitsDF$MaxValueSp[traiti]]) , prob=0.90 )[2]
#   traitsDF$MaxValue_lower[traiti] <- HPDI( as.vector(longMeans$traitMean[longMeans$speciesname == traitsDF$MaxValueSp[traiti]]) , prob=0.90 )[1]
#   
#   
#   
#   #Min trait value
#   traitsDF$MinValue[traiti] <- min(speciesMeans$traitMean)
#   #min trait value species id
#   traitsDF$MinValueSp[traiti] <- as.character(speciesMeans$speciesname[speciesMeans$traitMean == min(speciesMeans$traitMean)])
#   #Uncertainty around max species
#   traitsDF$MinValue_upper[traiti] <- HPDI( as.vector(longMeans$traitMean[longMeans$speciesname == traitsDF$MinValueSp[traiti]]) , prob=0.90 )[2]
#   traitsDF$MinValue_lower[traiti] <- HPDI( as.vector(longMeans$traitMean[longMeans$speciesname == traitsDF$MinValueSp[traiti]]) , prob=0.90 )[1]
#   
#   
#   
#   
#   
#   
#   
# }
# 
# 
# 
# #this plotting code needs the patchwork library 
# png("figures/FourTraitFit_37spp_wp.png", width = 14, height = 15, units = "in", res = 72)
# combined <- traitPlotList[[1]] + traitPlotList[[4]] + traitPlotList[[3]]+ traitPlotList[[2]]  & theme(legend.position = "bottom") # combien plots and put legend at teh bottom
# combined[[2]] <- combined[[2]] + theme(axis.title.y = element_blank() )#Remove y labels from plots 2 and 4
# combined[[4]] <- combined[[4]] + theme(axis.title.y = element_blank() )
# combined + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a")#add letter annotation to plots 
# dev.off()