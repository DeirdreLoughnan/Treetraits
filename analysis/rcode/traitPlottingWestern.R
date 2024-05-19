# started March 18, 2023 by Deirdre

# aim of this code is to make the plots for BC traits projects

#1. Figure comapring the site differences for each trait
#2. Figure comparing the betatraitCues for each figure

rm(list=ls())
options(stringsAsFactors = FALSE)

library(stringr)
library(plyr)
library(rstan)
library(rshape2)
library(cowplot)
require(dplyr)


if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/Treetraits") 
}  else{
  setwd("/home/deirdre/Treetraits") # for midge
}

load("output/rda/ht_western_Mar4.Rda")
sumer3 <- summary(mdl.ht3)$summary
postHtW <- rstan::extract(mdl.ht3)

load("output/rda/lma_western_Mar4.Rda")
postLMAW <- rstan::extract(mdl.lma3)

load("output/rda/dbh_western_Mar4.Rda")
postDbhW <- rstan::extract(mdl.dbh3)

load("output/rda/cn_western_Mar4.Rda")
postCNW <- rstan::extract(mdl.cn3)

load("output/rda/ssd_western_Mar4.Rda")
postSsdW <- rstan::extract(mdl.ssd2)

###### Compare the spp and site level effects across traits ###########
col1 <- rgb(204 / 255, 102 / 255, 119 / 255, alpha = 0.8)
col2 <- rgb(68 / 255, 170 / 255, 153 / 255, alpha = 0.6)

pdf("figures/florum_sigmaSpSite.pdf", width = 12, height = 4)
par(mfrow = c(1, 4))
hist(postHtW$sigma_site, main = "Height", xlab = "Variation", ylab = "Frequency", col = col1, breaks = 30, ylim = c(0, 4000))
hist(postHtW$sigma_sp, add = T, col = col2)
text(0.2,4000, label = "a)", cex = 2)
     
hist(postLMAW$sigma_site, main = "Leaf Mass Area", xlab = "Variation", ylab = "Frequency", col = col1, breaks = 50, xlim = c(0,0.2), ylim = c(0, 6000))
hist(postLMAW$sigma_sp, add = T, col = col2)
text(0.002,6000, label = "b)", cex = 2)

hist(postDbhW$sigma_site, main = "Diameter at Breast Height", xlab = "Variation", ylab = "Frequency", col = col1, ylim = c(0, 6000))
hist(postDbhW$sigma_sp, add = T, col = col2)
text(0.2,6000, label = "c)", cex = 2)

hist(postCNW$sigma_site, main = "Carbon:Nitrogen", xlab = "Variation", ylab = "Frequency", col = col1, breaks = 30, ylim = c(0, 4000), xlim = c(0,12))
hist(postCNW$sigma_sp, add = T, col = col2)
text(0.2,4000, label = "d)", cex = 2)

hist(postSsdW$sigma_site, main = "Stem Specific Density", xlab = "Variation", ylab = "Frequency", col = col1, breaks = 30, ylim = c(0, 4000), xlim = c(0,6))
hist(postCNW$sigma_sp, add = T, col = col2, breaks = 50)
text(0.1,4000, label = "e)", cex = 2)

legend("topright",legend = c(expression("Species variation"),
                             expression("Study variation")),
       col = c(col2,col1),
       inset = 0.02, cex = 1, bty = "n", lwd = 3)

dev.off()

# Compare the musites for the 4 sites
muSiteHtW <- data.frame(postHtW$musite)
site <- c("Alex Fraser", "Kamloops","Manning Park", "Smithers")
names(muSiteHtW) <- site

longSiteHtW <- melt(muSiteHtW)
names(longSiteHtW) <- c("Population", "value")

htSiteW <- ggplot() + 
  stat_eye(data = longSiteHtW, aes(x = Population, y = value, fill = Population), .width = c(.90, .5), cex = 0.75) +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  theme_classic() +  
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, angle = 78,hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Populations", y = "Population effect", main = "Height") +
    scale_fill_manual(values = c("cyan4","#cc6a70ff", "#f9b641ff", "mediumpurple2")) + 
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 5, label = "a)    Height", cex = 5) 

# LMA
muSiteLMAW <- data.frame(postLMAW$musite)
site <- c("Alex Fraser", "Kamloops","Manning Park", "Smithers")
names(muSiteLMAW) <- site

longSiteLMAW <- melt(muSiteLMAW)
names(longSiteLMAW) <- c("Population", "value")

LMASiteW <- ggplot() + 
  stat_eye(data = longSiteLMAW, aes(x = Population, y = value, fill = Population), .width = c(.90, .5), cex = 0.75) +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  theme_classic() +  
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, angle = 78,hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Populations", y = "Population effect", main = "LMA") +
  scale_fill_manual(values = c("cyan4","#cc6a70ff", "#f9b641ff", "mediumpurple2")) + 
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 0.1, label = "b)    LMA", cex =5) 

# DBH
muSiteDBHW <- data.frame(postDbhW$musite)
site <- c("Alex Fraser", "Kamloops","Manning Park", "Smithers")
names(muSiteDBHW) <- site

longSiteDBHW <- melt(muSiteDBHW)
names(longSiteDBHW) <- c("Population", "value")

DBHSiteW <- ggplot() + 
  stat_eye(data = longSiteDBHW, aes(x = Population, y = value, fill = Population), .width = c(.90, .5), cex = 0.75) +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  theme_classic() +  
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, angle = 78,hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Populations", y = "Population effect", main = "DBH") +
  scale_fill_manual(values = c("cyan4","#cc6a70ff", "#f9b641ff", "mediumpurple2")) + 
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 7, label = "c)    DBH", cex =5) 

# C:N
muSiteCNW <- data.frame(postCNW$musite)
site <- c("Alex Fraser", "Kamloops","Manning Park", "Smithers")
names(muSiteCNW) <- site

longSiteCNW <- melt(muSiteCNW)
names(longSiteCNW) <- c("Population", "value")

CNSiteW <- ggplot() + 
  stat_eye(data = longSiteCNW, aes(x = Population, y = value, fill = Population), .width = c(.90, .5), cex = 0.75) +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  theme_classic() +  
  theme(axis.text.x = element_text(size=12, angle = 78, hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Populations", y = "Population effect", main = "Carbon:Nitrogen") +
  scale_fill_manual(values = c("cyan4","#cc6a70ff", "#f9b641ff", "mediumpurple2")) + 
  theme(legend.title = element_blank()) +  annotate("text", x = 2, y = 10, label = "d)    Carbon:Nitrogen", cex =5) 

pdf("figures/traitMuSiteW.pdf", width = 20, height =5)
plot_grid(htSiteW,LMASiteW,DBHSiteW,CNSiteW, ncol = 4, align = "v")
dev.off()

pdf("figures/sigmaSpSite.pdf", width = 12, height = 4)
par(mfrow = c(1, 5))
hist(postHtW$sigma_site, main = "Height", xlab = "Variation", ylab = "Frequency", col = col1, breaks = 30, ylim = c(0, 4000))
hist(postHtW$sigma_sp, add = T, col = col2)
text(0.2,4000, label = "a)", cex = 2)

hist(postLMAW$sigma_site, main = "Leaf Mass Area", xlab = "Variation", ylab = "Frequency", col = col1, breaks = 50, xlim = c(0,0.2), ylim = c(0, 6000))
hist(postLMAW$sigma_sp, add = T, col = col2)
text(0.002,6000, label = "b)", cex = 2)

hist(postDbhW$sigma_site, main = "Diameter at Breast Height", xlab = "Variation", ylab = "Frequency", col = col1, ylim = c(0, 6000))
hist(postDbhW$sigma_sp, add = T, col = col2)
text(0.2,6000, label = "c)", cex = 2)

hist(postCNW$sigma_site, main = "Carbon:Nitrogen", xlab = "Variation", ylab = "Frequency", col = col1, breaks = 30, ylim = c(0, 4000), xlim = c(0,12))
hist(postCNW$sigma_sp, add = T, col = col2)
text(0.2,4000, label = "d)", cex = 2)

hist(postSsdW$sigma_site, main = "Stem Specific Density", xlab = "Variation", ylab = "Frequency", col = col1, breaks = 30, ylim = c(0, 4000), xlim = c(0,6))
hist(postCNW$sigma_sp, add = T, col = col2, breaks = 50)
text(0.1,4000, label = "e)", cex = 2)

legend("topright",legend = c(expression("Species variation"),
                             expression("Study variation")),
       col = c(col2,col1),
       inset = 0.02, cex = 1, bty = "n", lwd = 3)

dev.off()



load("output/rda/ht_eastern_Feb21.Rda")
Model <- mdl.ht4

load("output/rda/lma_eastern_Feb21.Rda")
Model <- mdl.lma4

load("output/rda/dbh_eastern_Feb21.Rda")
Model <- mdl.dbh4

load("output/rda/cn_eastern_Feb21.Rda")
Model <- mdl.cn4

load("output/rda/ssd_eastern_Feb21.Rda")
Model <- mdl.ssd1