#Started March 2020

#aim: create preliminary figures for use in my committee proposal


rm(list=ls()) 
options(stringsAsFactors = FALSE)

require(ggplot2)
setwd("~/Documents/github/Treetraits")

df<-read.csv("data/eastern/raw.indiv.values.wtrait.csv")

require(tidyr)
names(df)

traitsummary <-
  ddply(df, c("sp"), summarise,
        meanlma = mean(lma, na.rm=TRUE),
        sdlma = sd(lma, na.rm=TRUE),
        selma = sd(lma, na.rm=TRUE)/sqrt(length(lma)),
        meanht = mean(ht, na.rm=TRUE),
        sdht = sd(ht, na.rm=TRUE),
        seht = sd(ht, na.rm=TRUE)/sqrt(length(ht)),
        meancn = mean(cn, na.rm=TRUE),
        sdcn = sd(cn, na.rm=TRUE),
        secn = sd(cn, na.rm=TRUE)/sqrt(length(cn)))

traitsummary

shrub<-traitsummary[c(21,22,38,42,44),]; shrub$group<-"shrub"
tree<-traitsummary[c(2,11,18,27,30),]; tree$group<-"tree"

sub<-rbind(shrub, tree)

require(ggplot2)
order<-c("KALANG","LONCAN","SPIALB","VACMYR","VIBLAN","ACERUB", "BETPAP","FRANIG","POPGRA","QUERUB")

library(wesanderson)
### Plot for lma ############################################################
lma<- ggplot(sub, aes(x=sp, y=meanlma, col=group)) +
  geom_errorbar(aes(ymin = meanlma-selma, ymax= meanlma+selma))+
  geom_point()+ 
  labs(title="A)", y=bquote('Leaf Mass Area'~(g/cm^2)), x="Species")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  scale_color_brewer(palette="Dark2")


### Plot for ht ############################################################
ht<- ggplot(sub, aes(x=sp, y=meanht, col=group)) +
  geom_errorbar(aes(ymin = meanht-seht, ymax= meanht+seht))+
  geom_point()+
  labs(title="B)",y=bquote('Height'~(m)), x="Species")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  scale_color_brewer(palette="Dark2")


### Plot for cn ############################################################
cn<- ggplot(sub, aes(x=sp, y=meancn, col=group)) +
  geom_errorbar(aes(ymin = meancn-secn, ymax= meancn+secn))+
  geom_point()+
  labs(title="C)",y="C:N", x="Species")+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  scale_color_brewer(palette="Dark2")



## Want all 3 figures to be in one long panel
#install.packages("gridExtra")
require(gridExtra)
grid.arrange(lma,ht, cn, nrow=1, widths= c(2,2,2.4), heights = 0.5)
