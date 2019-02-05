### Started January 29 2019 ###

## DL writing dummy data for 507 Project ##

#The aim of this code is to build test data for modeling the relationship between phenology and functional traits


rm(list=ls()) 
options(stringsAsFactors = FALSE)

#For this early model, I will have the following variables that need estimates:
#Traits: SLA, ht, m.dbh, wood density, C:N

setwd("~/Documents/github/Treetraits")
source('Rcode/source/Cleaning_compiling_data.R')


#Starting by looking for high correlations between traits

# Correlation panel
#there are many ways to do this, but this provides a busy but I think useful visual
library("PerformanceAnalytics")
chart.Correlation(comb[,c(20:21, 28:32)], histogram=TRUE, pch=19)


#To start with the simpler model, I am going to develop a linear model with just the five traits
#in the source file the mean phenology is taken for each indiviudal tree of for a given treatment factor

int<-11 # days into the experiment of bb 
sigma<-2

nsite = 2 #number of sites
nsp = 28 #number of species
rep = 50 # I think the maximum number is 10, but for many it is closer to 5
ntot<-nsite*nsp*rep #2800
ntot
#Building the required dataframe for all of the replication and traits
#Assigning a categorical value to site?
site=gl(nsite, rep, length=ntot)

#randomly generating numbers for each trait based on what I think a small tree would be
#sitev=rnorm(ntot,5,2) #Ultimatley I will be pooling across sites, so this is not necessary at this point in time
slav= rnorm(ntot, 5, 2)
htv=rnorm(ntot, 11, 3) #sigma should be really high (3m) bc it includes trees and woody shrubs
cnv=rnorm(ntot, 10, 2)
woodv=rnorm(ntot, 0.5, 0.05)
stomv=rnorm(ntot, 20, 5)


#effect sizes
site.diff=2 # I think site 2 will be delayed by 2 days due to the 5 degree diff in lat
sladiff= -0.5
htdiff=0.5
cndiff=-0.5
wooddiff=0.3
stomdiff=1

# #sigmas
# sla.sd=50
# ht.sd=5
# cn.sd=5
# wood.sd=1
# stom.sd=50
# 

#now need to calculate the phenology or y
mm <- model.matrix(~(site+slav+htv+cnv+woodv+stomv), data.frame(site,slav,htv,cnv,woodv,stomv))
mm

df<-data.frame(site,slav,htv,cnv,woodv,stomv)

phen<-vector()
# #stomatal density on its own
# for (i in 1:ntot){
#   phen[i]<-int+stomdiff*stomv[i]+rnorm(1,0, sigma)
# }
# phen
# range(phen)
# #
# #cn on its own
# for (i in 1:ntot){
#   phen[i]<-int+wooddiff*woodv[i]+rnorm(1,0, sigma)
# }
# range(phen)

#cn on its own
for (i in 1:ntot){
  phen[i]<-int+cndiff*cnv[i]+rnorm(1,0, sigma)
}
range(phen)

#SLA on its own
for (i in 1:ntot){
  phen[i]<-int+sladiff*slav[i]+rnorm(1,0, sigma)
}
range(phen)
summary(lm(phen~slav))

#ht on its own
# for (i in 1:ntot){
#   phen[i]<-int+htdiff*htv[i]+rnorm(1,0, sigma)
# }
# # phen
# range(phen)
#
phenfull<-vector()
for (i in 1:ntot){
  phenfull[i]<-int+sladiff*slav[i]+htdiff*htv[i]+cndiff*cnv[i]+wooddiff*woodv[i]+stomdiff*stomv[i]+rnorm(1, 0, sigma)
}
range(phenfull)
summary(lm(phen~slav+htv+cnv+woodv+stomv))
#These values look good, they are all positive (unlike using the method below...)

#SLA isn't beign fit well with the data
phennsla<-vector()
for (i in 1:ntot){
  phenfull[i]<-int+htdiff*htv[i]+cndiff*cnv[i]+wooddiff*woodv[i]+stomdiff*stomv[i]+rnorm(1, 0, sigma)
}
range(phenfull)
###############################################
head(mm)
mm.cent<-scale(mm[,3:7]);head(mm.cent)
pairs(mm.cent)
mm.cent<-cbind(mm[,1:2], mm.cent)
head(mm.cent)

# coeff <- c(1, sladiff, htdiff, cndiff, wooddiff, stomdiff)
# 
# #centered
# phencent <- rnorm(n = length(site), mean = mm.cent %*% coeff, sd = 1) # This code works but the values are HUGE
# phencent
# 
# #uncentered
# phen <- rnorm(n = length(site), mean = mm %*% coeff, sd = 1) # This code works but the values are HUGE
# range(phen) #why are some of these values negative?
##############################################

require(rethinking)
simplehist(phenfull)

#Now to do the same phen calculation, but with centered data


slav<-scale(slav)
htv<-scale(htv)
cnv<-scale(cnv)
woodv<-scale(woodv)
stomv<-scale(stomv)

phencent<-vector()
for (i in 1:ntot){
  phencent[i]<-int+sladiff*slav[i]+htdiff*htv[i]+cndiff*cnv[i]+wooddiff*woodv[i]+stomdiff*stomv[i]+rnorm(1,0, sigma)
}
head(phencent)

#hist of sample distribution
simplehist(phencent)


#Now creating the combined datasets of the fake data for use with map
fake<- data.frame(phenfull,mm.cent)
head(fake)

fakecent<- data.frame(phencent,mm.cent)
head(fakecent)

summary(lm(phenfull ~ (slav+htv+cnv+woodv+stomv), data = fake)) # sanity check 

summary(lm(phencent ~ (slav+htv+cnv+woodv+stomv), data = fakecent)) # sanity check 

####### Uncentered data #############
# phen <- rnorm(n = length(site), mean = mm %*% coeff, sd = 1) # This code works but the values are HUGE 
# phen

require(rethinking)
simplehist(phenfull)

simplehist(phencent)



###### MAP MODEL ############
# an.temp<- map(
#   alist(
#     phen~dnorm(mu, sigma),
#     mu<-intercept+sladiff*slav+htdiff*htv+cndiff*cnv+wooddiff*woodv+stomdiff*stomv,
#     intercept~dnorm(10, 1),
#     sladiff~dnorm(20, 1),
#     htdiff~dnorm(8, 5), #sigma should be really high bc it includes trees and woody shrubs
#     cndiff~dnorm(30, 5),
#     wooddiff~dnorm(1.5, 1),
#     stomdiff~dnorm(250, 50),
#     sigma~dunif(5,1)
#   ),
#   data=fake)
# 
# #Now trying a centered model
# an.temp<- map(
#   alist(
#     phen~dnorm(mu, sigma),
#     mu<-intercept+sladiff*slav+htdiff*htv+cndiff*cnv+wooddiff*woodv+stomdiff*stomv,
#     intercept~dnorm(0, 1),
#     sladiff~dnorm(0, 1),
#     htdiff~dnorm(0, 5), #sigma should be really high bc it includes trees and woody shrubs
#     cndiff~dnorm(0, 5),
#     wooddiff~dnorm(0, 1),
#     stomdiff~dnorm(0, 50),
#     sigma~dunif(0,1)
#   ),
#   data=fakecent)

##############################################################
#Full model wont work, what about a simplified one?
#Now trying a centered model
full.m<- map(
  alist(
    phencent~dnorm(mu, sigmahere),
    mu<-a+bsla*slav+bht*htv+bcn*cnv+bstom*stomv+bwood*woodv,
    a~dnorm(0, 10),
    bsla~dnorm(0, 10),
    bht~dnorm(0, 10),
    bcn~dnorm(0, 10),
    bstom~dnorm(0, 10),
    bwood~dnorm(0, 10),
    sigmahere~dunif(0,10)
  ),
  data=fakecent)
precis(full.m)

#Now using map2stan
full.stan<- map2stan(
  alist(
    phencent~dnorm(mu, sigmahere),
    mu<-a+bsla*slav+bht*htv+bcn*cnv+bstom*stomv+bwood*woodv,
    a~dnorm(0, 10),
    bsla~dnorm(0, 10),
    bht~dnorm(0, 10),
    bcn~dnorm(0, 10),
    bstom~dnorm(0, 10),
    bwood~dnorm(0, 10),
    sigmahere~dunif(0,10)
  ),
  data=fakecent)
precis(full.stan)


