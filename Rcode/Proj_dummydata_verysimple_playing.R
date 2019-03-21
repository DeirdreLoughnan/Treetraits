### Started January 29 2019 ###

## DL writing dummy data for 507 Project ##

#The aim of this code is to build test data for modeling the relationship between phenology and functional traits


rm(list=ls()) 
options(stringsAsFactors = FALSE)

#For this early model, I will have the following variables that need estimates:
#Traits: SLA, ht, m.dbh, wood density, C:N

unique(comb)
setwd("~/Documents/github/Treetraits")
source('Rcode/source/Cleaning_compiling_data.R')


#Starting by looking for high correlations between traits

# Correlation panel
#there are many ways to do this, but this provides a busy but I think useful visual
library("PerformanceAnalytics")
chart.Correlation(comb[,c(20:21, 28:32)], histogram=TRUE, pch=19)


#To start with the simpler model, I am going to develop a linear model with the five traits
#in the source file the mean phenology is taken for each indiviudal tree of for a given treatment factor
int<-11 # days into the experiment of bb 
sigma<-2

nsite = 2 #number of sites
nsp = 28 #number of species
rep = 10 # I think the maximum number is 10, but for many it is closer to 5
ntot<-nsite*nsp*rep #560

#Building the required dataframe for all of the replication and traits
#Assigning a categorical value to site?
site=gl(nsite, rep, length=ntot); site

#randomly generating numbers for each trait based on what I think a small tree would be
site.v=rnorm(ntot,5,2) #Ultimatley I will be pooling across sites, so this might not be necessary
sla.v= rnorm(ntot, 150, 50)
ht.v=rnorm(ntot, 10, 5) #sigma should be really high bc it includes trees and woody shrubs
cn.v=rnorm(ntot, 30, 5)
wood.v=rnorm(ntot, 0.5, 0.05)
stom.v=rnorm(ntot, 400, 50)

#effect sizes
site.diff=2 # I think site 2 will be delayed by 2 days due to the 5 degree diff in lat
sla.diff= -5
ht.diff=3
cn.diff=-5
wood.diff=0.3
stom.diff=25

# #sigmas
# sla.sd=50
# ht.sd=5
# cn.sd=5
# wood.sd=1
# stom.sd=50
# 

#now need to calculate the phenology or y
mm <- model.matrix(~(site+sla.v+ht.v+cn.v+wood.v+stom.v), data.frame(site,sla.v,ht.v,cn.v,wood.v,stom.v))
mm
coeff <- c(1, site.diff, sla.diff, ht.diff, cn.diff, wood.diff, stom.diff)

phen <- rnorm(n = length(site), mean = mm %*% coeff, sd = 1) # This code works but the values are HUGE 
phen

## Here i am trying to do the same calculation a different way, but still getting gigantic values
# phen<-vector(length=560)
 for (i in 1:ntot){
   phen[i]<-int+sla.diff*sla.v[i]+ht.diff*ht.v[i]+cn.diff*cn.v[i]+wood.diff*wood.v[i]+stom.diff*stom.v[i]+rnorm(1,15, sigma)
 }
phen
warnings()
# phen

#Cheating so I can move on:
x<-seq(1:45)
m.day<-15
sigma<-3

phen<-rnorm(ntot, m.day,sigma);phen

fake <- data.frame(phen,site, sla.v,ht.v,cn.v,wood.v,stom.v)
fake

summary(lm(phen ~ (site+site+sla.v+ht.v+cn.v+wood.v+stom.v), data = fake)) # sanity check 

require(rethinking)
simplehist(phen)

###### MAP MODEL ############
an.temp<- map(
  alist(
    phen~dnorm(mu, sigma),
    #mu<-intercept+sla.diff*sla.v+ht.diff*ht.v+cn.diff*cn.v+wood.diff*wood.v+stom.diff*stom.v,
    mu<-intercept+stom.diff*stom.v,
    intercept~dnorm(20, 5),
    #sla.v~dnorm(200, 50),
   # ht.v~dnorm(8, 5), #sigma should be really high bc it includes trees and woody shrubs
    #cn.v~dnorm(30, 5),
    #wood.v~dnorm(1.5, 1),
    stom.v~dnorm(250, 50),
    sigma~dunif(5,1)
  ),
  data=fake)

precis(an.temp)
