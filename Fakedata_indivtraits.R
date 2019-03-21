### Started January 29 2019 ###

## DL writing dummy data for 507 Project ##

#The aim of this code is to build test data for modeling the relationship between phenology and functional traits
#Update: Added partial pooling across species & also changed data to not exclude sites with no phenology data

rm(list=ls()) 
options(stringsAsFactors = FALSE)

#For this early model, I will have the following variables:
#Traits: SLA, height (ht), m.dbh, wood density, C:N, & stomatal density
#Note: After looking for collinearity in the x varibales, I decided to exclude DBH and just use ht, it is biologically more relevant as well

library(rethinking)

setwd("~/Documents/github/Treetraits")
##################################################
# Making test data
###################################################
#Now developing a model with all traits and partial pooling across species

int<-11 # days into the experiment of bb 
sigma<-2
nsite = 2 #number of sites
nsp = 28 #number of species
rep = 50 #
ntot<-nsite*nsp*rep #2800

#Building the required dataframe for all of the replication and traits
#Assigning a variable to site **This will be returned to later when partially pooling across sites
site=gl(nsite, rep, length=ntot)

#randomly generating numbers for each trait based on what I think a small tree would be
#sitev=rnorm(ntot,5,2) #Ultimatley I will be pooling across sites, so this is not necessary at this point in time

# Trait values (trait name + v for value)
# Scaled here so all data is centered
slav=rnorm(ntot, 5, 1)
htv=rnorm(ntot, 11, 3) #sigma should be really high (3m) bc it includes trees and woody shrubs
cnv=rnorm(ntot, 10, 2)
woodv=rnorm(ntot, 0.5, 0.05)
stomv=rnorm(ntot, 200, 50)


#effect sizes
site.diff=2 # I think site 2 will be delayed by 2 days due to the 5 degree diff in lat
sladiff= -0.5
htdiff=0.5
cndiff=-0.5
wooddiff=0.3
stomdiff=1

#####################################################################
# Stomatal density alone
#####################################################################

#NOW working with CENTERED DATA
#Centering stomatal density
stomc=scale(rnorm(ntot, 200, 50))
range(stomc)

fake_cent<-vector()
sp_mean <- rnorm(nsp, mean=int, sd=5)
for (j in c(1:length(sp_mean))){
  phen_stom<- sp_mean[j]+stomdiff*stomc+rnorm(rep*nsite, 0, sigma)
  temp<-data.frame(phen_stom)
  fake_stom<- rbind(fake_cent, temp)
}

head(fake_stom)

#Add species - there is probably a better way to do this
sp=gl(nsp, rep*nsite, length= ntot)

#Combine the phenology data with the species identity 
phen_stom<-cbind(fake_stom,sp)
head(fake_stom)
###############################################
# Now combining the phenology and the trait data

fake_stom<- data.frame(phen_stom,stomc)
head(fake_stom)

stom_m <- map2stan(
  alist(
    phen_stom~ dnorm( mu , sigma_sp ) ,
    mu<-a_sp[sp]+bstom*stomc, 
    a_sp[sp] ~ dnorm( a,sigma_sp) ,
    a~dnorm(0,20),
    bstom~dnorm(0, 10),
    sigma_sp ~ dnorm(0,10)),
  data=fake_stom, iter=4000 , chains=4 
)

sp_mean
stomdiff

#Still not getting hte same output!
precis(stom_m, depth=2)

plot(precis(stom_m, depth=2))

#####################################################################
# SLA alone
#####################################################################

#NOW working with scaled values of sla
slac=scale(rnorm(ntot, 5, 1))

fake_cent<-vector()
sp_mean <- rnorm(nsp, mean=int, sd=5)
for (j in c(1:length(sp_mean))){
  #phen_spint<- int+sladiff*slav+rnorm(rep*nsite, sp_mean[j], sigma)
  #phen_spint<- sladiff*slav+htdiff*htv+cndiff*cnv+wooddiff*woodv+sladiff*slav+rnorm(rep*nsite, sp_mean[j], sigma)
  phen_sla<- sp_mean[j]+sladiff*slac+rnorm(rep*nsite, 0, sigma)
  temp<-data.frame(phen_sla)
  fake_sla<- rbind(fake_cent, temp)
}

head(fake_sla)

#Add species - there is probably a better way to do this
sp=gl(nsp, rep*nsite, length= ntot)

#Combine the phenology data with the species identity 
fake_sla<-cbind(fake_sla,sp)
head(fake_sla)

###############################################
# Now combining the phenology and the trait data

fake_sla<- data.frame(fake_sla,slac)
head(fake_sla)

sla_m <- map2stan(
  alist(
    phen_sla~ dnorm( mu , sigma_sp) ,
    mu<-a_sp[sp]+bsla*slac, 
    a_sp[sp] ~ dnorm( a,sigma_sp) ,
    a~dnorm(0,1),
    bsla~dnorm(0, 10),
    sigma_sp ~ dnorm(0,10)),
  data=fake_sla, iter=4000 , chains=4 
)

sp_mean
sladiff

#Still not getting hte same output!
precis(sla_m, depth=2)

plot(precis(sla_m, depth=2))

#####################################################################
# Height alone
#####################################################################

#NOW working with scaled values of ht
htc=scale(rnorm(ntot, 5, 1))

fake_cent<-vector()
sp_mean <- rnorm(nsp, mean=int, sd=5)
for (j in c(1:length(sp_mean))){
  #phen_spint<- int+htdiff*htv+rnorm(rep*nsite, sp_mean[j], sigma)
  #phen_spint<- htdiff*htv+htdiff*htv+cndiff*cnv+wooddiff*woodv+htdiff*htv+rnorm(rep*nsite, sp_mean[j], sigma)
  phen_ht<- sp_mean[j]+htdiff*htc+rnorm(rep*nsite, 0, sigma)
  temp<-data.frame(phen_ht)
  fake_ht<- rbind(fake_cent, temp)
}

head(fake_ht)

#Add species - there is probably a better way to do this
sp=gl(nsp, rep*nsite, length= ntot)

#Combine the phenology data with the species identity 
fake_ht<-cbind(fake_ht,sp)
head(fake_ht)

###############################################
# Now combining the phenology and the trait data

fake_ht<- data.frame(fake_ht,htc)
head(fake_ht)

ht_m <- map2stan(
  alist(
    phen_ht~ dnorm( mu , sigma_sp) ,
    mu<-a_sp[sp]+bht*htc, 
    a_sp[sp] ~ dnorm( a,sigma_sp) ,
    a~dnorm(0,1),
    bht~dnorm(0, 10),
    sigma_sp ~ dnorm(0,10)),
  data=fake_ht, iter=4000 , chains=4 
)

sp_mean
htdiff

#Still not getting hte same output!
precis(ht_m, depth=2)

plot(precis(ht_m, depth=2))