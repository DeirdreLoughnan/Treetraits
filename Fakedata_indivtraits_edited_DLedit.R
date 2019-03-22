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
sigma<- 0.1
nsite = 2 #number of sites
nsp = 28 #number of species
rep = 50 #
ntot<-nsite*nsp*rep #2800

#Building the required dataframe for all of the replication and traits
#Assigning a variable to site **This will be returned to later when partially pooling across sites
#site=gl(nsite, rep, length=ntot)

#randomly generating numbers for each trait based on what I think a small tree would be
#sitev=rnorm(ntot,5,2) #Ultimatley I will be pooling across sites, so this is not necessary at this point in time

# Trait values (trait name + v for value)

slav=rnorm(ntot, 5, 1)
htv=rnorm(ntot, 11, 3) 
cnv=rnorm(ntot, 10, 2)
woodv=rnorm(ntot, 0.7, 0.05)
stomv=rnorm(ntot, 200, 50)


#effect sizes
site.diff=2 
sladiff= -0.5
htdiff=0.5
cndiff=-0.5
wooddiff=1
stomdiff=1

fake_cent<-vector()
sp_mean <- rnorm(nsp, mean=int, sd=5);
sp_mean<-sp_mean[order=T]
for (j in c(1:length(sp_mean))){
  phen_stom<- sp_mean[j]+ rnorm(rep*nsite, 0, sigma)+stomdiff*stomv+sladiff*slav+htdiff*htv+cndiff*cnv+wooddiff*woodv
  temp<-data.frame(phen_stom)
  fake_cent<- rbind(fake_cent, temp) # changed so not overwriting fake_stom every time ... 
}

range(fake_cent)

#####################################################################
# Building up to the full model
#####################################################################

#NOW working with CENTERED DATA
#Centering stomatal density
stomc=scale(rnorm(rep*nsite, 200, 50))
slac= scale(rnorm(rep*nsite, 5, 1)) 
htc=scale(rnorm(rep*nsite, 11, 3)) 
cnc=scale(rnorm(rep*nsite, 10, 2))
woodc=scale(rnorm(rep*nsite, 0.7, 0.05))

# Here's what I think the loop should be doing.... 
fake_cent<-vector()
sp_mean <- rnorm(nsp, mean=int, sd=5);
sp_mean<-sp_mean[order=T]
for (j in c(1:length(sp_mean))){
  phen_c<- sp_mean[j]+ rnorm(rep*nsite, 0, sigma)+stomdiff*stomc+sladiff*slac+htdiff*htc+cndiff*cnc+wooddiff*woodc
  temp<-data.frame(phen_c)
  fake_cent<- rbind(fake_cent, temp) # changed so not overwriting fake_stom every time ... 
}

range(fake_cent)

dim(fake_cent)


#Add species - there is probably a better way to do this
sp=gl(nsp, rep*nsite, length= ntot, order=TRUE)

#sp<-(c(1:28, rep*nsite)); sp

#Combine the phenology data with the species identity 
phen<-cbind(fake_cent,sp)

hist(phen$phen_c)
#Based on R.code 5.40: this seems to work and keep the order of the intercepts consistent
phen$sp_id<-as.integer(phen$sp)

#This line is taking the mean phen_stom for each species
#checkcode <- aggregate(phen_stom["phen_stom"], phen_stom["sp"], FUN=mean) ; checkcode

###############################################
# Now combining the phenology and the trait data
head(phen)
fake_data<- data.frame(phen,stomc,slac,htc,cnc,woodc)
head(fake_data)

#fake_stom$sp_id <- coerce_index(fake_stom$sp)
fake<-fake_data[,c(1,3:8)]
head(fake)

full_m <- map2stan(
  alist(
    phen_c ~ dnorm(mu, sigma) , # you have two sigmas in your fake data, should have two here, but coded only one (so I changed)
    mu <- a_sp[sp_id]+bstom*stomc+bsla*slac+bht*htc+bcn*cnc+bwood*woodc, 
    a_sp[sp_id] ~ dnorm(a, sigma_sp) , # line 65 gives sigma for sp as 5 
    a~dnorm(0, 10),
    bstom~dnorm(0, 10),
    bsla~dnorm(0, 10),
    bht~dnorm(0,10),
    bcn~dnorm(0,10),
    bwood~dnorm(0,10),
    sigma ~ dnorm(0,1),
    sigma_sp ~ dnorm(0,5)),
  data=fake, iter=4000 , chains=4 
)

plot(full_m)

sp_mean
stomdiff
sladiff
htdiff
cndiff
wooddiff
sigma
int

#Still not getting hte same output!
precis(full_m, depth=2)

par(mfrow=c(1,1))
plot(precis(stom_m, depth=2))

#####################################################################
#####################################################################
# Stomatal density alone
#####################################################################

#NOW working with CENTERED DATA
#Centering stomatal density
stomc=scale(rnorm(rep*nsite, 200, 50))
range(stomc)

if(FALSE){ # Here's a check I did to look at what the loop was doing ... 
fake_cent<-vector()
sp_mean <- rnorm(nsp, mean=int, sd=5);sp_mean
for (j in c(1:length(sp_mean))){
  phen_stom<- sp_mean[j]+stomdiff*stomc+rnorm(rep*nsite, 0, sigma)
  temp<-data.frame(phen_stom)
  fake_stom<- rbind(fake_cent, temp) 
}
fake_stom
    }


# Here's what I think the loop should be doing.... 
fake_cent<-vector()
sp_mean <- rnorm(nsp, mean=int, sd=5);
sp_mean<-sp_mean[order=T]
for (j in c(1:length(sp_mean))){
  phen_stom<- sp_mean[j]+ rnorm(rep*nsite, 0, sigma)+stomdiff*stomc # stomc is 2800 rows already so it's too big, should be rep*nsite long ....  
  temp<-data.frame(phen_stom)
  fake_cent<- rbind(fake_cent, temp) # changed so not overwriting fake_stom every time ... 
}

range(fake_cent$phen_stom)

fake_cent
dim(fake_cent)
fake_stom <- fake_cent

#Add species - there is probably a better way to do this
sp=gl(nsp, rep*nsite, length= ntot, order=TRUE)

#sp<-(c(1:28, rep*nsite)); sp

#Combine the phenology data with the species identity 
phen_stom<-cbind(fake_stom,sp)
head(phen_stom)

#Based on R.code 5.40: this seems to work and keep the order of the intercepts consistent
phen_stom$sp_id<-as.integer(phen_stom$sp)

#This line is taking the mean phen_stom for each species
#checkcode <- aggregate(phen_stom["phen_stom"], phen_stom["sp"], FUN=mean) ; checkcode

###############################################
# Now combining the phenology and the trait data
head(phen_stom)
fake_stom<- data.frame(phen_stom,stomc)
head(fake_stom)

#fake_stom$sp_id <- coerce_index(fake_stom$sp)
fake_stom<-fake_stom[,c(1,3:4)]

stom_m <- map2stan(
  alist(
    phen_stom ~ dnorm(mu, sigma) , # you have two sigmas in your fake data, should have two here, but coded only one (so I changed)
    mu <- a_sp[sp_id]+bstom*stomc, # add back in when fixed fake data and model for intercept
    a_sp[sp_id] ~ dnorm(a, sigma_sp) , # line 65 gives sigma for sp as 5 
    a~dnorm(0, 10),
    bstom~dnorm(0, 10),
    sigma ~ dnorm(0,1),
    sigma_sp ~ dnorm(0,5)),
  data=fake_stom, iter=4000 , chains=2 
)

# The above code looks right to me (can compare with winepooling_p2.R test data and code) but the a_sp values returned do not match up ... .except I think the values appear to match, just not in the right order. I turned your sigma down (you can turn it back up) and see that map2stan is ordering them this way I believe: c(1, 12, 22, 23, 24, 25, 26, 27, 28 ...) .... would be good to figure out how to get the right order (read up on map2stan I suggest, and re-read ordering in index ... ordering in R is painful), but I could not. I would work on that alongside trying to get your slopes added back in. 

order(as.factor(unique(fake_stom$sp)))
order(as.character(unique(fake_stom$sp)))

sp_mean
stomdiff
sigma
int


#Still not getting hte same output!
precis(stom_m, depth=2)

plot(precis(stom_m, depth=2))

#####################################################################
#####################################################################
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

check <- map2stan(
  alist(
    phen_sla~ dnorm( mu , sigma) , # you have two sigmas in your fake data, should have two here, but coded only one (so I changed)
    mu<-a_sp[sp]+bsla*slac, 
    a_sp[sp] ~ dnorm(a, sigma_sp) , # line 65 gives sigma for sp as 5 
    a~dnorm(0,1),
    bsla~dnorm(0, 10),
    sigma ~ dnorm(0,10), 
    sigma_sp ~ dnorm(0,10)),
  data=fake_stom, iter=4000 , chains=4 
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
