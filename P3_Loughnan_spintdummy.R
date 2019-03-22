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

int<-10 # days into the experiment of bb 
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
  phen_stom<- sp_mean[j]+ rnorm(rep*nsite, 0, sigma)+stomdiff*stomc+sladiff*slac+htdiff*htc+cndiff*cnc+wooddiff*woodc
  temp<-data.frame(phen_stom)
  fake_cent<- rbind(fake_cent, temp) # changed so not overwriting fake_stom every time ... 
}

range(fake_cent)

dim(fake_cent)

#Add species - there is probably a better way to do this
sp=gl(nsp, rep*nsite, length= ntot, order=TRUE)

#sp<-(c(1:28, rep*nsite)); sp

#Combine the phenology data with the species identity 
phen_cent<-cbind(fake_cent,sp)
head(phen_cent)
hist(phen_cent$phen_stom)
#Based on R.code 5.40: this seems to work and keep the order of the intercepts consistent
phen_cent$sp_id<-as.integer(phen_cent$sp)

#This line is taking the mean phen_stom for each species
#checkcode <- aggregate(phen_stom["phen_stom"], phen_stom["sp"], FUN=mean) ; checkcode

###############################################
# Now combining the phenology and the trait data
head(phen_cent)
fake_cent<- data.frame(phen_cent,stomc,slac,htc,cnc,woodc)
head(fake_cent)

#fake_stom$sp_id <- coerce_index(fake_stom$sp)
fake_cent<-fake_cent[,c(1,3:8)]
head(fake_cent)

full_m <- map2stan(
  alist(
    phen_cent ~ dnorm(mu, sigma) , # you have two sigmas in your fake data, should have two here, but coded only one (so I changed)
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
  data=fake_cent, iter=4000 , chains=4 
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

#a and sigma_sp tend to be a little off, but not terribly 
precis(full_m, depth=2)

####################################################################
# working with the real data
####################################################################

source('Rcode/source/Cleaning_compiling_data.R') #Here I combined the budburst and trait data, averaging over individual trees for a given treatment


#plot the posterior
par(mfrow=c(1,1))
plot( phen_spint ~ slav , data=fake_spint)
#abline( a=coef(slam)["a"] , b=coef(slam)["bsla"] )

#To get the uncertainty plot many of these lines
post <- extract.samples( slam, n=1000)

# display raw data and sample size
plot( fake_spint$phen_spint , fake_spint$slav ,
      xlim=range(fake_spint$phen_spint) , ylim=range(fake_spint$slav) ,
      col=rangi2 , xlab="SLA" , ylab="Budburst day" )

# plot the lines, with transparency
for ( i in 1:20 )
  abline( a=post$a[i] , b=post$b[i] , col=col.alpha("black",0.3) )

#########################
# Plots
#########################
par(mfrow=c(1,1))

mu_ppc <- link( slam)
mu.mean <- apply( mu_ppc , 2 , mean )
mu.PI <- apply( mu_ppc , 2 , PI )

#simulate data
phen.sim <- sim( slam , n=1000 )
phen.PI <- apply( phen.sim , 2 , PI )

plot( mu.mean ~ fake_spint$phen_spint , col="darkgreen" , ylim=range(mu.PI) ,
      xlab="Observed phen" , ylab="Predicted phen" )
abline( a=0 , b=1 , lty=2, lwd=4)
for ( i in 1:nrow(fake_spint) )
  lines( rep(fake_spint$phen_spint[i],2) , c(mu.PI[1,i],mu.PI[2,i]) ,
         col="darkgreen" )

###############

# compute residuals
phen.resid <- fake_spint$phen_spint - mu.mean
# get ordering by divorce rate
ord <- order(phen.resid)
# make the plot
dotchart( phen.resid[ord] , labels=fake_spint$sp[ord]  , cex=0.6 )
abline( v=0 , col=col.alpha("black",0.2) )
for ( i in 1:nrow(fake_spint) ) {
  j <- ord[i] # which State in order
  lines( fake_spint$phen_spint[j]-c(mu.PI[1,j],mu.PI[2,j]) , rep(i,2) )
  points( fake_spint$phen_spint[j]-c(phen.PI[1,j],phen.PI[2,j]) , rep(i,2),
          pch=3 , cex=0.6 , col="gray" )
}


# Counterfactual plot
sla.avg <- mean( fake_spint$slav)
R.seq <- seq( from=-5 , to=5 , length.out=2800 )
pred.data <- data.frame(
  htv=R.seq,
  wood=R.seq,
  stomv=R.seq,
  cnv=R.seq,
  slav=A.avg
)

mu <- link( slam , data=fakecent )
mu.mean <- apply( mu , 2 , mean ); length(mu.mean)
mu.PI <- apply( mu , 2 , PI )
A.sim <- sim( full_m
              , data=fakecent , n=1e4 )
A.PI <- apply( A.sim , 2 , PI )
plot( phencent ~ slav , data=fakecent , type="n" )

length(R.seq); length(mu.mean)
lines( R.seq , mu.mean )
shade( mu.PI , R.seq )
shade( A.PI , R.seq )

##########################################
#FULL model
##########################################
# Developing better priors

# You know the lowest score is 7 and the max is 20.
# Your model needs parameters: sigma, mean.
names(comb)
range(comb$cn, na.rm=TRUE)
hist(rnorm(1000, 20, 10)) # sla values in data range from 11 to to 54
hist(rnorm(1000, 15, 3)) # hts range from 0.20 to 20.75
hist(rnorm(1000, 400, 100)) # stom_d values range from 30 to 770
hist(rnorm(1000, 2.5, 0.6)) # wood_den range from 0.69 to 4.2
hist(rnorm(1000, 20, 10)) # cn values in data range from 10.59 to 41.47


full<- map2stan(
  alist(
    phen_spint~ dnorm( mu , sigma ) ,
    mu<-a[sp]+bsla*slav+bht*htv+bcn*cnv+bstom*stomv+bwood*woodv, 
    a[sp] ~ dnorm( 0,20) ,
    bsla~dnorm(0, 10),
    bht~dnorm(0, 50),
    bcn~dnorm(0, 15),
    bstom~dnorm(0, 50),
    bwood~dnorm(0, 10),
    sigma ~ dnorm(0,10)),
  data=fake_spint, iter=4000 , chains=4 
)

precis(full, depth=2)

plot(full)

par(mfrow=c(1,1))
plot(precis(full,depth=2))

#<><><><><><><><><><><><><><><><><>
# the estimates are mucg better EXCEPT for wood density and the intercepts are almost all the same! 
infprior<- map2stan(
  alist(
    phen_spint~ dnorm( mu , sigma ) ,
    mu<-a[sp]+bsla*slav+bht*htv+bcn*cnv+bwood*woodv
    +bstom*stomv
    , 
    a[sp] ~ dnorm( 0,20) ,
    bsla~dnorm(20,10),
    bht~dnorm(15,3),
    bcn~dnorm(20,10),
    bwood~dnorm(2.5,0.6),
    bstom~dnorm(400, 100),# the intercepts all become uniform when this variable is added
    sigma ~ dnorm(0,10)),
  data=fake_spint, iter=4000 , chains=4 
)

precis(infprior, depth=2)
plot(precis(infprior, depth=2))



#########################
# Plots: based on code from CH 5 box 5.11, 5.12
#########################

par(mfrow=c(1,1))

mu_ppc <- link( full)
mu.mean <- apply( mu_ppc , 2 , mean )
mu.PI <- apply( mu_ppc , 2 , PI )

#simulate data
phen.sim <- sim( full , n=1000 )
phen.PI <- apply( phen.sim , 2 , PI )

plot( mu.mean ~ fake_spint$phen_spint , col="darkgreen" , ylim=range(mu.PI) ,
      xlab="Observed phen" , ylab="Predicted phen" )
abline( a=0 , b=1 , lty=2, lwd=4)
for ( i in 1:nrow(fake_spint) )
  lines( rep(fake_spint$phen_spint[i],2) , c(mu.PI[1,i],mu.PI[2,i]) ,
         col="darkgreen" )
# This looks pretty good, the data is fairly aligned with the perfect fit line#

##########################

# compute residuals
phen.resid <- fake_spint$phen_spint - mu.mean
# get ordering by divorce rate
ord <- order(phen.resid)
# make the plot
dotchart( phen.resid[ord] , labels=fake_spint$sp[ord]  , cex=0.6 )
abline( v=0 , col=col.alpha("black",0.2) )
for ( i in 1:nrow(fake_spint) ) {
  j <- ord[i] # which State in order
  lines( fake_spint$phen_spint[j]-c(mu.PI[1,j],mu.PI[2,j]) , rep(i,2) )
  points( fake_spint$phen_spint[j]-c(phen.PI[1,j],phen.PI[2,j]) , rep(i,2),
          pch=3 , cex=0.6 , col="gray" )
}

##########################

# Counterfactual plot
sla.avg <- mean( fake_spint$slav)
R.seq <- seq( from=-10 , to=10 , length.out=2800 )
pred.data <- data.frame(
  htv=R.seq,
  wood=R.seq,
  stomv=R.seq,
  cnv=R.seq,
  slav=sla.avg
)

mu <- link(full , data=fake_9spint )
mu.mean <- apply( mu , 2 , mean ); length(mu.mean)
mu.PI <- apply( mu , 2 , PI )
A.sim <- sim( full
              , data=fake_spint , n=1000 )
A.PI <- apply( A.sim , 2 , PI )
plot( phen_spint ~ slav , data=fake_spint , type="n" )

length(R.seq); length(mu.mean)
lines( R.seq , mu.mean )
shade( mu.PI , R.seq )
shade( A.PI , R.seq )

##########################

#PPC based on code from:
# https://cran.r-project.org/web/packages/bayesplot/vignettes/graphical-ppcs.html

library("bayesplot")
library("ggplot2")
library("rstanarm") #rethinking masks "se" and "loo"

#When I attempted to use code from this vignette, several functions did not work for map2stan objects

#1st I attempted to simply run a stan model using the code created from the stancode function, but I was unable to easily do so
stancode(par_pool_r)
head(ws_trim)

wine.model<-stan("wine_stan.stan",data=c("N","N_judge_in","score", "judge_in"),  iter=8000, chains=4)
