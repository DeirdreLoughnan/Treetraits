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
#source('Rcode/source/Cleaning_compiling_data.R') #Here I combined the budburst and trait data, averaging over individual trees for a given treatment
##################################################
# Making test data
###################################################
#Now developing a model with all traits and partial pooling across species

int<-11 # days into the experiment of bb 
sigma<-2

nsite = 2 #number of sites
nsp = 28 #number of species
rep = 50 # I think the maximum number is 10, but to increase the size of the dataset I have made it 50
ntot<-nsite*nsp*rep #2800

#Building the required dataframe for all of the replication and traits
#Assigning a variable to site **This will be returned to later when partially pooling across sites
site=gl(nsite, rep, length=ntot)

#randomly generating numbers for each trait based on what I think a small tree would be
#sitev=rnorm(ntot,5,2) #Ultimatley I will be pooling across sites, so this is not necessary at this point in time

# Trait values (trait name + v for value)
slav= rnorm(ntot, 5, 1); 
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

#making a matrix of the x dummy variables
 mm <- model.matrix(~(site+slav+htv+cnv+woodv+stomv), data.frame(site,slav,htv,cnv,woodv,stomv))
 mm

#making a simple data frame of the x dummy variables
df<-data.frame(site,slav,htv,cnv,woodv,stomv)
range(df$stomv)
#####################################################################
# Now adding separate intercepts for each species ** Starting very simple with just SLA

fake_spint<-vector()
sp_mean <- rnorm(nsp, mean=int, sd=5)
for (j in c(1:length(sp_mean))){
  #phen_spint<- int+sladiff*slav+rnorm(rep*nsite, sp_mean[j], sigma)
  phen_spint<- sladiff*slav+htdiff*htv+cndiff*cnv+wooddiff*woodv+rnorm(rep*nsite, sp_mean[j], sigma)
  #phen_spint<- sladiff*slav+htdiff*htv+cndiff*cnv+wooddiff*woodv+stomdiff*stomv+rnorm(rep*nsite, sp_mean[j], sigma)
  #phen_spint<- sp_mean[j]+sladiff*slav+htdiff*htv+cndiff*cnv+wooddiff*woodv+stomdiff*stomv+rnorm(rep*nsite, 0, sigma)
  
  temp<-data.frame(phen_spint)
  phen_spint_cent<- rbind(fake_spint, temp)
}

#Centered phenology values

## Stomatal density excluded

fake_spint<-vector()
sp_mean <- rnorm(nsp, mean=int, sd=sigma)
for (j in c(1:length(sp_mean))){
  phen_spint<- sladiff*slav+htdiff*htv+cndiff*cnv+wooddiff*woodv+rnorm(rep*nsite, sp_mean[j], sigma)
  #phen_spint<- int+sladiff*slav+htdiff*htv+cndiff*cnv+wooddiff*woodv+stomdiff*stomv+rnorm(rep*nsite, sp_mean[j], sigma)
  temp<-data.frame(phen_spint)
  phen_spint_cent<- rbind(fake_spint, temp)
}

head(phen_spint)
range(phen_spint)

#Add species - there is probably a better way to do this
sp=gl(nsp, rep*nsite, length= ntot)


#Combine the phenology data with the species identity 
phen_spint<-cbind(phen_spint,sp)
head(phen_spint)

phen_spint_cent<-cbind(phen_spint_cent,sp)
###############################################
# Now combining the phenology and the trait data

fake_spint<- data.frame(phen_spint,mm)
#Note Stan doesn't like dots so changing X.Intercept name
colnames(fake_spint)[colnames(fake_spint)=="X.Intercept."] <- "Intercept"

mm_cent<-scale(mm[,3:7]);head(mm.cent)
fake_spint_cent<-data.frame(phen_spint_cent,mm_cent)
head(fake_spint_cent)


##########################################
#Model with all traits except stomatal den
##########################################
# Developing better priors

names(comb)
range(comb$cn, na.rm=TRUE)
hist(rnorm(1000, 20, 10)) # sla values in data range from 11 to to 54
hist(rnorm(1000, 15, 3)) # hts range from 0.20 to 20.75
hist(rnorm(1000, 400, 100)) # stom_d values range from 30 to 770
hist(rnorm(1000, 2.5, 0.6)) # wood_den range from 0.69 to 4.2
hist(rnorm(1000, 20, 10)) # cn values in data range from 10.59 to 41.47

#trim stomatal den
fake_spint_cent<-fake_spint_cent[,1:6]

full<- map2stan(
  alist(
    phen_spint~ dnorm( mu , sigma ) ,
    mu<-a[sp]+bsla*slav+bht*htv+bcn*cnv
    #+bstom*stomv
    +bwood*woodv, 
    a[sp] ~ dnorm( 0,20) ,
    bsla~dnorm(0, 10),
    bht~dnorm(0, 50),
    bcn~dnorm(0, 15),
    #bstom~dnorm(0, 50),
    bwood~dnorm(0, 10),
    sigma ~ dnorm(0,10)),
  data=fake_spint_cent, iter=4000 , chains=4 
)

sp_mean
sigma
site.diff=2 
sladiff= -0.5
htdiff=0.5
cndiff=-0.5
wooddiff=0.3
stomdiff=1

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
