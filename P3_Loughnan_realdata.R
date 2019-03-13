### Started March 11 2019 ###

## DL running models using DF trait+phenology data ##

#The aim of this code is use the model developed using test data to model the relationship between phenology and functional traits

rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(rethinking)

setwd("~/Documents/github/Treetraits")
source('Rcode/source/Cleaning_compiling_data.R') #Here I combined the budburst and trait data, averaging over individual trees for a given treatment. 

#Source file creates the datafile in which the phenology and trait data is combined
head(comb)

#Need to trim unneeded columns:
d<-comb[,c(2,12:13, 15:18)]
head(d)


(d$stom_d)
d <- d[ complete.cases(d), ] #there are a lot of missing cases for height

# Again starting with the simpler model:
slam<- map2stan(
  alist(
    bday~ dnorm( mu , sigma ) ,
    mu<-a[sp]+bsla*sla, 
    a[sp] ~ dnorm( 0,20) ,
    bsla~dnorm(0, 10),
    sigma ~ dnorm(0,10)),
  data=d, iter=4000 , chains=4 
)

precis(slam, depth=2)

plot(slam)

par(mfrow=c(1,1))

#plot the posterior distribution estimates
plot(precis(slam,depth=2))

#the etimate of sla is essentially zero 
####################################################################
#Now building, adding one factor at a time until full model:

full<- map2stan(
  alist(
    bday~ dnorm( mu , sigma ) ,
    mu<-a[sp]+bsla*sla+bht*Height+bcn*cn+bstom*stom_d+bwood*wood_den, 
    a[sp] ~ dnorm( 0,20) ,
    bsla~dnorm(0, 10),
    bht~dnorm(0, 10),
    bcn~dnorm(0, 10),
    bstom~dnorm(0, 10),
    bwood~dnorm(0,10),
    sigma ~ dnorm(0,10)),
  data=d, iter=4000 , chains=4 
)

plot(full)
#all chains appear to converge, although the caterpillars for the intercepts are fuzzy, while it is too hard to see the ones for the actual traits
precis(full, depth=2)

#plot the posterior distribution estimates
par(mfrow=c(1,1))
plot(precis(full,depth=2))

#Interestingly we see negative intercepts for several species, suggesting there are issues with the model. 

#########################
# Plots
#Adapting code taken from box 5.11 & 5.12 Rethinking text
#########################
par(mfrow=c(1,1))
mean
mu_ppc <- link( full)
mu.mean <- apply( mu_ppc , 2 , mean )
mu.PI <- apply( mu_ppc , 2 , PI )

#simulate data
phen.sim <- sim( full , n=1000 )
phen.PI <- apply( phen.sim , 2 , PI )

plot( mu.mean ~ d$bday , col="darkgreen" , ylim=range(mu.PI) ,
      xlab="Observed phen" , ylab="Predicted phen" )
abline( a=0 , b=1 , lty=2, lwd=4)
for ( i in 1:nrow(d) )
  lines( rep(d$bday[i],2) , c(mu.PI[1,i],mu.PI[2,i]) ,
         col="darkgreen" )
#It appears there are some points that are really poorly predicted by the model
###############

# compute residuals
phen.resid <- d$bday- mu.mean
# get ordering by divorce rate
ord <- order(phen.resid); ord
# make the plot
dotchart( phen.resid[ord] , labels=d$sp[ord]  , cex=0.6 )
abline( v=0 , col=col.alpha("black",0.2) )
for ( i in 1:nrow(d) ) {
  j <- ord[i] # which State in order
  lines( d$bday[j]-c(mu.PI[1,j],mu.PI[2,j]) , rep(i,2) )
  points( d$bday[j]-c(phen.PI[1,j],phen.PI[2,j]) , rep(i,2),
          pch=3 , cex=0.6 , col="gray" )
}
#There is something critically wrong with this figure,
############################################################
#This code does not work, there is an issue with calculating mu that I still need to figure out
# Counterfactual plot
sla.avg <- mean( d$sla)
R.seq <- seq( from=-100 , to=100 , length.out=2800 )
pred.data <- data.frame(
  Height=R.seq,
  wood_den=R.seq,
  stom_d=R.seq,
  cn=R.seq,
  sla=sla.avg
)
#Mu is just NA's, unclear why
mu <- link(full , data=d ); mu
mu.mean <- apply( mu , 2 , mean ); length(mu.mean)
mu.PI <- apply( mu , 2 , PI)
A.sim <- sim( full
              , data=d , n=1000 )
A.PI <- apply( A.sim , 2 , PI )
plot( bday ~ sla , data=d , type="n" )

length(R.seq); length(mu.mean)
lines( R.seq , mu.mean )
shade( mu.PI , R.seq )
shade( A.PI , R.seq )

############################################################
