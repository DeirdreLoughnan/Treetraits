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
#sigma<-0.3
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
sp_mean <- rnorm(nsp, mean=int, sd=1);
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
    a~dnorm(0, 50),
    bstom~dnorm(0, 50),
    bsla~dnorm(0,50 ),
    bht~dnorm(0,50),
    bcn~dnorm(0,50),
    bwood~dnorm(0,50),
    sigma ~ dnorm(0,0.1),
    sigma_sp ~ dnorm(0,1)),
  data=fake, iter=4000 , chains=4 
)

plot(full_m)

sp_mean
stomdiff # 1
sladiff #-0.5
htdiff #0.5
cndiff #-0.5
wooddiff #1
sigma #0.1
int #11

precis(full_m, depth=2)
# Mean StdDev lower 0.89 upper 0.89 n_eff Rhat
# a_sp[1]  11.68   0.01      11.67      11.70 10337    1
# a_sp[2]  10.86   0.01      10.84      10.87 11137    1
# a_sp[3]  10.43   0.01      10.41      10.44 10392    1
# a_sp[4]  12.94   0.01      12.92      12.95 10431    1
# a_sp[5]   9.75   0.01       9.74       9.77 12011    1
# a_sp[6]  11.75   0.01      11.73      11.76 11609    1
# a_sp[7]  11.00   0.01      10.98      11.01  9752    1
# a_sp[8]   9.85   0.01       9.83       9.86 10920    1
# a_sp[9]  11.59   0.01      11.57      11.60 10543    1
# a_sp[10] 11.78   0.01      11.76      11.79 10780    1
# a_sp[11]  7.98   0.01       7.96       7.99  9775    1
# a_sp[12] 11.41   0.01      11.39      11.43 10579    1
# a_sp[13] 10.94   0.01      10.92      10.96 11169    1
# a_sp[14] 13.14   0.01      13.13      13.16 11915    1
# a_sp[15] 10.51   0.01      10.49      10.52 10754    1
# a_sp[16]  9.96   0.01       9.95       9.98 11341    1
# a_sp[17] 12.96   0.01      12.94      12.97 11817    1
# a_sp[18]  9.07   0.01       9.06       9.09  9646    1
# a_sp[19] 11.30   0.01      11.29      11.32 10947    1
# a_sp[20] 11.07   0.01      11.06      11.09 11036    1
# a_sp[21] 11.88   0.01      11.87      11.90 10761    1
# a_sp[22] 11.41   0.01      11.40      11.43 10040    1
# a_sp[23]  8.90   0.01       8.88       8.91 11691    1
# a_sp[24] 11.08   0.01      11.06      11.09  9999    1
# a_sp[25]  9.95   0.01       9.93       9.97  9308    1
# a_sp[26] 10.28   0.01      10.27      10.30  9664    1
# a_sp[27] 10.18   0.01      10.16      10.19 11303    1
# a_sp[28] 10.66   0.01      10.65      10.68 10463    1
# a        10.86   0.23      10.51      11.25  9907    1
# bstom     1.00   0.00       1.00       1.00 16984    1
# bsla     -0.50   0.00      -0.50      -0.50 22187    1
# bht       0.50   0.00       0.50       0.50 18583    1
# bcn      -0.50   0.00      -0.51      -0.50 16383    1
# bwood     1.00   0.00       1.00       1.01 19504    1
# sigma     0.10   0.00       0.10       0.10 10677    1
# sigma_sp  1.23   0.17       0.96       1.47  8170    1

par(mfrow=c(1,1))
plot(precis(full_m, depth=2))

#<><><><><><><><><><><><><><><><><><><><>
# Posterior Predictive Checks:

pppost <- extract.samples(full_m)
totaln <- nrow(fake)

# Make one new dataset from model parameters and means from posterior
newsp <- rnorm(nsp, mean(pppost$a[1:1000]), mean(pppost$sigma_sp[1:1000]))
newspdat <- rnorm(totaln, rep(newsp, each=nsite), mean(pppost$sigma[1:1000]))

par(mfrow=c(1,2))
hist(fake$phen_c, main="real data")
hist(newspdat, main="PPC")

# Now do it 10 times (still using means from posterior)
newspmat <- matrix(nrow=10, ncol=totaln)
for (i in 1:10){
  newsp <- rnorm(nsp, mean(pppost$a[1:1000]), mean(pppost$sigma_sp[1:1000]))
  newspdat <- rnorm(totaln, rep(newsp, each=nsite), mean(pppost$sigma[1:1000]))
  newspmat[i,] <- newspdat
}


str(pppost)
par(mfrow=c(2,6))
hist(fake$phen_c, main="real data")
for (i in 1:10){
  hist(newspmat[i,], main="PPC")
}

# Similarly ...
library(bayesplot)
color_scheme_set("brightblue")
ppc_dens_overlay(fake$phen_c, newspmat[1:5, ])
newspmat
dim(newspmat) #10, 2800
# Visualize the partial pooling ...
# See https://www.tjmahr.com/plotting-partial-pooling-in-mixed-effects-models/
##
coerce_index(unique(fake$sp)); head(fake)
d_sort <- fake[order(fake$sp_id),]

post <- extract.samples(full_m) # grab the posterior
str(post)
unique(d_sort$sp_id)
nsp=28
J <- nsp 

com_pool_mod <- lm(phen_c ~ 1, d_sort)
species <- unique(d_sort$sp_id)
no_pool <- data.frame(sp_id=unique(d_sort$sp_id),
                      intercept=rep(NA, J), 
                      slope=rep(NA, J))
head(d_sort)
with_pool <- data.frame(sp_id=unique(d_sort$sp_id),
                        intercept=rep(NA, J), 
                        slope=rep(NA, J))
df.gravity <- no_pool[1:2,2:3]
with_pool$model <- "partial pooling"
no_pool$model <- "no pooling"

data<-subset(d_sort, sp_id==unique(d_sort$sp_id)[1])


#data=subset(vino.formodel.alt.sort, judge==unique(vino.formodel.alt.sort$judge)[j]))
for (j in 1:J){
  modhere <- lm(phen_c~1,data=subset(d_sort, sp_id==unique(d_sort$sp_id)[j]))
  no_pool$intercept[j] <- coef(modhere)[1]
  no_pool$slope[j] <- j # cheat for no-slope model
}

for (j in 1:J){
  with_pool$intercept[j] <- mean(post$a_sp[1:1000, j])
  with_pool$slope[j] <- j # cheat for no-slope model
}

# head(no_pool)
# head(with_pool)

df.pulled <- rbind(no_pool, with_pool)
str(post)
df.gravity$model <- NA
df.gravity$intercept[1] <- coef(com_pool_mod)["(Intercept)"]
df.gravity$slope[1] <- 4.5 # cheat for no-slope model: middle of plot
df.gravity$model[1] <- "complete pooling"
df.gravity$intercept[2] <- mean(post$a[1:1000])
df.gravity$slope[2] <- 4.5 # cheat for no-slope model: middle of plot
df.gravity$model[2] <- "partial pooling (mu)"

library(ggplot2)
#pdf("modelscompare_pp.pdf", width = 9, height = 6)
ggplot(df.pulled) + 
  aes(x = intercept, y = slope, color = model) + 
  geom_point(size = 2) + 
  geom_point(data = df.gravity, size = 5) + 
  # Draw an arrow connecting the observations between models
  geom_path(aes(group = as.character(sp_id), color = NULL), 
            arrow = arrow(length = unit(.02, "npc"))) + 
  # Use ggrepel to jitter the labels away from the points
  ggrepel::geom_text_repel(
    aes(label = sp_id, color = NULL), 
    data = no_pool, size=2) + 
  theme(legend.position = "bottom") + 
  ggtitle("Pooling of regression parameters") + 
  xlab("Intercept estimate") + 
  ylab("") + 
  scale_color_brewer(palette = "Dark2") 
#dev.off()

##################
#Attempting to better predict the shape of the observed data

# full_in <- map2stan(
#   alist(
#     phen_c ~ dnorm(mu, sigma) , # you have two sigmas in your fake data, should have two here, but coded only one (so I changed)
#     mu <- a_sp[sp_id]+bstom*stomc+bsla*slac+bht*htc+bcn*cnc+bwood*woodc, 
#     a_sp[sp_id] ~ dnorm(mu_a, sigma_sp) , # line 65 gives sigma for sp as 5 
#     mu_a~dnorm(0, 2),
#     bsla~dnorm(0,0.5),
#         bht~dnorm(0,0.5),
#         bcn~dnorm(0,0.5),
#         bwood~dnorm(0,0.5),
#         bstom~dnorm(0,0.5),
#     sigma ~ dnorm(0,0.5),
#     sigma_sp ~ dnorm(0,5)),
#   data=fake, iter=4000 , chains=4 
# )
# 
# # full_in <- map2stan(
# #   alist(
# #     phen_c ~ dnorm(mu, sigma) , # you have two sigmas in your fake data, should have two here, but coded only one (so I changed)
# #     mu <- a_sp[sp_id]+bstom*stomc+bsla*slac+bht*htc+bcn*cnc+bwood*woodc, 
# #     a_sp[sp_id] ~ dnorm(mu_a, sigma_sp) , # line 65 gives sigma for sp as 5 
# #     mu_a~dnorm(0, 2),
# #     bsla~dnorm(20,10),
# #     bht~dnorm(15,3),
# #     bcn~dnorm(20,10),
# #     bwood~dnorm(2.5,0.6),
# #     bstom~dnorm(400, 100),
# #     sigma ~ dnorm(0,0.5),
# #     sigma_sp ~ dnorm(0,5)),
# #   data=fake, iter=4000 , chains=4 
# # )
# 
# precis(full_in, depth=2)
# sp_mean
# 
# # Posterior Predictive Checks:
# 
# pppost <- extract.samples(full_in)
# totaln <- nrow(fake)
# 
# # Make one new dataset from model parameters and means from posterior
# newsp <- rnorm(nsp, mean(pppost$a[1:1000]), mean(pppost$sigma_sp[1:1000]))
# newspdat <- rnorm(totaln, rep(newsp, each=nsite), mean(pppost$sigma[1:1000]))
# 
# par(mfrow=c(1,2))
# hist(fake$phen_c, main="real data")
# hist(newspdat, main="PPC")
# 
# # Now do it 10 times (still using means from posterior)
# newspmat <- matrix(nrow=10, ncol=totaln)
# for (i in 1:10){
#   newsp <- rnorm(nsp, mean(pppost$a[1:1000]), mean(pppost$sigma_sp[1:1000]))
#   newspdat <- rnorm(totaln, rep(newsp, each=nsite), mean(pppost$sigma[1:1000]))
#   newspmat[i,] <- newspdat
# }
# 
# par(mfrow=c(2,6))
# hist(fake$phen_c, main="real data")
# for (i in 1:10){
#   hist(newspmat[i,], main="PPC")
# }
# 
# # Similarly ...
# library(bayesplot)
# color_scheme_set("brightblue")
# ppc_dens_overlay(fake$phen_c, newspmat[1:5, ])
# 


####################################################################
# working with the real data
####################################################################

source('Rcode/source/Cleaning_compiling_data.R') #Here I combined the budburst and trait data, averaging over individual trees for a given treatment
head(comb)
length(unique(comb$sp))
d<-comb[,c(2,12:13, 15:18)]
d <- d[ complete.cases(d), ] #there are a lot of missing cases for height
length(unique(d$sp))

d$sp_in<-coerce_index(d$sp)

head(d)
d_in<-d[,c(2:8)]
head(d_in)
length(unique(d_in$sp_in))

# range(d$cn, na.rm=TRUE)
# hist(rnorm(1000, 1, 10)) # sla values in data range from 11 to to 54
# hist(rnorm(1000, 15, 3)) # hts range from 0.20 to 20.75
# hist(rnorm(1000, 400, 100)) # stom_d values range from 30 to 770
# hist(rnorm(1000, 2.5, 0.6)) # wood_den range from 0.69 to 4.2
# hist(rnorm(1000, 20, 10)) # cn values in data range from 10.59 to 41.47

unique(d_in$sp_in) 

# full_m <- map2stan(
#   alist(
#     bday ~ dnorm(mu, sigma) , # you have two sigmas in your fake data, should have two here, but coded only one (so I changed)
#     mu <- a_sp[sp]+bstom*stom_d+bsla*sla+bht*Height+bcn*cn+bwood*wood_den, 
#     a_sp[sp] ~ dnorm(a, sigma_sp) , # line 65 gives sigma for sp as 5 
#     a~dnorm(0, 10),
#     bstom~dnorm(0, 10),
#     bsla~dnorm(0, 10),
#     bht~dnorm(0,10),
#     bcn~dnorm(0,10),
#     bwood~dnorm(0,10),
#     sigma ~ dnorm(0,1),
#     sigma_sp ~ dnorm(0,5)),
#   data=d, iter=4000 , chains=4 
# )

#issue had something to do with making the index (w/ 28) and then removing the ones with no heights, left only 26st
full_in <- map2stan(
  alist(
    bday ~ dnorm(mu, sigma) , # you have two sigmas in your fake data, should have two here, but coded only one (so I changed)
    mu <- a[sp_in]+bstom*stom_d+bsla*sla+bht*Height+bcn*cn+bwood*wood_den,
    a[sp_in] ~ dnorm(mu_a, sigma_a) , # line 65 gives sigma for sp as 5
    mu_a~dnorm(0, 50),
    bstom~dnorm(0, 50),
    bsla~dnorm(0, 50),
    bht~dnorm(0,50),
    bcn~dnorm(0,50),
    bwood~dnorm(0,50),
    sigma ~ dnorm(0,0.1),
    sigma_a ~ dnorm(0,5)),
  data=d_in, iter=4000 , chains=4
)
precis(full_in, depth=2)
# Mean StdDev lower 0.89 upper 0.89 n_eff Rhat
# a[1]     8.12   5.45      -0.39      17.01   756 1.01
# a[2]    16.99   5.18       8.88      25.25   747 1.01
# a[3]    28.77   5.13      20.78      37.09   751 1.01
# a[4]    23.17   4.72      15.76      30.67   821 1.01
# a[5]     8.78   5.21       0.34      16.81  1112 1.01
# a[6]    16.89   4.79       9.04      24.17   836 1.01
# a[7]    24.50   5.69      15.48      33.64   890 1.01
# a[8]    11.22   4.35       4.33      18.09   663 1.01
# a[9]    15.78   5.25       7.42      24.15   755 1.01
# a[10]   34.63   5.05      26.37      42.38   797 1.01
# a[11]   29.61   4.90      21.96      37.40   838 1.01
# a[12]   36.15   5.54      26.72      44.36   813 1.01
# a[13]    4.42   5.23      -3.92      12.54   742 1.01
# a[14]   20.92   5.79      11.51      30.05  1009 1.01
# a[15]    6.05   6.59      -4.29      16.63   875 1.01
# a[16]   22.68   5.83      12.93      31.51   879 1.01
# a[17]   24.41   4.60      16.86      31.47   726 1.01
# a[18]    8.33   4.95       0.20      15.93   828 1.01
# a[19]   28.38   6.88      16.96      38.78  1229 1.01
# a[20]   23.98   4.80      16.47      31.58   809 1.01
# a[21]   23.60   5.63      14.55      32.48   927 1.01
# a[22]   20.01   5.39      10.99      28.09   855 1.01
# a[23]   10.13   5.94       0.73      19.46   924 1.01
# a[24]    8.84   5.27       0.68      17.37   839 1.01
# a[25]    9.73   4.91       2.07      17.52   724 1.01
# a[26]   22.12   5.05      13.93      29.96   724 1.01
# mu_a    18.17   4.94      10.16      25.86   672 1.01
# bstom    0.01   0.00       0.00       0.02  3734 1.00
# bsla    -0.02   0.06      -0.12       0.07  1613 1.00
# bht     -0.13   0.13      -0.34       0.06  5085 1.00
# bcn      0.17   0.11      -0.01       0.35  1453 1.00
# bwood    1.79   1.54      -0.70       4.16  1477 1.00
# sigma   10.21   0.23       9.85      10.56  5425 1.00
# sigma_a  9.17   1.28       7.13      11.09  4793 1.00
plot(full_in)

par(mfrow=c(1,1))
plot(precis(full_in,depth=2))

#<><><><><><><><><><><><><><><><><>
# # the estimates are mucg better EXCEPT for wood density and the intercepts are almost all the same! 
# infprior<- map2stan(
#   alist(
#     phen_spint~ dnorm( mu , sigma ) ,
#     mu<-a[sp]+bsla*slav+bht*htv+bcn*cnv+bwood*woodv
#     +bstom*stomv
#     , 
#     a[sp] ~ dnorm( 0,20) ,
#     bsla~dnorm(20,10),
#     bht~dnorm(15,3),
#     bcn~dnorm(20,10),
#     bwood~dnorm(2.5,0.6),
#     bstom~dnorm(400, 100),# the intercepts all become uniform when this variable is added
#     sigma ~ dnorm(0,10)),
#   data=fake_spint, iter=4000 , chains=4 
# )
# 
# precis(infprior, depth=2)
# plot(precis(infprior, depth=2))
# 

#########################
# Plots: based on code from CH 5 box 5.11, 5.12
#########################

par(mfrow=c(1,1))

mu_ppc <- link( full_in)
mu.mean <- apply( mu_ppc , 2 , mean )
mu.PI <- apply( mu_ppc , 2 , PI )

#simulate data
phen.sim <- sim( full_in , n=1000 )
phen.PI <- apply( phen.sim , 2 , PI )

plot( mu.mean ~ d$bday , col="darkgreen" , ylim=range(mu.PI) ,
      xlab="Observed phen" , ylab="Predicted phen" )
for ( i in 1:nrow(d) )
  lines( rep(d$bday[i],2) , c(mu.PI[1,i],mu.PI[2,i]) ,
         col="darkgreen" )
abline( a=0 , b=1 , lty=2, lwd=4)
# This looks pretty good, the data is fairly aligned with the perfect fit line#

##########################

# # compute residuals
# phen.resid <- d$bday - mu.mean
# # get ordering by divorce rate
# ord <- order(phen.resid)
# # make the plot
# dotchart( phen.resid[ord] , labels=d$sp[ord]  , cex=0.6 )
# abline( v=0 , col=col.alpha("black",0.2) )
# for ( i in 1:nrow(fake_spint) ) {
#   j <- ord[i] # which State in order
#   lines( fake_spint$phen_spint[j]-c(mu.PI[1,j],mu.PI[2,j]) , rep(i,2) )
#   points( fake_spint$phen_spint[j]-c(phen.PI[1,j],phen.PI[2,j]) , rep(i,2),
#           pch=3 , cex=0.6 , col="gray" )
# }

####################################################
# Visualize the partial pooling ...
# See https://www.tjmahr.com/plotting-partial-pooling-in-mixed-effects-models/
##
coerce_index(unique(d$sp))
d_sort <- d_in[order(d_in$sp_in),]
#d_sort
post <- extract.samples(full_in) # grab the posterior

nsp=26
J <- nsp #once NA's are removed there are only 26

com_pool_mod <- lm(bday ~ 1, d_sort)
species <- unique(d_sort$sp_id)
no_pool <- data.frame(sp_in=unique(d_sort$sp_in),
                      intercept=rep(NA, J), 
                      slope=rep(NA, J))
head(d_sort)
with_pool <- data.frame(sp_in=unique(d_sort$sp_in),
                        intercept=rep(NA, J), 
                        slope=rep(NA, J))
df.gravity <- no_pool[1:2,2:3]
with_pool$model <- "partial pooling"
no_pool$model <- "no pooling"

#d_sort
data<-subset(d_sort, sp_in==unique(d_sort$sp_in)[1])
#data

#data=subset(vino.formodel.alt.sort, judge==unique(vino.formodel.alt.sort$judge)[j]))
for (j in 1:J){
  modhere <- lm(bday~1,data=subset(d_sort, sp_in==unique(d_sort$sp_in)[j]))
  no_pool$intercept[j] <- coef(modhere)[1]
  no_pool$slope[j] <- j # cheat for no-slope model
}

for (j in 1:J){
  with_pool$intercept[j] <- mean(post$a[1:1000, j])
  with_pool$slope[j] <- j # cheat for no-slope model
}

# head(no_pool)
# head(with_pool)

df.pulled <- rbind(no_pool, with_pool)

df.gravity$model <- NA
df.gravity$intercept[1] <- coef(com_pool_mod)["(Intercept)"]
df.gravity$slope[1] <- 4.5 # cheat for no-slope model: middle of plot
df.gravity$model[1] <- "complete pooling"
df.gravity$intercept[2] <- mean(post$mu_a[1:1000])
df.gravity$slope[2] <- 4.5 # cheat for no-slope model: middle of plot
df.gravity$model[2] <- "partial pooling (mu)"

library(ggplot2)
#pdf("modelscompare_pp.pdf", width = 9, height = 6)
ggplot(df.pulled) + 
  aes(x = intercept, y = slope, color = model) + 
  geom_point(size = 2) + 
  geom_point(data = df.gravity, size = 5) + 
  # Draw an arrow connecting the observations between models
  geom_path(aes(group = as.character(sp_in), color = NULL), 
            arrow = arrow(length = unit(.02, "npc"))) + 
  # Use ggrepel to jitter the labels away from the points
  ggrepel::geom_text_repel(
    aes(label = sp_in, color = NULL), 
    data = no_pool, size=2) + 
  theme(legend.position = "bottom") + 
  ggtitle("Pooling of regression parameters") + 
  xlab("Intercept estimate") + 
  ylab("") + 
  scale_color_brewer(palette = "Dark2") 
#dev.off()

######################################################################
######################################################################
# Posterior Predictive Checks:

pppost <- extract.samples(full_in)
totaln <- nrow(d_sort)
nsite=2
# Make one new dataset from model parameters and means from posterior
newsp <- rnorm(nsp, mean(pppost$mu_a[1:1000]), mean(pppost$sigma_a[1:1000]))
newspdat <- rnorm(totaln, rep(newsp, each=nsite), mean(pppost$sigma[1:1000]))

par(mfrow=c(1,2))
hist(d_sort$bday, main="real data")
hist(newspdat, main="PPC")

# Now do it 10 times (still using means from posterior)
newspmat <- matrix(nrow=10, ncol=totaln)
for (i in 1:10){
  newsp <- rnorm(nsp, mean(pppost$mu_a[1:1000]), mean(pppost$sigma_a[1:1000]))
  newspdat <- rnorm(totaln, rep(newsp, each=nsite), mean(pppost$sigma[1:1000]))
  newspmat[i,] <- newspdat
}
newspmat
par(mfrow=c(2,6))
hist(d_sort$bday, main="real data")
for (i in 1:10){
  hist(newspmat[i,], main="PPC")
}

# Similarly ...
library(bayesplot)
color_scheme_set("brightblue")
ppc_dens_overlay(d_sort$bday, newspmat[1:5, ])


