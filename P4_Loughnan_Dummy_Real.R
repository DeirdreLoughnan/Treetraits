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
phen$sp_in<-as.integer(phen$sp)

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
str(fake)


full_f <- map2stan(
  alist(
    phen_c ~ dnorm(mu, sigma) , # you have two sigmas in your fake data, should have two here, but coded only one (so I changed)
    mu <- a[sp_in]+bstom*stomc+bsla*slac+bht*htc+bcn*cnc+bwood*woodc,
    a[sp_in] ~ dnorm(mu_a, sigma_a) , 
    mu_a~dnorm(0, 50),
    bstom~dnorm(0, 50),
    bsla~dnorm(0, 50),
    bht~dnorm(0,50),
    bcn~dnorm(0,50),
    bwood~dnorm(0,50),
    sigma ~ dnorm(0,0.1),
    sigma_a ~ dnorm(0,5)),
  data=fake, iter=4000 , chains=4
)

#plot(full_f)

sp_mean
stomdiff # 1
sladiff #-0.5
htdiff #0.5
cndiff #-0.5
wooddiff #1
sigma #0.1
int #11

precis(full_f, depth=2)

par(mfrow=c(1,1))
plot(precis(full_f, depth=2))

######################################################################
######################################################################
# Posterior Predictive Checks:

pppost <- extract.samples(full_f)
totalnf <- nrow(fake)
nsp<-28
reps=50
# Make one new dataset from model parameters and means from posterior
newsp <- rnorm(nsp, mean(pppost$mu_a[1:1000]), mean(pppost$sigma_a[1:1000]))
newspdat <- rnorm(totalnf, rep(newsp, each=rep), mean(pppost$sigma[1:1000]))

par(mfrow=c(1,2))
hist(fake$phen_c, main="real data")
hist(newspdat, main="PPC")

# Now do it 10 times (still using means from posterior)
newmat <- matrix(nrow=10, ncol=totalnf)
for (i in 1:10){
  newsp <- rnorm(nsp, mean(pppost$mu_a[1:1000]), mean(pppost$sigma_a[1:1000]))
  newspdat <- rnorm(totalnf, rep(newsp, each=reps), mean(pppost$sigma[1:1000]))
  newmat[i,] <- newspdat
}

par(mfrow=c(2,6))
hist(fake$phen_c, main="real data")
for (i in 1:10){
  hist(newmat[1,], main="PPC")
}


# Similarly ...
library(bayesplot)
color_scheme_set("brightblue")
ppc_dens_overlay(fake$phen_c, newmat)

dim(newmat)
#############################################################
##############################################################

# Visualize the partial pooling ...
# See https://www.tjmahr.com/plotting-partial-pooling-in-mixed-effects-models/
##
coerce_index(unique(fake$sp)); str(fake)
#d_sort <- fake[order(fake$sp_id),]
post <- extract.samples(full_f) # grab the posterior
str(post)

unique(fake$sp_in)
nsp=28
J <- nsp 

com.pool.mod <- lm(phen_c ~ stomc+slac+htc+cnc+woodc, fake)
spp <- sort(unique(fake$sp_in)); spp
no.pool <- data.frame(sp_in=rep(NA, length(spp)),
                      intercept=rep(NA, length(spp)), 
                      stomc=rep(NA, length(spp)), 
                      slac=rep(NA, length(spp)), 
                      htc=rep(NA, length(spp)), 
                      cnc=rep(NA, length(spp)), 
                      woodc=rep(NA, length(spp))
                      )

with.pool <- no.pool
dim(no.pool)
df.gravity <- no.pool[1:2,2:7] #why is it 1:2?
with.pool$model <- "partial pooling"
no.pool$model <- "no pooling"

for (sp in c(1:length(spp))){
  no.pool$sp_in[sp] <- spp[sp]
  subby <- subset(fake,  sp_in==spp[sp])
  lmfit <- lm(phen_c ~ stomc+slac+htc+cnc+woodc, data=subby)
  no.pool$intercept[sp] <- coef(lmfit)["(Intercept)"]
  no.pool$stomc[sp] <- coef(lmfit)["stomc"]
  no.pool$slac[sp] <- coef(lmfit)["slac"]
  no.pool$htc[sp] <- coef(lmfit)["htc"]
  no.pool$cnc[sp] <- coef(lmfit)["cnc"]
  no.pool$woodc[sp] <- coef(lmfit)["woodc"]
}

# sumer.ni <- summary(full_f)$summary
# sumer.ni[grep("mu_a", rownames(sumer.ni)),]
# 
# require(lme4)
# 
# modhere <- lmer(phen_c~stomc+slac+htc+cnc+woodc+(1|sp_in),data=fake)
# 
# for (sp in c(1:length(spp))){
#   with.pool$sp_in[sp] <- spp[sp]
#   with.pool$intercept[sp] <- coef(modhere)["Intercept"]; with.pool$intercept
#   with.pool$stomc[sp] <- coef(modhere)["stomc"]
#   with.pool$slac[sp] <- coef(modhere)["slac"]
#   with.pool$htc[sp] <- coef(modhere)["htc"]
#   with.pool$cnc[sp] <- coef(modhere)["cnc"]
#   with.pool$woodc[sp] <- coef(modhere)["woodc"]
# }


# for (sp in c(1:length(spp))){
#   with.pool$sp_in[sp] <- spp[sp]
#   subby <- subset(fake,  sp_in==spp[sp])
#   lmfit <- lmer(phen_c ~ stomc+slac+htc+cnc+woodc+(1|sp_in), data=fake)
#   with.pool$intercept[sp] <- fixef(lmfit)["(Intercept)"]
#   with.pool$stomc[sp] <- fixef(lmfit)["stomc"]
#   with.pool$slac[sp] <- fixef(lmfit)["slac"]
#   with.pool$htc[sp] <- fixef(lmfit)["htc"]
#   with.pool$cnc[sp] <- fixef(lmfit)["cnc"]
#   with.pool$woodc[sp] <- fixef(lmfit)["woodc"]
# }

#What should these values be, the coefficients or the fixed eff?
 modhere<-full_f
# coef(modhere)[1:28]
# for (sp in c(1:length(spp))){
#   with.pool$sp_in[sp] <- spp[sp]
#   with.pool$intercept[spp] <- coef(modhere)[1:28]
#   with.pool$stomc[spp] <- coef(modhere)["bstom"]
#   with.pool$slac[spp] <- coef(modhere)["bsla"]
#   with.pool$htc[spp] <-  coef(modhere)["bht"]
#   with.pool$cnc[spp] <- coef(modhere)["bcn"]
#   with.pool$woodc[spp] <- coef(modhere)["bwood"]
# }

no.pool$sp_no <- sort(unique(fake$sp_in))
with.pool$sp_no <- sort(unique(fake$sp_in))
df.pulled <- rbind(no.pool, with.pool)
#df.pulled <- bind_rows(no.pool, with.pool)

df.gravity$model <- NA
df.gravity$intercept[1] <-coef(com.pool.mod)["(Intercept)"]
df.gravity$stomc[1] <-coef(com.pool.mod)["stomc"]
df.gravity$slac[1] <-coef(com.pool.mod)["slac"]
df.gravity$htc[1] <-coef(com.pool.mod)["htc"]
df.gravity$cnc[1] <-coef(com.pool.mod)["cnc"]
df.gravity$woodc[1] <-coef(com.pool.mod)["woodc"]
df.gravity$model[1] <- "complete pooling"

# df.gravity$intercept[2] <- fixef(lmfit)["(Intercept)"]
# df.gravity$stomc[2] <- fixef(lmfit)["stomc"]
# df.gravity$slac[2] <- fixef(lmfit)["slac"]
# df.gravity$htc[2] <- fixef(lmfit)["htc"]
# df.gravity$cnc[2] <- fixef(lmfit)["ccc"]
# df.gravity$woodc[2] <- fixef(lmfit)["woodc"]
# df.gravity$model[2] <- "partial pooling (mu)"
modhere<-full_f
df.gravity$intercept[2] <- coef(modhere)["mu_a"]
df.gravity$stomc[2] <- coef(modhere)["bstom"]
df.gravity$slac[2] <- coef(modhere)["bsla"]
df.gravity$htc[2] <- coef(modhere)["bht"]
df.gravity$cnc[2] <- coef(modhere)["bcn"]
df.gravity$woodc[2] <- coef(modhere)["bwood"]
df.gravity$model[2] <- "partial pooling (mu)"

# ggplot(df.pulled) + 
#   aes(x = intercept, y = slope, color = model) + 
#   geom_point(size = 2) + 
#   geom_point(data = df.gravity, size = 5) + 
#   # Draw an arrow connecting the observations between models
#   geom_path(aes(group = as.character(sp_in), color = NULL), 
#             arrow = arrow(length = unit(.02, "npc"))) + 
#   # Use ggrepel to jitter the labels away from the points
#   ggrepel::geom_text_repel(
#     aes(label = sp_in, color = NULL), 
#     data = no_pool, size=2) + 
#   theme(legend.position = "bottom") + 
#   ggtitle("Pooling of regression parameters") + 
#   xlab("Intercept estimate") + 
#   ylab("") + 
#   scale_color_brewer(palette = "Dark2") 

ggplot(df.pulled) + 
  aes(x = intercept, y = slac, color = model)+ 
  geom_point(size = 2) + 
  geom_point(data = df.gravity, size = 5) + 
  # Draw an arrow connecting the observations between models
  geom_path(aes(group = as.character(sp_no), color = NULL), 
            arrow = arrow(length = unit(.02, "npc"))) + 
  # Use ggrepel to jitter the labels away from the points
  ggrepel::geom_text_repel(
    aes(label = sp_no, color = NULL), 
    data = no.pool, size=2) + 
  theme(legend.position = "bottom") + 
  ggtitle("Pooling of regression parameters: Int and SLA") + 
  xlab("Intercept estimate") + 
  ylab("Slope estimate") + 
  scale_color_brewer(palette = "Dark2") 


pdf(file.path(figpath, "modelscompare_pp_photo.pdf"), width = 9, height = 6)
ggplot(df.pulled) + 
  aes(x = intercept, y = photo, color = model) + 
  geom_point(size = 2) + 
  geom_point(data = df.gravity, size = 5) + 
  # Draw an arrow connecting the observations between models
  geom_path(aes(group = as.character(complex), color = NULL), 
            arrow = arrow(length = unit(.02, "npc"))) + 
  # Use ggrepel to jitter the labels away from the points
  ggrepel::geom_text_repel(
    aes(label = complex.wname, color = NULL), 
    data = no.pool, size=2) + 
  theme(legend.position = "bottom") + 
  ggtitle("Pooling of regression parameters: Int and photo") + 
  xlab("Intercept estimate") + 
  ylab("Slope estimate") + 
  scale_color_brewer(palette = "Dark2") 

#data<-subset(fake, sp_in==unique(fake$sp_in)[1])
#data=subset(vino.formodel.alt.sort, judge==unique(vino.formodel.alt.sort$judge)[j]))
for (j in 1:J){
  #no_pool$sp_in[sp]<-species[sp]
  modhere <- lm(phen_c~stomc+slac+htc+cnc+woodc,data=subset(fake, sp_in==unique(fake$sp_in)[j]))
  no_pool$intercept[j] <- coef(modhere)["Intercept"]
  no_pool$slope[j] <- coef(modhere)["slope"]
}

for (j in 1:J){
  with_pool$intercept[j] <- mean(post$a_sp[1:1000, j])
  with_pool$slope[j] <- coef(modhere)["slope"] 
}

# head(no_pool)
# head(with_pool)

df.pulled <- rbind(no_pool, with_pool)

df.gravity$model <- NA
df.gravity$intercept[1] <- coef(com_pool_mod)["(Intercept)"]
df.gravity$slope[1] <- 4.5 # cheat for no-slope model: middle of plot
df.gravity$model[1] <- "complete pooling"
df.gravity$intercept[2] <- mean(post$a[1:1000])
df.gravity$slope[2] <-mean(post$mu_a[1:1000])
df.gravity$model[2] <- "partial pooling (mu)"

# require(lme4)
# df_fixef <- data_frame(
#   Model = "Partial pooling (average)",
#   Intercept = fixef()[1],
#   Slope_Days = fixef(full_f)[2])

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

#issue had something to do with making the index (w/ 28) and then removing the ones with no heights, left only 26sp
full_in <- map2stan(
  alist(
    bday ~ dnorm(mu, sigma) , 
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

plot(full_in)

par(mfrow=c(1,1))
plot(precis(full_in,depth=2))


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


#########################################################################
# Visualize the partial pooling ...
# See https://www.tjmahr.com/plotting-partial-pooling-in-mixed-effects-models/
##
coerce_index(unique(d$sp))
#d_sort <- fake[order(fake$sp_id),]
post <- extract.samples(full_in) # grab the posterior
str(post)

unique(d$sp)
nsp=26
J <- nsp 

com.pool.mod <- lm(bday ~ stom_d+sla+Height+cn+wood_den, d_in)
spp <- sort(unique(d$sp_in)); spp
no.pool <- data.frame(sp_in=rep(NA, length(spp)),
                      intercept=rep(NA, length(spp)), 
                      stom=rep(NA, length(spp)), 
                      sla=rep(NA, length(spp)), 
                      ht=rep(NA, length(spp)), 
                      cn=rep(NA, length(spp)), 
                      wood=rep(NA, length(spp))
)
no.pool
with.pool <- no.pool
dim(no.pool)
df.gravity <- no.pool[1:2,2:7] #why is it 1:2?
with.pool$model <- "partial pooling"
no.pool$model <- "no pooling"

for (sp in c(1:length(spp))){
  no.pool$sp_in[sp] <- spp[sp]
  subby <- subset(d_in,  sp_in==spp[sp])
  lmfit <- lm(bday~ stom_d+sla+Height+cn+wood_den, data=subby)
  no.pool$intercept[sp] <- coef(lmfit)["(Intercept)"]
  no.pool$stom[sp] <- coef(lmfit)["stom_d"]
  no.pool$sla[sp] <- coef(lmfit)["sla"]
  no.pool$ht[sp] <- coef(lmfit)["Height"]
  no.pool$cn[sp] <- coef(lmfit)["cn"]
  no.pool$wood[sp] <- coef(lmfit)["wood_den"]
}
no.pool
# sumer.ni <- summary(full_f)$summary
# sumer.ni[grep("mu_a", rownames(sumer.ni)),]
# 
# require(lme4)
# 
# modhere <- lmer(phen_c~stomc+slac+htc+cnc+woodc+(1|sp_in),data=fake)
# 
# for (sp in c(1:length(spp))){
#   with.pool$sp_in[sp] <- spp[sp]
#   with.pool$intercept[sp] <- coef(modhere)["Intercept"]; with.pool$intercept
#   with.pool$stomc[sp] <- coef(modhere)["stomc"]
#   with.pool$slac[sp] <- coef(modhere)["slac"]
#   with.pool$htc[sp] <- coef(modhere)["htc"]
#   with.pool$cnc[sp] <- coef(modhere)["cnc"]
#   with.pool$woodc[sp] <- coef(modhere)["woodc"]
# }


# for (sp in c(1:length(spp))){
#   with.pool$sp_in[sp] <- spp[sp]
#   subby <- subset(fake,  sp_in==spp[sp])
#   lmfit <- lmer(phen_c ~ stomc+slac+htc+cnc+woodc+(1|sp_in), data=fake)
#   with.pool$intercept[sp] <- fixef(lmfit)["(Intercept)"]
#   with.pool$stomc[sp] <- fixef(lmfit)["stomc"]
#   with.pool$slac[sp] <- fixef(lmfit)["slac"]
#   with.pool$htc[sp] <- fixef(lmfit)["htc"]
#   with.pool$cnc[sp] <- fixef(lmfit)["cnc"]
#   with.pool$woodc[sp] <- fixef(lmfit)["woodc"]
# }

#What should these values be, the coefficients or the fixed eff?
modhere<-full_f
coef(modhere)[1:28]
for (sp in c(1:length(spp))){
  with.pool$sp_in[sp] <- spp[sp]
  with.pool$intercept[spp] <- coef(modhere)[1:28]
  with.pool$stomc[spp] <- coef(modhere)["bstom"]
  with.pool$slac[spp] <- coef(modhere)["bsla"]
  with.pool$htc[spp] <-  coef(modhere)["bht"]
  with.pool$cnc[spp] <- coef(modhere)["bcn"]
  with.pool$woodc[spp] <- coef(modhere)["bwood"]
}

no.pool$sp_no <- sort(unique(d$sp_in))
with.pool$sp_no <- sort(unique(d$sp_in))
df.pulled <- rbind(no.pool, with.pool)
#df.pulled <- bind_rows(no.pool, with.pool)

df.gravity$model <- NA
df.gravity$intercept[1] <-coef(com.pool.mod)["(Intercept)"]
df.gravity$stom[1] <-coef(com.pool.mod)["stom_d"]
df.gravity$sla[1] <-coef(com.pool.mod)["sla"]
df.gravity$ht[1] <-coef(com.pool.mod)["Height"]
df.gravity$cn[1] <-coef(com.pool.mod)["cn"]
df.gravity$wood[1] <-coef(com.pool.mod)["wood_den"]
df.gravity$model[1] <- "complete pooling"

# df.gravity$intercept[2] <- fixef(lmfit)["(Intercept)"]
# df.gravity$stomc[2] <- fixef(lmfit)["stomc"]
# df.gravity$slac[2] <- fixef(lmfit)["slac"]
# df.gravity$htc[2] <- fixef(lmfit)["htc"]
# df.gravity$cnc[2] <- fixef(lmfit)["ccc"]
# df.gravity$woodc[2] <- fixef(lmfit)["woodc"]
# df.gravity$model[2] <- "partial pooling (mu)"
modhere<-full_in
df.gravity$intercept[2] <- coef(modhere)["mu_a"]
df.gravity$stomc[2] <- coef(modhere)["bstom"]
df.gravity$slac[2] <- coef(modhere)["bsla"]
df.gravity$htc[2] <- coef(modhere)["bht"]
df.gravity$cnc[2] <- coef(modhere)["bcn"]
df.gravity$woodc[2] <- coef(modhere)["bwood"]
df.gravity$model[2] <- "partial pooling (mu)"

# ggplot(df.pulled) + 
#   aes(x = intercept, y = slope, color = model) + 
#   geom_point(size = 2) + 
#   geom_point(data = df.gravity, size = 5) + 
#   # Draw an arrow connecting the observations between models
#   geom_path(aes(group = as.character(sp_in), color = NULL), 
#             arrow = arrow(length = unit(.02, "npc"))) + 
#   # Use ggrepel to jitter the labels away from the points
#   ggrepel::geom_text_repel(
#     aes(label = sp_in, color = NULL), 
#     data = no_pool, size=2) + 
#   theme(legend.position = "bottom") + 
#   ggtitle("Pooling of regression parameters") + 
#   xlab("Intercept estimate") + 
#   ylab("") + 
#   scale_color_brewer(palette = "Dark2") 

ggplot(df.pulled) + 
  aes(x = intercept, y = sla, color = model)+ 
  geom_point(size = 2) + 
  geom_point(data = df.gravity, size = 5) + 
  # Draw an arrow connecting the observations between models
  geom_path(aes(group = as.character(sp_in), color = NULL), 
            arrow = arrow(length = unit(.02, "npc"))) + 
  # Use ggrepel to jitter the labels away from the points
  ggrepel::geom_text_repel(
    aes(label = sp_in, color = NULL), 
    data = no.pool, size=2) + 
  theme(legend.position = "bottom") + 
  ggtitle("Pooling of regression parameters: Int and forcing") + 
  xlab("Intercept estimate") + 
  ylab("Slope estimate") + 
  scale_color_brewer(palette = "Dark2") 
