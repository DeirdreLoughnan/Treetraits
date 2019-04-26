### Started March 30 2019 ###

## DL writing final code for creation and testing of fake data for 507 Project IV ##

#The aim of this code is to build test data for modeling the relationship between phenology and functional traits
#Update: Added partial pooling across species and multiple PPC

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

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Slide 3 content
int<-11 # days into the experiment of bb 
sigma<- 0.1

nsite = 2 #number of sites
nsp = 28 #number of species
rep = 50 #
ntot<-nsite*nsp*rep #2800

#Building the required dataframe for all of the replication and traits
# Trait values (trait name + v for value)
slav=rnorm(ntot, 5, 1)
htv=rnorm(ntot, 11, 3) 
cnv=rnorm(ntot, 10, 2)
woodv=rnorm(ntot, 0.7, 0.05)
stomv=rnorm(ntot, 200, 50)


#effect sizes
#site.diff=2 
sladiff= -0.5
htdiff=0.5
cndiff=-0.5
wooddiff=1
stomdiff=1
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

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

# Generatng test phenology data
fake_cent<-vector()
sp_mean <- rnorm(nsp, mean=int, sd=1);
sp_mean<-sp_mean[order=T]
for (j in c(1:length(sp_mean))){
  phen_c<- sp_mean[j]+ rnorm(rep*nsite, 0, sigma)+stomdiff*stomc+sladiff*slac+htdiff*htc+cndiff*cnc+wooddiff*woodc
  temp<-data.frame(phen_c)
  fake_cent<- rbind(fake_cent, temp) # changed so not overwriting fake_stom every time ... 
}

#Add species - there is probably a better way to do this
sp=gl(nsp, rep*nsite, length= ntot, order=TRUE)

#Combine the phenology data with the species identity 
phen<-cbind(fake_cent,sp)
hist(phen$phen_c)
#Based on R.code 5.40: this seems to work and keep the order of the intercepts consistent
phen$sp_in<-as.integer(phen$sp)

###############################################
# Now combining the phenology and the trait data
head(phen)
fake_data<- data.frame(phen,stomc,slac,htc,cnc,woodc)
head(fake_data)

#fake_stom$sp_id <- coerce_index(fake_stom$sp)
fake<-fake_data[,c(1,3:8)]
head(fake)
#Data looks good!

#Test Data: Partial pooling across species
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#Table 1: Model output for fake data 
full_f <- map2stan(
  alist(
    phen_c ~ dnorm(mu, sigma) , 
    mu <- a[sp_in]+bstom*stomc+bsla*slac+bht*htc+bcn*cnc+bwood*woodc,
    a[sp_in] ~ dnorm(mu_a, sigma_a) , 
    mu_a~dnorm(0, 50),
    bstom~dnorm(0, 50),
    bsla~dnorm(0, 50),
    bht~dnorm(0,50),
    bcn~dnorm(0,50),
    bwood~dnorm(0,50),
    sigma ~ dnorm(0,5),
    sigma_a ~ dnorm(0,5)),
  data=fake, iter=4000 , chains=4
)

plot(full_f)

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
post <- extract.samples(full_in) # grab the posterior
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


df.gravity$intercept[2] <- coef(modhere)["mu_a"]
df.gravity$stomc[2] <- coef(modhere)["bstom"]
df.gravity$slac[2] <- coef(modhere)["bsla"]
df.gravity$htc[2] <- coef(modhere)["bht"]
df.gravity$cnc[2] <- coef(modhere)["bcn"]
df.gravity$woodc[2] <- coef(modhere)["bwood"]
df.gravity$model[2] <- "partial pooling (mu)"

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
