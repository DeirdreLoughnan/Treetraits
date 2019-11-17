rm(list=ls()) 
options(stringsAsFactors = FALSE)

require(doBy)

setwd("~/Documents/github/Treetraits/input/western")
#leaf area
la<-read.csv("leafarea_Oct162019.csv")

#leaf mass
lm<-read.csv("traits_mass.csv")
unique(la$site)
summaryBy(area~ , data=dframe, FUN=c(mean, sd))

#dbh, ht, etc
etc<-read.csv("traits_ht_dbh_master.csv")
unique(etc$site)
names(etc);names(la)
?merge
test<-merge(etc,la, by=c("site","species","indiv.no"))
head(test)

write.csv(test, file="complete_data_Nov11.csv")

