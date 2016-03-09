###############################################################################
###############################################################################
#Effect of site and treatment on resistance of Pietin-verse
###############################################################################
###############################################################################

#set the right working directory
setwd("~/work/Rfichiers/Githuber/PietinVerse_data")

datapietin<-read.table("data_pietinverse.txt",header=TRUE,sep="\t",
                       colClasses = c("numeric","character","numeric",
                                      "factor","factor","factor","factor",
                                      "factor","factor","factor")
                         )

library(lme4)
library()
library(MASS)
library(car)
library(visreg)

#this code contains the details of the models used for the analysis of the 
#impact of the site and the treatment on the resistance of the strains


###############################################################################
#
###############################################################################


prochlo_dat<-datapietin[datapietin$Espece!="Oa" & datapietin$Modalite %in% 
                          c("Temoin","Prochloraze"),]
modProchlo<-glm(Prochloraze~ESSAI+Modalite,family=binomial,prochlo_dat)
summary(modProchlo)
table(prochlo_dat$Prochloraze,prochlo_dat$Modalite)
table(prochlo_dat$Prochloraze,prochlo_dat$ESSAI)




