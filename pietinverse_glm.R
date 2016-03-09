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
library(gdata)
library(visreg)

#this code contains the details of the models used for the analysis of the 
#impact of the site and the treatment on the resistance of the strains


###############################################################################
#
###############################################################################


prochlo_dat<-datapietin[datapietin$Espece!="Oa" & datapietin$Modalite %in% 
                          c("Temoin","Prochloraze"),]
prochlo_dat<-drop.levels(prochlo_dat)
modProchlo<-glm(Prochloraze~ESSAI+Modalite,family=binomial,prochlo_dat)
summary(modProchlo)
barplot(table(prochlo_dat$Prochloraze,prochlo_dat$Modalite),beside=TRUE)
barplot(table(prochlo_dat$Prochloraze,prochlo_dat$ESSAI))

barplot(table(prochlo_dat$Modalite,prochlo_dat$Prochloraze),beside=TRUE)
barplot(table(prochlo_dat$ESSAI,prochlo_dat$Prochloraze),beside=TRUE)



