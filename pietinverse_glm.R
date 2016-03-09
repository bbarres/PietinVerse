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

#this code contains the details of the models used for the analysis of the 
#impact of the site and the treatment on the resistance of the strains


###############################################################################
#Analysis of the effect of Prochloraze on Prochloraze resistance
###############################################################################

#restricted dataset for Prochloraze effect analysis
prochlo_dat<-datapietin[datapietin$Espece!="Oa" & datapietin$Modalite %in% 
                          c("Temoin","Prochloraze"),]
prochlo_dat<-drop.levels(prochlo_dat)

#the different models considered
modProchlo1<-glm(Prochloraze~ESSAI,family=binomial,prochlo_dat)
modProchlo2<-glm(Prochloraze~Modalite,family=binomial,prochlo_dat)
modProchlo<-glm(Prochloraze~ESSAI+Modalite,family=binomial,prochlo_dat)
anova(modProchlo1,modProchlo,test="Chisq")
anova(modProchlo2,modProchlo,test="Chisq")
drop1(modProchlo,test="Chisq")
summary(modProchlo)

#barplot by Traitement and Essai
barplot(table(prochlo_dat$Prochloraze,prochlo_dat$Modalite),beside=TRUE,
        col=c("green3","red2"),main="Par traitement")
barplot(table(prochlo_dat$Prochloraze,prochlo_dat$ESSAI),beside=TRUE,
        col=c("green3","red2"),main="Par essai")

#number of Prochloraze sensitive/resistant strains
table(prochlo_dat$Prochloraze,prochlo_dat$Modalite)

#barplot "temoin" against "treated
barplot(table(prochlo_dat[prochlo_dat$Prochloraze==1,]$ESSAI,
              prochlo_dat[prochlo_dat$Prochloraze==1,]$Modalite)[,c(2,1)]/
          table(prochlo_dat$ESSAI,prochlo_dat$Modalite)[,c(2,1)]*100,
        beside=TRUE,ylim=c(0,100))


###############################################################################
#Analysis of the effect of Cyprodinil on Cyprodinil resistance
###############################################################################

#restricted dataset for Prochloraze effect analysis
cyprod_dat<-datapietin[datapietin$Espece!="Oa" & datapietin$Modalite %in% 
                          c("Temoin","Cyprodinil"),]
cyprod_dat<-drop.levels(cyprod_dat)

#the different models considered
modCyprod1<-glm(Cyprodinil~ESSAI,family=binomial,cyprod_dat)
modCyprod2<-glm(Cyprodinil~Modalite,family=binomial,cyprod_dat)
modCyprod<-glm(Cyprodinil~ESSAI+Modalite,family=binomial,cyprod_dat)
anova(modCyprod1,modCyprod,test="Chisq")
anova(modCyprod2,modCyprod,test="Chisq")
drop1(modCyprod,test="Chisq")
summary(modCyprod)

#barplot by Traitement and Essai
barplot(table(cyprod_dat$Cyprodinil,cyprod_dat$Modalite),beside=TRUE,
        col=c("green3","red2"),main="Par traitement")
barplot(table(cyprod_dat$Cyprodinil,cyprod_dat$ESSAI),beside=TRUE,
        col=c("green3","red2"),main="Par essai")

#number of Cyprodinil sensitive/resistant strains
table(cyprod_dat$Cyprodinil,cyprod_dat$Modalite)

#barplot "temoin" against "treated
barplot(table(cyprod_dat[cyprod_dat$Cyprodinil==1,]$ESSAI,
              cyprod_dat[cyprod_dat$Cyprodinil==1,]$Modalite)[,c(2,1)]/
          table(cyprod_dat$ESSAI,cyprod_dat$Modalite)[,c(2,1)]*100,
        beside=TRUE,ylim=c(0,100))


###############################################################################
#Analysis of the effect of Boscalid_epoxiconazole on MDR resistance
###############################################################################

#restricted dataset for Prochloraze effect analysis
bosca_dat<-datapietin[datapietin$Espece!="Oa" & datapietin$Modalite %in% 
                         c("Temoin","Boscalid_epoxiconazole") & 
                        datapietin$ESSAI!="ARVALIS17",]
bosca_dat<-drop.levels(bosca_dat)

#the different models considered
modBosca1<-glm(MDR~ESSAI,family=binomial,bosca_dat)
modBosca2<-glm(MDR~Modalite,family=binomial,bosca_dat)
modBosca<-glm(MDR~ESSAI+Modalite,family=binomial,bosca_dat)
anova(modBosca1,modBosca,test="Chisq")
anova(modBosca2,modBosca,test="Chisq")
drop1(modBosca,test="Chisq")
summary(modBosca)

#barplot by Traitement and Essai
barplot(table(bosca_dat$MDR,bosca_dat$Modalite),beside=TRUE,
        col=c("green3","red2"),main="Par traitement")
barplot(table(bosca_dat$MDR,bosca_dat$ESSAI),beside=TRUE,
        col=c("green3","red2"),main="Par essai")

#number of MDR sensitive/resistant strains
table(bosca_dat$MDR,bosca_dat$Modalite)

#barplot "temoin" against "treated"
barplot(table(bosca_dat[bosca_dat$MDR==1,]$ESSAI,
              bosca_dat[bosca_dat$MDR==1,]$Modalite)[,c(2,1)]/
          table(bosca_dat$ESSAI,bosca_dat$Modalite)[,c(2,1)]*100,
        beside=TRUE,ylim=c(0,100))


###############################################################################
#Analysis of the effect of Prothioconazole on MDR resistance
###############################################################################

#restricted dataset for Prochloraze effect analysis
prothio_dat<-datapietin[datapietin$Espece!="Oa" & datapietin$Modalite %in% 
                        c("Temoin","Prothioconazole") & 
                          datapietin$ESSAI!="ARVALIS17",]
prothio_dat<-drop.levels(prothio_dat)

#the different models considered
modprothio1<-glm(MDR~ESSAI,family=binomial,prothio_dat)
modprothio2<-glm(MDR~Modalite,family=binomial,prothio_dat)
modprothio<-glm(MDR~ESSAI+Modalite,family=binomial,prothio_dat)
anova(modprothio1,modprothio,test="Chisq")
anova(modprothio2,modprothio,test="Chisq")
drop1(modprothio,test="Chisq")
summary(modprothio)

#barplot by Traitement and Essai
barplot(table(prothio_dat$MDR,prothio_dat$Modalite),beside=TRUE,
        col=c("green3","red2"),main="Par traitement")
barplot(table(prothio_dat$MDR,prothio_dat$ESSAI),beside=TRUE,
        col=c("green3","red2"),main="Par essai")

#number of MDR sensitive/resistant strains
table(prothio_dat$MDR,prothio_dat$Modalite)

#barplot "temoin" against "treated"
barplot(table(prothio_dat[prothio_dat$MDR==1,]$ESSAI,
              prothio_dat[prothio_dat$MDR==1,]$Modalite)[,c(2,1)]/
          table(prothio_dat$ESSAI,prothio_dat$Modalite)[,c(2,1)]*100,
        beside=TRUE,ylim=c(0,100))


###############################################################################
#Analysis of the effect of Prothioconazole on Prochloraze resistance
###############################################################################

#restricted dataset for Prochloraze effect analysis
prothio_dat<-datapietin[datapietin$Espece!="Oa" & datapietin$Modalite %in% 
                          c("Temoin","Prothioconazole") & 
                          datapietin$ESSAI!="ARVALIS17",]
prothio_dat<-drop.levels(prothio_dat)

#the different models considered
modProPro1<-glm(Prochloraze~ESSAI,family=binomial,prothio_dat)
modProPro2<-glm(Prochloraze~Modalite,family=binomial,prothio_dat)
modProPro<-glm(Prochloraze~ESSAI+Modalite,family=binomial,prothio_dat)
anova(modProPro1,modProPro,test="Chisq")
anova(modProPro2,modProPro,test="Chisq")
drop1(modProPro,test="Chisq")
summary(modProPro)

#barplot by Traitement and Essai
barplot(table(prothio_dat$Prochloraze,prothio_dat$Modalite),beside=TRUE,
        col=c("green3","red2"),main="Par traitement")
barplot(table(prothio_dat$Prochloraze,prothio_dat$ESSAI),beside=TRUE,
        col=c("green3","red2"),main="Par essai")

#number of Prochloraze sensitive/resistant strains
table(prothio_dat$Prochloraze,prothio_dat$Modalite)

#barplot "temoin" against "treated"
barplot(table(prothio_dat[prothio_dat$Prochloraze==1,]$ESSAI,
              prothio_dat[prothio_dat$Prochloraze==1,]$Modalite)[,c(2,1)]/
          table(prothio_dat$ESSAI,prothio_dat$Modalite)[,c(2,1)]*100,
        beside=TRUE,ylim=c(0,100))


###############################################################################
#END
###############################################################################