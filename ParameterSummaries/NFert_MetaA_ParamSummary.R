#Open necessary libraries
library(plyr)
library(dplyr)
library(Matrix)
library(Hmisc)
library(ggplot2)
library(nlme)
library(multcomp)
library(MuMIn)


#Download data file. Users will need to update link to reflect their own file path. 
urlfile = "https://raw.githubusercontent.com/gill20a/LitterDecompMetaAnalysis/master/Gill_NDeposition_MetaAnalysis_ParameterSummary.csv"
data<-read.csv(url(urlfile))
data[data == "<NA>"] <- NA

#Remove one outlier
data[, "weibull.mrt"][data[, "weibull.mrt"] > 3000] <- NA

#Evaluate range of study length, N addition rate, final mass remaining (Results paragraph 1)
min(data$Years, na.rm=T)
max(data$Years, na.rm=T)
min(data$N_addition_kgN_ha_yr, na.rm=T)
max(data$N_addition_kgN_ha_yr, na.rm=T)

min(data$FinalMass, na.rm=T)
max(data$FinalMass, na.rm=T)
mean(data$FinalMass, na.rm=T)
median(data$FinalMass, na.rm=T)
fm1<-lme((FinalMass)~(Years),random=~1|Paper, data=data, na.action=na.omit)
summary(fm1)

########################Best Model Fits
#####################Stacked Barplot Model Dist'n by Treatment

##Proportion best model, 3 way comparison. ALG updated 2020/02/18. Results paragrpah 2
sum(data$Single_3)/length(data$Single)
sum(data$Double_3)/length(data$Double)
sum(data$Asymp_3)/length(data$Asymp)

#Proportion if studies <1 year. Final discussion paragraph
less1<-data$Years < 1
sum(less1)/length(data$Years)


#Chi square of three way model. Results paragraph 2
TbN<-table(data$Treatment, data$Best_Model_3)
TbN
chisq.test(TbN)

###Proportion best model, four way comparison. Results paragraph 2
sum(data$Single_4)/length(data$Single)
sum(data$Double_4)/length(data$Double)
sum(data$Asymp_4)/length(data$Asymp)
sum(data$Weibull_4)/length(data$Weibull)


########### Figure 1
tib2<-data %>% group_by(Treatment) %>% count(Best_Model_3)
tib2<-as.data.frame(tib2)
tib2$perc<-tib2$n/c(rep(sum(tib2$n[1:3]),3), rep(sum(tib2$n[4:6]),3))

Control<-c(tib2$perc[3],tib2$perc[2],tib2$perc[1])
Nitrogen<-c(tib2$perc[6],tib2$perc[5],tib2$perc[4])
ModelDistn<-cbind(Control,Nitrogen)
rownames(ModelDistn)<-c("Single",	"Double",	 "Asymp")

par(mfrow=c(1,2), oma=c(0,0,0,0),mar = c(5, 4, 1.4, 0.2), xpd=TRUE,mgp=c(2.4,0.8,0), cex.lab=1.25, cex.main=1.25, cex.axis=1.15)
barplot(ModelDistn, col=c("#76b7b2", "#b07aa1", "#edc948"), border=NA, xlab="Treatment", 
        ylab="Proportion Best Model")
par(mar=c(1,2.25,1,1))
plot(NULL)
legend("topleft", inset=c(-0.15,0),c("Single",	"Double",	 "Asymp"), col=c("#76b7b2", "#b07aa1", "#edc948"),
       pch=c(16,16,16), bty="n")



