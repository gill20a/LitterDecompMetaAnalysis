# Open necessary libraries
library(boot)
library(plyr)
library(dplyr)
library(Matrix)
library(Hmisc)
library(ggplot2)
library(nlme)
library(multcomp)
library(MuMIn)

#Download data file. Users will need to update link to reflect their own file path. 
urlfile = "https://raw.githubusercontent.com/gill20a/LitterDecompMetaAnalysis/master/ResponseRatios/NFert_MetaA_RRSummary.csv"
data<-read.csv(url(urlfile))
data[data == "<NA>"] <- NA
data$N_asymp.A<-data$N_asymp.A+0.000000000001
data$C_asymp.A<-data$C_asymp.A+0.000000000001

My_Theme = theme(
  axis.title.x = element_text(size = 14),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 14))

meanfun <- function(data, i){
  d <- data[i, ]
  return(mean(d, na.rm=T))   
}

#Calculate log response ratios
data$C_MRT_HL<-data$C_weibull.mrt/data$C_weibull.half.life
data$N_MRT_HL<-data$N_weibull.mrt/data$N_weibull.half.life
data$single.k.RR<-log(data$N_single.k/data$C_single.k)
data$double.k.RR<-log(data$N_double.k/data$C_double.k)
data$double.ks.RR<-log(data$N_double.ks/data$C_double.ks)
data$double.A.RR<-log(data$N_double.A/data$C_double.A)
data$asymp.k.RR<-log(data$N_asymp.k/data$C_asymp.k)
data$asymp.A.RR<-log(data$N_asymp.A/data$C_asymp.A)
data$weibull.alpha.RR<-log(data$N_weibull.alpha/data$C_weibull.alpha)
data$weibull.mrt.RR<-log(data$N_weibull.mrt/data$C_weibull.mrt)
data$weibull.half.life.RR<-log(data$N_weibull.half.life/data$C_weibull.half.life)
data$weibull.tenth.RR<-log(data$N_weibull.tenth/data$C_weibull.tenth)
data$weibull.quarter.RR<-log(data$N_weibull.quarter/data$C_weibull.quarter)
data$MRT_HL.RR<-log(data$N_MRT_HL/data$C_MRT_HL)

#Transform predictor variables to approximate normaltiy
data$Sqrt_N_addition_kgN_ha_yr <-sqrt(data$N_addition_kgN_ha_yr)
data$Abs_Latitude<-abs(data$Latitude)
data$Log_MAP<-log(data$MAP)
data$Log_Litter_Ca<-log(data$Litter_Ca)
data$Sqrt_Litter_N<-sqrt(data$Litter_N)

#Break out "best model" subsets
asymp <- subset(data, subset = (C_Asymp_3 == "TRUE" | N_Asymp_3 == "TRUE") )
double <- subset(data, subset = (C_Double_3=="TRUE"| N_Double_3 == "TRUE"))
single<-subset(data, subset = (data$C_Single_3 == "TRUE" | data$N_Single_3 == "TRUE") )
weibull<-subset(data, subset = (data$C_Weibull_4 == "TRUE" | data$N_Weibull_4 == "TRUE") )


####################################################
### N form table Table S6
model <- lme(single.k.RR~N_Form_Cat2,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
singlek_mean_all<-tapply(data$single.k.RR, data$N_Form_Cat2, mean, na.rm=T)
singlek_se_all<-tapply(data$single.k.RR, data$N_Form_Cat2, sd, na.rm=T)/tapply(data$single.k.RR, data$N_Form_Cat2, nnzero, na.counted=F)
singlek_n_all<-tapply(data$single.k.RR+0.0000000001, data$N_Form_Cat2, nnzero, na.counted=F)
singlek_p_all<-anova(model)[2,4]
#
model <- lme(single.k.RR~N_Form_Cat2,random=~1|Paper,data=single, na.action=na.omit)
anova(model)
singlek_mean_single<-tapply(single$single.k.RR, single$N_Form_Cat2, mean, na.rm=T)
singlek_se_single<-tapply(single$single.k.RR, single$N_Form_Cat2, sd, na.rm=T)/tapply(single$single.k.RR, single$N_Form_Cat2, nnzero, na.counted=F)
singlek_n_single<-tapply(single$single.k.RR+0.0000000001, single$N_Form_Cat2, nnzero, na.counted=F)
singlek_p_single<-anova(model)[2,4]
###
model <- lme(double.k.RR~N_Form_Cat2,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
doublek_mean_all<-tapply(data$double.k.RR, data$N_Form_Cat2, mean, na.rm=T)
doublek_se_all<-tapply(data$double.k.RR, data$N_Form_Cat2, sd, na.rm=T)/tapply(data$double.k.RR, data$N_Form_Cat2, nnzero, na.counted=F)
doublek_n_all<-tapply(data$double.k.RR+0.0000000001, data$N_Form_Cat2, nnzero, na.counted=F)
doublek_p_all<-anova(model)[2,4]
#
model <- lme(double.k.RR~N_Form_Cat2,random=~1|Paper,data=double, na.action=na.omit)
anova(model)
doublek_mean_double<-tapply(double$double.k.RR, double$N_Form_Cat2, mean, na.rm=T)
doublek_se_double<-tapply(double$double.k.RR, double$N_Form_Cat2, sd, na.rm=T)/tapply(double$double.k.RR, double$N_Form_Cat2, nnzero, na.counted=F)
doublek_n_double<-tapply(double$double.k.RR+0.0000000001, double$N_Form_Cat2, nnzero, na.counted=F)
doublek_p_double<-anova(model)[2,4]
###
model <- lme(double.ks.RR~N_Form_Cat2,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
doubleks_mean_all<-tapply(data$double.ks.RR, data$N_Form_Cat2, mean, na.rm=T)
doubleks_se_all<-tapply(data$double.ks.RR, data$N_Form_Cat2, sd, na.rm=T)/tapply(data$double.ks.RR, data$N_Form_Cat2, nnzero, na.counted=F)
doubleks_n_all<-tapply(data$double.ks.RR, data$N_Form_Cat2, sd, na.rm=T)
doubleks_n_all<-tapply(data$double.ks.RR+0.0000000001, data$N_Form_Cat2, nnzero, na.counted=F)
doubleks_p_all<-anova(model)[2,4]
#
model <- lme(double.ks.RR~N_Form_Cat2,random=~1|Paper,data=double, na.action=na.omit)
anova(model)
doubleks_mean_double<-tapply(double$double.ks.RR, double$N_Form_Cat2, mean, na.rm=T)
doubleks_se_double<-tapply(double$double.ks.RR, double$N_Form_Cat2, sd, na.rm=T)/tapply(double$double.ks.RR, double$N_Form_Cat2, nnzero, na.counted=F)
doubleks_n_double<-tapply(double$double.ks.RR+0.0000000001, double$N_Form_Cat2, nnzero, na.counted=F)
doubleks_p_double<-anova(model)[2,4]
###
model <- lme(double.A.RR~N_Form_Cat2,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
doubleA_mean_all<-tapply(data$double.A.RR, data$N_Form_Cat2, mean, na.rm=T)
doubleA_se_all<-tapply(data$double.A.RR, data$N_Form_Cat2, sd, na.rm=T)/tapply(data$double.A.RR, data$N_Form_Cat2, nnzero, na.counted=F)
doubleA_n_all<-tapply(data$double.A.RR+0.0000000001, data$N_Form_Cat2, nnzero, na.counted=F)
doubleA_p_all<-anova(model)[2,4]
#
model <- lme(double.A.RR~N_Form_Cat2,random=~1|Paper,data=double, na.action=na.omit)
anova(model)
doubleA_mean_double<-tapply(double$double.A.RR, double$N_Form_Cat2, mean, na.rm=T)
doubleA_se_double<-tapply(double$double.A.RR, double$N_Form_Cat2, sd, na.rm=T)/tapply(double$double.A.RR, double$N_Form_Cat2, nnzero, na.counted=F)
doubleA_n_double<-tapply(double$double.A.RR+0.0000000001, double$N_Form_Cat2, nnzero, na.counted=F)
doubleA_p_double<-anova(model)[2,4]
###
model <- lme(asymp.k.RR~N_Form_Cat2,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
asympk_mean_all<-tapply(data$asymp.k.RR, data$N_Form_Cat2, mean, na.rm=T)
asympk_se_all<-tapply(data$asymp.k.RR, data$N_Form_Cat2, sd, na.rm=T)/tapply(data$asymp.k.RR, data$N_Form_Cat2, nnzero, na.counted=F)
asympk_n_all<-tapply(data$asymp.k.RR+0.0000000001, data$N_Form_Cat2, nnzero, na.counted=F)
asympk_p_all<-anova(model)[2,4]
#
model <- lme(asymp.k.RR~N_Form_Cat2,random=~1|Paper,data=asymp, na.action=na.omit)
anova(model)
asympk_mean_asymp<-tapply(asymp$asymp.k.RR, asymp$N_Form_Cat2, mean, na.rm=T)
asympk_se_asymp<-tapply(asymp$asymp.k.RR, asymp$N_Form_Cat2, sd, na.rm=T)/tapply(asymp$asymp.k.RR, asymp$N_Form_Cat2, nnzero, na.counted=F)
asympk_n_asymp<-tapply(asymp$asymp.k.RR+0.0000000001, asymp$N_Form_Cat2,nnzero, na.counted=F)
asympk_p_asymp<-anova(model)[2,4]
###
model <- lme(asymp.A.RR~N_Form_Cat2,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
asympA_mean_all<-tapply(data$asymp.A.RR, data$N_Form_Cat2, mean, na.rm=T)
asympA_se_all<-tapply(data$asymp.A.RR, data$N_Form_Cat2, sd, na.rm=T)/tapply(data$asymp.A.RR, data$N_Form_Cat2, nnzero, na.counted=F)
asympA_n_all<-tapply((data$asymp.A.RR+0.0000000001), data$N_Form_Cat2,nnzero, na.counted=F)
asympA_p_all<-anova(model)[2,4]
#
model <- lme(asymp.A.RR~N_Form_Cat2,random=~1|Paper,data=asymp, na.action=na.omit)
anova(model)
asympA_mean_asymp<-tapply(asymp$asymp.A.RR, asymp$N_Form_Cat2, mean, na.rm=T)
asympA_se_asymp<-tapply(asymp$asymp.A.RR, asymp$N_Form_Cat2, sd, na.rm=T)/tapply(asymp$asymp.A.RR, asymp$N_Form_Cat2, nnzero, na.counted=F)
asympA_n_asymp<-tapply(asymp$asymp.A.RR+0.0000000001, asymp$N_Form_Cat2, nnzero, na.counted=F)
asympA_p_asymp<-anova(model)[2,4]
###
model <- lme(weibull.tenth.RR~N_Form_Cat2,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
weibulltenth_mean_all<-tapply(data$weibull.tenth.RR, data$N_Form_Cat2, mean, na.rm=T)
weibulltenth_se_all<-tapply(data$weibull.tenth.RR, data$N_Form_Cat2, sd, na.rm=T)/tapply(data$weibull.tenth.RR, data$N_Form_Cat2, nnzero, na.counted=F)
weibulltenth_n_all<-tapply(data$weibull.tenth.RR+0.0000000001, data$N_Form_Cat2, nnzero, na.counted=F)
weibulltenth_p_all<-anova(model)[2,4]
#
model <- lme(weibull.tenth.RR~N_Form_Cat2,random=~1|Paper,data=weibull, na.action=na.omit)
anova(model)
weibulltenth_mean_weibull<-tapply(weibull$weibull.tenth.RR, weibull$N_Form_Cat2, mean, na.rm=T)
weibulltenth_se_weibull<-tapply(weibull$weibull.tenth.RR, weibull$N_Form_Cat2, sd, na.rm=T)/tapply(weibull$weibull.tenth.RR, weibull$N_Form_Cat2, nnzero, na.counted=F)
weibulltenth_n_weibull<-tapply(weibull$weibull.tenth.RR+0.0000000001, weibull$N_Form_Cat2, nnzero, na.counted=F)
weibulltenth_p_weibull<-anova(model)[2,4]
###
model <- lme(weibull.quarter.RR~N_Form_Cat2,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
weibullquarter_mean_all<-tapply(data$weibull.quarter.RR, data$N_Form_Cat2, mean, na.rm=T)
weibullquarter_se_all<-tapply(data$weibull.quarter.RR, data$N_Form_Cat2, sd, na.rm=T)/tapply(data$weibull.quarter.RR, data$N_Form_Cat2, nnzero, na.counted=F)
weibullquarter_n_all<-tapply(data$weibull.quarter.RR+0.0000000001, data$N_Form_Cat2, nnzero, na.counted=F)
weibullquarter_p_all<-anova(model)[2,4]
#
model <- lme(weibull.quarter.RR~N_Form_Cat2,random=~1|Paper,data=weibull, na.action=na.omit)
anova(model)
weibullquarter_mean_weibull<-tapply(weibull$weibull.quarter.RR, weibull$N_Form_Cat2, mean, na.rm=T)
weibullquarter_se_weibull<-tapply(weibull$weibull.quarter.RR, weibull$N_Form_Cat2, sd, na.rm=T)/tapply(weibull$weibull.quarter.RR, weibull$N_Form_Cat2, nnzero, na.counted=F)
weibullquarter_n_weibull<-tapply(weibull$weibull.quarter.RR+0.0000000001, weibull$N_Form_Cat2, nnzero, na.counted=F)
weibullquarter_p_weibull<-anova(model)[2,4]
###
model <- lme(weibull.half.life.RR~N_Form_Cat2,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
weibullhalflife_mean_all<-tapply(data$weibull.half.life.RR, data$N_Form_Cat2, mean, na.rm=T)
weibullhalflife_se_all<-tapply(data$weibull.half.life.RR, data$N_Form_Cat2, sd, na.rm=T)/tapply(data$weibull.half.life.RR, data$N_Form_Cat2, nnzero, na.counted=F)
weibullhalflife_n_all<-tapply(data$weibull.half.life.RR+0.0000000001, data$N_Form_Cat2, nnzero, na.counted=F)
weibullhalflife_p_all<-anova(model)[2,4]
#
model <- lme(weibull.half.life.RR~N_Form_Cat2,random=~1|Paper,data=weibull, na.action=na.omit)
anova(model)
weibullhalflife_mean_weibull<-tapply(weibull$weibull.half.life.RR, weibull$N_Form_Cat2, mean, na.rm=T)
weibullhalflife_se_weibull<-tapply(weibull$weibull.half.life.RR, weibull$N_Form_Cat2, sd, na.rm=T)/tapply(weibull$weibull.half.life.RR, weibull$N_Form_Cat2, nnzero, na.counted=F)
weibullhalflife_n_weibull<-tapply(weibull$weibull.half.life.RR+0.0000000001, weibull$N_Form_Cat2, nnzero, na.counted=F)
weibullhalflife_p_weibull<-anova(model)[2,4]
###
model <- lme(weibull.mrt.RR~N_Form_Cat2,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
weibullmrt_mean_all<-tapply(data$weibull.mrt.RR, data$N_Form_Cat2, mean, na.rm=T)
weibullmrt_se_all<-tapply(data$weibull.mrt.RR, data$N_Form_Cat2, sd, na.rm=T)/tapply(data$weibull.mrt.RR, data$N_Form_Cat2, nnzero, na.counted=F)
weibullmrt_n_all<-tapply(data$weibull.mrt.RR+0.0000000001, data$N_Form_Cat2, nnzero, na.counted=F)
weibullmrt_p_all<-anova(model)[2,4]
#
model <- lme(weibull.mrt.RR~N_Form_Cat2,random=~1|Paper,data=weibull, na.action=na.omit)
anova(model)
weibullmrt_mean_weibull<-tapply(weibull$weibull.half.life.RR, weibull$N_Form_Cat2, mean, na.rm=T)
weibullmrt_se_weibull<-tapply(weibull$weibull.half.life.RR, weibull$N_Form_Cat2, sd, na.rm=T)/tapply(weibull$weibull.mrt.RR, weibull$N_Form_Cat2, nnzero, na.counted=F)
weibullmrt_n_weibull<-tapply(weibull$weibull.half.life.RR+0.0000000001, weibull$N_Form_Cat2, nnzero, na.counted=F)
weibullmrt_p_weibull<-anova(model)[2,4]
###
model <- lme(weibull.alpha.RR~N_Form_Cat2,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
weibullalpha_mean_all<-tapply(data$weibull.alpha.RR, data$N_Form_Cat2, mean, na.rm=T)
weibullalpha_se_all<-tapply(data$weibull.alpha.RR, data$N_Form_Cat2, sd, na.rm=T)/tapply(data$weibull.alpha.RR, data$N_Form_Cat2, nnzero, na.counted=F)
weibullalpha_n_all<-tapply(data$weibull.alpha.RR+0.0000000001, data$N_Form_Cat2, nnzero, na.counted=F)
weibullalpha_p_all<-anova(model)[2,4]
#
model <- lme(weibull.alpha.RR~N_Form_Cat2,random=~1|Paper,data=weibull, na.action=na.omit)
anova(model)
weibullalpha_mean_weibull<-tapply(weibull$weibull.alpha.RR, weibull$N_Form_Cat2, mean, na.rm=T)
weibullalpha_se_weibull<-tapply(weibull$weibull.alpha.RR, weibull$N_Form_Cat2, sd, na.rm=T)/tapply(weibull$weibull.alpha.RR, weibull$N_Form_Cat2, nnzero, na.counted=F)
weibullalpha_n_weibull<-tapply(weibull$weibull.alpha.RR+0.0000000001, weibull$N_Form_Cat2,nnzero, na.counted=F)
weibullalpha_p_weibull<-anova(model)[2,4]

###
model <- lme(MRT_HL.RR~N_Form_Cat2,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
MRT_HL_mean_all<-tapply(data$MRT_HL.RR, data$N_Form_Cat2, mean, na.rm=T)
MRT_HL_se_all<-tapply(data$MRT_HL.RR, data$N_Form_Cat2, sd, na.rm=T)/tapply(data$MRT_HL.RR, data$N_Form_Cat2, nnzero, na.counted=F)
MRT_HL_n_all<-tapply(data$MRT_HL.RR+0.0000000001, data$N_Form_Cat2, nnzero, na.counted=F)
MRT_HL_p_all<-anova(model)[2,4]
#
model <- lme(MRT_HL.RR~N_Form_Cat2,random=~1|Paper,data=weibull, na.action=na.omit)
anova(model)
MRT_HL_mean_weibull<-tapply(weibull$MRT_HL.RR, weibull$N_Form_Cat2, mean, na.rm=T)
MRT_HL_se_weibull<-tapply(weibull$MRT_HL.RR, weibull$N_Form_Cat2, sd, na.rm=T)/tapply(weibull$MRT_HL.RR, weibull$N_Form_Cat2, nnzero, na.counted=F)
MRT_HL_n_weibull<-tapply(weibull$MRT_HL.RR+0.0000000001, weibull$N_Form_Cat2, nnzero, na.counted=F)
MRT_HL_p_weibull<-anova(model)[2,4]

means<-rbind(
  singlek_mean_all,
  singlek_mean_single,
  doublek_mean_all,
  doublek_mean_double,
  doubleks_mean_all,
  doubleks_mean_double,
  doubleA_mean_all,
  doubleA_mean_double,
  asympk_mean_all,
  asympk_mean_asymp,
  asympA_mean_all,
  asympA_mean_asymp,
  weibulltenth_mean_all,
  weibulltenth_mean_weibull,
  weibullquarter_mean_all,
  weibullquarter_mean_weibull,
  weibullhalflife_mean_all,
  weibullhalflife_mean_weibull,
  weibullmrt_mean_all,
  weibullmrt_mean_weibull,
  weibullalpha_mean_all,
  weibullalpha_mean_weibull,
  MRT_HL_mean_all,
  MRT_HL_mean_weibull
)

ses<-rbind(
  singlek_se_all,
  singlek_se_single,
  doublek_se_all,
  doublek_se_double,
  doubleks_se_all,
  doubleks_se_double,
  doubleA_se_all,
  doubleA_se_double,
  asympk_se_all,
  asympk_se_asymp,
  asympA_se_all,
  asympA_se_asymp,
  weibulltenth_se_all,
  weibulltenth_se_weibull,
  weibullquarter_se_all,
  weibullquarter_se_weibull,
  weibullhalflife_se_all,
  weibullhalflife_se_weibull,
  weibullmrt_se_all,
  weibullmrt_se_weibull,
  weibullalpha_se_all,
  weibullalpha_se_weibull,
  MRT_HL_se_all,
  MRT_HL_se_weibull
)

ps<-rbind(
  singlek_p_all,
  singlek_p_single,
  doublek_p_all,
  doublek_p_double,
  doubleks_p_all,
  doubleks_p_double,
  doubleA_p_all,
  doubleA_p_double,
  asympk_p_all,
  asympk_p_asymp,
  asympA_p_all,
  asympA_p_asymp,
  weibulltenth_p_all,
  weibulltenth_p_weibull,
  weibullquarter_p_all,
  weibullquarter_p_weibull,
  weibullhalflife_p_all,
  weibullhalflife_p_weibull,
  weibullmrt_p_all,
  weibullmrt_p_weibull,
  weibullalpha_p_all,
  weibullalpha_p_weibull,
  MRT_HL_p_all,
  MRT_HL_p_weibull
)

ns<-rbind(
  singlek_n_all,
  singlek_n_single,
  doublek_n_all,
  doublek_n_double,
  doubleks_n_all,
  doubleks_n_double,
  doubleA_n_all,
  doubleA_n_double,
  asympk_n_all,
  asympk_n_asymp,
  asympA_n_all,
  asympA_n_asymp,
  weibulltenth_n_all,
  weibulltenth_n_weibull,
  weibullquarter_n_all,
  weibullquarter_n_weibull,
  weibullhalflife_n_all,
  weibullhalflife_n_weibull,
  weibullmrt_n_all,
  weibullmrt_n_weibull,
  weibullalpha_n_all,
  weibullalpha_n_weibull,
  MRT_HL_n_all,
  MRT_HL_n_weibull
)

######################################################
###Data for substrate table S7
model <- lme(single.k.RR~Substrate_Cat,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
singlek_mean_all<-tapply(data$single.k.RR, data$Substrate_Cat, mean, na.rm=T)
singlek_se_all<-tapply(data$single.k.RR, data$Substrate_Cat, sd, na.rm=T)/tapply(data$single.k.RR, data$Substrate_Cat, nnzero, na.counted=F)
singlek_n_all<-tapply(data$single.k.RR+0.0000000001, data$Substrate_Cat, nnzero, na.counted=F)
singlek_p_all<-anova(model)[2,4]
#
model <- lme(single.k.RR~Substrate_Cat,random=~1|Paper,data=single, na.action=na.omit)
anova(model)
singlek_mean_single<-tapply(single$single.k.RR, single$Substrate_Cat, mean, na.rm=T)
singlek_se_single<-tapply(single$single.k.RR, single$Substrate_Cat, sd, na.rm=T)/tapply(single$single.k.RR, single$Substrate_Cat, nnzero, na.counted=F)
singlek_n_single<-tapply(single$single.k.RR+0.0000000001, single$Substrate_Cat, nnzero, na.counted=F)
singlek_p_single<-anova(model)[2,4]
###
model <- lme(double.k.RR~Substrate_Cat,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
doublek_mean_all<-tapply(data$double.k.RR, data$Substrate_Cat, mean, na.rm=T)
doublek_se_all<-tapply(data$double.k.RR, data$Substrate_Cat, sd, na.rm=T)/tapply(data$double.k.RR, data$Substrate_Cat, nnzero, na.counted=F)
doublek_n_all<-tapply(data$double.k.RR+0.0000000001, data$Substrate_Cat, nnzero, na.counted=F)
doublek_p_all<-anova(model)[2,4]
#
model <- lme(double.k.RR~Substrate_Cat,random=~1|Paper,data=double, na.action=na.omit)
anova(model)
doublek_mean_double<-tapply(double$double.k.RR, double$Substrate_Cat, mean, na.rm=T)
doublek_se_double<-tapply(double$double.k.RR, double$Substrate_Cat, sd, na.rm=T)/tapply(double$double.k.RR, double$Substrate_Cat, nnzero, na.counted=F)
doublek_n_double<-tapply(double$double.k.RR+0.0000000001, double$Substrate_Cat, nnzero, na.counted=F)
doublek_p_double<-anova(model)[2,4]
###
model <- lme(double.ks.RR~Substrate_Cat,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
doubleks_mean_all<-tapply(data$double.ks.RR, data$Substrate_Cat, mean, na.rm=T)
doubleks_se_all<-tapply(data$double.ks.RR, data$Substrate_Cat, sd, na.rm=T)/tapply(data$double.ks.RR, data$Substrate_Cat, nnzero, na.counted=F)
doubleks_n_all<-tapply(data$double.ks.RR, data$Substrate_Cat, sd, na.rm=T)
doubleks_n_all<-tapply(data$double.ks.RR+0.0000000001, data$Substrate_Cat, nnzero, na.counted=F)
doubleks_p_all<-anova(model)[2,4]
#
model <- lme(double.ks.RR~Substrate_Cat,random=~1|Paper,data=double, na.action=na.omit)
anova(model)
doubleks_mean_double<-tapply(double$double.ks.RR, double$Substrate_Cat, mean, na.rm=T)
doubleks_se_double<-tapply(double$double.ks.RR, double$Substrate_Cat, sd, na.rm=T)/tapply(double$double.ks.RR, double$Substrate_Cat, nnzero, na.counted=F)
doubleks_n_double<-tapply(double$double.ks.RR+0.0000000001, double$Substrate_Cat, nnzero, na.counted=F)
doubleks_p_double<-anova(model)[2,4]
###
model <- lme(double.A.RR~Substrate_Cat,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
doubleA_mean_all<-tapply(data$double.A.RR, data$Substrate_Cat, mean, na.rm=T)
doubleA_se_all<-tapply(data$double.A.RR, data$Substrate_Cat, sd, na.rm=T)/tapply(data$double.A.RR, data$Substrate_Cat, nnzero, na.counted=F)
doubleA_n_all<-tapply(data$double.A.RR+0.0000000001, data$Substrate_Cat, nnzero, na.counted=F)
doubleA_p_all<-anova(model)[2,4]
#
model <- lme(double.A.RR~Substrate_Cat,random=~1|Paper,data=double, na.action=na.omit)
anova(model)
doubleA_mean_double<-tapply(double$double.A.RR, double$Substrate_Cat, mean, na.rm=T)
doubleA_se_double<-tapply(double$double.A.RR, double$Substrate_Cat, sd, na.rm=T)/tapply(double$double.A.RR, double$Substrate_Cat, nnzero, na.counted=F)
doubleA_n_double<-tapply(double$double.A.RR+0.0000000001, double$Substrate_Cat, nnzero, na.counted=F)
doubleA_p_double<-anova(model)[2,4]
###
model <- lme(asymp.k.RR~Substrate_Cat,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
asympk_mean_all<-tapply(data$asymp.k.RR, data$Substrate_Cat, mean, na.rm=T)
asympk_se_all<-tapply(data$asymp.k.RR, data$Substrate_Cat, sd, na.rm=T)/tapply(data$asymp.k.RR, data$Substrate_Cat, nnzero, na.counted=F)
asympk_n_all<-tapply(data$asymp.k.RR+0.0000000001, data$Substrate_Cat, nnzero, na.counted=F)
asympk_p_all<-anova(model)[2,4]
#
model <- lme(asymp.k.RR~Substrate_Cat,random=~1|Paper,data=asymp, na.action=na.omit)
anova(model)
asympk_mean_asymp<-tapply(asymp$asymp.k.RR, asymp$Substrate_Cat, mean, na.rm=T)
asympk_se_asymp<-tapply(asymp$asymp.k.RR, asymp$Substrate_Cat, sd, na.rm=T)/tapply(asymp$asymp.k.RR, asymp$Substrate_Cat, nnzero, na.counted=F)
asympk_n_asymp<-tapply(asymp$asymp.k.RR+0.0000000001, asymp$Substrate_Cat,nnzero, na.counted=F)
asympk_p_asymp<-anova(model)[2,4]
###
model <- lme(asymp.A.RR~Substrate_Cat,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
asympA_mean_all<-tapply(data$asymp.A.RR, data$Substrate_Cat, mean, na.rm=T)
asympA_se_all<-tapply(data$asymp.A.RR, data$Substrate_Cat, sd, na.rm=T)/tapply(data$asymp.A.RR, data$Substrate_Cat, nnzero, na.counted=F)
asympA_n_all<-tapply((data$asymp.A.RR+0.0000000001), data$Substrate_Cat,nnzero, na.counted=F)
asympA_p_all<-anova(model)[2,4]
#
model <- lme(asymp.A.RR~Substrate_Cat,random=~1|Paper,data=asymp, na.action=na.omit)
anova(model)
asympA_mean_asymp<-tapply(asymp$asymp.A.RR, asymp$Substrate_Cat, mean, na.rm=T)
asympA_se_asymp<-tapply(asymp$asymp.A.RR, asymp$Substrate_Cat, sd, na.rm=T)/tapply(asymp$asymp.A.RR, asymp$Substrate_Cat, nnzero, na.counted=F)
asympA_n_asymp<-tapply(asymp$asymp.A.RR+0.0000000001, asymp$Substrate_Cat, nnzero, na.counted=F)
asympA_p_asymp<-anova(model)[2,4]
###
model <- lme(weibull.tenth.RR~Substrate_Cat,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
weibulltenth_mean_all<-tapply(data$weibull.tenth.RR, data$Substrate_Cat, mean, na.rm=T)
weibulltenth_se_all<-tapply(data$weibull.tenth.RR, data$Substrate_Cat, sd, na.rm=T)/tapply(data$weibull.tenth.RR, data$Substrate_Cat, nnzero, na.counted=F)
weibulltenth_n_all<-tapply(data$weibull.tenth.RR+0.0000000001, data$Substrate_Cat, nnzero, na.counted=F)
weibulltenth_p_all<-anova(model)[2,4]
#
model <- lme(weibull.tenth.RR~Substrate_Cat,random=~1|Paper,data=weibull, na.action=na.omit)
anova(model)
weibulltenth_mean_weibull<-tapply(weibull$weibull.tenth.RR, weibull$Substrate_Cat, mean, na.rm=T)
weibulltenth_se_weibull<-tapply(weibull$weibull.tenth.RR, weibull$Substrate_Cat, sd, na.rm=T)/tapply(weibull$weibull.tenth.RR, weibull$Substrate_Cat, nnzero, na.counted=F)
weibulltenth_n_weibull<-tapply(weibull$weibull.tenth.RR+0.0000000001, weibull$Substrate_Cat, nnzero, na.counted=F)
weibulltenth_p_weibull<-anova(model)[2,4]
###
model <- lme(weibull.quarter.RR~Substrate_Cat,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
weibullquarter_mean_all<-tapply(data$weibull.quarter.RR, data$Substrate_Cat, mean, na.rm=T)
weibullquarter_se_all<-tapply(data$weibull.quarter.RR, data$Substrate_Cat, sd, na.rm=T)/tapply(data$weibull.quarter.RR, data$Substrate_Cat, nnzero, na.counted=F)
weibullquarter_n_all<-tapply(data$weibull.quarter.RR+0.0000000001, data$Substrate_Cat, nnzero, na.counted=F)
weibullquarter_p_all<-anova(model)[2,4]
#
model <- lme(weibull.quarter.RR~Substrate_Cat,random=~1|Paper,data=weibull, na.action=na.omit)
anova(model)
weibullquarter_mean_weibull<-tapply(weibull$weibull.quarter.RR, weibull$Substrate_Cat, mean, na.rm=T)
weibullquarter_se_weibull<-tapply(weibull$weibull.quarter.RR, weibull$Substrate_Cat, sd, na.rm=T)/tapply(weibull$weibull.quarter.RR, weibull$Substrate_Cat, nnzero, na.counted=F)
weibullquarter_n_weibull<-tapply(weibull$weibull.quarter.RR+0.0000000001, weibull$Substrate_Cat, nnzero, na.counted=F)
weibullquarter_p_weibull<-anova(model)[2,4]
###
model <- lme(weibull.half.life.RR~Substrate_Cat,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
weibullhalflife_mean_all<-tapply(data$weibull.half.life.RR, data$Substrate_Cat, mean, na.rm=T)
weibullhalflife_se_all<-tapply(data$weibull.half.life.RR, data$Substrate_Cat, sd, na.rm=T)/tapply(data$weibull.half.life.RR, data$Substrate_Cat, nnzero, na.counted=F)
weibullhalflife_n_all<-tapply(data$weibull.half.life.RR+0.0000000001, data$Substrate_Cat, nnzero, na.counted=F)
weibullhalflife_p_all<-anova(model)[2,4]
#
model <- lme(weibull.half.life.RR~Substrate_Cat,random=~1|Paper,data=weibull, na.action=na.omit)
anova(model)
weibullhalflife_mean_weibull<-tapply(weibull$weibull.half.life.RR, weibull$Substrate_Cat, mean, na.rm=T)
weibullhalflife_se_weibull<-tapply(weibull$weibull.half.life.RR, weibull$Substrate_Cat, sd, na.rm=T)/tapply(weibull$weibull.half.life.RR, weibull$Substrate_Cat, nnzero, na.counted=F)
weibullhalflife_n_weibull<-tapply(weibull$weibull.half.life.RR+0.0000000001, weibull$Substrate_Cat, nnzero, na.counted=F)
weibullhalflife_p_weibull<-anova(model)[2,4]
###
model <- lme(weibull.mrt.RR~Substrate_Cat,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
weibullmrt_mean_all<-tapply(data$weibull.mrt.RR, data$Substrate_Cat, mean, na.rm=T)
weibullmrt_se_all<-tapply(data$weibull.mrt.RR, data$Substrate_Cat, sd, na.rm=T)/tapply(data$weibull.mrt.RR, data$Substrate_Cat, nnzero, na.counted=F)
weibullmrt_n_all<-tapply(data$weibull.mrt.RR+0.0000000001, data$Substrate_Cat, nnzero, na.counted=F)
weibullmrt_p_all<-anova(model)[2,4]
#
model <- lme(weibull.mrt.RR~Substrate_Cat,random=~1|Paper,data=weibull, na.action=na.omit)
anova(model)
weibullmrt_mean_weibull<-tapply(weibull$weibull.half.life.RR, weibull$Substrate_Cat, mean, na.rm=T)
weibullmrt_se_weibull<-tapply(weibull$weibull.half.life.RR, weibull$Substrate_Cat, sd, na.rm=T)/tapply(weibull$weibull.mrt.RR, weibull$Substrate_Cat, nnzero, na.counted=F)
weibullmrt_n_weibull<-tapply(weibull$weibull.half.life.RR+0.0000000001, weibull$Substrate_Cat, nnzero, na.counted=F)
weibullmrt_p_weibull<-anova(model)[2,4]
###
model <- lme(weibull.alpha.RR~Substrate_Cat,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
weibullalpha_mean_all<-tapply(data$weibull.alpha.RR, data$Substrate_Cat, mean, na.rm=T)
weibullalpha_se_all<-tapply(data$weibull.alpha.RR, data$Substrate_Cat, sd, na.rm=T)/tapply(data$weibull.alpha.RR, data$Substrate_Cat, nnzero, na.counted=F)
weibullalpha_n_all<-tapply(data$weibull.alpha.RR+0.0000000001, data$Substrate_Cat, nnzero, na.counted=F)
weibullalpha_p_all<-anova(model)[2,4]
#
model <- lme(weibull.alpha.RR~Substrate_Cat,random=~1|Paper,data=weibull, na.action=na.omit)
anova(model)
weibullalpha_mean_weibull<-tapply(weibull$weibull.alpha.RR, weibull$Substrate_Cat, mean, na.rm=T)
weibullalpha_se_weibull<-tapply(weibull$weibull.alpha.RR, weibull$Substrate_Cat, sd, na.rm=T)/tapply(weibull$weibull.alpha.RR, weibull$Substrate_Cat, nnzero, na.counted=F)
weibullalpha_n_weibull<-tapply(weibull$weibull.alpha.RR+0.0000000001, weibull$Substrate_Cat,nnzero, na.counted=F)
weibullalpha_p_weibull<-anova(model)[2,4]

###
model <- lme(MRT_HL.RR~Substrate_Cat,random=~1|Paper,data=data, na.action=na.omit)
anova(model)
MRT_HL_mean_all<-tapply(data$MRT_HL.RR, data$Substrate_Cat, mean, na.rm=T)
MRT_HL_se_all<-tapply(data$MRT_HL.RR, data$Substrate_Cat, sd, na.rm=T)/tapply(data$MRT_HL.RR, data$Substrate_Cat, nnzero, na.counted=F)
MRT_HL_n_all<-tapply(data$MRT_HL.RR+0.0000000001, data$Substrate_Cat, nnzero, na.counted=F)
MRT_HL_p_all<-anova(model)[2,4]
#
model <- lme(MRT_HL.RR~Substrate_Cat,random=~1|Paper,data=weibull, na.action=na.omit)
anova(model)
MRT_HL_mean_weibull<-tapply(weibull$MRT_HL.RR, weibull$Substrate_Cat, mean, na.rm=T)
MRT_HL_se_weibull<-tapply(weibull$MRT_HL.RR, weibull$Substrate_Cat, sd, na.rm=T)/tapply(weibull$MRT_HL.RR, weibull$Substrate_Cat, nnzero, na.counted=F)
MRT_HL_n_weibull<-tapply(weibull$MRT_HL.RR+0.0000000001, weibull$Substrate_Cat, nnzero, na.counted=F)
MRT_HL_p_weibull<-anova(model)[2,4]

means<-rbind(
  singlek_mean_all,
  singlek_mean_single,
  doublek_mean_all,
  doublek_mean_double,
  doubleks_mean_all,
  doubleks_mean_double,
  doubleA_mean_all,
  doubleA_mean_double,
  asympk_mean_all,
  asympk_mean_asymp,
  asympA_mean_all,
  asympA_mean_asymp,
  weibulltenth_mean_all,
  weibulltenth_mean_weibull,
  weibullquarter_mean_all,
  weibullquarter_mean_weibull,
  weibullhalflife_mean_all,
  weibullhalflife_mean_weibull,
  weibullmrt_mean_all,
  weibullmrt_mean_weibull,
  weibullalpha_mean_all,
  weibullalpha_mean_weibull,
  MRT_HL_mean_all,
  MRT_HL_mean_weibull
)

ses<-rbind(
  singlek_se_all,
  singlek_se_single,
  doublek_se_all,
  doublek_se_double,
  doubleks_se_all,
  doubleks_se_double,
  doubleA_se_all,
  doubleA_se_double,
  asympk_se_all,
  asympk_se_asymp,
  asympA_se_all,
  asympA_se_asymp,
  weibulltenth_se_all,
  weibulltenth_se_weibull,
  weibullquarter_se_all,
  weibullquarter_se_weibull,
  weibullhalflife_se_all,
  weibullhalflife_se_weibull,
  weibullmrt_se_all,
  weibullmrt_se_weibull,
  weibullalpha_se_all,
  weibullalpha_se_weibull,
  MRT_HL_se_all,
  MRT_HL_se_weibull
)

ps<-rbind(
  singlek_p_all,
  singlek_p_single,
  doublek_p_all,
  doublek_p_double,
  doubleks_p_all,
  doubleks_p_double,
  doubleA_p_all,
  doubleA_p_double,
  asympk_p_all,
  asympk_p_asymp,
  asympA_p_all,
  asympA_p_asymp,
  weibulltenth_p_all,
  weibulltenth_p_weibull,
  weibullquarter_p_all,
  weibullquarter_p_weibull,
  weibullhalflife_p_all,
  weibullhalflife_p_weibull,
  weibullmrt_p_all,
  weibullmrt_p_weibull,
  weibullalpha_p_all,
  weibullalpha_p_weibull,
  MRT_HL_p_all,
  MRT_HL_p_weibull
)

ns<-rbind(
  singlek_n_all,
  singlek_n_single,
  doublek_n_all,
  doublek_n_double,
  doubleks_n_all,
  doubleks_n_double,
  doubleA_n_all,
  doubleA_n_double,
  asympk_n_all,
  asympk_n_asymp,
  asympA_n_all,
  asympA_n_asymp,
  weibulltenth_n_all,
  weibulltenth_n_weibull,
  weibullquarter_n_all,
  weibullquarter_n_weibull,
  weibullhalflife_n_all,
  weibullhalflife_n_weibull,
  weibullmrt_n_all,
  weibullmrt_n_weibull,
  weibullalpha_n_all,
  weibullalpha_n_weibull,
  MRT_HL_n_all,
  MRT_HL_n_weibull
)


#######################################################################
####################### Means and confidence intervals for figure 2
data_single_k_mean<-mean(data$single.k.RR, na.rm=T)
boot.single.k.RR <- boot(as.data.frame(data$single.k.RR), statistic=meanfun, R=5000)
ci.single.k.RR <-boot.ci(boot.single.k.RR, conf=0.95)
norm.ci.single.k.RR<-ci.single.k.RR[4]

data_double_k_mean<-mean(data$double.k.RR, na.rm=T)
boot.double.k.RR <- boot(as.data.frame(data$double.k.RR), statistic=meanfun, R=5000)
ci.double.k.RR <-boot.ci(boot.double.k.RR, conf=0.95)
norm.ci.double.k.RR<-ci.double.k.RR[4]

data_double_ks_mean<-mean(data$double.ks.RR, na.rm=T)
boot.double.ks.RR <- boot(as.data.frame(data$double.ks.RR), statistic=meanfun, R=5000)
ci.double.ks.RR <-boot.ci(boot.double.ks.RR, conf=0.95)
norm.ci.double.ks.RR<-ci.double.ks.RR[4]

data_double_A_mean<-mean(data$double.A.RR, na.rm=T)
boot.double.A.RR <- boot(as.data.frame(data$double.A.RR), statistic=meanfun, R=5000)
ci.double.A.RR <-boot.ci(boot.double.A.RR, conf=0.95)
norm.ci.double.A.RR<-ci.double.A.RR[4]

data_asymp_k_mean<-mean(data$asymp.k.RR, na.rm=T)
boot.asymp.k.RR <- boot(as.data.frame(data$asymp.k.RR), statistic=meanfun, R=5000)
ci.asymp.k.RR <-boot.ci(boot.asymp.k.RR, conf=0.95)
norm.ci.asymp.k.RR<-ci.asymp.k.RR[4]

data_asymp_A_mean<-mean(data$asymp.A.RR, na.rm=T)
boot.asymp.A.RR <- boot(as.data.frame(data$asymp.A.RR), statistic=meanfun, R=5000)
ci.asymp.A.RR <-boot.ci(boot.asymp.A.RR, conf=0.95)
norm.ci.asymp.A.RR<-ci.asymp.A.RR[4]

data_tenth_mean<-mean(data$weibull.tenth.RR, na.rm=T)
boot.weibull.tenth.RR <- boot(as.data.frame(data$weibull.tenth.RR), statistic=meanfun, R=5000)
ci.weibull.tenth.RR <-boot.ci(boot.weibull.tenth.RR, conf=0.95)
norm.ci.weibull.tenth.RR<-ci.weibull.tenth.RR[4]

data_quarter_mean<-mean(data$weibull.quarter.RR, na.rm=T)
boot.weibull.quarter.RR <- boot(as.data.frame(data$weibull.quarter.RR), statistic=meanfun, R=5000)
ci.weibull.quarter.RR <-boot.ci(boot.weibull.quarter.RR, conf=0.95)
norm.ci.weibull.quarter.RR<-ci.weibull.quarter.RR[4]

data_half_life_mean<-mean(data$weibull.half.life.RR, na.rm=T)
boot.weibull.half.life.RR <- boot(as.data.frame(data$weibull.half.life.RR), statistic=meanfun, R=5000)
ci.weibull.half.life.RR <-boot.ci(boot.weibull.half.life.RR, conf=0.95)
norm.ci.weibull.half.life.RR<-ci.weibull.half.life.RR[4]

data_mrt_mean<-mean(data$weibull.mrt.RR, na.rm=T)
boot.weibull.mrt.RR <- boot(as.data.frame(data$weibull.mrt.RR), statistic=meanfun, R=5000)
ci.weibull.mrt.RR <-boot.ci(boot.weibull.mrt.RR, conf=0.95)
norm.ci.weibull.mrt.RR<-ci.weibull.mrt.RR[4]

data_alpha_mean<-mean(data$weibull.alpha.RR, na.rm=T)
boot.weibull.alpha.RR <- boot(as.data.frame(data$weibull.alpha.RR), statistic=meanfun, R=5000)
ci.weibull.alpha.RR <-boot.ci(boot.weibull.alpha.RR, conf=0.95)
norm.ci.weibull.alpha.RR<-ci.weibull.alpha.RR[4]

data_MRT_HL_mean<-mean(data$MRT_HL.RR, na.rm=T)
boot.MRT_HL.RR <- boot(as.data.frame(data$MRT_HL.RR), statistic=meanfun, R=5000)
ci.MRT_HL.RR <-boot.ci(boot.MRT_HL.RR, conf=0.95)
norm.ci.MRT_HL.RR<-ci.MRT_HL.RR[4]

####################################################################################
single_single_k_mean<-mean(single$single.k.RR, na.rm=T)
boot.single.k.RR <- boot(as.data.frame(single$single.k.RR), statistic=meanfun, R=5000)
ci.single.k.RR <-boot.ci(boot.single.k.RR, conf=0.95)
specific.norm.ci.single.k.RR<-ci.single.k.RR[4]

double_double_k_mean<-mean(double$double.k.RR, na.rm=T)
boot.double.k.RR <- boot(as.data.frame(double$double.k.RR), statistic=meanfun, R=5000)
ci.double.k.RR <-boot.ci(boot.double.k.RR, conf=0.95)
specific.norm.ci.double.k.RR<-ci.double.k.RR[4]

double_double_ks_mean<-mean(double$double.ks.RR, na.rm=T)
boot.double.ks.RR <- boot(as.data.frame(double$double.ks.RR), statistic=meanfun, R=5000)
ci.double.ks.RR <-boot.ci(boot.double.ks.RR, conf=0.95)
specific.norm.ci.double.ks.RR<-ci.double.ks.RR[4]

double_double_A_mean<-mean(double$double.A.RR, na.rm=T)
boot.double.A.RR <- boot(as.data.frame(double$double.A.RR), statistic=meanfun, R=5000)
ci.double.A.RR <-boot.ci(boot.double.A.RR, conf=0.95)
specific.norm.ci.double.A.RR<-ci.double.A.RR[4]

asymp_asymp_k_mean<-mean(asymp$asymp.k.RR, na.rm=T)
boot.asymp.k.RR <- boot(as.data.frame(asymp$asymp.k.RR), statistic=meanfun, R=5000)
ci.asymp.k.RR <-boot.ci(boot.asymp.k.RR, conf=0.95)
specific.norm.ci.asymp.k.RR<-ci.asymp.k.RR[4]

asymp_asymp_A_mean<-mean(asymp$asymp.A.RR, na.rm=T)
boot.asymp.A.RR <- boot(as.data.frame(asymp$asymp.A.RR), statistic=meanfun, R=5000)
ci.asymp.A.RR <-boot.ci(boot.asymp.A.RR, conf=0.95)
specific.norm.ci.asymp.A.RR<-ci.asymp.A.RR[4]

weibull_tenth_mean<-mean(weibull$weibull.tenth.RR, na.rm=T)
boot.weibull.tenth.RR <- boot(as.data.frame(weibull$weibull.tenth.RR), statistic=meanfun, R=5000)
ci.weibull.tenth.RR <-boot.ci(boot.weibull.tenth.RR, conf=0.95)
specific.norm.ci.weibull.tenth.RR<-ci.weibull.tenth.RR[4]

weibull_quarter_mean<-mean(weibull$weibull.quarter.RR, na.rm=T)
boot.weibull.quarter.RR <- boot(as.data.frame(weibull$weibull.quarter.RR), statistic=meanfun, R=5000)
ci.weibull.quarter.RR <-boot.ci(boot.weibull.quarter.RR, conf=0.95)
specific.norm.ci.weibull.quarter.RR<-ci.weibull.quarter.RR[4]

weibull_half_life_mean<-mean(weibull$weibull.half.life.RR, na.rm=T)
boot.weibull.half.life.RR <- boot(as.data.frame(weibull$weibull.half.life.RR), statistic=meanfun, R=5000)
ci.weibull.half.life.RR <-boot.ci(boot.weibull.half.life.RR, conf=0.95)
specific.norm.ci.weibull.half.life.RR<-ci.weibull.half.life.RR[4]

weibull_mrt_mean<-mean(weibull$weibull.mrt.RR, na.rm=T)
boot.weibull.mrt.RR <- boot(as.data.frame(weibull$weibull.mrt.RR), statistic=meanfun, R=5000)
ci.weibull.mrt.RR <-boot.ci(boot.weibull.mrt.RR, conf=0.95)
specific.norm.ci.weibull.mrt.RR<-ci.weibull.mrt.RR[4]

weibull_alpha_mean<-mean(weibull$weibull.alpha.RR, na.rm=T)
boot.weibull.alpha.RR <- boot(as.data.frame(weibull$weibull.alpha.RR), statistic=meanfun, R=5000)
ci.weibull.alpha.RR <-boot.ci(boot.weibull.alpha.RR, conf=0.95)
specific.norm.ci.weibull.alpha.RR<-ci.weibull.alpha.RR[4]

weibull_MRT_HL_mean<-mean(weibull$MRT_HL.RR, na.rm=T)
boot.MRT_HL.RR <- boot(as.data.frame(weibull$MRT_HL.RR), statistic=meanfun, R=5000)
ci.MRT_HL.RR <-boot.ci(boot.MRT_HL.RR, conf=0.95)
specific.norm.ci.MRT_HL.RR<-ci.MRT_HL.RR[4]
as.data.frame(specific.norm.ci.MRT_HL.RR)[1,2]


#################################################################################################
#################################################### Figure 2 
both_means<-c(single_single_k_mean,data_single_k_mean, double_double_k_mean,data_double_k_mean, double_double_ks_mean, data_double_ks_mean,double_double_A_mean, data_double_A_mean, 
                  asymp_asymp_k_mean, data_asymp_k_mean,asymp_asymp_A_mean, data_asymp_A_mean,weibull_tenth_mean, data_tenth_mean,weibull_quarter_mean, data_quarter_mean, 
                  weibull_half_life_mean, data_half_life_mean, weibull_mrt_mean, data_mrt_mean, weibull_alpha_mean,data_alpha_mean,weibull_MRT_HL_mean,data_MRT_HL_mean)
both_upperCI<-c(as.data.frame(specific.norm.ci.single.k.RR)[1,3], as.data.frame(norm.ci.single.k.RR)[1,3],as.data.frame(specific.norm.ci.double.k.RR)[1,3], as.data.frame(norm.ci.double.k.RR)[1,3],as.data.frame(specific.norm.ci.double.ks.RR)[1,3], as.data.frame(norm.ci.double.ks.RR)[1,3],as.data.frame(specific.norm.ci.double.A.RR)[1,3], as.data.frame(norm.ci.double.A.RR)[1,3],
                    as.data.frame(specific.norm.ci.asymp.k.RR)[1,3], as.data.frame(norm.ci.asymp.k.RR)[1,3],as.data.frame(specific.norm.ci.asymp.A.RR)[1,3], as.data.frame(norm.ci.asymp.A.RR)[1,3],as.data.frame(specific.norm.ci.weibull.tenth.RR)[1,3], as.data.frame(norm.ci.weibull.tenth.RR)[1,3],as.data.frame(specific.norm.ci.weibull.quarter.RR)[1,3],as.data.frame(norm.ci.weibull.quarter.RR)[1,3],
                    as.data.frame(specific.norm.ci.weibull.half.life.RR)[1,3], as.data.frame(norm.ci.weibull.half.life.RR)[1,3],as.data.frame(specific.norm.ci.weibull.mrt.RR)[1,3], as.data.frame(norm.ci.weibull.mrt.RR)[1,3],as.data.frame(specific.norm.ci.weibull.alpha.RR)[1,3],as.data.frame(norm.ci.weibull.alpha.RR)[1,3],as.data.frame(specific.norm.ci.MRT_HL.RR)[1,3],as.data.frame(norm.ci.MRT_HL.RR)[1,3])
both_lowerCI<-c(as.data.frame(specific.norm.ci.single.k.RR)[1,2], as.data.frame(norm.ci.single.k.RR)[1,2],as.data.frame(specific.norm.ci.double.k.RR)[1,2], as.data.frame(norm.ci.double.k.RR)[1,2],as.data.frame(specific.norm.ci.double.ks.RR)[1,2], as.data.frame(norm.ci.double.ks.RR)[1,2],as.data.frame(specific.norm.ci.double.A.RR)[1,2], as.data.frame(norm.ci.double.A.RR)[1,2],
                    as.data.frame(specific.norm.ci.asymp.k.RR)[1,2], as.data.frame(norm.ci.asymp.k.RR)[1,2],as.data.frame(specific.norm.ci.asymp.A.RR)[1,2], as.data.frame(norm.ci.asymp.A.RR)[1,2],as.data.frame(specific.norm.ci.weibull.tenth.RR)[1,2], as.data.frame(norm.ci.weibull.tenth.RR)[1,2],as.data.frame(specific.norm.ci.weibull.quarter.RR)[1,2],as.data.frame(norm.ci.weibull.quarter.RR)[1,2],
                    as.data.frame(specific.norm.ci.weibull.half.life.RR)[1,2], as.data.frame(norm.ci.weibull.half.life.RR)[1,2],as.data.frame(specific.norm.ci.weibull.mrt.RR)[1,2], as.data.frame(norm.ci.weibull.mrt.RR)[1,2],as.data.frame(specific.norm.ci.weibull.alpha.RR)[1,2],as.data.frame(norm.ci.weibull.alpha.RR)[1,2],as.data.frame(specific.norm.ci.MRT_HL.RR)[1,2],as.data.frame(norm.ci.MRT_HL.RR)[1,2])

xs<-c(1.5,3.5,5.5,7.5,9.5,11.5,13.5,15.5,17.5,19.5,21.5,23.5)
new_xs<-c(1.25,1.75,3.25,3.75,5.25,5.75,7.25,7.75,9.25,9.75,11.25,11.75,13.25,13.75,15.25,15.75,17.25,17.75,19.25,19.75,21.25,21.75,23.25,23.75)
length(both_means)
circle<-rep(.3,24)

# tiff(file = "20200206_RR_2.tiff", width = 11, height = 11, units = "cm", res = 300)
# png(filename="20201206_RR.png", width = 11, height = 11, units = "cm", res=300)
par(mfrow=c(1,1),mar = c(8, 3.25, 2, 3.25), mgp=c(2,1,0), xpd=FALSE,cex.lab=1, cex.main=1, cex.axis=1.)
plot(new_xs, both_means, ylim=c(-1.2,4), xaxt="n", xlab="",ylab="Ln(Fertilized/Reference)", pch=16,  col=c("#76b7b2","dark grey"))
abline(h=0, col="#DD6F6A", lwd=2)
circle<-rep(.3,24)
symbols(x=new_xs, both_means, circles=circle,inches=1/14,ann=F, bg=c("#76b7b2","dark grey"), fg=NULL, add=T)
# errbar(new_xs, both_means, both_upperCI, both_lowerCI, cap=0.015, col=c("dark grey", "dodgerblue"), xaxt="n", add=T, lwd=3)
arrows(new_xs, both_lowerCI, new_xs, both_upperCI, angle = 90,length=0, lwd=3, col=c("#76b7b2","dark grey"),code=3)
axis(1, las=2,  at=xs,
     labels=c(expression('Single '*italic(k[s])*''), expression('Double '*italic(k[1])*''), expression('Double '*italic(k[2])*''), expression('Double '*italic(C)*''), expression('Asymp. '*italic(k[a])*''), expression('Asymp. '*italic(A)*'') ,expression('Weibull '*italic(t[1/10])*''), expression('Weibull '*italic(t[1/4])*''), expression('Weibull '*italic(t[1/2])*''),"Weibull MRT","Weibull alpha",
              expression('Weibull MRT/'*italic(t[1/2])*'')), cex.axis=1.20, tick=FALSE, pos=-1.28)
abline(v=2.5, lty=2, col="light grey")
abline(v=4.5, lty=2, col="light grey")
abline(v=6.5, lty=2, col="light grey")
abline(v=8.5, lty=2, col="light grey")
abline(v=10.5, lty=2, col="light grey")
abline(v=12.5, lty=2, col="light grey")
abline(v=14.5, lty=2, col="light grey")
abline(v=16.5, lty=2, col="light grey")
abline(v=18.5, lty=2, col="light grey")
abline(v=20.5, lty=2, col="light grey")
abline(v=22.5, lty=2, col="light grey")
legend("topright", c("Best Model","All Data") , lty=c(1, 1), 
       lwd=c(4, 4), col=c("#76b7b2","dark grey"), bty="n")
letters=c("","","","","","*","", "*","*","*","*", "*","*","*","*","*","","","","*","*","*","*","*")
text(new_xs,-1.1, pos=1, label=letters, cex=1.5, col=c("#76b7b2","dark grey"))
dev.off()

#Stats for figure 2 
t.test(single$single.k.RR, na.rm=T)
t.test(data$single.k.RR, na.rm=T)
t.test(double$double.k.RR, na.rm=T)
t.test(data$double.k.RR, na.rm=T)
t.test(double$double.ks.RR, na.rm=T)
t.test(data$double.ks.RR, na.rm=T)
t.test(double$double.A.RR, na.rm=T)
t.test(data$double.A.RR, na.rm=T)
t.test(asymp$asymp.k.RR, na.rm=T)
t.test(data$asymp.k.RR, na.rm=T)
t.test(asymp$asymp.A.RR, na.rm=T)
t.test(data$asymp.A.RR, na.rm=T)
t.test(weibull$weibull.tenth.RR, na.rm=T)
t.test(data$weibull.tenth.RR, na.rm=T)
t.test(weibull$weibull.quarter.RR, na.rm=T)
t.test(data$weibull.quarter.RR, na.rm=T)
t.test(weibull$weibull.half.life.RR, na.rm=T)
t.test(data$weibull.half.life.RR, na.rm=T)
t.test(weibull$weibull.mrt.RR, na.rm=T)
t.test(data$weibull.mrt.RR, na.rm=T)
t.test(weibull$weibull.alpha.RR, na.rm=T)
t.test(data$weibull.alpha.RR, na.rm=T)
t.test(weibull$MRT_HL.RR, na.rm=T)
t.test(data$MRT_HL.RR, na.rm=T)

#####################################################
################ Mixed effects models in Table S2
model <- lme(MRT_HL.RR~Sqrt_N_addition_kgN_ha_yr,random=~1|Paper,data=data, na.action=na.omit)
modela <- lme(MRT_HL.RR~Sqrt_N_addition_kgN_ha_yr+Years,random=~1|Paper,data=data, na.action=na.omit)
modelb <- lme(MRT_HL.RR~Sqrt_N_addition_kgN_ha_yr*Years,random=~1|Paper,data=data, na.action=na.omit)
modelc <- lme(MRT_HL.RR~Sqrt_N_addition_kgN_ha_yr,random=~1|Paper,data=weibull, na.action=na.omit)
modeld <- lme(MRT_HL.RR~Sqrt_N_addition_kgN_ha_yr+Years,random=~1|Paper,data=weibull, na.action=na.omit)
modele <- lme(MRT_HL.RR~Sqrt_N_addition_kgN_ha_yr*Years,random=~1|Paper,data=weibull, na.action=na.omit)

cbind(rbind(AIC(model)[1], df2[1,], AIC(modela)[1], df2[1,],df2[1,],AIC(modelb)[1], df2[1,], df2[1,],df2[1,],
            AIC(modelc)[1], df2[1,], AIC(modeld)[1], df2[1,],df2[1,],AIC(modele)[1], df2[1,]),
      rbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"][2])[c(1),1])), df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modela)$tTable[,"Value"][c(2:3)])[c(1:2),1])), 
            df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modelb)$tTable[,"Value"][c(2:4)])[c(1:3),1])), as.data.frame(as.matrix(as.data.frame(summary(modelc)$tTable[,"Value"][2])[c(1),1])),
            df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modeld)$tTable[,"Value"][c(2:3)])[c(1:2),1])), 
            df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modele)$tTable[,"Value"][c(2:4)])[c(1:3),1]))),
      rbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"p-value"][2])[c(1),1])), df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modela)$tTable[,"p-value"][c(2:3)])[c(1:2),1])), 
            df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modelb)$tTable[,"p-value"][c(2:4)])[c(1:3),1])), as.data.frame(as.matrix(as.data.frame(summary(modelc)$tTable[,"p-value"][2])[c(1),1])), df2[1,],
            as.data.frame(as.matrix(as.data.frame(summary(modeld)$tTable[,"p-value"][c(2:3)])[c(1:2),1])), 
            df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modele)$tTable[,"p-value"][c(2:4)])[c(1:3),1]))))


#################################
## Hierarchical mixed effects models reported in Table S3. Table reports parameters associated with most parsimonius
### model structure (dAIC < 3). Significant relationships are summarized in Table S4. The code in the following section
### should be updated to reflect each parameter (singl.k, double.k., etc.) and response (Latitude, MAT) etc. 
### Response variables were transformed to approximate normality as follows: Abs_Latitude, Sqrt_N_addition_kgN_ha_yr,
### Log_Litter_Ca, Sqrt_Litter_N, Log_MAP.
model <- lme(asymp.k.RR~Abs_Latitude,random=~1|Paper,data=data, na.action=na.omit)
modela <- lme(asymp.k.RR~Abs_Latitude+Sqrt_N_addition_kgN_ha_yr,random=~1|Paper,data=data, na.action=na.omit)
modelb <- lme(asymp.k.RR~Abs_Latitude*Sqrt_N_addition_kgN_ha_yr,random=~1|Paper,data=data, na.action=na.omit)
modelc <- lme(asymp.k.RR~Abs_Latitude+Duration_years,random=~1|Paper,data=data, na.action=na.omit)
modeld <- lme(asymp.k.RR~Abs_Latitude*Duration_years,random=~1|Paper,data=data, na.action=na.omit)
modele <- lme(asymp.k.RR~Abs_Latitude+Sqrt_N_addition_kgN_ha_yr+Duration_years,random=~1|Paper,data=data, na.action=na.omit)
modelf <- lme(asymp.k.RR~Abs_Latitude*Sqrt_N_addition_kgN_ha_yr*Duration_years,random=~1|Paper,data=data, na.action=na.omit)
modelg <- lme(asymp.k.RR~Abs_Latitude,random=~1|Paper,data=asymp, na.action=na.omit)
modelh <- lme(asymp.k.RR~Abs_Latitude+Sqrt_N_addition_kgN_ha_yr,random=~1|Paper,data=asymp, na.action=na.omit)
modeli <- lme(asymp.k.RR~Abs_Latitude*Sqrt_N_addition_kgN_ha_yr,random=~1|Paper,data=asymp, na.action=na.omit)
modelj <- lme(asymp.k.RR~Abs_Latitude+Duration_years,random=~1|Paper,data=asymp, na.action=na.omit)
modelk <- lme(asymp.k.RR~Abs_Latitude*Duration_years,random=~1|Paper,data=asymp, na.action=na.omit)
modell <- lme(asymp.k.RR~Abs_Latitude+Sqrt_N_addition_kgN_ha_yr+Duration_years,random=~1|Paper,data=asymp, na.action=na.omit)
modelm <- lme(asymp.k.RR~Abs_Latitude*Sqrt_N_addition_kgN_ha_yr*Duration_years,random=~1|Paper,data=asymp, na.action=na.omit)

cbind(rbind(AIC(model)[1], df2[1,], AIC(modela)[1], df2[1,],df2[1,],AIC(modelb)[1], df2[1,], df2[1,],df2[1,],
            AIC(modelc)[1], df2[1,], df2[1,],AIC(modeld)[1], df2[1,], df2[1,],df2[1,],AIC(modele)[1],df2[1,], df2[1,],df2[1,],
            AIC(modelf)[1],df2[1,], df2[1,],df2[1,],df2[1,], df2[1,],df2[1,],AIC(model)[1], df2[1,], AIC(modela)[1], df2[1,],df2[1,],AIC(modelb)[1], df2[1,], df2[1,],df2[1,],
            AIC(modelc)[1], df2[1,], df2[1,],AIC(modeld)[1], df2[1,], df2[1,],df2[1,],AIC(modele)[1],df2[1,], df2[1,],df2[1,],
            AIC(modelf)[1],df2[1,], df2[1,],df2[1,],df2[1,], df2[1,],df2[1,]),
      rbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"][2])[c(1),1])), df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modela)$tTable[,"Value"][c(2:3)])[c(1:2),1])), 
            df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modelb)$tTable[,"Value"][c(2:4)])[c(1:3),1])), df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modelc)$tTable[,"Value"][c(2:3)])[c(1:2),1])),
            df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modeld)$tTable[,"Value"][c(2:4)])[c(1:3),1])),df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modele)$tTable[,"Value"][c(2:4)])[c(1:3),1])),
            df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modelf)$tTable[,"Value"][c(2:8)])[c(1:7),1])),as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"][2])[c(1),1])), df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modela)$tTable[,"Value"][c(2:3)])[c(1:2),1])), 
            df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modelb)$tTable[,"Value"][c(2:4)])[c(1:3),1])), df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modelc)$tTable[,"Value"][c(2:3)])[c(1:2),1])),
            df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modeld)$tTable[,"Value"][c(2:4)])[c(1:3),1])),df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modele)$tTable[,"Value"][c(2:4)])[c(1:3),1])),
            df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modelf)$tTable[,"Value"][c(2:8)])[c(1:7),1]))),
      rbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"p-value"][2])[c(1),1])), df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modela)$tTable[,"p-value"][c(2:3)])[c(1:2),1])), 
            df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modelb)$tTable[,"p-value"][c(2:4)])[c(1:3),1])), df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modelc)$tTable[,"p-value"][c(2:3)])[c(1:2),1])),
            df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modeld)$tTable[,"p-value"][c(2:4)])[c(1:3),1])),df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modele)$tTable[,"p-value"][c(2:4)])[c(1:3),1])),
            df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modelf)$tTable[,"p-value"][c(2:8)])[c(1:7),1])),as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"p-value"][2])[c(1),1])), df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modela)$tTable[,"p-value"][c(2:3)])[c(1:2),1])), 
            df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modelb)$tTable[,"p-value"][c(2:4)])[c(1:3),1])), df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modelc)$tTable[,"p-value"][c(2:3)])[c(1:2),1])),
            df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modeld)$tTable[,"p-value"][c(2:4)])[c(1:3),1])),df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modele)$tTable[,"p-value"][c(2:4)])[c(1:3),1])),
            df2[1,],as.data.frame(as.matrix(as.data.frame(summary(modelf)$tTable[,"p-value"][c(2:8)])[c(1:7),1]))))

#############################################################
## Code associated with Figure 3

model <-lme(asymp.k.RR~Log_Litter_Ca,random=~1|Paper,data=data, na.action=na.omit)
r.squaredGLMM(model)
anova(model)
P<-ggpredict(model, "Log_Litter_Ca") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2, fill="grey") +
  geom_line(colour = "grey44", size=1, linetype="solid")

Q<-P +  xlab("Log(% Litter Ca)") + ylab(expression('Asymptotic '*italic(k[a])*' RR'))  + geom_hline(yintercept=0, linetype="dotted", 
                                                                                                    color = "#DD6F6A", size=1)+ ylim(-1, 1)+xlim(-2.25,1.5)
Asympk_LitterCa_Data<-Q + theme( panel.border = element_rect(colour = "black", fill=NA, size=1),panel.grid.major = element_blank(), axis.title.y=element_text(size=10), axis.title.x=element_text(size=12),panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                 plot.margin = margin(0.25, 0.5, 0.25, 0.25, "cm"),)+
  annotate("text",colour="black", x = -Inf, y = -Inf, hjust=-0.05, vjust=-2, label = expression('All: p=0.046; '*R^2*'-M=0.050'), cex=4)+
  annotate("text",colour="#76b7b2", x = -Inf, y = -Inf, hjust=-0.05, vjust=-1, label = expression('Best: p=0.17; '*R^2*'-M=0.22'), cex=4)+
  annotate("text",colour="black",x = Inf, y = Inf, hjust=2, vjust=1.75, label = expression(bold('A')), cex=5)

Asympk_LitterCa_Data


##
model <-lme(weibull.tenth.RR~Log_Litter_Ca,random=~1|Paper,data=data, na.action=na.omit)
# model <-lme(weibull.tenth.RR~Log_Litter_Ca,random=~1|Paper,data=weibull, na.action=na.omit)
r.squaredGLMM(model)
anova(model)

P<-ggpredict(model, "Log_Litter_Ca") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2, fill="grey") +
  geom_line(colour = "grey44", size=1, linetype="dashed")

Q<-P +  xlab("Log(% Litter Ca)") + ylab(expression('Weibull '*italic(t[1/10])*' RR'))  + geom_hline(yintercept=0, linetype="dotted", 
                                                                                                    color = "#DD6F6A", size=1)+ ylim(-2.5, 2.5)+xlim(-2.25,1.5)
Weib10_LitterCa_Data<-Q + theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.grid.major = element_blank(), axis.title.y=element_text(size=10), axis.title.x=element_text(size=12),panel.grid.minor = element_blank(),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                plot.margin = margin(0.25, 0.5, 0.25, 0.5, "cm"),)+
  annotate("text",colour="black", x = -Inf, y = -Inf, hjust=-0.05, vjust=-2, label = expression('All: p=0.079; '*R^2*'-M=0.048'), cex=4)+
  annotate("text",colour="#76b7b2",x = -Inf, y = -Inf, hjust=-0.05, vjust=-1, label = expression('Best: p=0.074; '*R^2*'-M=0.065'), cex=4)+
  annotate("text",colour="black",x = Inf, y = Inf, hjust=2, vjust=1.75, label = expression(bold('B')), cex=5)
Weib10_LitterCa_Data

##
model <-lme(weibull.quarter.RR~Log_Litter_Ca,random=~1|Paper,data=data, na.action=na.omit)
# model <-lme(weibull.quarter.RR~Log_Litter_Ca,random=~1|Paper,data=weibull, na.action=na.omit)
r.squaredGLMM(model)
anova(model)

P<-ggpredict(model, "Log_Litter_Ca") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2, fill="grey") +
  geom_line(colour = "grey44", size=1, linetype="dashed")

Q<-P +  xlab("Log(% Litter Ca)") + ylab(expression('Weibull '*italic(t[1/4])*' RR'))  + geom_hline(yintercept=0, linetype="dotted", 
                                                                                                   color = "#DD6F6A", size=1)+ ylim(-1.25, 1.25)+xlim(-2.25,1.5)
Weib25_LitterCa_Data<-Q + theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.grid.major = element_blank(), axis.title.y=element_text(size=10), axis.title.x=element_text(size=12),panel.grid.minor = element_blank(),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                plot.margin = margin(0.25, 0.5, 0.25, 0.25, "cm"),)+annotate("text",colour="black", x = -Inf, y = -Inf, hjust=-0.05, vjust=-2, label = expression('All: p=0.065; '*R^2*'-M=0.053'), cex=4)+annotate("text",colour="#76b7b2",
                                                                                                                                                                                                                                    x = -Inf, y = -Inf, hjust=-0.05, vjust=-1, label = expression('Best: p=0.081; '*R^2*'-M=0.051'), cex=4)+
  annotate("text",colour="black",x = Inf, y = Inf, hjust=2, vjust=1.75, label = expression(bold('C')), cex=5)
Weib25_LitterCa_Data

##
model <-lme(single.k.RR~pH,random=~1|Paper,data=data, na.action=na.omit)
# model <-lme(single.k.RR~pH,random=~1|Paper,data=single, na.action=na.omit)
r.squaredGLMM(model)
anova(model)

P<-ggpredict(model, "pH") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2, fill="grey") +
  geom_line(colour = "grey44", size=1, linetype="solid")

Q<-P +  xlab("Soil pH") + ylab(expression('Single '*italic(k[s])*' RR'))   + geom_hline(yintercept=0, linetype="dotted", 
                                                                                        color = "#DD6F6A", size=1)+ ylim(-.15, .325)+xlim(4.5,8.5)
Singlek_pH_Data<-Q + theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.grid.major = element_blank(), axis.title.y=element_text(size=10), axis.title.x=element_text(size=12),panel.grid.minor = element_blank(),
                           panel.background = element_blank(), axis.line = element_line(colour = "black"),
                           plot.margin = margin(0.25, 0.5, 0.25, 0.25, "cm"),)+annotate("text",colour="black", x = -Inf, y = -Inf, hjust=-0.05, vjust=-2, label = expression('All: p=0.043; '*R^2*'-M=0.05'), cex=4)+annotate("text",colour="#76b7b2",
                                                                                                                                                                                                                              x = -Inf, y = -Inf, hjust=-0.05, vjust=-1, label = expression('Best: p=0.015; '*R^2*'-M=0.10'), cex=4)+
  annotate("text",colour="black",x = Inf, y = Inf, hjust=2, vjust=1.75, label = expression(bold('D')), cex=5)
Singlek_pH_Data

##
model <-lme(double.A.RR~Litter_Lignin*Years,random=~1|Paper,data=data, na.action=na.omit)
# model <-lme(double.A.RR~Litter_Lignin+Sqrt_N_addition_kgN_ha_yr,random=~1|Paper,data=double, na.action=na.omit)
r.squaredGLMM(model)
summary(model)

P<-ggpredict(model, "Litter_Lignin") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2, fill="grey") +
  geom_line(colour = "grey44", size=1, linetype="solid")

Q<-P +  xlab("% Litter Lignin") + ylab(expression('Double '*italic(C)*' RR'))   + geom_hline(yintercept=0, linetype="dotted", 
                                                                                             color = "#DD6F6A", size=1)+ ylim(-2.5, 2)+xlim(0,50)
DoubleC_Lignin_Data<-Q + theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.grid.major = element_blank(), axis.title.y=element_text(size=10), axis.title.x=element_text(size=12),panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"),
                               plot.margin = margin(0.25, 0.5, 0.25, 0.5, "cm"),)+annotate("text",colour="black", x = -Inf, y = -Inf, hjust=-0.05, vjust=-2, label = expression('All: p < 0.001; '*R^2*'-M=0.21'), cex=4)+annotate("text",colour="#76b7b2",
                                                                                                                                                                                                                                   x = -Inf, y = -Inf, hjust=-0.05, vjust=-1, label = expression('Best: p=0.14; '*R^2*'-M=0.043'), cex=4)+
  annotate("text",colour="black",x = Inf, y = Inf, hjust=2, vjust=1.75, label = expression(bold('E')), cex=5)
DoubleC_Lignin_Data


##
model <-lme(asymp.A.RR~Litter_Lignin,random=~1|Paper,data=data, na.action=na.omit)
# model <-lme(asymp.A.RR~Litter_Lignin,random=~1|Paper,data=asymp, na.action=na.omit)
r.squaredGLMM(model)
anova(model)

P<-ggpredict(model, "Litter_Lignin") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2, fill="grey") +
  geom_line(colour = "grey44", size=1, linetype="solid")

Q<-P +  xlab("% Litter Lignin") + ylab(expression('Asymptotic '*italic(A)*' RR'))   + geom_hline(yintercept=0, linetype="dotted", 
                                                                                                 color = "#DD6F6A", size=1)+ ylim(-18, 25)+xlim(0,50)
AsympA_Lignin_Data<-Q + theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.grid.major = element_blank(), axis.title.y=element_text(size=10), axis.title.x=element_text(size=12),panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"),
                              plot.margin = margin(0.25, 0.5, 0.25, 0.25, "cm"),)+annotate("text",colour="black", x = -Inf, y = -Inf, hjust=-0.05, vjust=-2, label = expression('All: p=0.0001; '*R^2*'-M=0.09'), cex=4)+annotate("text",colour="#76b7b2",
                                                                                                                                                                                                                                  x = -Inf, y = -Inf, hjust=-0.05, vjust=-1, label = expression('Best: p<0.0001; '*R^2*'-M=0.19'), cex=4)+
  annotate("text",colour="black",x = Inf, y = Inf, hjust=2, vjust=1.75, label = expression(bold('F')), cex=5)
AsympA_Lignin_Data
##



SixPanel<-ggarrange(Asympk_LitterCa_Data, Singlek_pH_Data,
                    Weib10_LitterCa_Data,DoubleC_Lignin_Data,
                    Weib25_LitterCa_Data,AsympA_Lignin_Data,
                    ncol = 2, nrow = 3)

ggsave(filename = "ContRR_AllData3.pdf",
       plot = SixPanel,
       dpi=300,
       width = 18, height = 22, units = "cm")

##############################
##Figure 3 Code - Best model panel
##
model1 <-lme(asymp.k.RR~Log_Litter_Ca,random=~1|Paper,data=asymp, na.action=na.omit)
r.squaredGLMM(model1)
anova(model1)

P<-ggpredict(model1, "Log_Litter_Ca") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2, fill="#76b7b2") +
  geom_line(colour = "#76b7b2", size=1, linetype="dotted")

Q<-P + xlab("Log(% Litter Ca)") + ylab(expression('Asymptotic '*italic(k[a])*' RR'))  + ylim(-1, 1)+xlim(-2.25,1.5)
Asympk_LitterCa_Best<-Q + theme( panel.grid.major = element_blank(), axis.title.y=element_text(size=10), axis.title.x=element_text(size=12),panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                 plot.margin =  margin(0.25, 0.5, 0.25, 0.25, "cm"),)+theme(
                                   panel.background = element_rect(fill = "transparent"), # bg of the panel
                                   plot.background = element_rect(fill = "transparent", color = NA), # bg of the plo
                                 )
Asympk_LitterCa_Best


##
model1 <-lme(weibull.tenth.RR~Log_Litter_Ca,random=~1|Paper,data=weibull, na.action=na.omit)
r.squaredGLMM(model1)
anova(model1)
P<-ggpredict(model1, "Log_Litter_Ca") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2, fill="#76b7b2") +
  geom_line(colour = "#76b7b2", size=1, linetype="dashed")

Q<-P +  xlab("Log(% Litter Ca)") + ylab(expression('Weibull '*italic(t[1/10])*' RR')) + ylim(-2.5, 2.5)+xlim(-2.25,1.5)
Weib10_LitterCa_Best<-Q + theme( panel.grid.major = element_blank(), axis.title.y=element_text(size=10), axis.title.x=element_text(size=12),panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                 plot.margin = margin(0.25, 0.5, 0.25, 0.5, "cm"),)+theme(
                                   panel.background = element_rect(fill = "transparent"), # bg of the panel
                                   plot.background = element_rect(fill = "transparent", color = NA), # bg of the plo
                                 )
Weib10_LitterCa_Best

##
model1 <-lme(weibull.quarter.RR~Log_Litter_Ca,random=~1|Paper,data=weibull, na.action=na.omit)
r.squaredGLMM(model1)
anova(model1)
P<-ggpredict(model1, "Log_Litter_Ca") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2, fill="#76b7b2") +
  geom_line(colour = "#76b7b2", size=1, linetype="dashed")

Q<-P +   xlab("Log(% Litter Ca)") + ylab(expression('Weibull '*italic(t[1/4])*' RR'))  + ylim(-1.25, 1.25)+xlim(-2.25,1.5)
Weib25_LitterCa_Best<-Q + theme( panel.grid.major = element_blank(), axis.title.y=element_text(size=10), axis.title.x=element_text(size=12),panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                 plot.margin = margin(0.25, 0.5, 0.25, 0.25, "cm"),)+theme(
                                   panel.background = element_rect(fill = "transparent"), # bg of the panel
                                   plot.background = element_rect(fill = "transparent", color = NA), # bg of the plo
                                 )
Weib25_LitterCa_Best


##
model1 <-lme(single.k.RR~pH,random=~1|Paper,data=single, na.action=na.omit)
r.squaredGLMM(model1)
anova(model1)
P<-ggpredict(model1, "pH") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2, fill="#76b7b2") +
  geom_line(colour = "#76b7b2", size=1, linetype="solid")

Q<-P +   xlab("Soil pH") + ylab(expression('Single '*italic(k[s])*' RR'))  + ylim(-.15, .325)+xlim(4.5,8.5)
Singlek_pH_Best<-Q + theme( panel.grid.major = element_blank(), axis.title.y=element_text(size=10), axis.title.x=element_text(size=12),panel.grid.minor = element_blank(),
                            panel.background = element_blank(), axis.line = element_line(colour = "black"),
                            plot.margin = margin(0.25, 0.5, 0.25, 0.25, "cm"),)+theme(
                              panel.background = element_rect(fill = "transparent"), # bg of the panel
                              plot.background = element_rect(fill = "transparent", color = NA), # bg of the plo
                            )
Singlek_pH_Best

##
# model <-lme(double.A.RR~Litter_Lignin*Years,random=~1|Paper,data=data, na.action=na.omit)
model <-lme(double.A.RR~Litter_Lignin+Sqrt_N_addition_kgN_ha_yr,random=~1|Paper,data=double, na.action=na.omit)
anova(model)

P<-ggpredict(model, "Litter_Lignin") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2, fill="#76b7b2") +
  geom_line(colour = "#76b7b2", size=1, linetype="dotted")

Q<-P +  xlab("% Litter Lignin") + ylab(expression('Double '*italic(C)*' RR'))   +  ylim(-2.5, 2)+xlim(0,50)
DoubleC_Lignin_Best<-Q + theme( panel.grid.major = element_blank(), axis.title.y=element_text(size=10), axis.title.x=element_text(size=12),panel.grid.minor = element_blank(),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                plot.margin = margin(0.25, 0.5, 0.25, 0.5, "cm"),)+theme(
                                  panel.background = element_rect(fill = "transparent"), # bg of the panel
                                  plot.background = element_rect(fill = "transparent", color = NA), # bg of the plo
                                )
DoubleC_Lignin_Best


##
model1 <-lme(asymp.A.RR~Litter_Lignin,random=~1|Paper,data=single, na.action=na.omit)
r.squaredGLMM(model1)
anova(model1)
P<-ggpredict(model1, "Litter_Lignin") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2, fill="#76b7b2") +
  geom_line(colour = "#76b7b2", size=1, linetype="solid")

Q<-P +   xlab("% Litter Lignin") + ylab(expression('Asymptotic '*italic(A)*' RR'))  + ylim(-18, 25)+xlim(0,50)
AsympA_Lignin_Best<-Q + theme( panel.grid.major = element_blank(), axis.title.y=element_text(size=10), axis.title.x=element_text(size=12),panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"),
                               plot.margin = margin(0.25, 0.5, 0.25, 0.25, "cm"),)+theme(
                                 panel.background = element_rect(fill = "transparent"), # bg of the panel
                                 plot.background = element_rect(fill = "transparent", color = NA), # bg of the plo
                               )
AsympA_Lignin_Best


SixPanel2<-ggarrange( Asympk_LitterCa_Best, Singlek_pH_Best,
                      Weib10_LitterCa_Best,DoubleC_Lignin_Best,
                      Weib25_LitterCa_Best,AsympA_Lignin_Best,
                      ncol = 2, nrow = 3)

ggsave(filename = "ContRR_Best2.pdf",
       plot = SixPanel2,
       bg = "transparent", dpi=300,
       width = 18, height = 22, units = "cm")

#############################################################
#Figure 4 All Data
##
model <- lme(asymp.A.RR~asymp.k.RR,random=~1|Paper,data=data, na.action=na.omit)
r.squaredGLMM(model)
summary(model)
# model <-lme(asymp.A.RR~asymp.k.RR,random=~1|Paper,data=asymp, na.action=na.omit)
# r.squaredGLMM(model)
# summary(model)

P<-ggpredict(model, "asymp.k.RR") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2, fill="grey") +
  geom_line(colour = "grey44", size=1, linetype="solid")

Q<-P +  xlab(expression('Asymptotic '*italic(k[a])*' RR')) + ylab(expression('Asymptotic '*italic(A)*' RR'))   + geom_hline(yintercept=0, linetype="dotted", 
                                                                                                                            color = "#DD6F6A", size=1)+ geom_vline(xintercept=0, linetype="dotted", color = "#DD6F6A", size=1)+ ylim(-25,35)+xlim(-2,3)
AsympARR_AsympkRR_Data<-Q + theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.grid.major = element_blank(), axis.title.y=element_text(size=10), axis.title.x=element_text(size=12),panel.grid.minor = element_blank(),
                                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                  plot.margin = margin(0.25, 0.5, 0.25, 0.5, "cm"),)+annotate("text",colour="black", x = -Inf, y = -Inf, hjust=-0.05, vjust=-2, label = expression('All: p<0.0001; '*R^2*'-M=0.11'), cex=3)+annotate("text",colour="#76b7b2",
                                                                                                                                                                                                                                     x = -Inf, y = -Inf, hjust=-0.05, vjust=-1, label = expression('Best: p<0.0001; '*R^2*'-M=0.08'), cex=3)+
  annotate("text",colour="black",x = Inf, y = Inf, hjust=2, vjust=1.75, label = expression(bold('A')), cex=5)
AsympARR_AsympkRR_Data

##
model <- lme(weibull.mrt.RR~weibull.tenth.RR,random=~1|Paper,data=data, na.action=na.omit)
r.squaredGLMM(model)
summary(model)
# model <-lme(weibull.mrt.RR~weibull.tenth.RR,random=~1|Paper,data=weibull, na.action=na.omit)
# r.squaredGLMM(model)
# summary(model)

P<-ggpredict(model, "weibull.tenth.RR") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2, fill="grey") +
  geom_line(colour = "grey44", size=1, linetype="solid")

Q<-P +  xlab(expression('Weibull '*italic(t[1/10])*' RR')) + ylab(expression('Weibull MRT RR'))   + geom_hline(yintercept=0, linetype="dotted", 
                                                                                                               color = "#DD6F6A", size=1)+ geom_vline(xintercept=0, linetype="dotted", color = "#DD6F6A", size=1)+ylim(-5,12)+xlim(-10,5)
WeibMRT_Weib10ML_Data<-Q + theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.grid.major = element_blank(), axis.title.y=element_text(size=10), axis.title.x=element_text(size=12),panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                 plot.margin = margin(0.25, 0.5, 0.25, 0.5, "cm"),)+annotate("text",colour="black", x = -Inf, y = -Inf, hjust=-0.05, vjust=-2, label = expression('All: p=<0.0001; '*R^2*'-M=0.40'), cex=3)+annotate("text",colour="#76b7b2",
                                                                                                                                                                                                                                     x = -Inf, y = -Inf, hjust=-0.05, vjust=-1, label = expression('Best: p<0.0001; '*R^2*'-M=0.31'), cex=3)+
  annotate("text",colour="black",x = Inf, y = Inf, hjust=2, vjust=1.75, label = expression(bold('B')), cex=5)
WeibMRT_Weib10ML_Data

##
model <- lme(weibull.mrt.RR~weibull.quarter.RR,random=~1|Paper,data=data, na.action=na.omit)
r.squaredGLMM(model)
summary(model)
# model <-lme(weibull.mrt.RR~weibull.quarter.RR,random=~1|Paper,data=weibull, na.action=na.omit)
# r.squaredGLMM(model)
# summary(model)

P<-ggpredict(model, "weibull.quarter.RR") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2, fill="grey") +
  geom_line(colour = "grey44", size=1, linetype="solid")

Q<-P +  xlab(expression('Weibull '*italic(t[1/4])*' RR')) + ylab(expression('Weibull MRT RR'))   + geom_hline(yintercept=0, linetype="dotted", 
                                                                                                              color = "#DD6F6A", size=1)+ geom_vline(xintercept=0, linetype="dotted", color = "#DD6F6A", size=1)+ ylim(-2,5.5)+xlim(-5,2)
WeibMRT_Weib25ML_Data<-Q + theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.grid.major = element_blank(), axis.title.y=element_text(size=10), axis.title.x=element_text(size=12),panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                 plot.margin = margin(0.25, 0.5, 0.25, 0.5, "cm"),)+annotate("text",colour="black", x = -Inf, y = -Inf, hjust=-0.05, vjust=-2, label = expression('All: p=<0.0001; '*R^2*'-M=0.05'), cex=3)+annotate("text",colour="#76b7b2",
                                                                                                                                                                                                                                     x = -Inf, y = -Inf, hjust=-0.05, vjust=-1, label = expression('Best: p=0.021; '*R^2*'-M=0.013'), cex=3)+
  annotate("text",colour="black",x = Inf, y = Inf, hjust=2, vjust=1.75, label = expression(bold('C')), cex=5)
WeibMRT_Weib25ML_Data



ThreePanel<-ggarrange( AsympARR_AsympkRR_Data,NA,
                       WeibMRT_Weib10ML_Data,NA,
                       WeibMRT_Weib25ML_Data,NA,
                       ncol = 2, nrow = 3)

ggsave(filename = "BothRR_Data2.pdf",
       plot = ThreePanel,
       bg = "transparent", dpi=300,
       width = 18, height = 22, units = "cm")

##########################################
#Figure 4 Best Model

model <-lme(asymp.A.RR~asymp.k.RR,random=~1|Paper,data=asymp, na.action=na.omit)
r.squaredGLMM(model)
summary(model)

P<-ggpredict(model, "asymp.k.RR") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2, fill="#76b7b2") +
  geom_line(colour = "#76b7b2", size=1, linetype="solid")

Q<-P + xlab(expression('Asymptotic '*italic(k[a])*' RR')) + ylab(expression('Asymptotic '*italic(A)*' RR'))   + geom_hline(yintercept=0, linetype="dotted", 
                                                                                                                           color = "#DD6F6A", size=1)+ ylim(-25,35)+xlim(-2,3)
AsympARR_AsympkRR_Best<-Q + theme( panel.grid.major = element_blank(), axis.title.y=element_text(size=10), axis.title.x=element_text(size=12),panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                   plot.margin = margin(0.25, 0.5, 0.25, 0.5, "cm"),)+theme(
                                     panel.background = element_rect(fill = "transparent"), # bg of the panel
                                     plot.background = element_rect(fill = "transparent", color = NA), # bg of the plo
                                   )
AsympARR_AsympkRR_Best


##
model <-lme(weibull.mrt.RR~weibull.tenth.RR,random=~1|Paper,data=weibull, na.action=na.omit)
r.squaredGLMM(model)
summary(model)

P<-ggpredict(model, "weibull.tenth.RR") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2, fill="#76b7b2") +
  geom_line(colour = "#76b7b2", size=1, linetype="solid")

Q<-P +  xlab(expression('Weibull '*italic(t[1/10])*' RR')) + ylab(expression('Weibull MRT RR'))   + geom_hline(yintercept=0, linetype="dotted", 
                                                                                                               color = "#DD6F6A", size=1)+ ylim(-5,12)+xlim(-10,5)
WeibMRT_Weib10ML_Best<-Q + theme( panel.grid.major = element_blank(), axis.title.y=element_text(size=10), axis.title.x=element_text(size=12),panel.grid.minor = element_blank(),
                                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                  plot.margin = margin(0.25, 0.5, 0.25, 0.5, "cm"),)+theme(
                                    panel.background = element_rect(fill = "transparent"), # bg of the panel
                                    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plo
                                  )                                                                                                                                                                                                                   
WeibMRT_Weib10ML_Best

##
model <-lme(weibull.mrt.RR~weibull.quarter.RR,random=~1|Paper,data=weibull, na.action=na.omit)
# r.squaredGLMM(model)
# summary(model)

P<-ggpredict(model, "weibull.quarter.RR") %>%
  ggplot(aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2, fill="#76b7b2") +
  geom_line(colour = "#76b7b2", size=1, linetype="solid")

Q<-P +  xlab(expression('Weibull '*italic(t[1/4])*' RR')) + ylab(expression('Weibull MRT RR'))   + geom_hline(yintercept=0, linetype="dotted", 
                                                                                                              color = "#DD6F6A", size=1)+ ylim(-2,5.5)+xlim(-5,2)
WeibMRT_Weib25ML_Best<-Q +  theme( panel.grid.major = element_blank(), axis.title.y=element_text(size=10), axis.title.x=element_text(size=12),panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                   plot.margin = margin(0.25, 0.5, 0.25, 0.5, "cm"),)+theme(
                                     panel.background = element_rect(fill = "transparent"), # bg of the panel
                                     plot.background = element_rect(fill = "transparent", color = NA), # bg of the plo
                                   )        
WeibMRT_Weib25ML_Best



ThreePanel2<-ggarrange( AsympARR_AsympkRR_Best,NA,
                        WeibMRT_Weib10ML_Best,NA,
                        WeibMRT_Weib25ML_Best,NA,
                        ncol = 2, nrow = 3)

ggsave(filename = "BothRR_Best.pdf",
       plot = ThreePanel2,
       bg = "transparent", dpi=300,
       width = 18, height = 22, units = "cm")

#####################################

#################################
