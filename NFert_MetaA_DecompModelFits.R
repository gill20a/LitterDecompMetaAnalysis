#Template to fit single, double, asymptotic, and weibull litter decomposition models. Revised by ALG 2020-02-20

#Open necessary libraries
library(readr)
library(stats4)
library(bbmle)

#Download data file. Users will need to update link to reflect their own file path. 
urlfile = "https://raw.githubusercontent.com/gill20a/LitterDecompMetaAnalysis/master/Gill_NFertilization_MetaAnalysis_LitterHarvestData.csv"
data<-read.csv(url(urlfile))
data$propinit<-as.numeric(as.character(data$propinit))

#Create empty vectors to store fitted model parameters
single.k=array(0,dim=c(max(data$DecompID, na.rm=T),1))
double.k=array(0,dim=c(max(data$DecompID, na.rm=T),3))
asymptotic.k=array(0,dim=c(max(data$DecompID, na.rm=T),2))
weibull=array(0,dim=c(max(data$DecompID, na.rm=T),7))

#Create Empty array for AIC values
a.aicc=array(0,dim=c(max(data$DecompID, na.rm=T),4)) #3aicc; 3weights; single, asym, double

##Functions for fitting single, double, and asymptotic exponential models
LLS = function(y,k){
  Mhat=1*exp(-k*xNA$t) 
  ssq = sum((Mhat - xNA$Mt)^2)
  sigma = sqrt(ssq/length(xNA$Mt))
  return(-sum(dnorm(xNA$Mt,mean=Mhat,sd=sigma,log = TRUE)))
}

LLA = function(y,k,A){
  Mhat=A+(1-A)*exp(-k*xNA$t) 
  ssq = sum((Mhat - xNA$Mt)^2)
  sigma = sqrt(ssq/length(xNA$Mt))
  return(-sum(dnorm(xNA$Mt,mean=Mhat,sd=sigma,log = TRUE)))
}

LLD = function(y,ks,k,A){
  Mhat=A*exp(-ks*xNA$t)+(1-A)*exp(-k*xNA$t)
  ssq = sum((Mhat - xNA$Mt)^2)
  #browser()
  sigma = sqrt(ssq/length(xNA$Mt))
  return(-sum(dnorm(xNA$Mt,mean=Mhat,sd=sigma,log = TRUE)))
}

##Function to calculate litter tenth, quarter, and half life associated with weibull model
quarter.life.calc=function(nls.mod){
  pars= coef(nls.mod)
  hl=pars[1] * (log(1/(3/4)))^(1/pars[2])
  names(hl)="half.life"
  return(hl)
}

tenth.life.calc=function(nls.mod){
  pars= coef(nls.mod)
  hl=pars[1] * (log(1/(9/10)))^(1/pars[2])
  names(hl)="half.life"
  return(hl)
}
half.life.calc=function(nls.mod){
  pars= coef(nls.mod)
  hl=pars[1] * (log(2))^(1/pars[2])
  names(hl)="half.life"
  return(hl)
}

##Function to calculate litter mean residence time associated with weibull model
mrt.calc=function(nls.mod){
  pars= coef(nls.mod)
  mrt=pars[1] * gamma(1+(1/pars[2]))
  names(mrt)="mrt"
  return(mrt)
}





######################################################################################################
i=1 ##This should be updated to reflect the DecompID for each decay series you would like to ft. It will also reflect the
#### row in which the model parameters are stored in summary file

s1=subset(data,data$DecompID==i)
t=s1$Years
t1<-seq(0,5,length.out=10000)
Mt=s1$propinit
x <- data.frame(t, Mt)
xNA <-na.exclude(x) #omit lines with missing data (NAs)


#Weibull 
fit<- nls(Mt ~ exp(-(t/beta)^alpha), start =list(beta=1, alpha = 1), algorithm="port", lower=c(0.0001,0.0001))
weibull[i,1]=coef(fit)[2] #weibull alpha value
weibull[i,2]=coef(fit)[1] #weibull beta value
weibull[i,3]<-mrt.calc(fit)
weibull[i,4]<-half.life.calc(fit)
weibull[i,5]<-sum(resid(fit)^2)
weibull[i,6]<-quarter.life.calc(fit)
weibull[i,7]<-tenth.life.calc(fit)

weibull.fit<-exp(-(t1/weibull[i,2])^weibull[i,1]) #For figure
weibull.AIC<-AIC(fit)


#Single pool 
singleLL = mle2(minuslogl = LLS, start = list(k = 0.5), data = list(y=xNA$Mt),method="L-BFGS-B", lower=c(0),upper=c(1000))
summary(singleLL)
single.k[i,1]=coef(singleLL)[1] #Extract k value
single.exp.fit<-exp(-single.k[i,1]*t1) #For figure


#Asymptotic: M(t)=A+(1-A)*exp(-k*t)
asymLL = mle2(minuslogl = LLA, start = list(k = .9,A=0.2), data = list(y=xNA$Mt),method="L-BFGS-B", lower=c(0.00001,0.00001),upper=c(1000,1000))#,control=list(maxint=10e6))
summary(asymLL)

#Extract k & A value
asymptotic.k[i,1]=coef(asymLL)[1] #ka value
asymptotic.k[i,2]=coef(asymLL)[2] #A value
asymptotic.k.fit<-asymptotic.k[i,2]+((1-asymptotic.k[i,2])*exp(-asymptotic.k[i,1]*t1)) #For figure


#Double pool: M(t) = A*exp(-ks*t)+(1-A)*exp(-k*t);
#Liklihood function	
start=list(k=1,ks=.3,A =.4)
doubleLL = mle2(minuslogl = LLD, start = start, data = list(y=xNA$Mt),method="L-BFGS-B", lower=c(0.00001,0.00001,0.00001),upper=c(1000,1000,1000),control=list(maxit=1000, trace=TRUE, parscale=abs(unlist(start))))
summary(doubleLL)

#Extract k value; ks&A, k&(1-A)
double.k[i,1]=coef(doubleLL)[1]#ks
double.k[i,2]=coef(doubleLL)[2]#k
double.k[i,3]=coef(doubleLL)[3]#C
double.exp.fit<-(double.k[i,3]*exp(-double.k[i,1]*t1))+((1-double.k[i,3])*exp(-double.k[i,2]*t1))


#SummaryData
Sample1<-AICctab(singleLL,asymLL,doubleLL, nobs=nrow(x), sort=FALSE,base=TRUE, weights=TRUE)
Sample1

a.aicc[i,1]=Sample1$AICc[1]
a.aicc[i,2]=Sample1$AICc[2]
a.aicc[i,3]=Sample1$AICc[3]
a.aicc[i,4]=weibull.AIC

par(mfrow=c(3,4), oma=c(0,0,0,0), mar=c(4,4,1,1), mgp=c(2.4,0.8,0), cex.lab=1.25, cex.main=1.25, cex.axis=1.15)
plot(t,Mt, pch=16, col="dark grey", xlab="Years",ylab="Prop. Remaining", xlim=c(0,5), ylim=c(0,1), main="Sample1")
lines(t1,exp(-single.k[i,1]*t1), col="grey")
lines(t1, asymptotic.k.fit, col="mediumpurple", lty=2)
lines(t1, double.exp.fit, col="orange", lty=3)
lines(t1, weibull.fit, col="dodgerblue1", lty=1)
#legend("topright", c("single", "asymp.", "double", "weibull"), col=c("grey","mediumpurple","orange","dodgerblue1"), lty=c(1,2,3,1), bty="n")


##############################
###Bind parameter output together after fitting models to all relevant decomp sequences

colnames(a.aicc)<-c("Single_AICc", "Asym_AICc", "Double_AICc", "weibull.AICc")
head(a.aicc)

ParameterSummary<-cbind(
single.k,
double.k,
asymptotic.k,
weibull, a.aicc)

colnames(ParameterSummary)<-c("single.k", "double.ks", "double.k","double.A", "asymp.k", "asymp.A", "weibull.alpha","weibull.beta",  "weibull.mrt", 
              "weibull.half.life", "weibull.RSS", "weibull.quarter", "weibull.tenth", "Single_AICc", "Asym_AICc", "Double_AICc", "weibull.AICc")

#Export data
write.csv(a.aicc, "a.aicc_table.csv")
write.csv(k_vals, "20190801_k_vals_table.csv")

