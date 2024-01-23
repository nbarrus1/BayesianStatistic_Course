#Lecture 3, estimating a mean parameter

library(R2jags)
library(tidyverse)
library(ggmcmc)
theme_set(theme_bw())
#Set correct working directory
#setwd("~/Dropbox/Babcock Bayesian/lectures/lecture 2. Estimate a mean")

#Example data with frequentist calculations
Y<-c(42, 43, 58, 70, 47, 51, 85, 63, 58, 46)
meanval<-mean(Y)
sdval<-sd(Y)
standard.error<-function(x) sd(x)/sqrt(length(x))
standard.error(Y)
CI<-c(mean(Y)-1.96*standard.error(Y),mean(Y)+1.96*standard.error(Y))
CI

#Likelihood grid search
meanGrid<-data.frame(m=seq(40,70,0.1))
for(i in 1:nrow(meanGrid)) {
  meanGrid$Likelihood[i]<-prod(dnorm(Y,meanGrid$m[i],sd(Y)))
  meanGrid$LL[i]<-sum(dnorm(Y,meanGrid$m[i],sd(Y)),log=TRUE)
}
meanGrid$LL<-log(meanGrid$Likelihood)

ggplot(meanGrid,aes(x=m,y=Likelihood*1E18))+
  geom_line(linewidth=1.5)+
  theme_bw()+
  ylab(bquote("Likelihood " ~x10^-18))
meanGrid$m[meanGrid$Likelihood==max(meanGrid$Likelihood)]

## Bayesian grid search
priormean<-53
priorsd<-5
for(i in 1:nrow(meanGrid)) {
  meanGrid$logPrior[i]<-dnorm(meanGrid$m[i],priormean,priorsd,log=TRUE)
  meanGrid$logPost[i]<-meanGrid$LL[i]+meanGrid$logPrior[i]
}  
ggplot(meanGrid,aes(x=m))+
  geom_line(aes(y=exp(logPrior)/sum(exp(logPrior))), linewidth=1.5,color="red")+
  geom_line(aes(y=exp(logPost)/sum(exp(logPost))), linewidth=1.5, color="blue")+
  geom_line(aes(y=exp(LL)/sum(exp(LL))), linewidth=1.5)+
  theme_bw()+
  ylab("Probablity")
meanGrid$m[meanGrid$logPost==max(meanGrid$logPost)]

# Bayesian analytic solution with varying prior sd
meanval<-mean(Y)
sdval<-sd(Y)
BayesMean<-data.frame(prior.sd=c(5,10,25,100,1000))
BayesMean$posterior.mean<-(priormean/BayesMean$prior.sd^2+meanval*10/sdval^2)/(1/BayesMean$prior.sd^2+10/sdval^2)
BayesMean$posterior.sd<-1/sqrt(10/sdval^2+1/BayesMean$prior.sd^2)
round(BayesMean,2)

#Minimization to find mode of posterior

getLogPosterior<-function(par,data) {
  x<-data$x
  m<-par[1]
  LL<-dnorm(x,m,sdval,log=TRUE)
  prior<-dnorm(m,53,5,log=TRUE)
  -sum(LL+prior)
}
posteriorMode<-optim(par=50,
                     fn=getLogPosterior,
                     data=data.frame(x=Y),
                     method="Brent",
                     lower=-100,upper=200)
posteriorMode

## Run the Bayesian model with MCMC in JAGS ###

write("model{
  m ~ dnorm(53, 0.04)	# prior for mean, with sd=5, so precision=1/5^2
#  m ~ dnorm(0, 1.0E-6)	# use this for an uninformative prior
  prec <- 1 / ( stdev*stdev)	# precision of the data = 1/variance
  for (i in 1:10)		# for each of the ten trees
  {
    Y[i] ~ dnorm(m, prec) # diameter drawn from normal (likelihood)
  }
}", file="Box1.8NormalMeanKnownVariance.txt")

meanResult<-jags(data=list(Y=Y,stdev=sd(Y)),
  parameters.to.save = c("m"),
  model.file = "Box1.8NormalMeanKnownVariance.txt",
  n.iter=100000,
  n.burnin=1000
)

meanResult$BUGSoutput$summary  

##ggmcmc for looking at whether the model converged
gg1<-ggs(as.mcmc(meanResult))  #Make the MCMC object a ggs object for plotting
head(gg1)
table(gg1$Parameter)
summary(gg1)
gg1<-filter(gg1,Parameter=="m")
ggs_traceplot(gg1)  # Trace plot
ggs_density(gg1)  #Density plot

#Summary statistics
round(meanResult$BUGSoutput$summary  ,3)
print(meanResult)


