#Script Title: Bayesian Statistics Homework2-Normal, bionomial, neg-binomial
#Date first created: 2/1/2024
#purpose: The purpose of this script is to work through Homework 2 which 
#contains problems with normal, binomial, and negative binomial. The problems 
#related to lectures 3 and 4 

#----------------------------
##general start to scripting##
#----------------------------

rm(list = ls())

######github#####
#note, only needed after 90 days from 1/16/2024

#  usethis::create_github_token()  
#  gitcreds::gitcreds_set()

#####check for r updates#####
#note, updateing may take some time so plan accordingly

require(installr)

check.for.updates.R()

#updateR() #only if needed

#######check for package updates#####
#note, updateing may take some time so plan accordingly

old.packages()

update.packages()

#----------------------------
####load packages####
#----------------------------

library(tidyverse)
library(R2jags)
library(ggmcmc)
library(here)
theme_set(theme_bw(base_size=15))


#-------------------------------
####1. Mean fish weight####
#-------------------------------

#A) Calculate the frequentist mean and standard deviation of the mean 
#(i.e. standard error) for the following 20 data points of fish weight in 
#kg: weightVals<-c(46, 51, 35, 46, 26, 37, 49, 48, 30, 46, 21, 26, 37, 36, 49,
#42, 50, 33, 52, 37). 

fish_kg <- c(46, 51, 35, 46, 26, 37, 49, 48, 30, 46, 21, 26, 37, 36, 49, 42, 
             50, 33, 52, 37)

fish_kg.ave <- mean(fish_kg)
fish_kg.sd <- sd(fish_kg)

fish_kg.ave
fish_kg.sd

#B) b. Calculate the mean and standard deviation of the Bayesian posterior 
#distribution of the mean fish weight for the same data with an informative 
#Normal prior with mean of 35 and standard deviation of 2, using the analytic 
#solutions for the Normal distribution (assuming the standard deviation is 
#known and equal to the standard deviation calculated from the data). 
#Use this posterior distribution to calculate the 95% credible interval of 
#the posterior distribution, using 1.96 as the critical value of the normal 
#distribution. Show your work. 

###done on paper

#C) Use a grid search to find the maximum posterior probability estimate for 
#the mean using the same prior as in part b, assuming the standard deviation is
#known and equal to the standard deviation calculated from the data 
#(ie sd(weigthVals)). 

#Likelihood grid search
meanGrid<-data.frame(m=seq(min(fish_kg)-20,max(fish_kg)+20,0.1))
meanGrid$Likelihood <- 0
meanGrid$LL <- 0

for(i in 1:nrow(meanGrid)) {
  meanGrid$Likelihood[i]<-prod(dnorm(fish_kg,meanGrid$m[i],sd(fish_kg)))
  meanGrid$LL[i]<-sum(dnorm(fish_kg,meanGrid$m[i],sd(fish_kg)),log=TRUE)
}
meanGrid$LL<-log(meanGrid$Likelihood)

ggplot(meanGrid,aes(x=m,y=Likelihood*1E18))+
  geom_line(linewidth=1.5)+
  theme_bw()+
  ylab(bquote("Likelihood " ~x10^-18))

meanGrid$m[meanGrid$Likelihood==max(meanGrid$Likelihood)]

## Bayesian grid search


priormean<-35
priorsd<-2

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


#D)Modify the JAGS code in Box 1.8 to calculate the same posterior distribution.
#With a burn in of 10000, and a sample size of 100000 iterations, what were 
#the mean and its standard deviation and 95% C.I.?

write("model{
  m ~ dnorm(35, 0.25)	# prior for mean, with sd=2, so precision=1/2^2
  prec <- 1 / ( stdev*stdev)	# precision of the data = 1/variance
  for (i in 1:20)		# for each of the ten trees
  {
    Y[i] ~ dnorm(m, prec) # diameter drawn from normal (likelihood)
  }
}", file=here("JAGS_mods","HW2_1d_mean.txt"))

meanResult_d<-jags(data=list(Y=fish_kg,stdev=sd(fish_kg)),
                 parameters.to.save = c("m"),
                 model.file = here("JAGS_mods","HW2_1d_mean.txt"),
                 n.iter=100000,
                 n.burnin=10000
)

meanResult_d$BUGSoutput$summary  

#E) . Did you get the same estimates for the mean in part b, c and d?  If not, 
#why are they different and which (if any) is correct? 

### see document

#F) Now run the JAGS code again with prior mean of 35 and standard deviation of 
#1000. Give the mean and 95% C.I. with the same burn in and number of iterations.
#How and why do the results differ from part d?

write("model{
  m ~ dnorm(35, 1.0E-06)	# prior for mean, with sd=1000, so precision=1/1000^2
  prec <- 1 / ( stdev*stdev)	# precision of the data = 1/variance
  for (i in 1:20)		# for each of the ten trees
  {
    Y[i] ~ dnorm(m, prec) # diameter drawn from normal (likelihood)
  }
}", file=here("JAGS_mods","HW2_1f_mean.txt"))

meanResult_f<-jags(data=list(Y=fish_kg,stdev=sd(fish_kg)),
                 parameters.to.save = c("m"),
                 model.file = here("JAGS_mods","HW2_1f_mean.txt"),
                 n.iter=100000,
                 n.burnin=10000
)

meanResult_f$BUGSoutput$summary  


#G) Based on Rhat and n.eff from the summary, and the traceplot and density plot,
#have the models in part d and f converged adequately? 

gg1<-ggs(as.mcmc(meanResult_d))  #Make the MCMC object a ggs object for plotting
head(gg1)
table(gg1$Parameter)
summary(gg1)
gg1<-filter(gg1,Parameter=="m")
ggs_traceplot(gg1)  # Trace plot
ggs_density(gg1)  
#Summary statistics
round(meanResult_d$BUGSoutput$summary  ,3)


gg2<-ggs(as.mcmc(meanResult_f))  #Make the MCMC object a ggs object for plotting
head(gg1)
table(gg1$Parameter)
summary(gg1)
gg1<-filter(gg1,Parameter=="m")
ggs_traceplot(gg1)  # Trace plot
ggs_density(gg1)  

round(meanResult_f$BUGSoutput$summary  ,3)
#------------------------
###2: Binomial and negative binomial###
#-------------------------

#Researchers recorded the number of mortalities and survivals in a group of 
#animals that were diagnoses with a virus.  There were 12 survivors and 5 
#mortalities out of 17 animals. We want to estimate the mortality rate. 

#A) Write out how you would set up a binomial model for this analysis, including 
#the data, parameter, likelihood, and prior.  

write("model{
  x ~ dbin(r, 17)	# data sampled binomially with n=12
  r ~ dunif(0, 1)	# uninformative prior for the sex ratio of pouch young
}
",file=here("JAGS_mods", "HW2_2a_binomial.txt"))


#B) b. Fit this model in JAGS using 100000 iterations and give the posterior 
#statistics. Judging by Rhat and n.eff in the summary, and the traceplot and 
#density plot, has it converged?  What is the mortality rate and its 95% credible
#interval? 

mortality_b<-jags(data=list(x=5),
                    parameters.to.save="r",
                    model.file = here("JAGS_mods", "HW2_2a_binomial.txt"),
                    n.iter=100000,
                    n.burnin=10000,
                    n.thin=1)

print(mortality_b)
ggsBinom<-ggs(as.mcmc(mortality_b)) %>% filter(Parameter!="deviance")
ggs_traceplot(ggsBinom)
ggs_density(ggsBinom)

#C) Write out how you would set up a negative binomial model for this analysis, 
#including the data, parameter, likelihood, and prior.  

write("model{
  x ~ dnegbin(r, 5)	# data sampled binomially with n=12
  r ~ dunif(0, 1)	# uninformative prior for the sex ratio of pouch young
}
",file=here("JAGS_mods", "HW2_2c_negbinomial.txt"))

#D) Fit this model in JAGS using 100000 iterations and give the posterior 
#statistics. Judging by Rhat and n.eff in the summary, and the traceplot and 
#density plot, has it converged?  What is the mortality rate and its 95% credible 
#interval? 

mortality_d<-jags(data=list(x=12),
                  parameters.to.save="r",
                  model.file = here("JAGS_mods", "HW2_2c_negbinomial.txt"),
                  n.iter=100000,
                  n.burnin=10000,
                  n.thin=1)

print(mortality_d)
ggsNegBinom<-ggs(as.mcmc(mortality_d)) %>% filter(Parameter!="deviance")
ggs_traceplot(ggsBinom)
ggs_density(ggsBinom)


#E) How similar are your results for the binomial and negative binomial? 


