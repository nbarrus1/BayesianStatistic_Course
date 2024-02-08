#Script Title: Bayesian Statistics Homework2-Normal, bionomial, neg-binomial
#Date first created: 2/6/2024
#purpose: The purpose of this script is to work through Homework 3 which 
#contains problems with normal. Specifically,
#looking at setting up different priors, and models the problems related to lectures 5 & 6?

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

#-----------------------------
####Problem 1###
#-----------------------------

#problem
#The following data are measured shark lengths, in mm, assumed to be from a
#Normal distribution:
#  y<- c(793, 1105, 1268, 678, 1363, 981, 990, 879, 739, 876, 1393, 
#        774, 1037, 1183, 1101, 744, 1407, 683, 999, 940, 1431, 1270, 
#        969, 917, 1547, 1193, 965, 1068, 1042, 1085, 970, 776, 872, 1302, 
#        625, 780, 1297, 894, 1273, 1140, 1349, 709, 657, 714, 780, 874, 
#        885, 803, 838, 841)



#load data 

shark.length <- list(Y = c(793, 1105, 1268, 678, 1363, 981, 990, 879, 739, 876, 1393, 
                           774, 1037, 1183, 1101, 744, 1407, 683, 999, 940, 1431, 1270, 
                           969, 917, 1547, 1193, 965, 1068, 1042, 1085, 970, 776, 872, 1302, 
                           625, 780, 1297, 894, 1273, 1140, 1349, 709, 657, 714, 780, 874, 
                           885, 803, 838, 841))


#A) Using an uninformative Normal prior for the mean, and an uninformative gamma
#prior for the precision, use JAGS to calculate the mean, variance, standard 
#deviation and precision.  Show the summary statistics and posterior marginal 
#density and verify that the model is converged by Rhat and n.eff. Which #
#parameters are estimated vs. derived?

### write the model

write("model
{
  mean ~ dnorm(100, 1E-6)  #uninformative prior for mean	

  #uniformative prior for variance
  prec ~ dgamma(0.01,0.001)   
  Var <- 1/prec
  sd <- sqrt(Var)
  sd.diff <- sd-250
  sd.diff.test <- step(sd.diff)

  #likelihood
  for (i in 1:10)  		 	
  {
    Y[i] ~ dnorm(mean, prec)
  }
}
", file=here("JAGS_mods","HW3_sharklength_sd-gamma.txt"))

####run the model

#If you initialize, must be on estimated parameter
model_HW3_sharklength_gam<-jags(data=shark.length,
                      inits=list(list(prec=1),list(prec=2)),
                      parameters.to.save = c("mean","Var","prec","sd","sd.diff.test"),
                      n.chains = 2,
                      n.burnin = 10000,
                      n.iter = 200000,
                      n.thin = 1,
                      model.file = here("JAGS_mods","HW3_sharklength_sd-gamma.txt"))

print(model_HW3_sharklength_gam)
ggs_density(ggs(as.mcmc(model_HW3_sharklength_gam)))


#B) To make sure you understand what the model is estimating, fill in the 
#following numbers from the summary statistics:


###posterior standard deviation of the mean

#79.570

###posterior mean of the standard deviation

#243.896

###posterior standard deviation of the standard deviation

#66.064

###Which of these is most analogous to the standard error of the mean? 

###the posterior standard deviation of the mean is most analogous to the 
#the standard error of the mean

#C) Repeat part a using an uninformative uniform prior for the standard deviation
#rather than a gamma prior on precision.  Show the posterior marginal densities 
#of each parameter, as well as the summary statistics.  Did you use a wide 
#enough range in the standard deviation? If not, use a wider range
 
#write the model

write("model
{
  mean ~ dnorm(100, 1E-6)  #prior for mean
  

  #Prior on sd
  sd~dunif(0,400)   
  Var <- sd *sd
  prec <- 1/Var
  sd.diff <- sd-250
  sd.diff.test <- step(sd.diff)
  
  #likelihood
  for (i in 1:10)  		 	
  {
    Y[i] ~ dnorm(mean, prec)
  }
}
", file=here("JAGS_mods","HW3_sharklength_sd-uniform.txt"))

#If you initialize, must be on estimated parameter
model_HW3_sharklength_unif<-jags(data=shark.length,
                                 inits=list(list(sd=2),list(sd=1)),  #for sd as the estimated parameter
                                 parameters.to.save = c("mean","Var","prec","sd","sd.diff.test"),
                                 n.chains = 2,
                                 n.burnin = 10000,
                                 n.iter = 200000,
                                 n.thin = 1,
                                 model.file = here("JAGS_mods","HW3_sharklength_sd-uniform.txt"))

print(model_HW3_sharklength_unif)
ggs_density(ggs(as.mcmc(model_HW3_sharklength_unif)))

#D)  Repeat part a using an exponential prior for the standard deviation.  Show 
#the posterior marginal densities of each parameter, as well as the summary 
#statistics. Make sure that the lambda parameter is consistent with the range of
#the data (where lambda=1/mean).

#write the model

write("model
{
  mean ~ dnorm(100, 1E-6)  #prior for mean	

  #Prior on sd
  sd~dexp(0.1) 
  Var <- sd *sd
  prec <- 1/Var
  sd.diff <- sd-250
  sd.diff.test <- step(sd.diff)

  
  #likelihood
  for (i in 1:10)  		 	
  {
    Y[i] ~ dnorm(mean, prec)
  }
}
", file=here("JAGS_mods","HW3_sharklength_sd-exp.txt"))

#If you initialize, must be on estimated parameter
model_HW3_sharklength_exp<-jags(data=shark.length,
                                 inits=list(list(sd=2),list(sd=1)),  #for sd as the estimated parameter
                                 parameters.to.save = c("mean","Var","prec","sd","sd.diff.test"),
                                 n.chains = 2,
                                 n.burnin = 10000,
                                 n.iter = 200000,
                                 n.thin = 1,
                                 model.file = here("JAGS_mods","HW3_sharklength_sd-exp.txt"))

print(model_HW3_sharklength_exp)
ggs_density(ggs(as.mcmc(model_HW3_sharklength_unif)))


#E)How do the results differ between the models with different priors in parts 
#a and b? Did the choice of prior parameterization matter?

round(model_HW3_sharklength_gam$BUGSoutput$summary,4)
round(model_HW3_sharklength_unif$BUGSoutput$summary,4)
round(model_HW3_sharklength_exp$BUGSoutput$summary,4)

#F) For both all three models, use the step() function to calculate the 
#probability that the standard deviation is larger than 250.  (The researchers 
#want to know if their sharks are more variable in size than those from a 
#related population). Do you get similar answers from each parameterization?

round(model_HW3_sharklength_gam$BUGSoutput$summary[5:6,],4)
round(model_HW3_sharklength_unif$BUGSoutput$summary[5:6,],4)
round(model_HW3_sharklength_exp$BUGSoutput$summary[5:6,],4)


#-------------------
###Problem 2####
#-------------------

#data problem

#We have prior information from a previous study that the fraction of streams 
#where kingfishers are present is 0.25, with a CV of 0.2. (CV, or coefficient 
#of variation, is the standard deviation divided by the mean).  Our survey 
#observed kingfishers at 6 streams out of the 19 that we visited. 
####Assume presence or absence at a stream is a binomial process. 

#A) Using an uninformative beta prior of dbeta(1,1) estimate the fraction of 
#streams with kingfishers. Show the summary statistics and posterior density 
#of the fraction. 

write("model{
  r ~ dbeta(1, 1)	# uninformative prior for propotion present
  r.prior ~ dbeta(1, 1)
  x ~ dbin(r, n)	# data sampled binomially with total samples = 19, success = 6
  x.prior.sim ~ dbin(r.prior,n)
  x.sim ~ dbin(r,n)
  }
",file= here("JAGS_mods","HW3_kingfish_beta-unifor.txt"))

kf.dat.uninfor <- list(x = 6, n = 19)

presMod_HW2_2a <- jags(data=kf.dat.uninfor,
     parameters.to.save=c("r","r.prior","x.prior.sim","x.sim"),
     model.file = here("JAGS_mods","HW3_kingfish_beta-unifor.txt"),
     n.iter=100000,
     n.burnin=10000,
     n.thin=1)

print(presMod_HW2_2a)
ggsBinom_2a<-ggs(as.mcmc(presMod_HW2_2a)) %>% filter(Parameter=="r")
ggs_traceplot(ggsBinom_2a)
ggs_density(ggsBinom_2a)

#B) Using an informative beta prior based on the previous study, estimate the 
#fraction of streams with kingfishers. What are your values for a and b? Show 
#the summary statistics and posterior of the fraction. 

mu = 0.25
cv = 0.2

#cv = sd/mu
#sd = CV*mu
sd = cv*mu

afunc<-function(mu,sig) mu*(mu*(1-mu)/sig^2-1)
bfunc<-function(mu,sig) (1-mu)*(mu*(1-mu)/sig^2-1)

a <- afunc(mu,sd)
b <- bfunc(mu,sd)


write("model{
  r ~ dbeta(a, b)	# informative prior for propotion present
  r.prior ~ dbeta(a, b)
  x ~ dbin(r, n)	# data sampled binomially with total samples = 19, success = 6
  x.prior.sim ~ dbin(r.prior,n)
  x.sim~dbin(r,n)
  }
",file= here("JAGS_mods","HW3_kingfish_beta-infor.txt"))

kf.dat.infor <- list(x = 6,n = 19,a = a, b = b)

presMod_HW2_2b <- jags(data=kf.dat.infor,
                       parameters.to.save=c("r","r.prior","x.prior.sim","x.sim"),
                       model.file = here("JAGS_mods","HW3_kingfish_beta-infor.txt"),
                       n.iter=100000,
                       n.burnin=10000,
                       n.thin=1)

print(presMod_HW2_2b)
ggsBinom_2b<-ggs(as.mcmc(presMod_HW2_2b)) %>% filter(Parameter=="r")
ggs_traceplot(ggsBinom_2b)
ggs_density(ggsBinom_2b)


#C) How do your results differ between a and b?

print(presMod_HW2_2a)
print(presMod_HW2_2b)

#-----------------------------
###Problem 3####
#----------------------------

##Using the two models from question 2, do the following checks:

#A) Add lines of code in the models to calculate the prior for the fraction 
#parameter and plot both the prior and the posterior using ggs_density.  How 
#did the data update the priors?

ggsBinom_2a<-ggs(as.mcmc(presMod_HW2_2a)) %>% filter(Parameter=="r"|Parameter == "r.prior")
ggs_density(ggsBinom_2a)


ggsBinom_2b<-ggs(as.mcmc(presMod_HW2_2b)) %>% filter(Parameter=="r"|Parameter == "r.prior")
ggs_density(ggsBinom_2b)


#B) Now do a post model pre data run, by replacing the values of the number of 
#streams with kingfishers with NA. and re-running the model. Look at the 
#posterior summary. How do the results differ from what you got in part b and 
#from each other?

kf.pred.uninfor <- list(x = NA, n = 19)



presMod_HW2_2a_pred <- jags(data=kf.pred.uninfor,
                       parameters.to.save=c("r","r.prior","x", "x.sim","x.prior.sim"),
                       model.file = here("JAGS_mods","HW3_kingfish_beta-unifor.txt"),
                       n.iter=100000,
                       n.burnin=10000,
                       n.thin=1,DIC = FALSE)

print(presMod_HW2_2a_pred)

kf.pred.infor <- list(x = NA,n = 19,a = a, b = b)



presMod_HW2_2b_pred <- jags(data=kf.pred.infor,inits =list(init1,init2,init3),
                       parameters.to.save=c("r","r.prior","x","x.sim","x.prior.sim"),
                       model.file = here("JAGS_mods","HW3_kingfish_beta-infor.txt"),
                       n.iter=100000,
                       n.burnin=10000,
                       n.thin=1,DIC = FALSE)

print(presMod_HW2_2b_pred)


#C) Now look at the posterior and prior predictive distributions. Add lines to the 
#code making simulated data from both the posterior and the prior from the 
#model with data and plot the distributions of the simulated data under both 
#the prior and the posterior for both models. Does the binomial distribution 
#seem appropriate for this data?

sim.df <-bind_rows(list(`Uninformative` = data.frame(presMod_HW2_2a_pred$BUGSoutput$sims.matrix),
                        `Informative` = data.frame(presMod_HW2_2b_pred$BUGSoutput$sims.matrix)),
                        .id = "Prior")

ggplot(sim.df)+
  geom_bar(aes(x = x.sim, fill = "Posterior"), alpha = 0.5)+
  geom_bar(aes(x=x.prior.sim,fill="Prior"),alpha=0.5)+
  facet_wrap(Prior~.,ncol=1)+
  scale_fill_manual(values = c("Posterior" = "blue","Prior" = "red"))+
  ggtitle("Prior and posterior predictive distributions")
