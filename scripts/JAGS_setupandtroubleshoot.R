#Script Title: Bayesian Statistics Course Setup
#Date first created: 1/16/2024
#purpose: The purpose of this script is to run some software set-up at the beginning
#of the the Bayesian Statistics Course, and to see if the software was installed 
#correctly


#----------------------------
##general start to scripting##
#----------------------------



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

library(R2jags)
library(tidyverse)
library(ggmcmc)
library(here)

#-----------------------------
###check the software###
#-----------------------------

#test JAGS model 
write("model {
 mu~dnorm(0,0.01)
 for(i in 1:10) {
    x[i]~dnorm(mu,1)
 }
}",
  file=here("JAGS_mods","temp.txt"))

#run jags and save to list called model1

model1 <- jags(data = list(x=rnorm(10)),
               model.file = here("JAGS_mods","temp.txt"),
               parameters.to.save = "mu")

#look at the model summary

model1$BUGSoutput$summary

#make a figure

gg1 <- ggs(as.mcmc(model1))
ggs_density(gg1)
