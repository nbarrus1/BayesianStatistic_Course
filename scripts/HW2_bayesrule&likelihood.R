#Script Title: Bayesian Statistics Homework2-Bayes Rule & Likelihoods
#Date first created: 1/18/2024
#purpose: The purpose of this script is to work through Homework 2 which 
#which contains our first problem set. The problems are related to lectures
#1 and 2 which are over bayes rule and likelihoods, respectively. 

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
theme_set(theme_bw(base_size=15))

#------------------------
####Homework 2: Code####
#---------------------------

####problem 1: Bayes Rule ####
#problem one will be done analytically and does not have code.

####problem 2: Likelihood using the binomial distribution####

#data

females<-c(8, 10, 8, 8, 7, 8, 9, 10, 6, 7, 7, 4, 8, 6, 8, 
           7, 8, 10, 6, 7)

###Use a grid search to find the maximum likelihood estimate (MLE)
#of the propotion female.

sex.ratio <- tibble(p = seq(0,1, by = 0.01)) |>  #create sequence of p
  mutate(LL = 0)                                 #set up column to save loop results

##for loop to obtain LL across the sequences of ps

for(i in 1:length(sex.ratio$p)) {
 sex.ratio$LL[i] <- sum(dbinom(females,20,sex.ratio$p[i], log = TRUE))
}

#now find the maximum log likelihood
sex.ratio$p [sex.ratio$LL == max(sex.ratio$LL)]

###Make a plot of the log-likelihood profile for propotion female. 

##loglikelhood plot

sex.ratio |> 
  ggplot(aes(x = p, y = LL)) +
  geom_line()
  

##likelihood plot

sex.ratio |> 
  ggplot(aes(x = p, y = exp(LL))) +
  geom_line()

###Use the optim() fuction to get the MLE fo the propotion female.
###is it the same value you got in part a?

getNegLogLike.binom <- function(p,nfemale){
  x = nfemale
  LL <- dbinom(x,20,p,log = TRUE)
  -sum(LL)
}


sex.ratio.model <- optim(par = 0.5,
                          fn = getNegLogLike.binom,
                          nfemale = females,
                          method = "Brent",
                          lower = 0,
                          upper = 1
                          )

sex.ratio.model

#do the two methods give the same result
sex.ratio.model$par                                
sex.ratio$p [sex.ratio$LL == max(sex.ratio$LL)]
#yes!!!

####problem 3: Likelihood with normal mean and variance####

#data

lengthDat<-c(26.4, 16.1, 26.8, 26.5, 29, 22.4, 21.4, 25.2, 26.1, 16.3, 24.2, 
             20.9, 11.5, 22.2, 22.1, 19.1, 24.4, 24.1, 21.5, 20.2)

sd(lengthDat)
mean(lengthDat)

##### Use a grid search to find the maximum likelihood estimate (MLE) of
####the mean and standard deviation. What are the values?

meanSDVal <-as_tibble(expand.grid(mean = seq(min(lengthDat),max(lengthDat), length = 200),
                                  sd = seq(0.1,6, length = 200))) |> 
  mutate(LL = 0)

#do for loop to get LL 

for (i in 1:length(meanSDVal$mean)){
  meanSDVal$LL[i] <- sum(dnorm(lengthDat,meanSDVal$mean[i],meanSDVal$sd[i],log=TRUE))
}

head(meanSDVal)

meanSDVal[meanSDVal$LL==max(meanSDVal$LL),]

####Use the optim() function to get the MLE of the mean and sd. Is it the
####same value you got in part a?

getNegLogLikeNorm2<-function(pars,data) {
  x=data
  mu=pars[1]
  sd=pars[2]
  LL<-dnorm(x,mu,sd,log=TRUE)
  -sum(LL)
}

meanSdModel<-optim(par=c(2,2),
                   fn=getNegLogLikeNorm2,
                   data=lengthDat,
                   hessian = TRUE)

meanSdModel

#do the two methods give the same result
meanSdModel$par
meanSDVal[meanSDVal$LL==max(meanSDVal$LL),]

#yes!! sd wasn't quite the same but it was really close

#####Calculate the variance/covariance matrix of mean and sd. Are mean 
#####and sd correlated? 

varCov <- solve(meanSdModel$hessian)
varCov

#no they are not related because the covariance is small

#####For the MLE of both the mean and the sd, what is the pointwise 
#####log-likelihood of each of the 20 data points?  Plot the pointwise LL 
#####against the length data. What pattern do you see? 

MLE_mean = meanSdModel$par[1]
MLE_sd = meanSdModel$par[2]
point.LL <- dnorm(x = lengthDat, mean = MLE_mean, sd = MLE_sd, log = TRUE)
point.LL


tibble (x = lengthDat, y = exp(point.LL)) |> 
ggplot(aes(x = x, y = y)) +
  geom_point()+
  geom_line()+
  labs(x = "Lengths", y = "Point Log-likelihoods")
