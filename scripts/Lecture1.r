##Frog presence/absence example from McCarthy, chapter 1

#Calculations from Bayes' theorem
frog<-data.frame(Presence=c("Present","Absent"), Prob.detect=c(0.8,0), Prob.not.detect=c(0.2,1))
frog$prior<-c(0.5,0.5)
#frog$prior<-c(0.75,0.25)
frog
prob.present<-frog$prior[1]*frog$Prob.not.detect[1]/sum(frog$prior*frog$Prob.not.detect)
prob.present


#To run JAGS from R, you must download and install JAGS from https://sourceforge.net/projects/mcmc-jags/
#Then install R2jags in RStudio
#install.packages("R2jags")  #Install package if needed
library(R2jags)

#To look at help file for jags
#?jags

#Write model to a text file so that JAGS can access it.
write("model{
  prior <- 0.5	# the prior probability of presence
#  prior<-0.75
  present ~ dbern(prior)	# actual presence drawn from a Bernoulli dist'n
  prob_detect <- 0.8*present	# Pr(detection) depends on presence/absence
  detected ~ dbern(prob_detect)	# actual detection occurs with random variation
}
",file="Box1.5frog.txt")

#Run JAGS using jags function
frogJags<-jags(data=list(detected=0), #list with data inputs
                 parameters.to.save = "present",  #character vector of output parameters
                 model.file="Box1.5frog.txt" , #name of text file with model
                 n.iter=1000) #Number of MCMC iterations to run
#Look at output summary
print(frogJags)
frogJags$BUGSoutput$summary

