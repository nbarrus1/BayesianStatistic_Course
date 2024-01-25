library(tidyverse)
library(R2jags)
library(ggmcmc)
library(here)
theme_set(theme_bw())

# Sea turtle binomial null distribution
df1<-data.frame(Females=0:20,Probability=dbinom(0:20,20,0.5))
df1$`In P` = ifelse(df1$Females %in% 6:14, "no","yes")
ggplot(df1,aes(x=Females,y=Probability))+geom_col()
ggplot(df1,aes(x=Females,y=Probability,fill=`In P`))+geom_col()+
  scale_fill_manual(values=c("darkgrey","red"))
sum(df1$Probability[df1$`In P`=="yes"])
binom.test(15,20)
#likelihood profile
LL<-data.frame(p=seq(0.3,0.99,0.01)) %>%
  mutate(Likelihood=dbinom(15,20,p)) %>%
  mutate(LL=log(Likelihood))
ggplot(LL,aes(x=p,y=Likelihood))+geom_line(linewidth=2)
ggplot(LL,aes(x=p,y=LL))+geom_line(linewidth=2)

#P values and likelihood for koalas

#P value for binomial is the probability of getting 3 or fewer males out of 12
sum(dbinom(c(0:3),12,0.5))
#or
pbinom(3,12,0.5)

#P value for negative binomial is the probability of getting 9 or more females before the 3rd male
sum(dnbinom(9:100000,3,0.5))

#Likelihood for binomial
dfbin<-data.frame(p=seq(0.01,.99,by=0.001))
dfbin$Likelihood<-dbinom(3,12,dfbin$p)
dfbin$LL<-log(dfbin$Likelihood)
ggplot(dfbin,aes(x=p,y=Likelihood))+geom_line(linewidth=2)+ggtitle("Binomial")

#Likelihood for negative binomial
dfnb<-data.frame(p=seq(0.01,.99,by=0.001))
dfnb$Likelihood<-dnbinom(9,3,dfnb$p)
dfnb$LL<-log(dfnb$Likelihood)
ggplot(dfnb,aes(x=p,y=Likelihood))+geom_line(linewidth=2)+ggtitle("Negative Binomial")

### Box 2.2 Bayesian binomial and negative binomial

write("model{
  x ~ dbin(r, 12)	# data sampled binomially with n=12
  r ~ dunif(0, 1)	# uninformative prior for the sex ratio of pouch young
}
",file=here("JAGS_mods", "box2.2Binomial.txt"))

KoalaBinomial<-jags(data=list(x=3),
  parameters.to.save="r",
  model.file = here("JAGS_mods", "box2.2Binomial.txt"),
  n.iter=100000,
  n.burnin=10000,
  n.thin=1)

print(KoalaBinomial)
ggsBinom<-ggs(as.mcmc(KoalaBinomial)) %>% filter(Parameter!="deviance")
ggs_density(ggsBinom)

write("model{
  x ~ dnegbin(r, 3)	# number of females sampled neg. binomially
  r ~ dunif(0, 1)		# uninformative prior for the sex ratio of pouch young
}
",file=here("JAGS_mods","box2.2NegBin.txt"))

KoalaNegBin<-jags(data=list(x=9),
  parameters.to.save="r",
  model.file = here("JAGS_mods","box2.2NegBin.txt"),
  n.iter=100000,
  n.burnin=10000,
  n.thin=1)

print(KoalaNegBin)
ggsNegBin<-ggs(as.mcmc(KoalaNegBin)) %>% filter(Parameter!="deviance")
ggs_density(ggsNegBin)


