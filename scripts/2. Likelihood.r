# Lecture 2. Likelihood
library(tidyverse)
theme_set(theme_bw(base_size=15))

#Sex ratio example

# Grid search
sexRatio<-data.frame(p=seq(0,1,by=0.01))
sexRatio
sexRatio$Likelihood<-dbinom(15,20,sexRatio$p)
ggplot(sexRatio,aes(x=p,y=Likelihood))+geom_line()+
  ggtitle("Likelihood")
ggplot(sexRatio,aes(x=p,y=log(Likelihood)))+geom_line()+
  ggtitle("Log likelihood")
sexRatio$p[sexRatio$Likelihood==max(sexRatio$Likelihood)]   #Maximum likelihood parameter
sexRatio$LL<-log(sexRatio$Likelihood)
sexRatio$p[sexRatio$LL==max(sexRatio$LL)]   #Maximum likelihood parameter

# Minimization
?optim
getNegLogLike<-function(p) {
  -dbinom(15,20,p,log=TRUE)
}

sexRatioModel<-optim(par=0.5,
                      fn=getNegLogLike,
                     method="Brent",
                     lower=0,
                     upper=1
                     )
sexRatioModel

## Example of a mean, with normal data and known variance

data<-data.frame(x=c(6.4, 3.3, 3.7, 2.2, 6.2, 6.3, 1.3, 5.1, 4.9, 0.9),
                 index=1:10)
ggplot(data,aes(x=index,y=x))+geom_point()+
  geom_hline(yintercept=5)

# grid search
meanVal<-data.frame(mean=seq(2,7,0.01))
for(i in 1:nrow(meanVal)) {
  meanVal$LL[i]<-sum(dnorm(data$x,meanVal$mean[i],2,log=TRUE))
}

#Note that R is summing over all the data points in one step. 
# We could have done a nested loope


ggplot(meanVal,aes(x=mean,y=LL))+geom_line()
meanVal$mean[meanVal$LL==max(meanVal$LL)]

mean(data$x)

#minimization

getNegLogLikeNorm<-function(mu,data) {
  x=data$x
  LL<-dnorm(x,mu,2,log=TRUE)
  -sum(LL)
}
meanModel<-optim(par=10,
                  fn=getNegLogLikeNorm,
                  data=data,
                 method="Brent",
                 lower=-100,upper=100)
meanModel

## Mean and standard deviation

# grid search 
meanSdVal<-expand.grid(mean=seq(1,6,length=200),sd=seq(0.1,6,length=200))
for(i in 1:nrow(meanSdVal)) {
  meanSdVal$LL[i]<-sum(dnorm(data$x,meanSdVal$mean[i],meanSdVal$sd[i],log=TRUE))
}
head(meanSdVal)
meanSdVal[meanSdVal$LL==max(meanSdVal$LL),]
bestmean<-meanSdVal$mean[meanSdVal$LL==max(meanSdVal$LL)]
bestsd<-meanSdVal$sd[meanSdVal$LL==max(meanSdVal$LL)]

#joint likelihood
ggplot(mutate(meanSdVal,LL=ifelse(LL>-30,LL,-30)),aes(x=mean,y=sd,fill=LL))+
  geom_tile()+
  scale_fill_viridis_c()

#for mean and best sd
ggplot(filter(meanSdVal,sd==bestsd),aes(x=mean,y=LL))+
  geom_line()+
  ylim(-30,-20)

# for sd at best mean
ggplot(filter(meanSdVal,mean==bestmean),aes(x=sd,y=LL))+
  geom_line()+
  ylim(-30,-20)

#in likelihoods not Log likelihoods
ggplot(filter(meanSdVal,mean==bestmean),aes(x=sd,y=exp(LL)))+
  geom_line()

# Minimizer function
getNegLogLikeNorm2<-function(pars,data) {
  x=data$x
  mu=pars[1]
  sd=pars[2]
  LL<-dnorm(x,mu,sd,log=TRUE)
  -sum(LL)
}
meanSdModel<-optim(par=c(2,2),
                  fn=getNegLogLikeNorm2,
                  data=data)
meanSdModel$par

#analytic values
mean(data$x)
sd(data$x)
sqrt(sum((data$x-mean(data$x))^2/10))  #population sd

# Adding Hessian

meanSdModel<-optim(par=c(2,2),
                     fn=getNegLogLikeNorm2,
                     data=data,
                     hessian=TRUE)
meanSdModel

#This is the inversion to get the variance covariance model
solve(meanSdModel$hessian)

# Confidence intervals

round(meanSdModel$hessian,4)
varCov<-solve(meanSdModel$hessian)
round(varCov,6)

meanCI<-c(meanSdModel$par[1]-2*sqrt(varCov[1,1]),
          meanSdModel$par[1]+2*sqrt(varCov[1,1]))
sdCI<-c(meanSdModel$par[2]-2*sqrt(varCov[2,2]),
          meanSdModel$par[2]+2*sqrt(varCov[2,2]))

#for mean and best sd
ggplot(filter(meanSdVal,sd==bestsd),aes(x=mean,y=LL))+
  geom_line() +
  geom_line(data=data.frame(CI=meanCI,y=rep(-23,2)),
               aes(x=CI,y=y),color="red",lwd=2)

#for  sd and best mean
ggplot(filter(meanSdVal,mean==bestmean,LL>-30),aes(x=sd,y=LL))+
  geom_line()+
  geom_line(data=data.frame(CI=sdCI,y=rep(-23,2)),
            aes(x=CI,y=y),color="red",lwd=2)

