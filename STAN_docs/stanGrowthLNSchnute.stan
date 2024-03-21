
data {
  int<lower=0> N;
  vector[N] ObsL;
  vector[N] Age;
  real  Age1;
  real  Age2;
}
parameters {
  real<lower=10, upper=2000> L1;
  real<lower=10, upper=3000> L2;
  real<lower=0, upper=4> K;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] PredL;
  vector[N] logPredL;
  for(i in 1:N) {			
    PredL[i] = L1 +(L2-L1)* (1 - exp(- K* (Age[i] -Age1)))/(1 - exp(- K* (Age2 -Age1))); 
		logPredL[i] = log(PredL[i]);  
  }
}
model {
  K ~ uniform(0,4);
  L1 ~ uniform(10,3000); 	
  L2 ~ uniform(10,2000); 	
  sigma ~ uniform(0.01,1);
  for(i in 1:N) {
    ObsL[i] ~ lognormal(logPredL[i], sigma);
  }
}
generated quantities {
  vector[N] LL;
  vector[N] residual;
  vector[N] repL;
  real      Linf;
  Linf= (L2-L1*exp(-K*(Age2-Age1)))/(1-exp(-K*(Age2-Age1)));
  for(i in 1:N) {
  		LL[i] = lognormal_lpdf(ObsL[i]|logPredL[i],sigma);
  		residual[i] = log(ObsL[i])-logPredL[i];
      repL[i] = lognormal_rng(logPredL[i],sigma);
  }
}

