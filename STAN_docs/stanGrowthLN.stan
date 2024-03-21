data {
  int<lower=0> N;
  vector[N] ObsL;
  vector[N] Age;
}
parameters {
  real<lower=10, upper=2000> Linf;
  real<lower=0, upper=4> K;
  real<lower=-15,upper=0> Age0;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] PredL;
  vector[N] logPredL;
  for(i in 1:N) {
		PredL[i] = Linf * (1 - exp(- K* (Age[i] -Age0))); 
		logPredL[i] = log(PredL[i]);  
  }
}
model {
  K ~ uniform(0,4);
  Age0~ uniform(-15,0);
  Linf ~ uniform(10,3000); 	
  sigma ~ uniform(0.01,1);
  for(i in 1:N) {
    ObsL[i] ~ lognormal(logPredL[i], sigma);
  }
}
generated quantities {
  vector[N] LL;
  vector[N] residual;
  vector[N] repL;
  for(i in 1:N) {
  		LL[i] = lognormal_lpdf(ObsL[i]|logPredL[i],sigma);
  		residual[i] = log(ObsL[i])-logPredL[i];
      repL[i] = lognormal_rng(logPredL[i],sigma);
  }
}

