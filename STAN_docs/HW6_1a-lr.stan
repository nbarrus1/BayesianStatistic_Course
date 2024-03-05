
data {
  int<lower=0> N;
  int<lower=0> y[N];
}

// The parameters accepted by the model. 
parameters {
  real a; //intercept term
  real b; //slope term
}

//model to be estimated
model {
  a ~ normal(0,1000000);
  b ~ normal(0,1000000);
  for(i in 1:N) {
  y[i] ~ bernoulli_logit(a + b * x[i]);
  }
}

//values we wante to generate

generated quantities {
  real lr[N];
  real resid[N];
  real ypred[N];
  real loglik[N];
  
  for (i in 1:N) {
    ypred = bernoulli_logit_rng(a + b * x[i]);
    loglik[i] = bernoulli_logit_lpdf(y[i]|(a+b*x[i]));
  }
  
}