
data {
  int<lower=0> N;
  int<lower=0> Nsubject;
  vector[N] y;
  int subject[N];
  vector[N] x;
}
parameters {
  real<lower=0> interceptSigma;
  real<lower=0> slopeSigma;
  real<lower=0> errorSigma;
  real interceptMu;
  real slopeMu;
  real alpha[Nsubject];
  real beta[Nsubject];
}
model {
  errorSigma ~ exponential(0.01);
  interceptSigma ~ exponential(0.01);
  slopeSigma ~ exponential(0.01);
  interceptMu ~ normal(0,100);
  slopeMu ~ normal(0,100);
  alpha~normal(interceptMu,interceptSigma);
  beta~normal(slopeMu,slopeSigma);
  for(i in 1:N) {
   y[i] ~ normal(alpha[subject[i]]+x[i]*beta[subject[i]], errorSigma);
  }
}
generated quantities {
  real LL[N];
  for(i in 1:N) {
   LL[i] = normal_lpdf(y[i] | alpha[subject[i]]+x[i]*beta[subject[i]], errorSigma);
  }
}

  
