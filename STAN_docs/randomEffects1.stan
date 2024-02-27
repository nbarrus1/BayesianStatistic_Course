
data {
  int<lower=0> N;
  int<lower=0> Nsubject;
  vector[N] y;
  int subject[N];
  vector[N] x;
}
parameters {
  real interceptSigma;
  real slopeSigma;
  real<lower=0> errorSigma;
  real alpha[Nsubject];
  real beta[Nsubject];
}
model {
  errorSigma ~ lognormal(0,10);
  interceptSigma ~ lognormal(0,10);
  slopeSigma ~ lognormal(0,10);
  alpha~normal(0,interceptSigma);
  beta~normal(0,slopeSigma);
  for(i in 1:N) {
   y ~ normal(alpha[subject[i]]+x[i]*beta[subject[i]], errorSigma);
  }
}

  
