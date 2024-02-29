
data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
}
parameters {
  real<lower=0> errorSigma;
  real alpha;
  real beta;
}
model {
  errorSigma ~ exponential(0.01);
  alpha ~ normal(0,100);
  beta~normal(0,100);
  for(i in 1:N) {
   y[i] ~ normal(alpha+x[i]*beta, errorSigma);
  }
}
generated quantities {
  real LL[N];
  for(i in 1:N) {
   LL[i] = normal_lpdf(y[i] | alpha+x[i]*beta, errorSigma);
  }
}


