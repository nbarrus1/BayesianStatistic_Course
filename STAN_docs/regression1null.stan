data {
  int N;
  vector[N] y;
}
parameters {
  real a;
  real<lower=0> sigma;
}
model {
a ~ normal(0,100);
sigma ~ exponential(0.01);
for(i in 1:N) {
 y[i] ~ normal(a , sigma);
}
}
generated quantities {
  real loglik[N];
  for (i in 1:N)  {
    loglik[i] = normal_lpdf(y[i]|a,sigma);
  }
}
