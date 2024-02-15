data {
  int N;
  vector[N] y;
  vector[N] x;
}
parameters {
  real a;
  real b;
  real<lower=0> sigma;
}
model {
a ~ normal(0,100);
b ~ normal(0,100);
sigma ~ exponential(0.01);
for(i in 1:N) {
 y[i] ~ normal(a + b * x[i], sigma);
}
}
