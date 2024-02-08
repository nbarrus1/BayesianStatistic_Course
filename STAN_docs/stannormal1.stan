data {
  vector[10] Y;  
}
parameters {
  real mu;
  real<lower=0> sigma;
}
model {
  mu ~ normal(100,1000); //prior for mean
  sigma ~ exponential(0.1);  //prior for standard deviation
  Y ~ normal(mu, sigma); //likelihood vectorized
}

