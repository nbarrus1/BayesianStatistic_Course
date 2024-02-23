// The input data 
data {
  int<lower=0> N;
  int<lower=0> Y[N];
}

// The parameters accepted by the model. lamda
parameters {
  real<lower=0> r;
  real<lower=0, upper=1> p;
}

transformed parameters{
  real<lower=0> m;
  m = r*(1-p)/p;
}

// The model to be estimated. We model the output
// 'y' to be negative binomial distributed with propotion p
// and number of success r
model {
  p ~ uniform(0,1);
  r ~ lognormal(0.0,1000000);
  Y ~ neg_binomial_2(m,r);
}

generated quantities{
  real v;
  v=m+m*m/r;
  real LL[N];
  for(i in 1:N) {
    LL[i] = neg_binomial_2_lpmf(Y[i] | m, r);
  }
}