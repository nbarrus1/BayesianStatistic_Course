// The input data 
data {
  int<lower=0> N;
  int<lower=0> Y[N];
}

// The parameters accepted by the model. lamda
parameters {
  real<lower=0, upper=10> r;
  real<lower=0> m;
}


// The model to be estimated. We model the output
// 'y' to be negative binomial distributed with propotion p
// and number of success r
model {
  m ~ lognormal(0.0,1000);
  r ~ lognormal(0.0,1000);
  Y ~ neg_binomial_2(m,r);
}

generated quantities{
  real v;
  real p;
  v=m+m*m/r;
  p = r/(r+m);
}