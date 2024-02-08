data {
  int<lower=0> N;
  int<lower=0> y[N];
}
parameters {
  real<lower=0,upper=1> p;
  real<lower=0> r;
}
transformed parameters {
  real<lower=0> m;
  m = r *(1-p)/p;
}
model{
  p~uniform(0,1);
  r~gamma(0.01,0.01);
  y ~ neg_binomial_2(m,r);			
}
generated quantities{
  real v;
  v=m+m*m/r;
}

