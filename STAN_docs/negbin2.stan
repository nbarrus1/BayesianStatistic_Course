data {
  int<lower=0> N;
  int<lower=0> y[N];
}
parameters {
  real<lower=0,upper=100> r;
  real<lower=0> m;
}
model{
  r ~lognormal(0.0, 10);
  m ~lognormal(0.0, 10);  
  y ~ neg_binomial_2(m,r);			
}
generated quantities{
  real p;
  real v;
  real dispersion;
  real stepProb;
  p=r/(r+m);
  v=m+m*m/r;
  dispersion=1/p;
  stepProb=if_else(dispersion>1,1,0);
}

