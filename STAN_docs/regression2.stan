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
generated quantities {
  real chi2dat;
  real chi2sim;
  real pvalue;
  real resid[N];
  real ypred[N];
  real pearson2dat[N];
  real pearson2sim[N];
  real loglik[N];
  
  
  for (i in 1:N)  {
    resid[i] = y[i]-(a + b * x[i]);
    ypred[i] = normal_rng(a + b *x[i], sigma); 
    pearson2dat[i] = (y[i]-(a + b * x[i]))/sigma^2;
    pearson2sim[i] = (ypred[i]-(a + b * x[i]))/sigma^2;
    loglik[i] = normal_lpdf(y[i]|(a + b * x[i]),sigma);
    
  }
  chi2dat = sum(pearson2dat);
  chi2sim = sum(pearson2sim);
  pvalue = chi2dat>chi2sim;
}

