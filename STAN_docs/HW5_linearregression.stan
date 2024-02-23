//set up the data sturcture

data {
  int N;
  vector[N] y;
  vector[N] x;
  real volDose;
}

//set up the data for the parameter estimates

parameters {
  real a;
  real b;
  real <lower = 0 > sigma;
  
}

//set up the model

model {
  //Priors
  a ~ normal(0,100);
  b ~ normal(0,100);
  sigma ~ exponential(0.01);
  
  //model
 for(i in 1:N) {
   y[i] ~ normal(a + b * x[i], sigma);
 }
}

//calculations from the model to assess fit, model selection, etc.

generated quantities {
  real chi2dat;
  real chi2sim;
  real pvalue;
  real probPos;
  real resid[N];
  real ypred[N];
  real pearson2dat[N];
  real pearson2sim[N];
  real loglik[N];
  real ymean[N];
  real volPred;
  real volMean;
  
  for (i in 1:N)  {
    resid[i] = y[i]-(a + b * x[i]);
    ypred[i] = normal_rng(a + b *x[i], sigma); 
    pearson2dat[i] = (y[i]-(a + b * x[i]))/sigma;
    pearson2sim[i] = (ypred[i]-(a + b * x[i]))/sigma;
    loglik[i] = normal_lpdf(y[i]|(a + b * x[i]),sigma);
    ymean[i] = a + b *x[i];
  }
  
  volMean = a+b*volDose;
  volPred = normal_rng(a + b*volDose, sigma);
  probPos = if_else(b>0,1,0);
  chi2dat = sum(pearson2dat^2);
  chi2sim = sum(pearson2sim^2);
  pvalue = chi2dat>chi2sim;
}