
// The input data 
data {
  int<lower=0> N;
  int<lower=0> Y[N];
}

// The parameters accepted by the model. lamda
parameters {
  real<lower=0> lamda;
}

// The model to be estimated. We model the output
// 'y' to be poisson distributed with mean 'Lambda'
model {
  lamda ~ lognormal(0.0,1000000);
  Y ~ poisson(lamda);
}
//things we want to generate (i.e., like the LL profiles for model selection etc)
generated quantities{
  real LL[N];
  for (i in 1:N) {
    LL[i] = poisson_lpmf(Y[i] | lamda);
  }
}