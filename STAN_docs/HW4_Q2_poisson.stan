
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
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  lamda ~ lognormal(0.0,1000000);
  Y ~ poisson(lamda);
}

