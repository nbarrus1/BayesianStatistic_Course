data{
  int<lower=0> Offspring[35];
  int<lower=1,upper=2> Intact[35];
}
parameters{
  real<lower=0> lambda[2];
}
model
{
  lambda ~ gamma(0.001, 0.001);  // broad prior for mean productivity
  for(i in 1:35) {
    Offspring[i] ~ poisson(lambda[Intact[i]]);  
  }
}
generated quantities{
  real LL[35];
  for(i in 1:35) {
   LL[i]= poisson_lpmf(Offspring[i] | lambda[Intact[i]]);
 }
}
