data{
  int<lower=0> Offspring[35];
}
parameters{
  real<lower=0> lambda;
}
model
{
  lambda ~ gamma(0.001, 0.001);  // broad prior for mean productivity
  Offspring ~ poisson(lambda);  // productivity drawn from a Poisson dist'n
}
generated quantities{
  real LL[35];
  for(i in 1:35) {
   LL[i]= poisson_lpmf(Offspring[i] | lambda);
 }
}
