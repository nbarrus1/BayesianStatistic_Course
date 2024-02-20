
data{
 int N;
 int Ncoef;
 vector[N] Y;
 matrix[N,Ncoef] xMatrix;
}
parameters{
 vector[Ncoef] b;
 real<lower=0> Sigma;
}
model{
  b~normal(0,10);
  Sigma~exponential(1);
  Y~normal(xMatrix*b,Sigma);
}

