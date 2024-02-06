####script title: stan setup
####datae fist created: 2/6/2024
###purpose: the purpose of this script is to set up stan in r, in order to run
###the different MCMC that stan has to offer for bayesian statistics.  For this
####to have worked you had to have the rtools43 downloaded which enables script
####to be compiled in C++ and other computer coding languages that deal directly
#### with binary. it can help the code run faster

install.packages("rstan")
library(rstan)


example(stan_model, package = "rstan", run.dontrun = TRUE)
