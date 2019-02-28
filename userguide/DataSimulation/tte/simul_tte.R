# library("rJava")
# library("rCMA")
library("mlxR")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)


model2 <- inlineModel("

  [LONGITUDINAL]
  input = {beta,lambda}  

  EQUATION:
  h=(beta/lambda)*(t/lambda)^(beta-1)

  DEFINITION:
  e = {type               = event, 
       rightCensoringTime = 6,  
       hazard             = h}
  [INDIVIDUAL]
  input={lambda_pop, o_lambda,beta_pop, o_beta}
                        
  DEFINITION:
  lambda  ={distribution=lognormal, prediction=lambda_pop,  sd=o_lambda}
  beta  ={distribution=lognormal, prediction=beta_pop,  sd=o_beta}
       ")

lambda_true <- 5
o_lambda_true <- 0.3
beta_true <- 3
o_beta_true <- 0.3


  p <- c(lambda_pop=lambda_true, o_lambda= o_lambda_true,
         beta_pop = beta_true,o_beta = o_beta_true)
  h <- list(name='h', time=seq(0, 6, by=1))
  e <- list(name='e', time=0)

  N <- 10
  res <- simulx(model     = model2, 
                settings  = list(seed=13),
                parameter = p, 
                output    = list(h,e), 
                 group     = list(size = N))
  

writeDatamlx(res, result.file = "data_tte.csv")



