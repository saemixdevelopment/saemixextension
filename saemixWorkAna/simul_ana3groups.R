library("mlxR")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)

catModel <- inlineModel("
[LONGITUDINAL]
input =  {beta0,gamma0,delta0, dose}
dose = {use=regressor}
EQUATION:
lm0 = beta0+gamma0*t + delta0*dose

D = exp(lm0)+1
p0 = exp(lm0)/D
p1 = 1/D

DEFINITION:
y = {type=categorical, categories={0, 1}, 
     P(y=0)=p0,
     P(y=1)=p1}

[INDIVIDUAL]
input={beta0_pop, o_beta0,
      gamma0_pop, o_gamma0,
      delta0_pop, o_delta0}

DEFINITION:
beta0  ={distribution=normal, prediction=beta0_pop,  sd=o_beta0}
gamma0  ={distribution=lognormal, prediction=gamma0_pop,  sd=o_gamma0}
delta0  ={distribution=lognormal, prediction=delta0_pop,  sd=o_delta0}
")

# Design
nobs <- 15
tobs <- seq(-10, 60, by=nobs) # nobs understood as the time between 2 measurements, not as number of observations, which would be seq(-10, 60, length.out=nobs)

popAna <- c(beta0_pop=1, o_beta0=0.3, 
            gamma0_pop= 0.5, o_gamma0=0.3, 
            delta0_pop=1, o_delta0=0.2)

N  <- 500
reg1 <- list(name='dose',time=tobs,value=10*(tobs>0))
reg2 <- list(name='dose',time=tobs,value=20*(tobs>0))
reg3 <- list(name='dose',time=tobs,value=30*(tobs>0))
reg<-list(reg1,reg2,reg3)
out  <- list(name='y', time=tobs)
g <- list(size=N, level='individual')

res1 <- simulx(model = catModel, parameter = popAna, regressor = reg1, group=g, output    = out)
res2 <- simulx(model = catModel, parameter = popAna, regressor = reg2, group=g, output    = out)
res3 <- simulx(model = catModel, parameter = popAna, regressor = reg3, group=g, output    = out)
res2$y$id <- res2$y$id+N
res3$y$id <- res3$y$id+2*N

simRes<-rbind(cbind(res1$y,dose=rep(reg1$value,N)), cbind(res2$y,dose=rep(reg2$value,N)), cbind(res3$y,dose=rep(reg3$value,N)))

ggplot(simRes, aes(x = time, fill=y)) + geom_bar(position='dodge') + facet_wrap(.~dose)
