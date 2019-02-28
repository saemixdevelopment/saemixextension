# library("rJava")
# library("rCMA")
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
gamma0  ={distribution=normal, prediction=gamma0_pop,  sd=o_gamma0}
delta0  ={distribution=normal, prediction=delta0_pop,  sd=o_delta0}
")


nobs = 15
tobs<- seq(-20, 50, by=nobs)

reg1 <- list(name='dose',
            time=tobs,
            value=10*(tobs>0))

reg2 <- list(name='dose',
            time=tobs,
            value=20*(tobs>0))

reg3 <- list(name='dose',
            time=tobs,
            value=30*(tobs>0))

out  <- list(name='y', time=tobs)
N  <- 100
p <- c(beta0_pop=-4, o_beta0=0.3, 
       gamma0_pop= -0.5, o_gamma0=0.2,
       delta0_pop=1, o_delta0=0.2)

g1 <- list(size=N,regressor = reg1)
g2 <- list(size=N,regressor = reg2)
g3 <- list(size=N,regressor = reg3)
g <- list(g1,g2,g3)
res <- simulx(model=catModel,output=out, group=g,parameter=p)
plot1 <- catplotmlx(res$y)
print(plot1)

writeDatamlx(res, result.file = "data_ord.csv")



