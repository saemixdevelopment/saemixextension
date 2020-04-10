saemixDir<-"/home/eco/work/saemix/saemixextension"
progDir<-file.path(saemixDir,"R")
anaDir<-file.path(saemixDir,"saemixWorkAna")
data.ana<-read.table(file.path(anaDir,"data","data_cat.csv"), header=T, sep=",")

compute.prob<-function(psi,id,xidep) {
  dose<-xidep[,1]
  time<-xidep[,2]
  beta0 <- psi[id,1]
  gamma0 <- psi[id,2]
  delta0 <- psi[id,3]
  
  lm0 <- beta0+gamma0*time + delta0*dose
  #  p0<-exp(lm0)/(exp(lm0)+1)
  p1<-1/(exp(lm0)+1)
  return(p1)
}

simcat<-function(prob) {
  if(min(prob)<0 | max(prob)>1) {
    cat("Input a vector of probabilities\n")
    return(NULL)
  }
  x<-runif(length(prob))
  return(as.integer(x>=prob))
}

nsuj<-length(unique(data.ana$id))

id1<-data.ana$id
xidep1<-data.ana[,c(5,3)]

# population predictions
psi0<-c(1,0.5,1)
sd.psi0<-c(0.3,0.3,0.2)
psi1<-do.call(rbind,rep(list(psi0),nsuj))
p1<-compute.prob(psi1,id1,xidep1)
cat("Using the population parameters\n")
cat("Probability of an event at t=(-10)\n")
print(summary(p1[xidep1[,2]==(-10)]))
cat("Probability of an event at t>0\n")
print(summary(p1[xidep1[,2]>0]))

# Checking data characteristics
times<-unique(data.ana$time)
doses<-unique(data.ana$dose)
nobs<-tapply(data.ana$id,data.ana$id,length)
table(nobs)

# Predicted probabilities for the different groups with the population parameters
levs<-expand.grid(dose=doses[doses>0],times=times[times>0])
levs<-rbind(levs,c(0,-10))
levs<-levs[order(levs[,2]),]
psi1.pop<-do.call(rbind,rep(list(psi0),dim(levs)[1]))
id1.pop<-1:dim(levs)[1]
p1.pop<-compute.prob(psi1.pop,id1.pop,levs)
print(cbind(levs,"P(Y=1)"=p1.pop))

# individual predictions
psi2<-psi1
i<-1
psi1[,i]<-psi1[,i]*exp(rnorm(nsuj,0,sd=sd.psi0[i]))
for(i in 2:3) psi1[,i]<-psi1[,i]*exp(rnorm(nsuj,0,sd=sd.psi0[i]))
head(psi)

cat("Using individual parameters\n")
pind1<-compute.prob(psi1,id1,xidep1)
ysim<-as.integer(runif(length(pind1))<pind1)
sum(ysim)
cat("Frequency of simulated events at t=(-10):",sum(ysim[data.ana$time==(-10)])/length(ysim[data.ana$time==(-10)]),"\n")
cat("Frequency of simulated events at t>0:",sum(ysim[data.ana$time>0])/length(ysim[data.ana$time>0]),"\n")

#### Simulx version
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
gamma0  ={distribution=lognormal, prediction=gamma0_pop,  sd=o_gamma0}
delta0  ={distribution=lognormal, prediction=delta0_pop,  sd=o_delta0}
")


nobs = 15
tobs<- seq(-10, 60, by=nobs) # nobs understood as the time between 2 measurements, not as number of observations, which would be seq(-10, 60, length.out=nobs)

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
N  <- 500
p <- c(beta0_pop=1, o_beta0=0.3, 
       gamma0_pop= 0.5, o_gamma0=0.3,
       delta0_pop=1, o_delta0=0.2)

g1 <- list(size=N,regressor = reg1)
g2 <- list(size=N,regressor = reg2)
g3 <- list(size=N,regressor = reg3)
g <- list(g1,g2,g3)
res <- simulx(model=catModel,output=out, group=g,parameter=p)
plot1 <- catplotmlx(res$y)
print(plot1)

if(FALSE) {
  writeDatamlx(res, result.file = "data/test2.csv")
}

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
N  <- 500
popAna <- c(beta0_pop=1, o_beta0=0.3, gamma0_pop= 0.5, o_gamma0=0.3, delta0_pop=1, o_delta0=0.2)
out  <- list(name='y', time=tobs)

reg1 <- list(name='dose',time=tobs,value=10*(tobs>0))
reg2 <- list(name='dose',time=tobs,value=20*(tobs>0))
reg3 <- list(name='dose',time=tobs,value=30*(tobs>0))
reg<-list(reg1,reg2,reg3)
g <- list(size=N, level='individual')

res1 <- simulx(model = catModel, parameter = popAna, regressor = reg1, group=g, output    = out)
res2 <- simulx(model = catModel, parameter = popAna, regressor = reg2, group=g, output    = out)
res3 <- simulx(model = catModel, parameter = popAna, regressor = reg3, group=g, output    = out)
res2$y$id <- res2$y$id+N
res3$y$id <- res3$y$id+2*N

simRes<-rbind(cbind(res1$y,dose=rep(reg1$value,N)), cbind(res2$y,dose=rep(reg2$value,N)), cbind(res3$y,dose=rep(reg3$value,N)))

ggplot(simRes, aes(x = time, fill=y)) + geom_bar(position='dodge') + facet_wrap(.~dose)

if(FALSE) {
  writeDatamlx(simRes, result.file = "data/test2.csv")
}
