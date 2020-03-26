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
tobs<- seq(-10, 60, by=nobs)

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
g <- list(g1)
res <- simulx(model=catModel,output=out, group=g,parameter=p)
plot1 <- catplotmlx(res$y)
print(plot1)


writeDatamlx(res, result.file = "data/test.csv")




# library(saemixextension)
source('newR/aaa_generics.R') 
source('newR/compute_LL.R') 
source('newR/func_aux.R') 
source('newR/func_distcond.R') 
source('newR/func_FIM.R')
source('newR/func_plots.R') 
source('newR/func_simulations.R') 

source('newR/main.R')
source('newR/main_estep.R')
source('newR/main_initialiseMainAlgo.R') 
source('newR/main_mstep.R') 
source('newR/SaemixData.R')
source('newR/SaemixModel.R') 
source('newR/SaemixRes.R') 
# source('newR/SaemixRes_c.R') 
source('newR/SaemixObject.R') 
source('newR/zzz.R') 
library("mlxR")

cat_data.saemix<-read.table("data/test.csv", header=T, sep=",")
saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"), name.response=c("y"), name.predictors=c("dose","time"))



cat_data.model<-function(psi,id,xidep) {
level<-xidep[,1]
dose<-xidep[,2]
time<-xidep[,3]

beta0 <- psi[id,1]
gamma0 <- psi[id,2]
delta0 <- psi[id,3]

lm0 <- beta0+gamma0*time + delta0*dose

D <- exp(lm0)+1
P0 <- exp(lm0)/D
P1 <- 1/D

P.obs = (level==0)*P0+(level==1)*P1
# logpdf <- log(P.obs)
return(P.obs)
# return(logpdf)
}


saemix.model<-saemixModel(model=cat_data.model,description="cat model",type="likelihood",   
  psi0=matrix(c(1,1,1),ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("th1","th2","th3"))), 
  transform.par=c(0,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),error.model="constant")


K1 = 300
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2


options<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=TRUE,save.graphs=FALSE,nbiter.burn =0)
cat.fit<-saemix(saemix.model,saemix.data,options)


