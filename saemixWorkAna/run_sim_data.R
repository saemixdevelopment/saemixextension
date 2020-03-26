# library(saemixextension)
source('R/aaa_generics.R') 
source('R/compute_LL.R') 
source('R/func_aux.R') 
source('R/func_distcond.R') 
source('R/func_FIM.R')
source('R/func_plots.R') 
source('R/func_simulations.R') 

source('R/main.R')
source('R/main_estep.R')
source('R/main_initialiseMainAlgo.R') 
source('R/main_mstep.R') 
source('R/SaemixData.R')
source('R/SaemixModel.R') 
source('R/SaemixRes.R') 
# source('R/SaemixRes_c.R') 
source('R/SaemixObject.R') 
source('R/zzz.R') 
library("mlxR")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)


catModel <- inlineModel("
[LONGITUDINAL]
input =  {gamma0,delta0, dose}
dose = {use=regressor}
EQUATION:
lm0 = gamma0*t + delta0*dose

D = exp(lm0)+1
p0 = exp(lm0)/D
p1 = 1/D

DEFINITION:
y = {type=categorical, categories={0, 1}, 
     P(y=0)=p0,
     P(y=1)=p1}

[INDIVIDUAL]
input={gamma0_pop, o_gamma0,
      delta0_pop, o_delta0}


DEFINITION:
gamma0  ={distribution=normal, prediction=gamma0_pop,  sd=o_gamma0}
delta0  ={distribution=normal, prediction=delta0_pop,  sd=o_delta0}
")


nobs = 15
tobs<- seq(-10, 60, by=nobs)

reg1 <- list(name='dose',
            time=tobs,
            value=10*(tobs>0))

out  <- list(name='y', time=tobs)
N  <- 500
p <- c(gamma0_pop= -0.5, o_gamma0=0.2,
       delta0_pop=1, o_delta0=0.2)

g1 <- list(size=N,regressor = reg1)
g <- list(g1)
res <- simulx(model=catModel,output=out, group=g,parameter=p)

plot1 <- catplotmlx(res$y)
print(plot1)

head(res$y)
writeDatamlx(res, result.file = "data/data_cat_easy.csv")




cat_data.saemix<-read.table("data/data_cat_easy.csv", header=T, sep=",")
saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"), name.response=c("y"), name.predictors=c("dose","time"))



cat_data.model<-function(psi,id,xidep) {
level<-xidep[,1]
dose<-xidep[,2]
time<-xidep[,3]

gamma0 <- psi[id,1]
delta0 <- psi[id,2]

lm0 <- gamma0*time + delta0*dose

D <- exp(lm0)+1
P0 <- exp(lm0)/D
P1 <- 1/D

P.obs = (level==0)*P0+(level==1)*P1
logpdf <- log(P.obs)
return(P.obs)
# return(logpdf)
}


saemix.model<-saemixModel(model=cat_data.model,description="cat model",type="likelihood",   
  psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("th1","th2"))), 
  transform.par=c(0,0),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE),omega.init=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE),error.model="constant")


K1 = 300
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2


options<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=TRUE,save.graphs=FALSE,nbiter.burn =0)
cat.fit<-saemix(saemix.model,saemix.data,options)
