# library("rJava")
# library("rCMA")
library("mlxR")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)



catModel <- inlineModel("
[LONGITUDINAL]
input = {a,b}

EQUATION:
lp0 = a-b*t
lp1 = a-b*t/2
p0  = 1/(1+exp(-lp0))
p1  = 1/(1+exp(-lp1)) -p0
p2  = 1-p0-p1

DEFINITION:
y = {type       = categorical, 
     categories = {0, 1, 2},
     P(y=0)     = p0,
     P(y=1)     = p1}


[INDIVIDUAL]
input={a_pop, o_a,
      b_pop, o_b}


DEFINITION:
a  ={distribution=normal, prediction=a_pop,  sd=o_a}
b  ={distribution=normal, prediction=b_pop,  sd=o_b}
")

seed <- 12345 
p    <- c(a_pop=40, o_a=0.5,
          b_pop=1,o_b=0.6)

pr   <- list(name=c('p0','p1','p2'), time=0:300)
y    <- list(name='y', time=seq(0, 300, by=2))
res  <- simulx(model     = catModel, 
               parameter = p, 
               output    = list(pr, y),
               settings  = list(seed=seed))

plot1 <- catplotmlx(res$y)
print(plot1)

writeDatamlx(res, result.file = "data/cat_data_test.csv")

res$y
head(res$y)


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


###read and prepare data for saemixData Object
cat_data.saemix<-read.table("data/cat_data_test.csv", header=T, sep=",")
head(cat_data.saemix)
data = cat_data.saemix[cat_data.saemix$ytype==1,]
nrow(data)
data$id = 1:nrow(data)

###Create saemixData Object
saemix.data<-saemixData(name.data=data,header=TRUE,sep=" ",
  na=NA, name.group=c("id"), name.response=c("y"), name.predictors=c("time"))



cat_data.model<-function(psi,id,xidep) {
level<-xidep[,1]
time<-xidep[,2]

a <- psi[id,1]
b <- psi[id,2]

lp0 = a-b*time
lp1 = a-b*time/2
p0  = 1/(1+exp(-lp0))
p1  = 1/(1+exp(-lp1)) -p0
p2  = 1-p0-p1


p.obs = (level==0)*p0+(level==1)*p1+(level==2)*p2
return(p.obs)
# logpdf <- log(p.obs)
# return(logpdf)
}


saemix.model<-saemixModel(model=cat_data.model,description="cat model",type="likelihood",   
  psi0=matrix(c(35,2),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("a","b"))), 
  transform.par=c(0,0),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE),omega.init=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE),error.model="constant")


K1 = 700
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2


options<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=TRUE,save.graphs=FALSE,nbiter.burn =0)
cat.fit<-saemix(saemix.model,saemix.data,options)

