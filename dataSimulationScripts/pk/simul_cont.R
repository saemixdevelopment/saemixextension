library("mlxR")
library(saemix)

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  Cl<-psi[id,3]
  k <- Cl/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}

model <- inlineModel("


[INDIVIDUAL]
input = {ka_pop, V_pop, Cl_pop, omega_ka, omega_V, omega_Cl}
DEFINITION:
ka = {distribution=lognormal, reference=ka_pop, sd=omega_ka}
V  = {distribution=lognormal, reference=V_pop,  sd=omega_V }
Cl = {distribution=lognormal, reference=Cl_pop, sd=omega_Cl}


[LONGITUDINAL]
input = {ka, V, Cl,a}
EQUATION:
C = pkmodel(ka,V,Cl)
DEFINITION:
y = {distribution=normal, prediction=C, sd=a}
")

N=5

param   <- c(
  ka_pop  = 1,    omega_ka  = 0.5,
  V_pop   = 10,   omega_V   = 0.4,
  Cl_pop  = 1,    omega_Cl  = 0.3, a =1)
  
res <- simulx(model     = model,
              parameter = param,
              treatment = list(time=0, amount=100),
              group     = list(size=N, level='individual'),
              output    = list(name='y', time=seq(1,4,by=1)))

 
res$y[2,3]
typeof(res$y[2,3])

writeDatamlx(res, result.file = "data_cont.txt", sep="\t")
