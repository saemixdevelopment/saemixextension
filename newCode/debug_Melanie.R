rm(list=ls())

library(saemix)
library(numDeriv)

saemixPath = "C:/Users/MelanieGUHL/Documents/Packages/saemixextension_06092023/newCode/"

source(file.path(saemixPath,"compute_derivFIMlin.R"))

###################################### Theophylline example (N=150, n=10)

# Initial parameters
ka = 1.5
Cl = 0.04
V = 0.5

# PK model ----
model1cpt_firstorder=function(psi,id,xidep) { 
  dose=4
  tim=xidep[,1]  
  ka=psi[id,1]
  CL=psi[id,2]
  V = psi[id,3]
  k = CL/V
  pred = dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(pred)
}

set.seed(123456)

simu = read.table(paste0(saemixPath,"dataset1_N150_n10.txt"),header=T)
simu = simu[simu$Time>0,]
saemix_simu=saemixData(name.data=simu,
                       header=TRUE,
                       sep=" ",
                       na=NA, 
                       name.group=c("Id"),
                       name.predictors=c("Time"),
                       name.response=c("Concentration"),
                       name.covariate="Tr",
                       name.X="Time",
                       units=list(x="day",y="mg/mL"))

saemix_model=saemixModel(model=model1cpt_firstorder,
                         psi0=matrix(c(ka,Cl,V),ncol=3, 
                                     dimnames=list(NULL, c("ka","CL","V"))),
                         transform.par=c(1,1,1),
                         covariance.model=matrix(c(1,0,0, 0,1,0, 0,0,1),ncol=3,byrow=TRUE),
                         covariate.model=c(1,1,1),
                         error.model="proportional")

saemix_options=list(save=FALSE,save.graphs=FALSE,nb.chains = 10, nbiter.saemix=c(300,100))
saemix_fit=saemix(saemix_model,saemix_simu,saemix_options)

# Comparing resulting confidence intervals
saemix_fit.fim <- computeFIM.lin(saemix_fit)
#saemix_fit.fim@results@conf.int
#saemix_fit@results@conf.int

# Comparing LL and SE
invfim1<-computeFIMinv.lin(saemix_fit)
ll1<-computeLL.lin(saemix_fit)
cat("LL, saemix 3.1=",saemix_fit@results@ll.lin, "\n")
cat("LL, new code=",ll1, "\n")

# LL, saemix 3.1= -1117.219
#   LL, new code= -1117.219

cat("SE, saemix 3.1=",sqrt(diag(solve(saemix_fit@results@fim))),"\n")
cat("SE, new code=",sqrt(diag(invfim1)),"\n")

# Inversion by bloc versus the whole matrix
cat("bloc inversion", sqrt(diag(solve(solve(invfim1[1:4,1:4])))), sqrt(diag(solve(solve(invfim1[5:10,5:10])))),"\n")

# SE, saemix 3.1= 0.03990953 0.03763673 0.0005889104 0.02045829 0.01147454 0.03373173 0.006094592 0.001811402 0.004873111 0.002230607
#   SE, new code= 0.04052455 0.03786594 0.0005867822 0.02040922 0.01150242 0.03377142 0.006231433 0.001836635 0.004899144 0.002233921
#  bloc inversion 0.04052455 0.03786594 0.0005867822 0.02040922 0.01150242 0.03377142 0.006231433 0.001836635 0.004899144 0.002233921



###################################### Theophylline example (N=12, n=3)

set.seed(123456)

simu = read.table(paste0(saemixPath,"dataset1_N12_n3.txt"),header=T)
simu = simu[simu$Time>0,]
saemix_simu=saemixData(name.data=simu,
                       header=TRUE,
                       sep=" ",
                       na=NA, 
                       name.group=c("Id"),
                       name.predictors=c("Time"),
                       name.response=c("Concentration"),
                       name.covariate="Tr",
                       name.X="Time",
                       units=list(x="day",y="mg/mL"))

saemix_model=saemixModel(model=model1cpt_firstorder,
                         psi0=matrix(c(ka,Cl,V),ncol=3, 
                                     dimnames=list(NULL, c("ka","CL","V"))),
                         transform.par=c(1,1,1),
                         covariance.model=matrix(c(1,0,0, 0,1,0, 0,0,1),ncol=3,byrow=TRUE),
                         covariate.model=c(1,1,1),
                         error.model="proportional")

saemix_options=list(save=FALSE,save.graphs=FALSE,nb.chains = 10, nbiter.saemix=c(300,100))
saemix_fit=saemix(saemix_model,saemix_simu,saemix_options)

# Comparing resulting confidence intervals
saemix_fit.fim <- computeFIM.lin(saemix_fit)
#saemix_fit.fim@results@conf.int
#saemix_fit@results@conf.int

# Comparing LL and SE
invfim1<-computeFIMinv.lin(saemix_fit)
ll1<-computeLL.lin(saemix_fit)
cat("LL, saemix 3.1=",saemix_fit@results@ll.lin, "\n")
cat("LL, new code=",ll1, "\n")

# LL, saemix 3.1= -15.26239
#   LL, new code= -15.26239

cat("SE, saemix 3.1=",sqrt(diag(solve(saemix_fit@results@fim))),"\n")
cat("SE, new code=",sqrt(diag(invfim1)),"\n")

# Inversion by bloc versus the whole matrix
cat("bloc inversion", sqrt(diag(solve(solve(invfim1[1:4,1:4])))), sqrt(diag(solve(solve(invfim1[5:10,5:10])))),"\n")

# SE, saemix 3.1= 0.0907322 0.06960148 0.002067061 0.06369995 0.04312661 0.119726 0.01331816 0.005146355 0.01743589 0.04380425
#   SE, new code= 0.0913758 0.06979865 0.002065077 0.06369577 0.04315726 0.119797 0.01061585 0.004919629 0.01770119 0.03933122 
#  bloc inversion 0.0913758 0.06979865 0.002065077 0.06369577 0.04315726 0.119797 0.01061585 0.004919629 0.01770119 0.03933122


###################################### Theophylline example (N=12, n=3) with corr = 0.90,0.98,0.9 and omega2 = 1.2,1.2,1.2

set.seed(123456)


simu = read.table(paste0(saemixPath,"dataset1_N12_n3_highIIV_corr.txt"),header=T)
simu = simu[simu$Time>0,]
saemix_simu=saemixData(name.data=simu,
                       header=TRUE,
                       sep=" ",
                       na=NA, 
                       name.group=c("Id"),
                       name.predictors=c("Time"),
                       name.response=c("Concentration"),
                       name.covariate="Tr",
                       name.X="Time",
                       units=list(x="day",y="mg/mL"))

saemix_model=saemixModel(model=model1cpt_firstorder,
                         psi0=matrix(c(ka,Cl,V),ncol=3, 
                                     dimnames=list(NULL, c("ka","CL","V"))),
                         transform.par=c(1,1,1),
                         covariance.model=matrix(1,nrow=3,ncol=3),
                         covariate.model=c(1,1,1),
                         error.model="proportional")

saemix_options=list(save=FALSE,save.graphs=FALSE,nb.chains = 10, nbiter.saemix=c(300,100))

saemix_fit=saemix(saemix_model,saemix_simu,saemix_options)

# Comparing resulting confidence intervals
saemix_fit.fim <- computeFIM.lin(saemix_fit)
#saemix_fit.fim@results@conf.int
#saemix_fit@results@conf.int

# Comparing LL and SE
invfim1<-computeFIMinv.lin(saemix_fit)
ll1<-computeLL.lin(saemix_fit)
cat("LL, saemix 3.1=",saemix_fit@results@ll.lin, "\n")
cat("LL, new code=",ll1, "\n")

# LL, saemix 3.1= -60.71805
#   LL, new code= -60.71805

cat("SE, saemix 3.1=",sqrt(diag(solve(saemix_fit@results@fim))),"\n")
cat("SE, new code=",sqrt(diag(invfim1)),"\n")

# Inversion by bloc versus the whole matrix
cat("bloc inversion", sqrt(diag(solve(solve(invfim1[1:4,1:4])))), sqrt(diag(solve(solve(invfim1[5:10,5:10])))),"\n")

# SE, saemix 3.1= 0.4528102 0.5490196 0.01358747 0.5072418 0.1067797 0.4759535 0.3789116 0.287953 0.3184628 0.3150379 0.2617716 0.2795965 0.07560011
#   SE, new code= NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
#  bloc inversion NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN