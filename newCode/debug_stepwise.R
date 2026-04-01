library(saemix)

# doc de theo.saemix
data(theo.saemix)
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA,
                        name.group=c("Id"),name.predictors=c("Dose","Time"),
                        name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
                        units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
model1cpt<-function(psi,id,xidep) {
  dose<-xidep[,1]
  tim<-xidep[,2]
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  k<-CL/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}
# Default model, no covariate
saemix.model<-saemixModel(model=model1cpt,
                          description="One-compartment model with first-order absorption",
                          psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE,
                                      dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1))


saemix.options<-list(seed=39546,save=FALSE,save.graphs=FALSE,displayProgress=FALSE)
saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)


# ensuite on definit covariate.init et on fait tourner stepwise.procedure:
covariate1 <- matrix(data=0,ncol=3,nrow=2)
myfit<-stepwise.procedure(saemix.fit, covariate1)

print(myfit)

myfit@transform.par

# on fait tourner
myfit@transform.par <- saemix.fit@model@transform.par

saemix.fit.cov<-saemix(myfit,saemix.data,saemix.options)
print(saemix.fit.cov) 
