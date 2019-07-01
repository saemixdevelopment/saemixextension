###################################################################################
cat("Running example Theophylline\n")

theo.saemix<-read.table(file.path(datDir,"theo.saemix.tab"),header=T)
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

# Model with covariate Weight

saemix.model<-saemixModel(model=model1cpt,modeltype="structural",
                          description="One-compartment model with first-order absorption",
                          psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),
                          transform.par=c(1,1,1),covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)

# Model with 2 covariates
saemix.model2<-saemixModel(model=model1cpt,modeltype="structural",
                           description="One-compartment model with first-order absorption",
                           psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),
                           transform.par=c(1,1,1),covariate.model=matrix(c(0,0,1,0,1,0),ncol=3,byrow=TRUE))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
saemix.fit2<-saemix(saemix.model2,saemix.data,saemix.options)


# Model with 3 covariates
saemix.model2bis<-saemixModel(model=model1cpt,modeltype="structural",
                           description="One-compartment model with first-order absorption",
                           psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),
                           transform.par=c(1,1,1),covariate.model=matrix(c(1,0,1,0,1,0),ncol=3,byrow=TRUE))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
saemix.fit2bis<-saemix(saemix.model2bis,saemix.data,saemix.options)

# Model with two random effects, V not random

saemix.model3<-saemixModel(model=model1cpt,modeltype="structural",
                          description="One-compartment model with first-order absorption",
                          psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),
                          transform.par=c(1,1,1),covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE),
                          covariance.model = matrix(c(1,0,0,0,0,0,0,0,1),ncol=3,byrow=TRUE))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
saemix.fit3<-saemix(saemix.model3,saemix.data,saemix.options)

# Model with two random effects, Cl not random

saemix.model4<-saemixModel(model=model1cpt,modeltype="structural",
                           description="One-compartment model with first-order absorption",
                           psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),
                           transform.par=c(1,1,1),covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE),
                           covariance.model = matrix(c(1,0,0,0,1,0,0,0,0),ncol=3,byrow=TRUE))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
saemix.fit4<-saemix(saemix.model4,saemix.data,saemix.options)


theo.fit1<-saemix.fit
theo.fit2<-saemix.fit2
theo.fit3<-saemix.fit3