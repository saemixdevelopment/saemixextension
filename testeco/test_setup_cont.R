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

# Model with 2 covariates and a covariance model
saemix.model3<-saemixModel(model=model1cpt,modeltype="structural",
          description="One-compartment model with first-order absorption",
          psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),
          covariance.model=matrix(c(1,0,0,0,1,1,0,1,1),ncol=3,byrow=TRUE),
          covariate.model=matrix(c(0,0,1,0,1,0),ncol=3,byrow=TRUE),
          transform.par=c(1,1,1),error.model="proportional")

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
saemix.fit2<-saemix(saemix.model2,saemix.data,saemix.options)
saemix.fit3<-saemix(saemix.model3,saemix.data,saemix.options)

theo.fit1<-saemix.fit
theo.fit2<-saemix.fit2
theo.fit3<-saemix.fit3

###################################################################################
# New dataset
xtim<-seq(0,24,2)
nsuj<-5
xwei<-seq(50,90,length.out = nsuj)
xsex<-rep(c("F","M"),length.out=nsuj)
xdose<-seq(280,320,length.out=nsuj)
theo.newdata<-data.frame(Id=rep(1:nsuj,each=length(xtim)),Time=rep(xtim,nsuj),Dose=rep(xdose,each=length(xtim)), Weight=rep(xwei,each=length(xtim)),Sex=rep(xsex,each=length(xtim)))

saemixObject<-saemix.fit2
psiM<-data.frame(ka=seq(1.6,2,0.1),V=seq(34,30),CL=c(2,2.5,2,2.5,2))
fpred<-saemixObject["model"]["model"](psiM, theo.newdata$Id, theo.newdata[,c("Dose","Time")])
theo.newdata$Concentration<-fpred+rnorm(length(fpred),sd=0.5)

theo.psiM<-psiM
