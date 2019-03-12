###################################################################################
cat("Running example Theophylline\n")

theo.saemix<-read.table(file.path(datDir,"theo.saemix.tab"),header=T)
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
#                         name.group=c("Id"),name.predictors=c("Dose","Time"),
#                         name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#                         units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
                        name.group=c("Id"),name.predictors=c("Dose","Time"),
                        name.response=c("Concentration"), name.X="Time")

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
                          transform.par=c(1,1,1))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)

# Model with 2 covariates
saemix.model2<-saemixModel(model=model1cpt,modeltype="structural",
                           description="One-compartment model with first-order absorption",
                           psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),
                           transform.par=c(1,1,1),covariate.model=matrix(c(0,0,1,0,1,0),ncol=3,byrow=TRUE))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
saemix.fit2<-saemix(saemix.model2,saemix.data,saemix.options)

theo.fit1<-saemix.fit
theo.fit2<-saemix.fit2

###################################################################################
# New dataset
xtim<-seq(0,24,1)
nsuj<-10
xwei<-seq(50,90,length.out = nsuj)
xsex<-rep(c("F","M"),length.out=nsuj)
xdose<-seq(280,320,length.out=nsuj)
theo.newdata<-data.frame(Id=rep(1:nsuj,each=length(xtim)),Time=rep(xtim,nsuj),Dose=rep(xdose,each=length(xtim)), Weight=rep(xwei,each=length(xtim)),Sex=rep(xsex,each=length(xtim)))

saemixObject<-saemix.fit
psiM<-data.frame(ka=seq(1.5,1.69,0.02),V=seq(31.5,31.69,0.02),CL=seq(2.5,2.69,0.02))
fpred<-saemixObject["model"]["model"](psiM, theo.newdata$Id, theo.newdata[,c("Dose","Time")])
theo.newdata$Concentration<-fpred+rnorm(length(fpred),sd=0.5)

# theo.newdata<-read.table(file.path(datDir,"theo.saemix.tab"),header=T)
theo.newdata <-theo.saemix[theo.saemix$Id<11,]
theo.psiM<-psiM
