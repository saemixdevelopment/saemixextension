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
# Note: remove the options save=FALSE and save.graphs=FALSE 
# to save the results and graphs
saveDir <- "myresults"
saemix.options<-list(seed=632545,save=TRUE,save.graphs=TRUE, displayProgress=FALSE, directory=saveDir)

# Not run (strict time constraints for CRAN)
saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
saemix.fit <- conddist.saemix(saemix.fit)
saemix.fit <- compute.sres(saemix.fit)

tab<-saemix.fit@data@data[,c(saemix.fit@data["name.group"], saemix.fit@data["name.predictors"],saemix.fit@data["name.response"],saemix.fit@data["name.covariates"], saemix.fit@data["name.mdv"])]
tab2<-data.frame(ipred=saemix.fit@results@ipred, icpred=saemix.fit@results@icpred, ypred=saemix.fit@results@ypred, ppred=saemix.fit@results@ppred)
tab <- cbind(tab, tab2)
write.table(tab, file.path(saveDir,"results.res"),row.names=FALSE, quote=F)

head(fitted(saemix.fit, type=c("ypred")))
head(fitted(saemix.fit, type=c("ppred")))
head(fitted(saemix.fit, type=c("ipred")))
head(fitted(saemix.fit, type=c("icpred")))

yfit <- saemix.predict(saemix.fit)
