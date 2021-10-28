###################################################################################
# Individual and population predictions for a new dataset

context("Testing saemixPredictNewdata for a structural model \n")

test_that("Comparing parameters", {
  test.newdata<-theo.newdata
  mylist<-saemixPredictNewdata(saemix.fit, theo.newdata, type=c("ipred", "ypred", "ppred", "icpred"))
  expect_is(mylist, "list") # tests for particular class
  param<-mylist$param$map.psi
  par(mfrow=c(2,2))
  for(i in 1:3) {
    plot(theo.psiM[,i],param[,i],main=colnames(psiM)[i],xlab="Simulated",ylab="Estimated")
    abline(0,1)
  }
#  expect_gte(cor(theo.psiM[,1],param[,1]),0.9)
#  expect_gte(cor(theo.psiM[,2],param[,2]),0.9)
})

test_that("Computing individual and population predictions for a new dataset with a valid structure and individual observations", {
  mylist<-saemixPredictNewdata(saemix.fit, theo.newdata, type=c("ipred", "ypred", "ppred", "icpred"))
  expect_is(mylist, "list") # tests for particular class
  apred<-mylist$predictions
  par(mfrow=c(2,2))
  plot(theo.newdata$Concentration,apred$ipred,pch=20,col="Blue")
  points(theo.newdata$Concentration,apred$icpred,pch=20,col="Red")
  abline(0,1)
  legend(0.5,8,pch=20,col=c("Blue","Red"),c("icpred","ipred"))
  plot(apred$icpred,apred$ipred,pch=20,col="Black")
  abline(0,1)
  plot(theo.newdata$Concentration,apred$ypred,pch=20,col="Black")
  points(theo.newdata$Concentration,apred$ppred,pch=20,col="Red")
  abline(0,1)
  legend(0.5,8,pch=20,col=c("Black","Red"),c("ypred","ppred"))
  plot(apred$ypred,apred$ppred,pch=20,col="Black")
  abline(0,1)
  expect_gte(cor(apred$ypred,apred$ppred),0.95)
  expect_gte(cor(apred$ipred,apred$icpred),0.95)
})

test_that("Computing MAP and population predictions for a new dataset with a valid structure and individual observations", {
  mylist<-saemixPredictNewdata(saemix.fit, theo.newdata, type=c("ipred", "ypred"))
  expect_is(mylist, "list") # tests for particular class
  apred<-mylist$predictions
  par(mfrow=c(1,1))
  plot(theo.newdata$Concentration,apred$ipred,pch=20,col="Blue")
  points(theo.newdata$Concentration,apred$ypred,pch=20,col="Red")
  abline(0,1)
  legend(0.5,9,pch=20,col=c("Blue","Red"),c("Individual predictions (MAP)","Population predictions"))
  expect_gte(cor(apred$ypred,apred$ipred),0.9)
})

test_that("Computing individual conditional and population predictions for a new dataset with a valid structure and individual observations", {
  mylist<-saemixPredictNewdata(saemix.fit, theo.newdata, type=c("icpred", "ypred"))
  expect_is(mylist, "list") # tests for particular class
  apred<-mylist$predictions
  par(mfrow=c(1,1))
  plot(theo.newdata$Concentration,apred$icpred,pch=20,col="Blue",xlab="Observed concentration in the new dataset",ylab="Predictions")
  points(theo.newdata$Concentration,apred$ypred,pch=20,col="Red")
  abline(0,1)
  legend(0.5,9,pch=20,col=c("Blue","Red"),c("Individual predictions (cond mean)","Population predictions"))
  expect_gte(cor(apred$ypred,apred$icpred),0.9)
})

test_that("Computing individual and population predictions for a new dataset with a valid structure, no individual observations", {
  theo2<-theo.newdata[,1:5]
  mylist<-saemixPredictNewdata(saemix.fit, theo2, type=c("ipred", "ypred", "ppred", "icpred"))
  expect_is(mylist, "list") # tests for particular class
  apred<-mylist$predictions
  par(mfrow=c(1,2))
  plot(theo.newdata$Concentration,apred$ypred,pch=20,col="Black")
  points(theo.newdata$Concentration,apred$ppred,pch=20,col="Red")
  abline(0,1)
  legend(0.5,8,pch=20,col=c("Black","Red"),c("ypred","ppred"))
  plot(apred$ypred,apred$ppred,pch=20,col="Black")
  abline(0,1)
  expect_gte(cor(apred$ypred,apred$ppred),0.95)
  expect_length(apred$ipred,0)
  expect_length(apred$icpred,0)
})

test_that("Computing individual and population predictions for a new dataset with missing covariates", {
  theo2<-theo.newdata[,-c(4)]
  mylist<-saemixPredictNewdata(saemix.fit, theo2, type=c("ipred", "ypred", "ppred", "icpred"))
  expect_is(mylist, "list") # tests for particular class
  apred<-mylist$predictions
  par(mfrow=c(1,2))
  plot(theo.newdata$Concentration,apred$ypred,pch=20,col="Black")
  points(theo.newdata$Concentration,apred$ppred,pch=20,col="Red")
  legend(0.5,8,pch=20,col=c("Black","Red"),c("ypred","ppred"))
  abline(0,1)
  plot(apred$ypred,apred$ppred,pch=20,col="Black")
  abline(0,1)
  expect_gte(cor(apred$ypred,apred$ppred),0.95)
})

###################################################################################
# Using predict to return a vector of predictions

context("Testing predict for a structural model \n")

test_that("Computing population predictions for a new dataset with a valid structure and individual observations using predict()", {
  fit.pred<-saemix.predict(saemix.fit)
  vec<-predict(saemix.fit)
  expect_equal(vec,fit.pred@results@predictions$ipred)
  expect_length(vec,saemix.fit@data@ntot.obs)
  expect_gte(cor(vec,saemix.fit@data@data[,saemix.fit@data@name.response]),0.8)
  vec<-predict(saemix.fit,type="ppred")
  expect_equal(vec,fit.pred@results@predictions$ppred)
  vec<-predict(saemix.fit,type="icpred")
  expect_equal(vec,fit.pred@results@predictions$icpred)
  vec<-predict(saemix.fit,type="ipred")
  expect_equal(vec,fit.pred@results@predictions$ipred)
})

test_that("Computing default (individual) predictions for a new dataset with a valid structure and individual observations using predict()", {
  vec<-predict(saemix.fit,theo.newdata)
  expect_gte(cor(vec,theo.newdata$Concentration),0.95)
  vec<-predict(saemix.fit,theo.newdata,type="ppred")
  expect_gte(cor(vec,theo.newdata$Concentration),0.95)
  vec<-predict(saemix.fit,theo.newdata,type="icpred")
  expect_gte(cor(vec,theo.newdata$Concentration),0.95)
})

test_that("Computing population predictions for a new dataset with a valid structure but no individual observations using predict()", {
  vec0<-predict(saemix.fit,theo.newdata,type="ppred")
  theo2<-theo.newdata[,1:5]
  vec<-predict(saemix.fit,theo2)
  expect_gte(cor(vec,theo.newdata$Concentration),0.95)
  par(mfrow=c(1,1))
  plot(vec0,vec,pch=20,col="Black")
  abline(0,1)
  expect_equal(vec,vec0)
})

###################################################################################

test_that("Estimating individual parameters using final estimates", {
  smx.data<-saemixData(name.data=file.path(datDir,"theo.saemix.tab"),header=T,na=".", name.group=c("Id"),name.predictors=c("Dose","Time"),name.covariates=c("Weight","Sex"), name.response=c("Concentration"),units=list(x="hr",y="mg/L"), name.X="Time",verbose=F)
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
  covmat.final<-matrix(data=0,nrow = 2,ncol=3)
  covmat.final[theo.fit2@model@covariate.model==1]<-theo.fit2@results@betaC
  psi0.final<-matrix(theo.fit2@results@fixed.effects[theo.fit2@results@indx.fix],ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL")))
  psi0.final<-rbind(psi0.final,covmat.final)
  smx.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", psi0=psi0.final, transform.par=c(1,1,1), covariate.model=matrix(c(0,0,1,0,1,0),ncol=3,byrow=TRUE), fixed.estim=c(1,1,1), covariance.model=matrix(c(1,0,0,0,1,1,0,1,1),ncol=3,byrow=TRUE), omega.init=theo.fit2@results@omega, error.init=theo.fit2@results@respar)
  smx.opt<-saemixControl(nb.chains=5,nbiter.saemix = c(500,300), ipar.lmcmc = 100)
  x<-createSaemixObject.initial(smx.model,smx.data,smx.opt)
  expect_equal(x@results@status,"initial")
  expect_equal(x@results@name.fixed,x@model@name.fixed)
  xpred<-saemix.predict(theo.fit2)
  x1<-predict(theo.fit2,type="ppred")
  expect_identical(x1,xpred@results@predictions$ppred)
  x2<-saemixPredictNewdata(x,smx.data@data)
  expect_identical(saemixObject["results"]["fixed.effects"],theo.fit2@results@fixed.effects)
  expect_identical(theo.fit2@results@omega,theo.fit2@results@omega)
  expect_identical(x2$predictions$ppred,xpred@results@predictions$ppred)
  expect_gt(cor(x2$predictions$ipred,xpred@results@predictions$ipred),0.99)
  expect_gt(cor(x2$predictions$icpred,xpred@results@predictions$icpred),0.99)
})

test_that("Estimating individual parameters using initial estimates", {
  smx.data<-saemixData(name.data=file.path(datDir,"theo.saemix.tab"),header=T,na=".", name.group=c("Id"),name.predictors=c("Dose","Time"),name.covariates=c("Weight","Sex"), name.response=c("Concentration"),units=list(x="hr",y="mg/L"), name.X="Time",verbose=F)
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
  covmat.final<-matrix(data=0,nrow = 2,ncol=3)
  covmat.final[theo.fit2@model@covariate.model==1]<-c(1,0.005)
  psi0.final<-matrix(c(1.5,20,2),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL")))
  psi0.final<-rbind(psi0.final,covmat.final)
  smx.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", psi0=psi0.final, transform.par=c(1,1,1), covariate.model=matrix(c(0,0,1,0,1,0),ncol=3,byrow=TRUE), fixed.estim=c(1,1,1), covariance.model=matrix(c(1,0,0,0,1,1,0,1,1),ncol=3,byrow=TRUE), omega.init=theo.fit2@results@omega, error.init=theo.fit2@results@respar)
  smx.opt<-saemixControl(nb.chains=5,nbiter.saemix = c(500,300), ipar.lmcmc = 100)
  x<-createSaemixObject.initial(smx.model,smx.data,smx.opt)
  expect_equal(x@results@status,"initial")
  expect_equal(x@results@name.fixed,x@model@name.fixed)
  xpred<-saemix.predict(theo.fit2)
  x2<-saemixPredictNewdata(x,smx.data@data)
  expect_gt(cor(x2$predictions$ppred,xpred@results@predictions$ppred),0.7)
  expect_gt(cor(x2$predictions$ipred,xpred@results@predictions$ipred),0.95)
  expect_gt(cor(x2$predictions$icpred,xpred@results@predictions$icpred),0.95)
})

if(FALSE) {
  saemixObject<-x
  print(head(x2$predictions))
  print(head(xpred@results@predictions))
  print(head(x2$param$population))
  print(head(xpred@results@mean.phi))
}
