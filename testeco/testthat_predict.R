###################################################################################
# Individual and population predictions for a new dataset

context("Testing predict.newdata for a structural model \n")

test_that("Comparing parameters", {
  test.newdata<-theo.newdata
  mylist<-predict.newdata(saemix.fit, theo.newdata, type=c("ipred", "ypred", "ppred", "icpred"))
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

context("Testing predict.newdata for a structural model \n")

test_that("Computing individual and population predictions for a new dataset with a valid structure and individual observations", {
  mylist<-predict.newdata(saemix.fit, theo.newdata, type=c("ipred", "ypred", "ppred", "icpred"))
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
  mylist<-predict.newdata(saemix.fit, theo.newdata, type=c("ipred", "ypred"))
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
  mylist<-predict.newdata(saemix.fit, theo.newdata, type=c("icpred", "ypred"))
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
  mylist<-predict.newdata(saemix.fit, theo2, type=c("ipred", "ypred", "ppred", "icpred"))
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
  mylist<-predict.newdata(saemix.fit, theo2, type=c("ipred", "ypred", "ppred", "icpred"))
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

context("Testing predict.newdata for a structural model \n")

test_that("Computing population predictions for a new dataset with a valid structure and individual observations using predict()", {
  fit.pred<-saemix.predict(saemix.fit)
  vec<-predict(saemix.fit)
  expect_equal(vec,fit.pred@results@predictions$ipred)
  expect_length(vec,saemix.fit@data@ntot.obs)
  expect_gte(cor(vec,saemix.fit@data@data[,saemix.fit@data@name.response]),0.8)
  vec<-predict(saemix.fit,type="ypred")
  expect_equal(vec,fit.pred@results@predictions$ypred)
  vec<-predict(saemix.fit,type="icpred")
  expect_equal(vec,fit.pred@results@predictions$icpred)
  vec<-predict(saemix.fit,type="ipred")
  expect_equal(vec,fit.pred@results@predictions$ipred)
})

test_that("Computing default (individual) predictions for a new dataset with a valid structure and individual observations using predict()", {
  vec<-predict(saemix.fit,theo.newdata)
  expect_gte(cor(vec,theo.newdata$Concentration),0.95)
  vec<-predict(saemix.fit,theo.newdata,type="ypred")
  expect_gte(cor(vec,theo.newdata$Concentration),0.95)
  vec<-predict(saemix.fit,theo.newdata,type="icpred")
  expect_gte(cor(vec,theo.newdata$Concentration),0.95)
})

test_that("Computing population predictions for a new dataset with a valid structure but no individual observations using predict()", {
  vec0<-predict(saemix.fit,theo.newdata,type="ypred")
  theo2<-theo.newdata[,1:5]
  vec<-predict(saemix.fit,theo2)
  expect_gte(cor(vec,theo.newdata$Concentration),0.95)
  par(mfrow=c(1,1))
  plot(vec0,vec,pch=20,col="Black")
  abline(0,1)
  expect_equal(vec,vec0)
})

  
###################################################################################
