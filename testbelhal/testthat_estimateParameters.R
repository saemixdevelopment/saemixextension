rootDir<-"/Users/karimimohammedbelhal/Desktop/Phd/R_Package/contributions/FinalsaemixExtension/ecomets"
datDir<-"/Users/karimimohammedbelhal/Desktop/Phd/R_Package/contributions/FinalsaemixExtension/ecomets/data"


## ORDINAL DATA
source(file.path(rootDir,"otherCode_noncontinuous","test_setup_ord.R"))

## TTE DATA
source(file.path(rootDir,"otherCode_noncontinuous","test_setup_tte.R"))

###################################################################################
# Individual and population predictions for a new dataset

context("Testing predict.newdata for a structural model \n")

test_that("Computing individual and population predictions for a new dataset with a valid structure and individual observations", {
  mylist<-predict.newdata(saemix.fit, test.newdata, type=c("ipred", "ypred", "ppred", "icpred"))
  expect_is(mylist, "list") # tests for particular class
  apred<-mylist$predictions
  par(mfrow=c(2,2))
  plot(test.newdata$LogProbs,apred$ipred,pch=20,col="Blue")
  points(test.newdata$LogProbs,apred$icpred,pch=20,col="Red")
  abline(0,1)
  plot(apred$icpred,apred$ipred,pch=20,col="Black")
  abline(0,1)
  plot(test.newdata$LogProbs,apred$ypred,pch=20,col="Black")
  points(test.newdata$LogProbs,apred$ppred,pch=20,col="Red")
  abline(0,1)
  plot(apred$ypred,apred$ppred,pch=20,col="Black")
  abline(0,1)
  expect_gte(cor(apred$ypred,apred$ppred),0.95)
  expect_gte(cor(apred$ipred,apred$icpred),0.95)
})


test_that("Computing MAP and population predictions for a new dataset with a valid structure and individual observations", {
  mylist<-predict.newdata(saemix.fit, test.newdata, type=c("ipred", "ypred"))
  expect_is(mylist, "list") # tests for particular class
  apred<-mylist$predictions
  par(mfrow=c(1,1))
  plot(test.newdata$LogProbs,apred$ipred,pch=20,col="Blue")
  points(test.newdata$LogProbs,apred$ypred,pch=20,col="Red")
  abline(0,1)
  expect_gte(cor(apred$ypred,apred$ipred),0.9)
})

test_that("Computing individual conditional and population predictions for a new dataset with a valid structure and individual observations", {
  mylist<-predict.newdata(saemix.fit, test.newdata, type=c("icpred", "ypred"))
  expect_is(mylist, "list") # tests for particular class
  apred<-mylist$predictions
  par(mfrow=c(1,1))
  plot(test.newdata$LogProbs,apred$icpred,pch=20,col="Blue")
  points(test.newdata$LogProbs,apred$ypred,pch=20,col="Red")
  abline(0,1)
  expect_gte(cor(apred$ypred,apred$icpred),0.9)
})


test_that("Computing individual and population predictions for a new dataset with a valid structure, no individual observations", {
  theo2<-test.newdata[,1:5]
  mylist<-predict.newdata(saemix.fit, theo2, type=c("ipred", "ypred", "ppred", "icpred"))
  expect_is(mylist, "list") # tests for particular class
  apred<-mylist$predictions
  par(mfrow=c(1,2))
  plot(test.newdata$LogProbs,apred$ypred,pch=20,col="Black")
  points(test.newdata$LogProbs,apred$ppred,pch=20,col="Red")
  abline(0,1)
  plot(apred$ypred,apred$ppred,pch=20,col="Black")
  abline(0,1)
  expect_gte(cor(apred$ypred,apred$ppred),0.95)
  expect_length(apred$ipred,0)
  expect_length(apred$icpred,0)
})


test_that("Computing individual and population predictions for a new dataset with missing covariates", {
  theo2<-test.newdata[,-c(4)]
  mylist<-predict.newdata(saemix.fit, theo2, type=c("ipred", "ypred", "ppred", "icpred"))
  expect_is(mylist, "list") # tests for particular class
  apred<-mylist$predictions
  par(mfrow=c(1,2))
  plot(test.newdata$LogProbs,apred$ypred,pch=20,col="Black")
  points(test.newdata$LogProbs,apred$ppred,pch=20,col="Red")
  abline(0,1)
  plot(apred$ypred,apred$ppred,pch=20,col="Black")
  abline(0,1)
  expect_gte(cor(apred$ypred,apred$ppred),0.95)
})

###################################################################################
# Using predict to return a vector of predictions

context("Testing predict.newdata for a structural model \n")

### FOR CAT DATA
psiM<-data.frame(id = seq(1,1000,1), alp1=seq(2.308,3,1),alp2 = seq(0.716,1,1),alp3 = seq(0.762,1,1))
logp<-saemixObject["model"]["model"](psiM[,2:4], saemix.fit@data@data[,c("ID")], saemix.fit@data@data[,c("TIME","Y")])

### FOR TTE
psiM<-data.frame(id = seq(1,10,1), lambda=seq(7.01,8,1),beta = seq(1.06,2,1))
logp<-saemixObject["model"]["model"](psiM[,2:3], saemix.fit@data@data[,c("id")], saemix.fit@data@data[,c("time","y")])


test_that("Computing population predictions for a new dataset with a valid structure and individual observations using predict()", {
  fit.pred<-saemix.predict(saemix.fit)
  vec<-predict(saemix.fit)
  expect_equal(vec,fit.pred@results@predictions$ipred)
  expect_length(vec,saemix.fit@data@ntot.obs)
  # expect_gte(cor(vec,saemix.fit@data@data[,saemix.fit@data@name.response]),0.8)
  #For ORD or TTE data we compare vec and the vector of log probabilities from the observed data
  expect_gte(cor(vec,logp),0.8)
  vec<-predict(saemix.fit,type="ypred")
  expect_equal(vec,fit.pred@results@predictions$ypred)
  vec<-predict(saemix.fit,type="icpred")
  expect_equal(vec,fit.pred@results@predictions$icpred)
  vec<-predict(saemix.fit,type="ipred")
  expect_equal(vec,fit.pred@results@predictions$ipred)
})


test_that("Computing population predictions for a new dataset with a valid structure and individual observations using predict()", {
  vec<-predict(saemix.fit,test.newdata)
  expect_gte(cor(vec,test.newdata$LogProbs),0.95)
  vec<-predict(saemix.fit,test.newdata,type="ypred")
  expect_gte(cor(vec,test.newdata$LogProbs),0.95)
  vec<-predict(saemix.fit,test.newdata,type="icpred")
  expect_gte(cor(vec,test.newdata$LogProbs),0.95)
})

  
###################################################################################
