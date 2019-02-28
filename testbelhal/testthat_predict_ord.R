###################################################################################
# Individual and population predictions for a new dataset

context("Checking ORD fit \n")

test_that("Plotting results from a ORD fit", {
  plot(ord.fit)
  vec<-summary(ord.fit)
})

context("Testing predict.newdata for a ORD model \n")

test_that("Computing individual and population predictions for a new dataset with a valid structure and individual observations", {
  test.newdata<-ord.newdata
  mylist<-predict.newdata(ord.fit, ord.newdata, type=c("ipred", "ypred", "ppred", "icpred"))
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
  test.newdata<-ord.newdata
  mylist<-predict.newdata(ord.fit, ord.newdata, type=c("ipred", "ypred"))
  expect_is(mylist, "list") # tests for particular class
  apred<-mylist$predictions
  par(mfrow=c(1,1))
  plot(test.newdata$LogProbs,apred$ipred,pch=20,col="Blue")
  points(test.newdata$LogProbs,apred$ypred,pch=20,col="Red")
  abline(0,1)
  expect_gte(cor(apred$ypred,apred$ipred),0.9)
})


###the response is a predictor so NA
# test_that("Computing individual and population predictions for a new dataset with a valid structure, no individual observations", {
#   missingY<-ord.newdata[,-c(3)]
#   mylist<-predict.newdata(ord.fit, missingY, type=c("ipred", "ypred", "ppred", "icpred"))
#   expect_null(mylist)
# })

test_that("Comparing parameters", {
  test.newdata<-ord.newdata
  mylist<-predict.newdata(ord.fit, ord.newdata, type=c("ipred", "ypred", "ppred", "icpred"))
  expect_is(mylist, "list") # tests for particular class
  param<-mylist$param$population
  par(mfrow=c(2,2))
  for(i in 1:3) {
    plot(ord.psiM[,i],param[,i],main=colnames(psiM)[i],xlab="Simulated",ylab="Estimated")
    abline(0,1)
  }
  param<-mylist$param$map.psi
  par(mfrow=c(2,2))
  for(i in 1:3) {
    plot(ord.psiM[,i],param[,i],main=colnames(psiM)[i],xlab="Simulated",ylab="Estimated")
    abline(0,1)
  }
  expect_gte(cor(ord.psiM[,1],param[,1]),0.9)
})



##################################################################################
##Using predict to return a vector of predictions

context("Testing predict.newdata for a likelihood model \n")

### FOR CAT DATA
psiM<-data.frame(id = seq(1,1000,1), alp1=seq(2.308,3,1),alp2 = seq(0.716,1,1),alp3 = seq(0.762,1,1))
logp<-saemixObject["model"]["model"](psiM[,2:4], saemix.fit@data@data[,c("ID")], saemix.fit@data@data[,c("TIME","Y")])


test_that("Computing individual (default) predictions using predict() for the original dataset", {
  vec<-predict(ord.fit)
  ord.fit<-saemix.predict(ord.fit)
  expect_equal(vec,ord.fit@results@predictions$ipred)
  expect_length(vec,saemix.fit@data@ntot.obs)
  # expect_gte(cor(vec,saemix.fit@data@data[,saemix.fit@data@name.response]),0.8)
  #For ORD or ord data we compare vec and the vector of log probabilities from the observed data
  expect_gte(cor(vec,ord.fit@results@predictions$ypred),0.8)
  vec<-predict(ord.fit,type="ypred")
  expect_equal(vec,ord.fit@results@predictions$ypred)
  vec<-predict(ord.fit,type="icpred")
  expect_equal(vec,ord.fit@results@predictions$icpred)
  vec<-predict(ord.fit,type="ipred")
  expect_equal(vec,ord.fit@results@predictions$ipred)
})



test_that("Computing population predictions for a new dataset with a valid structure and individual observations using predict()", {
  vec<-predict(ord.fit,ord.newdata)
  par(mfrow=c(1,1))  
  plot(ord.newdata$LogProbs, vec,xlab="Individual log-prob", ylab="Predicted log-prob",pch=20)
  abline(0,1)
  expect_gte(cor(vec,test.newdata$LogProbs),0.8)
  vec<-predict(saemix.fit,test.newdata,type="ypred")
  expect_gte(cor(vec,test.newdata$LogProbs),0.8)
  vec<-predict(saemix.fit,test.newdata,type="icpred")
  expect_gte(cor(vec,test.newdata$LogProbs),0.8)
})

  
###################################################################################
