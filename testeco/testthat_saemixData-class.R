context("Testing saemixData \n")

test_that("Using saemixData with dataframe", {
  theo.saemix<-read.table(file.path(datDir,"theo.saemix.tab"),header=T,na=".")
  expect_error(x<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L", covariates=c("kg","-")), name.X="Time"))
})

test_that("Using saemixData with file on disk", {
  x<-try(saemixData(name.data=file.path(datDir,"PD1.saemix.tab"),header=T,na=".", name.group="subject",name.predictors="dose",name.response="response", name.covariates="gender",units=list(x="mg",y="-",covariates="-"),verbose=TRUE))
})

context("Environment issue")

test_that("Environment problem with saemixData... looks for a dataframe in the global environment", {
  theo.saemix<-read.table(file.path(datDir,"theo.saemix.tab"),header=T,na=".")
  expect_error(x<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L", covariates=c("kg","-")), name.X="Time"))
})

context("Testing creating a new element of class")

test_that("Error in creating a SaemixData object", {
  expect_error(new(Class="SaemixData"));    # did not test for particular error message (not sure if good practice), but could be changed if desired
})

test_that("Empty SaemixData object, no automatic recognition", {
  x<-new(Class="SaemixData", name.data="mydata",automatic=FALSE);
  expect_equal(x@name.data,"mydata");
  expect_length(x@name.group,0)
  expect_error(validObject(x))
})

test_that("Empty SaemixData object, valid because automatic recognition sets name group/predictor/response to 1/2/3", {
  x<-new(Class="SaemixData", name.data="mydata");
  expect_equal(x@name.data,"mydata");
  expect_equal(x@name.group,'1')
  expect_equal(x@name.predictors,'2')
  expect_equal(x@name.response,'3')
  expect_true(validObject(x))
})


context("Testing creation of SaemixSimxData and SaemixRepData object\n")

test_that("Creating an object of class SaemixRepData", {
  x<-saemixData(name.data=file.path(datDir,"PD1.saemix.tab"),header=T,na=".", name.group="subject",name.predictors="dose",name.response="response", name.covariates="gender",units=list(x="mg",y="-",covariates="-"),verbose=FALSE)
  xrep<-new(Class="SaemixRepData",data=x)
  expect_is(xrep, "SaemixRepData") # tests for particular class
  expect_equal(xrep@N,100)
  expect_equal(xrep@NM,100)
  expect_equal(xrep@nrep,1)
})

test_that("Creating an object of class SaemixSimData", {
  x<-saemixData(name.data=file.path(datDir,"PD1.saemix.tab"),header=T,na=".", name.group="subject",name.predictors="dose",name.response="response", name.covariates="gender",units=list(x="mg",y="-",covariates="-"),verbose=FALSE)
  xsim<-new(Class="SaemixSimData",data=x)
  expect_is(xsim, "SaemixSimData") # tests for particular class
  expect_equal(x@N,100)
  expect_equal(xsim@nsim,0)
  expect_equal(x@name.group,xsim@name.group)
})
