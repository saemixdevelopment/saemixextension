context("Data provided in a dataframe")

test_that("Reading a file with binary observation data", {
  tte.saemix<-read.table(file.path(datDir,"rttellis.csv"),header=T,sep=",",na=".")
  tte.saemix<-tte.saemix[tte.saemix$ytype==2,]
#  expect_error(x<-saemixData(name.data=tte.saemix,header=TRUE,sep=" ",na=NA,
#    name.group=c("id"),name.predictors=c("time"), name.response="y"))
  x<-saemixData(name.data=tte.saemix,header=TRUE,sep=" ",na=NA,
                name.group=c("id"),name.predictors=c("time"), name.response="y")
  expect_equal(x@N,10)
  expect_equal(x@name.group,'id')
  expect_equal(x@name.predictors,'time')
  expect_equal(x@name.response,'y')
  expect_equal(sum(x@data$y),23)
})

test_that("Reading a file with TTE observation data (not sure what this data represents, infinite y)", {
  tte.saemix<-read.table(file.path(datDir,"rttellis.csv"),header=T,sep=",",na=".")
  tte.saemix<-tte.saemix[tte.saemix$ytype==1,]
  #  expect_error(x<-saemixData(name.data=tte.saemix,header=TRUE,sep=" ",na=NA,
  #    name.group=c("id"),name.predictors=c("time"), name.response="y"))
  x<-saemixData(name.data=tte.saemix,header=TRUE,sep=" ",na=NA,
                name.group=c("id"),name.predictors=c("time"), name.response="y")
  expect_equal(x@N,10)
  expect_equal(x@name.group,'id')
  expect_equal(x@name.predictors,'time')
  expect_equal(x@name.response,'y')
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
  x<-saemixData(name.data=file.path(datDir,"rtte1.csv"),header=T, sep=",",na=".", 
    name.group=c("id"),name.predictors=c("time"), name.response="y",verbose=FALSE)
  xrep<-new(Class="SaemixRepData",data=x)
  expect_is(xrep, "SaemixRepData") # tests for particular class
  expect_equal(xrep@N,10)
  expect_equal(xrep@NM,10)
  expect_equal(xrep@nrep,1)
})

test_that("Creating an object of class SaemixSimData", {
  x<-saemixData(name.data=file.path(datDir,"rtte1.csv"),header=T, sep=",",na=".",
   name.group=c("id"),name.predictors=c("time"), name.response="y",verbose=FALSE)
  xsim<-new(Class="SaemixSimData",data=x)
  expect_is(xsim, "SaemixSimData") # tests for particular class
  expect_equal(x@N,10)
  expect_equal(xsim@nsim,0)
  expect_equal(x@name.group,xsim@name.group)
})
