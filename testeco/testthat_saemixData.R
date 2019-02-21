context("Creating a SaemixData object through saemixData")

test_that("Errors in creating saemixData - no data", {
  x<-saemixData();    # did not test for particular error message (not sure if good practice), but could be changed if desired
#  expect_error(x<-saemixData(),"Error in saemixData: please provide the name of the datafile or dataframe (between quotes)")
  expect_is(x,"character")
  expect_equal(x,"Creation of SaemixData object failed")
})

test_that("Errors in creating saemixData - Wrong group name", {
  x<-saemixData(name.data=file.path(datDir,"theo.saemix.tab"),header=T,na=".", name.group="wrongid",name.predictors=c("Dose","Time"), name.response="Concentration", automatic=FALSE); 
  expect_is(x,"character")
  expect_equal(x,"Creation of SaemixData object failed")
})

test_that("Errors in creating saemixData - Wrong predictor name", {
  x<-saemixData(name.data=file.path(datDir,"theo.saemix.tab"),header=T,na=".", name.predictors="wrongpred", automatic=FALSE); 
  expect_is(x,"character")
  expect_equal(x,"Creation of SaemixData object failed")
})

test_that("Errors in creating saemixData - Wrong response name", {
  x<-saemixData(name.data=file.path(datDir,"theo.saemix.tab"),header=T,na=".", name.response="wrongresp", automatic=FALSE); 
  expect_is(x,"character")
  expect_equal(x,"Creation of SaemixData object failed")
})

test_that("Errors in creating saemixData - Wrong separator", {
  x<-saemixData(name.data=file.path(datDir,"theo.saemix.tab"),header=T,na=".", sep=",", name.predictors=c("Dose","Time"), automatic=FALSE); 
  expect_is(x,"character")
  expect_equal(x,"Creation of SaemixData object failed")
})

test_that("Successful creation of a SaemixData object from data on disk, automatic recognition", {
  x<-saemixData(name.data=file.path(datDir,"PD1.saemix.tab"),header=T,na=".", verbose=FALSE)
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("dose"))
  expect_equal(x@name.group,c("subject"))
  expect_equal(x@name.response,"response")
  expect_equal(length(x@name.covariates),0)
  expect_equal(x@header,TRUE)
  expect_equal(x@sep,"")
  expect_equal(x@na,".")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"dose")
  expect_equal(x@N,100)
  expect_equal(x@ntot.obs,300)
})

test_that("Successful creation of a SaemixData object, column given as numbers", {
  x<-saemixData(name.data=file.path(datDir,"PD1.saemix.tab"),header=T,na=".", name.group=1,name.predictors=2,name.response=3, name.covariates=4,units=list(x="mg",y="-",covariates="-"),verbose=F)
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("dose"))
  expect_equal(x@name.group,c("subject"))
  expect_equal(x@name.response,"response")
  expect_equal(x@name.covariates,"gender")
  expect_equal(x@header,TRUE)
  expect_equal(x@sep,"")
  expect_equal(x@na,".")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"dose")
  expect_equal(x@N,100)
  expect_equal(x@ntot.obs,300)
  expect_equal(sort(unique(x@data$gender)),c(0,1))
})

test_that("Successful creation of a SaemixData object, column given as names", {
  x<-saemixData(name.data=file.path(datDir,"PD1.saemix.tab"),header=T,na=".", name.group=c("subject"),name.predictors=c("dose"),name.response=c("response"), name.covariates=c("gender"),units=list(x="mg",y="-",covariates="-"),verbose=F)
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("dose"))
  expect_equal(x@name.group,c("subject"))
  expect_equal(x@name.response,"response")
  expect_equal(x@name.covariates,"gender")
  expect_equal(x@header,TRUE)
  expect_equal(x@sep,"")
  expect_equal(x@na,".")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"dose")
  expect_equal(x@N,100)
  expect_equal(x@ntot.obs,300)
  expect_equal(sort(unique(x@data$gender)),c(0,1))
})

test_that("Successful creation of a SaemixData object, full specification, theophylline example", {
  x<-saemixData(name.data=file.path(datDir,"theo.saemix.tab"),header=T,na=".", name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),units=list(x="hr",y="mg/L"), name.X="Time",verbose=F)
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("Dose","Time"))
  expect_equal(x@name.group,c("Id"))
  expect_equal(x@name.response,"Concentration")
  expect_equal(x@header,TRUE)
  expect_equal(x@na,".")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"Time")
  expect_equal(x@N,12)
  expect_equal(x@ntot.obs,120)
})

# From a data frame instead of a file on disk (doesn't work within test_that)
theo.saemix<-read.table(file.path(datDir,"theo.saemix.tab"),header=T,na=".")

test_that("Successful creation of a SaemixData object from dataframe object, full specification", {
  x<-saemixData(name.data=theo.saemix, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),units=list(x="hr",y="mg/L"), name.X="Time", verbose=FALSE)
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("Dose","Time"))
  expect_equal(x@name.group,c("Id"))
  expect_equal(x@name.response,"Concentration")
  expect_equal(x@header,TRUE)
  expect_equal(x@na,"NA")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"Time")
  expect_equal(x@N,12)
  expect_equal(x@ntot.obs,120)
})

# Automatic recognition on
test_that("Errors in creating saemixData - Wrong group name, but automatic recognition on", {
  x<-saemixData(name.data=theo.saemix, name.group="wrongid",name.predictors=c("Dose","Time"), name.response="Concentration"); 
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("Dose","Time"))
  expect_equal(x@name.group,c("Id"))
  expect_equal(x@name.response,"Concentration")
  expect_equal(x@header,TRUE)
  expect_equal(x@na,"NA")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"Dose")
  expect_equal(x@N,12)
  expect_equal(x@ntot.obs,120)
})

