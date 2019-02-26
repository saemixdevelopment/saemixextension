context("Creating a SaemixData object through saemixData")

test_that("Errors in creating saemixData - no data", {
  x<-saemixData();    # did not test for particular error message (not sure if good practice), but could be changed if desired
#  expect_error(x<-saemixData(),"Error in saemixData: please provide the name of the datafile or dataframe (between quotes)")
  expect_is(x,"character")
  expect_equal(x,"Creation of SaemixData object failed")
})

test_that("Errors in creating saemixData - Wrong group name", {
  tte.saemix<-read.table(file.path(datDir,"rttellis.csv"),header=T, sep=",",na=".")
  x<-saemixData(name.data=tte.saemix,header=TRUE,na=".", name.group="wrongid",name.predictors=c("time"), name.response="y", automatic=FALSE)
  expect_is(x,"character")
  expect_equal(x,"Creation of SaemixData object failed")
})

test_that("Errors in creating saemixData - Wrong predictor name", {
  x<-saemixData(name.data=file.path(datDir,"rttellis.csv"),header=T,na=".", name.predictors="wrongpred", automatic=FALSE); 
  expect_is(x,"character")
  expect_equal(x,"Creation of SaemixData object failed")
})

test_that("Errors in creating saemixData - Wrong response name", {
  x<-saemixData(name.data=file.path(datDir,"rttellis.csv"),header=T,na=".", name.response="wrongresp", automatic=FALSE); 
  expect_is(x,"character")
  expect_equal(x,"Creation of SaemixData object failed")
})

test_that("Errors in creating saemixData - Wrong separator", {
  x<-saemixData(name.data=file.path(datDir,"rttellis.csv"),header=T,na=".", sep=" ", name.group=c("id"),name.predictors=c("time"), name.response="y", automatic=FALSE); 
  expect_is(x,"character")
  expect_equal(x,"Creation of SaemixData object failed")
})

test_that("Successful creation of a SaemixData object from data on disk, automatic recognition", {
  x<-saemixData(name.data=file.path(datDir,"rtte1.csv"),header=TRUE,sep=",",na=".",verbose=FALSE)
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("time"))
  expect_equal(x@name.group,c("id"))
  expect_equal(x@name.response,"y")
  expect_equal(length(x@name.covariates),0)
  expect_equal(x@header,TRUE)
  expect_equal(x@sep,",")
  expect_equal(x@na,".")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"time")
  expect_equal(x@N,10)
  expect_equal(x@ntot.obs,1605)
})

test_that("Successful creation of a SaemixData object, column given as numbers", {
  x<-saemixData(name.data=file.path(datDir,"rtte1.csv"),header=TRUE,sep=",",na=".",name.group=1,name.predictors=2,name.response=3,verbose=F)
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("time"))
  expect_equal(x@name.group,c("id"))
  expect_equal(x@name.response,"y")
  expect_equal(x@header,TRUE)
  expect_equal(x@sep,",")
  expect_equal(x@na,".")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"time")
  expect_equal(x@N,10)
  expect_equal(x@ntot.obs,1605)
})

test_that("Successful creation of a SaemixData object, column given as names", {
  x<-saemixData(name.data=file.path(datDir,"rtte1.csv"),header=TRUE,sep=",",na=".",name.group=c("id"),name.predictors=c("time"),name.response=c("y"),verbose=F)
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("time"))
  expect_equal(x@name.group,c("id"))
  expect_equal(x@name.response,"y")
  expect_equal(x@header,TRUE)
  expect_equal(x@sep,",")
  expect_equal(x@na,".")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"time")
  expect_equal(x@N,10)
  expect_equal(x@ntot.obs,1605)
})

# test_that("Successful creation of a SaemixData object, full specification, theophylline example", {
#   x<-saemixData(name.data=file.path(datDir,"rttellis.csv"),header=T,na=".", name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),units=list(x="hr",y="mg/L"), name.X="Time",verbose=F)
#   expect_is(x, "SaemixData") # tests for particular class
#   expect_equal(x@name.predictors,c("Dose","Time"))
#   expect_equal(x@name.group,c("Id"))
#   expect_equal(x@name.response,"Concentration")
#   expect_equal(x@header,TRUE)
#   expect_equal(x@na,".")
#   expect_equal(x@name.mdv,"mdv")
#   expect_equal(x@name.cens,"cens")
#   expect_equal(x@name.occ,"occ")
#   expect_equal(x@name.ytype,"ytype")
#   expect_equal(x@name.X,"Time")
#   expect_equal(x@N,12)
#   expect_equal(x@ntot.obs,120)
# })

# Automatic recognition on
test_that("Errors in creating saemixData - Wrong group name, but automatic recognition on", {
  x<-saemixData(name.data=file.path(datDir,"rtte1.csv"),header=TRUE,sep=",",na=".",name.group="wrongid",name.predictors=c("time"),name.response=c("y"),verbose=F);
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("time"))
  expect_equal(x@name.group,c("id"))
  expect_equal(x@name.response,"y")
  expect_equal(x@header,TRUE)
  expect_equal(x@sep,",")
  expect_equal(x@na,".")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"time")
  expect_equal(x@N,10)
  expect_equal(x@ntot.obs,1605)
})

context("Testing various discrepancies (no warnings as verbose+F)\n")

test_that("Successful creation of a SaemixData object, two predictors, one does not exist", {
  x<-saemixData(name.data=file.path(datDir,"rtte1.csv"),header=TRUE,sep=",",na=".",name.group="wrongid",name.predictors=c("time","wrongpred"),name.response=c("y"),verbose=F);
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("time"))
  expect_equal(x@name.group,c("id"))
  expect_equal(x@name.response,"y")
  expect_equal(x@header,TRUE)
  expect_equal(x@sep,",")
  expect_equal(x@na,".")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"time")
  expect_equal(x@N,10)
  expect_equal(x@ntot.obs,1605)
})


# test_that("Successful creation of a SaemixData object, wrong number of units for covariates", {
#   x<-saemixData(name.data=file.path(datDir,"PD1.saemix.tab"),header=T,na=".", name.covariates=4,units=list(x="mg",y="-",covariates=c("-","-")),verbose=F)
#   expect_is(x, "SaemixData") # tests for particular class
#   expect_equal(x@name.predictors,c("dose"))
#   expect_equal(x@name.group,c("subject"))
#   expect_equal(x@name.response,"response")
#   expect_equal(x@name.covariates,"gender")
#   expect_equal(x@header,TRUE)
#   expect_equal(x@sep,"")
#   expect_equal(x@na,".")
#   expect_equal(x@name.mdv,"mdv")
#   expect_equal(x@name.cens,"cens")
#   expect_equal(x@name.occ,"occ")
#   expect_equal(x@name.ytype,"ytype")
#   expect_equal(x@name.X,"dose")
#   expect_equal(x@N,100)
#   expect_equal(x@ntot.obs,300)
#   expect_equal(sort(unique(x@data$gender)),c(0,1))
# })

# test_that("Successful creation of a SaemixData object, no units for covariates", {
#   x<-saemixData(name.data=file.path(datDir,"PD1.saemix.tab"),header=T,na=".", name.covariates=4,units=list(x="mg",y="-"),verbose=F)
#   expect_is(x, "SaemixData") # tests for particular class
#   expect_equal(x@name.predictors,c("dose"))
#   expect_equal(x@name.group,c("subject"))
#   expect_equal(x@name.response,"response")
#   expect_equal(x@name.covariates,"gender")
#   expect_equal(x@header,TRUE)
#   expect_equal(x@sep,"")
#   expect_equal(x@na,".")
#   expect_equal(x@name.mdv,"mdv")
#   expect_equal(x@name.cens,"cens")
#   expect_equal(x@name.occ,"occ")
#   expect_equal(x@name.ytype,"ytype")
#   expect_equal(x@name.X,"dose")
#   expect_equal(x@N,100)
#   expect_equal(x@ntot.obs,300)
#   expect_equal(sort(unique(x@data$gender)),c(0,1))
# })

# test_that("Successful creation of a SaemixData object, no units for x", {
#   x<-saemixData(name.data=file.path(datDir,"PD1.saemix.tab"),header=T,na=".", name.covariates=4,units=list(y="-"),verbose=F)
#   expect_is(x, "SaemixData") # tests for particular class
#   expect_equal(x@name.predictors,c("dose"))
#   expect_equal(x@name.group,c("subject"))
#   expect_equal(x@name.response,"response")
#   expect_equal(x@name.covariates,"gender")
#   expect_equal(x@header,TRUE)
#   expect_equal(x@sep,"")
#   expect_equal(x@na,".")
#   expect_equal(x@name.mdv,"mdv")
#   expect_equal(x@name.cens,"cens")
#   expect_equal(x@name.occ,"occ")
#   expect_equal(x@name.ytype,"ytype")
#   expect_equal(x@name.X,"dose")
#   expect_equal(x@N,100)
#   expect_equal(x@ntot.obs,300)
#   expect_equal(sort(unique(x@data$gender)),c(0,1))
# })

# test_that("Successful creation of a SaemixData object, no units for x", {
#   expect_warning(x<-saemixData(name.data=file.path(datDir,"PD1.saemix.tab"),header=T,na=".",  name.response="gender",verbose=F))
#   expect_is(x, "SaemixData") # tests for particular class
#   expect_equal(x@name.predictors,c("dose"))
#   expect_equal(x@name.group,c("subject"))
#   expect_equal(x@name.response,"gender")
#   expect_equal(x@header,TRUE)
#   expect_equal(x@sep,"")
#   expect_equal(x@na,".")
#   expect_equal(x@name.mdv,"mdv")
#   expect_equal(x@name.cens,"cens")
#   expect_equal(x@name.occ,"occ")
#   expect_equal(x@name.ytype,"ytype")
#   expect_equal(x@name.X,"dose")
#   expect_equal(x@N,100)
#   expect_equal(x@ntot.obs,300)
#   expect_equal(sort(unique(x@data$gender)),c(0,1))
# })

# test_that("Successful creation of a SaemixData object, one covariate does not exist in the dataset", {
#   x<-saemixData(name.data=file.path(datDir,"rttellis.csv"),header=T,na=".",sep=" ", name.group="Id",name.predictors=c("Time","Dose"), name.response="Concentration",name.covariates=c("Height","Sex"), units=list(x="mg",y="-",covariates=c("kg","-")),verbose=F)
#   expect_is(x, "SaemixData") # tests for particular class
#   expect_equal(x@name.covariates,"Sex")
# })
