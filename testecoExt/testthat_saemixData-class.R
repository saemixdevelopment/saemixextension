context("Creating a SaemixData object through its class")

test_that("Creating an object via its class", {
  x<-new(Class="SaemixData",name.data=file.path(datDir,"theo.saemix.tab"), name.group="wrongid",name.predictors=c("Dose","Time"), name.response="Concentration", name.covariates=c("gender"),units=list(x="mg",y="-",covariates="-"),verbose=F)
  expect_equal(x@outcome[[1]]@name,"Concentration")
  expect_equal(x@outcome[[1]]@type,"continuous")
  expect_equal(x@outcome[[1]]@unit,"")
  expect_equal(x@name.group,"wrongid")
})

### TODO: creator functions, reading data, sanitising, etc...
if(FALSE) {
test_that("Testing readSaemix", {
  x<-new(Class="SaemixData",name.data=file.path(datDir,"theo.saemix.tab"), header=T,na=".", name.group="wrongid",name.predictors=c("Dose","Time"), name.response="Concentration", name.covariates=c("gender"),units=list(x="mg",y="-",covariates="-"),verbose=F)
  x1<-readSaemix(x)
  verbose<-TRUE
  x<-x1
  if(length(unique(x@data$ytype))==1) { # Only one response
    if(length(unique(x@data[,x@name.response]))>2) {
      x@outcome<-c(y1="continuous")
      names(x@outcome)<-x@name.response
      if(verbose) message("Assuming response is continuous\n")
    } else {
      x@outcome<-c(y1="categorical")
      names(x@outcome)<-x@name.response
      if(verbose) message("Only two modalities in the response, assuming the response is binary\n")
    }
  }
  expect_equal(x@outcome,c(Concentration="continuous"))
})


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

# Multiple responses
# From a data frame instead of a file on disk (doesn't work within test_that)
pkpd.saemix<-read.table(file.path(datDir,"../","data40","warfarinPKPD.tab"),header=T,na=".")
pkrtte.saemix<-read.table(file.path(datDir,"../","data40","pkRTTE.tab"),header=T,na=".")
pkcount.saemix<-read.table(file.path(datDir,"../","data40","jointCount.tab"),header=T,na=".")
pkcat.saemix<-read.table(file.path(datDir,"../","data40","warfarinCatPD.tab"),header=T,na=".")

context("Multiple responses")

test_that("Deriving outcomes automatically, no outcome given or types only", {
  x1<-saemixData(name.data=pkpd.saemix, header=T,na=".", name.group=c("id"), name.predictors=c("time","amt"), name.response=c("dv"), name.ytype = "dvid", name.covariates=c("sex", "wt", "age"), units=list(x="hr",y="mg/L"), verbose=TRUE)
  expect_equal(x1@outcome,c(y1="continuous", y2="continuous"))
  x2<-saemixData(name.data=pkpd.saemix, header=T,na=".", name.group=c("id"), name.predictors=c("time","amt"), name.response=c("dv"), name.ytype = "dvid", name.covariates=c("sex", "wt", "age"), units=list(x="hr",y="mg/L"), verbose=TRUE, outcome=c("continuous", "categorical"))
  expect_equal(x2@outcome,c(y1="continuous", y2="categorical"))
})

test_that("Using explicit outcomes, possibly incomplete or oversized", {
  x1<-saemixData(name.data=pkpd.saemix, header=T,na=".", name.group=c("id"), name.predictors=c("time","amt"), name.response=c("dv"), name.ytype = "dvid", name.covariates=c("sex", "wt", "age"), units=list(x="hr",y="mg/L"), verbose=TRUE, outcome=c(conc="continuous", PCA="continuous"))
  expect_equal(x1@outcome,c(conc="continuous", PCA="continuous"))
  x2<-saemixData(name.data=pkpd.saemix, header=T,na=".", name.group=c("id"), name.predictors=c("time","amt"), name.response=c("dv"), name.ytype = "dvid", name.covariates=c("sex", "wt", "age"), units=list(x="hr",y="mg/L"), verbose=TRUE, outcome=c(conc="continuous", "categorical"))
  expect_equal(x2@outcome,c(conc="continuous",y2="categorical"))
  x2<-saemixData(name.data=pkpd.saemix, header=T,na=".", name.group=c("id"), name.predictors=c("time","amt"), name.response=c("dv"), name.ytype = "dvid", name.covariates=c("sex", "wt", "age"), units=list(x="hr",y="mg/L"), verbose=TRUE, outcome=c(conc="continuous"))
  expect_equal(x2@outcome,c(conc="continuous",y2="continuous"))
  x1<-saemixData(name.data=pkpd.saemix, header=T,na=".", name.group=c("id"), name.predictors=c("time","amt"), name.response=c("dv"), name.ytype = "dvid", name.covariates=c("sex", "wt", "age"), units=list(x="hr",y="mg/L"), verbose=TRUE, outcome=c(conc="continuous", PCA="continuous", bleeding="event"))
  expect_equal(x1@outcome,c(conc="continuous", PCA="continuous"))
})

context("Multiple responses, categorical and event")

test_that("Using explicit outcomes", {
  x1<-saemixData(name.data=pkrtte.saemix, header=T,na=".", name.group=c("id"), name.predictors=c("time","dose","ii"), name.response=c("y"), name.ytype = "ytype", units=list(x="hr",y="mg/L"), verbose=TRUE , outcome=c(conc="continuous", event="event"))
  expect_equal(x1@outcome,c(conc="continuous", event="event"))
  expect_equal(length(unique(x1@data$ytype)), 2)
  
  x2<-saemixData(name.data=pkcount.saemix, header=T,na=".", name.group=c("id"), name.predictors=c("time","dose","ii","ndoses"), name.response=c("y"), name.ytype = "ytype", units=list(x="hr",y="mg/L"), verbose=TRUE , outcome=c(conc="continuous", counts="categorical"))
  expect_equal(x2@outcome,c(conc="continuous", counts="categorical"))
  expect_equal(length(unique(x2@data$ytype)), 2)
  
  x3<-saemixData(name.data=pkcat.saemix, header=T,na=".", name.group=c("id"), name.predictors=c("time","amt"), name.response=c("dv"), name.ytype = "dvid", name.covariates=c("sex", "wt", "age"), units=list(x="hr",y="mg/L"), verbose=TRUE , outcome=c(conc="continuous", catPCA="categorical"))
  expect_equal(x3@outcome,c(conc="continuous", catPCA="categorical"))
  expect_equal(length(unique(x3@data$ytype)), 2)
})

}
