###################################################################################

context("Testing valid replacements\n")

test_that("Newdata object with the same structure as the original data - 1 covariate", {
  saemixObject<-theo.fit1
  xtim<-seq(0,24,2)
  nsuj<-5
  xwei<-seq(50,90,length.out = nsuj)
  xsex<-rep(c("F","M"),length.out=nsuj)
  xdose<-seq(280,320,length.out=nsuj)
  theo.newdat<-data.frame(Id=rep(1:nsuj,each=length(xtim)),Time=rep(xtim,nsuj),Dose=rep(xdose,each=length(xtim)), Weight=rep(xwei,each=length(xtim)),Sex=rep(xsex,each=length(xtim)))
  psiM<-data.frame(ka=seq(1.6,2,0.1),V=seq(34,30),CL=c(2,2.5,2,2.5,2))
  fpred<-saemixObject["model"]["model"](psiM, theo.newdat$Id, theo.newdat[,c("Dose","Time")])
  theo.newdat$Concentration<-fpred+rnorm(length(fpred),sd=0.5)
  
  expect_warning(newObj<-replaceData.saemixObject(saemixObject,theo.newdat))
  expect_is(newObj, "SaemixObject") # tests for particular class
  expect_length(newObj@results@cond.mean.phi ,0)
  expect_length(newObj@results@phi.samp ,0)
  expect_equal(newObj@data@name.covariates,c("Weight"))
})

test_that("Newdata object with the same structure as the original data - 2 covariates", {
  saemixObject<-theo.fit3
  xtim<-seq(0,24,2)
  nsuj<-5
  xwei<-seq(50,90,length.out = nsuj)
  xsex<-rep(c("F","M"),length.out=nsuj)
  xdose<-seq(280,320,length.out=nsuj)
  theo.newdat<-data.frame(Id=rep(1:nsuj,each=length(xtim)),Time=rep(xtim,nsuj),Dose=rep(xdose,each=length(xtim)), Weight=rep(xwei,each=length(xtim)),Sex=rep(xsex,each=length(xtim)))
  psiM<-data.frame(ka=seq(1.6,2,0.1),V=seq(34,30),CL=c(2,2.5,2,2.5,2))
  fpred<-saemixObject["model"]["model"](psiM, theo.newdat$Id, theo.newdat[,c("Dose","Time")])
  theo.newdat$Concentration<-fpred+rnorm(length(fpred),sd=0.5)
  
  expect_warning(newObj<-replaceData.saemixObject(saemixObject,theo.newdat))
  expect_is(newObj, "SaemixObject") # tests for particular class
  expect_length(newObj@results@cond.mean.phi ,0)
  expect_length(newObj@results@phi.samp ,0)
  expect_equal(newObj@data@name.covariates,c("Weight","Sex"))
})

context("Testing invalid replacements\n")

test_that("Newdata object without group", {
  saemixObject<-theo.fit1
  xtim<-seq(0,24,2)
  nsuj<-5
  xwei<-seq(50,90,length.out = nsuj)
  xsex<-rep(c("F","M"),length.out=nsuj)
  xdose<-seq(280,320,length.out=nsuj)
  theo.newdat<-data.frame(Time=rep(xtim,nsuj),Dose=rep(xdose,each=length(xtim)), Weight=rep(xwei,each=length(xtim)),Sex=rep(xsex,each=length(xtim)))
  newObj<-replaceData.saemixObject(saemixObject,theo.newdat)
  expect_null(newObj)
})


test_that("Newdata object without predictors", {
  saemixObject<-theo.fit1
  xtim<-seq(0,24,2)
  nsuj<-5
  xwei<-seq(50,90,length.out = nsuj)
  xsex<-rep(c("F","M"),length.out=nsuj)
  xdose<-seq(280,320,length.out=nsuj)
  theo.newdat<-data.frame(Id=rep(1:nsuj,each=length(xtim)),Dose=rep(xdose,each=length(xtim)), Weight=rep(xwei,each=length(xtim)),Sex=rep(xsex,each=length(xtim)))
  newObj<-replaceData.saemixObject(saemixObject,theo.newdat)
  expect_null(newObj)
})

context("Testing replacement objects with missing items\n")

test_that("Newdata object without a response", {
  saemixObject<-theo.fit3
  xtim<-seq(0,24,2)
  nsuj<-5
  xwei<-seq(50,90,length.out = nsuj)
  xsex<-rep(c("F","M"),length.out=nsuj)
  xdose<-seq(280,320,length.out=nsuj)
  theo.newdat<-data.frame(Id=rep(1:nsuj,each=length(xtim)),Time=rep(xtim,nsuj),Dose=rep(xdose,each=length(xtim)), Weight=rep(xwei,each=length(xtim)),Sex=rep(xsex,each=length(xtim)))
  expect_warning(newObj<-replaceData.saemixObject(saemixObject,theo.newdat))
  expect_is(newObj, "SaemixObject") # tests for particular class
  expect_length(newObj@results@cond.mean.phi ,0)
  expect_length(newObj@results@phi.samp ,0)
  expect_equal(newObj@data@name.covariates,c("Weight","Sex"))
  expect_equal(sum(!is.na(newObj@data@data$Concentration)),0)
})

test_that("Newdata object without Weight", {
  saemixObject<-theo.fit3
  xtim<-seq(0,24,2)
  nsuj<-5
  xwei<-seq(50,90,length.out = nsuj)
  xsex<-rep(c("F","M"),length.out=nsuj)
  xdose<-seq(280,320,length.out=nsuj)
  theo.newdat<-data.frame(Id=rep(1:nsuj,each=length(xtim)),Time=rep(xtim,nsuj),Dose=rep(xdose,each=length(xtim)), Sex=rep(xsex,each=length(xtim)))

  expect_warning(newObj<-replaceData.saemixObject(saemixObject,theo.newdat))
  expect_is(newObj, "SaemixObject") # tests for particular class
  expect_length(newObj@results@cond.mean.phi ,0)
  expect_length(newObj@results@phi.samp ,0)
  expect_equal(newObj@data@name.covariates,c("Weight","Sex"))
  expect_equal(min(newObj@data@data$Weight),70.5)
  expect_equal(max(newObj@data@data$Weight),70.5)
})

test_that("Newdata object without covariates or response", {
  saemixObject<-theo.fit3
  xtim<-seq(0,24,2)
  nsuj<-5
  xwei<-seq(50,90,length.out = nsuj)
  xsex<-rep(c("F","M"),length.out=nsuj)
  xdose<-seq(280,320,length.out=nsuj)
  theo.newdat<-data.frame(Id=rep(1:nsuj,each=length(xtim)),Time=rep(xtim,nsuj),Dose=rep(xdose,each=length(xtim)))

  expect_warning(newObj<-replaceData.saemixObject(saemixObject,theo.newdat))
  expect_is(newObj, "SaemixObject") # tests for particular class
  expect_length(newObj@results@cond.mean.phi ,0)
  expect_length(newObj@results@phi.samp ,0)
  expect_equal(newObj@data@name.covariates,c("Weight","Sex"))
  expect_equal(min(newObj@data@data$Weight),70.5)
  expect_equal(max(newObj@data@data$Weight),70.5)
  expect_equal(min(newObj@data@data$Sex[1]),1)
})

