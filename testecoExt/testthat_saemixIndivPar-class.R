context("Testing SaemixVarPhi class")

datDir <- "/home/eco/work/saemix/saemixextension/data"
test_that("Creating an SaemixVarPhi object for IIV only - no covariates in model", {
  theo.saemix<-read.table(file.path(datDir,"theo.saemix.tab"),header=T,na=".")
  x<-saemixData(name.data=theo.saemix, name.group="Id",name.predictors=c("Dose","Time"), name.response="Concentration", name.covariates=c("Sex","Weight"))
  #   xstr <- createStructData(x)
  phiM<-data.frame(ka=log(rnorm(x@N,mean=1, sd=0.1)), V=log(rnorm(x@N,mean=20, sd=2)), CL=log(rnorm(x@N,mean=2, sd=0.3)))
  Mcovar<-matrix(rep(1,dim(phiM)[1]), ncol=1, dimnames=list(NULL,c("mu")))
  phi1 <- new(Class="SaemixVarPhi", phi=as.matrix(phiM), index=x@data$index, Mcovar=Mcovar)
  expect_equal(phi1@name,"iiv")
  expect_equal(nrow(phi1@phi),x@N)
  expect_equal(length(phi1@index),x@ntot.obs)
})

test_that("Creating SaemixVarPhi objects for IIV and IOV levels, no covariates", {
  theonocov<-read.csv(file.path(datDir,"../","data40","theoSimul_nocov.csv"),header=T,na=".")
  x<-saemixData(name.data=theonocov, name.predictors=c("tim","dose"), name.response="ysim",varlevel=c("id","occ")) 
  psinocov<-read.csv(file.path(datDir,"../","data40","theoSimul_psinocov.csv"),header=T,na=".")
  Mcovar<-matrix(rep(1,x@N), ncol=1, dimnames=list(NULL,c("mu")))
  phinocov<-psinocov
  for(icol in 2:4) phinocov[,icol]<-log(phinocov[,icol])
  phi.iiv <- new(Class="SaemixVarPhi", phi=as.matrix(phinocov[phinocov$occ==1,2:4]), index=x@var.index[[1]], Mcovar=Mcovar)
  expect_equal(dim(phi.iiv@Mcovar), c(x@N,1))
  expect_equal(sum(phi.iiv@Mcovar), x@N)
  phi.iov <- new(Class="SaemixVarPhi", phi=as.matrix(phinocov[,2:4]), index=x@var.index[[1]])
  expect_equal(phi.iov@Mcovar, matrix(nrow=0,ncol=0))
})


test_that("Creating SaemixVarPhi objects for IIV and IOV levels, with covariates", {
  theocov<-read.csv(file.path(datDir,"../","data40","theoSimul_cov.csv"),header=T,na=".")
  x<-saemixData(name.data=theocov, name.predictors=c("tim","dose"), name.response="ysim",varlevel=c("id","occ"), name.covariates = c("sex","wt","crcl","comed"), covariate.varlevel = c("id","id","occ","occ")) 
  psicov<-read.csv(file.path(datDir,"../","data40","theoSimul_psicov.csv"),header=T,na=".")
  Mcovar<-as.matrix(cbind(mu=rep(1,x@N), x@var.covmat$id))
  phinocov<-psinocov
  for(icol in 2:4) phinocov[,icol]<-log(phinocov[,icol])
  phi.iiv <- new(Class="SaemixVarPhi", phi=as.matrix(phinocov[phinocov$occ==1,2:4]), index=x@var.index[[1]], Mcovar=Mcovar)
  expect_equal(dim(phi.iiv@Mcovar), c(x@N,3))
  expect_equal(sum(phi.iiv@Mcovar[,"mu"]), x@N)
  expect_equal(sum(phi.iiv@Mcovar[,"sex"]), 50)
  phi.iov <- new(Class="SaemixVarPhi", phi=as.matrix(phinocov[,2:4]), index=x@var.index[[2]], Mcovar=as.matrix(x@var.covmat$occ))
  expect_equal(dim(phi.iov@Mcovar), c(dim(psicov)[1],2))
})

context("Testing SaemixIndivPar class")
