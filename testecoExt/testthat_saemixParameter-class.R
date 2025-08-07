context("Creating SaemixParameter objects through the class constructor")

test_that("Defining parameters with different distributions", {
  ka<-new(Class="SaemixParameter")
  vd<-new(Class="SaemixParameter", mu.init=10)
  biodisp<-new(Class="SaemixParameter", distribution="logit", mu.init=0.8)
  expect_equal(ka@transform,log)
  expect_equal(ka@transform(ka@mu.init),0)
  expect_equal(vd@transform,log)
  expect_equal(vd@transform(vd@mu.init),log(10))
  expect_equal(biodisp@transform(biodisp@mu.init),log(0.8/.2))
  expect_equal(biodisp@invtransform(0),0.5)
})

test_that("Defining parameters with covariate effects", {
  cl<-new(Class="SaemixParameter", mu.init=2, covariate=c("wt","sex","age"), covariate.init=c(0.75))
  vd<-new(Class="SaemixParameter", mu.init=10, covariate="wt", varlevel=c("id","iov"))
  biodisp<-new(Class="SaemixParameter", distribution="logit", mu.init=0.8, covariate=c("comed",'age'), covariate.init=c(0.5,-0.02))
  expect_equal(vd@covariate,"wt")
  expect_equal(vd@covariate.init,c(0))
  expect_equal(vd@varlevel,c("id","iov"))
  # ToDo: add check for variability levels but need to decide on format (list of varlevels or simply vector of grouping levels)
  expect_equal(biodisp@covariate,c("comed","age"))
  expect_equal(biodisp@covariate.init,c(0.5,-0.02))
  expect_equal(biodisp@transform(biodisp@mu.init),log(0.8/.2))
  expect_equal(biodisp@invtransform(0),0.5)
  expect_equal(cl@covariate.init,c(0.75,0,0))
})

test_that("Defining parameters without IIV", {
  ka<-new(Class="SaemixParameter", varlevel=c())
  vd<-new(Class="SaemixParameter", mu.init=10, varlevel=c(),covariate=c("comed",'age'), covariate.init=c(0.5,-0.02))
  expect_equal(ka@varlevel,"iiv")
  expect_equal(ka@transform(ka@mu.init),0)
  expect_equal(vd@transform,log)
  expect_equal(vd@transform(vd@mu.init),log(10))
  expect_equal(vd@varlevel,"iiv")
  expect_equal(vd@covariate.varlevel,rep("id",2))
})

context("Auxiliary functions to match names") 

test_that("Testing matching levels in list to a vector", {
  varlev2 <-c(iiv="iiv",iov="iov")
  listcov <- list(iov=c("ka","vd"), iiv=c("cl","ka"))
  xvec <- checkMatchingVectorList(varlev2, listcov, name.vec1="varlevels", name.list2 = "alist", defaultValue = c(), type="character")
  expect_equal(names(xvec$list),names(varlev2))
  expect_equal(xvec$logmsg,"")
})

test_that("Testing matching levels in list to a vector - order or number of levels different", {
  varlev2 <-c(iiv="iiv",iov="iov")
  listcov <- list(iov=c("ka","vd"), iiv=c("cl","ka"))
  xvec <- checkMatchingVectorList(varlev2, listcov, name.vec1="varlevels", name.list2 = "alist", defaultValue = c(), type="character")
  expect_equal(names(xvec$list),names(varlev2))
  expect_equal(xvec$logmsg,"")
  varlev2 <-c(iiv="iiv")
  xvec <- checkMatchingVectorList(varlev2, listcov, name.vec1="varlevels", name.list2 = "alist", defaultValue = c(), type="character")
  expect_equal(xvec$list,list(iiv=listcov$iiv))
})

test_that("Testing matching levels in list to a vector - mismatched levels", {
  varlev2 <-c(iiv="iiv",iov="iov")
  listcov <- list(iov=c("ka","vd"), iab=c("cl","ka"))
  xvec <- checkMatchingVectorList(varlev2, listcov, name.vec1="varlevels", name.list2 = "alist", defaultValue = c(), type="character")
  expect_equal(xvec$list,list(iov=listcov$iov))
})

# ToDo: same tests for list matching

context("Creating SaemixParameter objects through the saemixParam function") 

# ToDo
## Add parameters with covariates on different variability levels

test_that("Defining a default parameter", {
  ka<-saemixParam()
  expect_equal(ka@mu.estim,"estimated")
  expect_equal(ka@varlevel,"iiv")
  expect_equal(ka@transform,log)
  expect_equal(ka@transform(ka@mu.init),0)
  print(ka)
  showall(ka)
})

test_that("Defining parameters with covariate effects", {
  cl<-saemixParam(name="cl",mu.init=2, covariate=c("wt","sex","age"), covariate.init=c(0.75), covariate.estim=c("fixed"))
  expect_equal(cl@covariate.init,c(0.75,0,0))
  expect_equal(cl@covariate.estim,c("fixed","estimated","estimated"))
  print(cl)
  showall(cl)
})

test_that("Defining parameters with two levels of variability and correlations", {
  cl<-saemixParam(name="cl",mu.init=2, varlevel=c("iiv","iov"), sd.init=c(1,0.01),  corr = list(iiv=c("ka","V"), iov=c("V")), covariate=c("wt","sex","age"), covariate.init=c(0.75), covariate.estim=c("fixed"), corr.init=list(iiv=c(-0.5,0.7), iov=0.7))
  expect_equal(cl@covariate.init,c(0.75,0,0))
  expect_equal(cl@covariate.estim,c("fixed","estimated","estimated"))
  print(cl)
  showall(cl)
})



