context("Creating SaemixVarModel objects (representing a level of variability)")

test_that("Variability levels", {
  x1<-new(Class="SaemixVarModel", omega.model=diag(1,nrow = 2,ncol=2))
  expect_equal(x1@name.level,"iiv")
  x2<-new(Class="SaemixVarModel", omega.model=diag(1,nrow = 2,ncol=2), name.level="iov")
  expect_equal(x2@name.level,"iov")
})


context("Creating SaemixVarModel objects ")

test_that("Create variability structures from parameters - one parameter", {
  param1<-list(ka=saemixParam(mu.init=c(2), sd.init=c(0.5), covariate=c("age")))
  xvar <- extractVarModel(param1)
  expect_equal(xvar$iiv@param, c(0.25))
})

test_that("Create variability structures from parameters - one level, diagonal, all defaults", {
  param1<-list(ka=saemixParam(),saemixParam("vd"))  
  xvar <- extractVarModel(param1)
  expect_equal(xvar$iiv@param, c(1,1))
  expect_equal(colnames(xvar$iiv@omega.model),c("ka","vd"))
  expect_equal(length(xvar$iiv@idvec.cov),0)
  expect_equal(xvar$iiv@idvec.var,c(1,2))
})

test_that("Create variability structures from parameters - two levels, diagonal, no names given for parameters", {
  param2<-list(saemixParam(varlevel=c("iiv","iov")),saemixParam(),saemixParam())
  xvar <- extractVarModel(param2)
  expect_equal(xvar$iiv@param, c(1,1,1))
})

test_that("Create variability structures from parameters - two levels, correlations", {
  param3<-list(ka=saemixParam(varlevel=c("iiv","iov")),vd=saemixParam(),cl=saemixParam(corr = list(iiv=c("ka","vd"),iov=c("vd"))))
  xvar <- extractVarModel(param3)
  expect_equal(colnames(xvar$iiv@omega.model),c("ka","vd","cl"))
  expect_equal(names(xvar),c("iiv","iov"))
})

test_that("Create variability structures from parameters - two levels, correlations, initial values", {
  param4<-list(ka=saemixParam(sd.init=c(0.8,0.5), varlevel=c("iiv","iov")),vd=saemixParam(sd.init=0.7),  
               cl=saemixParam(name="cl",mu.init=2, varlevel=c("iiv","iov"), sd.init=c(0.6,0.3), corr = list(iiv=c("ka","vd"),iov=c("vd")), covariate=c("wt","sex","age"), covariate.init=c(0.75,0,0), covariate.estim=c("fixed","estimated","estimated"), corr.init=list(iiv=c(-0.5,0.7), iov=0.7)))
  xvar <- extractVarModel(param4)
  expect_equal(names(xvar),c("iiv","iov"))
  expect_equal(sum(xvar$iiv@omega.model),9)
  expect_equal(xvar$iov@param.names,c("Var.ka","Var.cl"))
  expect_equal(grep("Var",xvar$iiv@param.names), xvar$iiv@idvec.var)
})

context("Creating SaemixVarModelHat objects ")

test_that("Create variability structures from parameters - two levels, correlations, initial values", {
  param4<-list(ka=saemixParam(sd.init=c(0.8,0.5), varlevel=c("iiv","iov")),vd=saemixParam(sd.init=0.7),  
               cl=saemixParam(name="cl",mu.init=2, varlevel=c("iiv","iov"), sd.init=c(0.6,0.3), corr = list(iiv=c("ka","vd"),iov=c("vd")), covariate=c("wt","sex","age"), covariate.init=c(0.75,0,0), covariate.estim=c("fixed","estimated","estimated"), corr.init=list(iiv=c(-0.5,0.7), iov=0.7)))
  xvar <- extractVarModel(param4)
  xvar.iiv.hat<-new(Class="SaemixVarModelHat", saemix.varmodel=xvar$iiv)
  for(i in slotNames(xvar$iiv))
    expect_identical(slot(xvar$iiv, i), slot(xvar.iiv.hat,i))
})
