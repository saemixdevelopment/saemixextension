
context("Creating SaemixIndivModel objects via Class constructor")

test_that("Individual statistical model class SaemixIndivModel, all defaults", {
  varlevel1<-new(Class="SaemixVarLevel")
  xindiv <- new(Class="SaemixIndivModel", var.model=varlevel1)
  expect_equal(xindiv@nb.modpar,1)
  expect_equal(xindiv@nb.cov,0)
  expect_equal(sum(xindiv@omega.model.fix),0)
  expect_equal(sum(xindiv@omega.model),1)
})


test_that("Individual statistical model class SaemixIndivModel, specifying covariance matrix", {
  omega1<-vec2mat(c(1,0.9,0,1,0,0.5))
  varlevel1<-new(Class="SaemixVarLevel", variable="id", size=3, omega=omega1, omega.model=matrix(data=1, nrow=3, ncol=3))
  xindiv <- new(Class="SaemixIndivModel", var.model=varlevel1)
  expect_equal(xindiv@nb.modpar,3)
  expect_equal(xindiv@nb.cov,0)
  expect_equal(sum(xindiv@omega.model.fix),0)
  expect_equal(sum(xindiv@omega.model),9)
})

test_that("Individual statistical model class SaemixIndivModel, specifying covariance matrix and covariate model", {
  omega1<-vec2mat(c(1,0.9,0,1,0,0.5))
  omega.model<-matrix(as.integer(omega1>0),ncol=3)
  covmodel1 <- matrix(c(1,0,0,1,0,0),ncol=3)
  varlevel1<-new(Class="SaemixVarLevel", variable="id", size=3, omega=omega1, omega.model=omega.model)
  xindiv <- new(Class="SaemixIndivModel", var.model=varlevel1, covariate.model=covmodel1)
  expect_equal(xindiv@nb.modpar,3)
  expect_equal(sum(xindiv@omega.model),5)
  expect_equal(sum(xindiv@omega.model.fix),0)
  expect_equal(xindiv@nb.cov,2)
  expect_equal(sum(xindiv@covariate.model),2)
})

context("Creating SaemixIndivModel objects using a list of parameters defined through lpar")

test_that("List of parameters used to infer covariance matrix and covariate model, no mismatch", {
  lpar <- list(ka=saemixPar(mu.start=1, covariate=c(trt=catCov(ncat=3))),
               cl=saemixPar(mu.start=20, covariate=c(age=contCov(), sex=binCov(name="gender"))),
               vd=saemixPar(mu.start=10, covariate=c(wt=contCov(), gender=binCov(name="gender")), rho.param=c("cl"), rho=0.5))
  lpar<-removeDuplicateCovDef(lpar)
  varlevel1 <- getVarianceModel(lpar)
  x<-getCovariateModel(lpar)
  xindiv <- new(Class="SaemixIndivModel", var.model=varlevel1, covariate.model=x$covariate.model, covariate.model.fix=x$covariate.model.fix)
  expect_equal(xindiv@nb.modpar,3)
  expect_equal(sum(xindiv@omega.model),5)
  expect_equal(sum(xindiv@omega.model.fix),0)
  expect_equal(xindiv@index.eta,1:3)
  expect_equal(xindiv@index.omega,c(1,4:6))
  expect_equal(xindiv@nb.cov,4)
  expect_equal(xindiv@name.covariates,c("trt","age","sex","wt"))
  expect_equal(sum(xindiv@covariate.model),5)
})


context("Creating SaemixIndivModel objects using a list of parameters defined through lpar, covariate model given separately (NOT RECOMMENDED as no covariate model is associated)")

test_that("List of parameters used to infer covariance matrix and covariate model, no mismatch", {
  lpar <- list(ka=saemixPar(mu.start=1),
               cl=saemixPar(mu.start=20),
               vd=saemixPar(mu.start=10, rho.param=c("cl"), rho=0.5))
  varlevel1 <- getVarianceModel(lpar)
  covmodel1 <- matrix(c(1,0,0,1,0,0),ncol=3)
  
  xindiv <- new(Class="SaemixIndivModel", var.model=varlevel1, covariate.model=covmodel1)
  expect_equal(xindiv@nb.modpar,3)
  expect_equal(sum(xindiv@omega.model),5)
  expect_equal(sum(xindiv@omega.model.fix),0)
  expect_equal(xindiv@index.eta,1:3)
  expect_equal(xindiv@index.omega,c(1,4:6))
  expect_equal(xindiv@nb.cov,2)
  expect_equal(xindiv@name.covariates, character())
})


