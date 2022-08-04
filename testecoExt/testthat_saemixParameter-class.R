# Not sure I'm going this way...

context("Creating SaemixParameter objects via Class constructor")

test_that("Parameter class, all defaults", {
  vd<-new(Class="SaemixParameter")
  expect_equal(vd@mu, 1)
  expect_equal(vd@omega, 1)
  expect_equal(length(vd@covariate), 0)
})

test_that("Parameter with starting value and covariate model", {
  vd<-new(Class="SaemixParameter", mu.start=10, covariate=c(wt=contCov()))
  expect_equal(vd@mu, 10)
  expect_equal(vd@omega, 1)
  expect_equal(length(vd@covariate), 1)
  expect_equal(vd@covariate[[1]]@name, "wt")
  expect_equal(vd@covariate[[1]]@beta, 0)
  expect_equal(vd@covariate[[1]]@beta.fix, 0)
})

test_that("Default lognormal parameter", {
  cl<-lognormalPar()
  expect_equal(cl@transform(c(1,2)),exp(1:2))
  expect_equal(cl@inversetransform(c(1,2)), c(0,log(2)))
  expect_equal(cl@mu, 1)
  expect_equal(cl@mu.fix, 0)
})

test_that("Default logit-normal parameter", {
  imax<-logitPar()
  expect_equal(imax@transform(1),1/(1+exp(-1)))
  expect_equal(imax@inversetransform(0.5),0)
  expect_equal(imax@mu, 0.5)
  expect_equal(imax@mu.fix, 0)
})

test_that("Default probit parameter", {
  f1 <- probitPar()
  expect_equal(f1@mu, 0.5)
  expect_equal(f1@mu.fix, 0)
})


test_that("Parameter with starting value, covariate model and no IIV", {
  vd<-new(Class="SaemixParameter", mu.start=10, covariate=c(wt=contCov()), omega.level=c())
  expect_equal(vd@mu, 10)
  expect_equal(vd@mu.fix, 0)
  expect_equal(length(vd@omega), 0)
  expect_equal(length(vd@omega.level), 0)
  expect_equal(length(vd@omega.fix), 0)
  expect_equal(length(vd@rho), 0)
  expect_equal(vd@covariate[[1]]@name, "wt")
  expect_equal(vd@covariate[[1]]@beta, 0)
  expect_equal(vd@covariate[[1]]@beta.fix, 0)
})


test_that("Parameter with starting value, covariate model and two covariances", {
  vd<-new(Class="SaemixParameter", mu.start=10, covariate=c(wt=contCov()), rho.param=c("cl","ka"))
  expect_equal(vd@mu, 10)
  expect_equal(vd@mu.fix, 0)
  expect_equal(length(vd@rho.param[[1]]), 2)
  expect_equal(length(vd@rho.fix[[1]]), 2)
  expect_equal(vd@rho[[1]], rep(0.5,2))
  expect_equal(vd@covariate[[1]]@name, "wt")
  expect_equal(vd@covariate[[1]]@beta, 0)
  expect_equal(vd@covariate[[1]]@beta.fix, 0)
})


context("Creating SaemixParameter objects via constructor function")

test_that("Parameter class, all defaults", {
  vd<-saemixPar()
  expect_equal(vd@distribution, "lognormal")
  expect_equal(vd@mu, 1)
  expect_equal(vd@omega, 1)
  expect_equal(length(vd@covariate), 0)
})


test_that("List of parameters with different options", {
  lpar <- list(ka=saemixPar(mu.start=2),
               cl=saemixPar(mu.start=20, covariate=c(age=contCov(), sex=binCov())),
               vd=saemixPar(mu.start=10, covariate=c(wt=contCov()), rho.param=c("cl","ka")),
               imax=saemixPar(distribution="logit"))
  expect_equal(lpar[[1]]@distribution, "lognormal")
  expect_equal(lpar[[1]]@mu, 2)
  expect_equal(lpar[[1]]@omega, 1)
  expect_equal(length(lpar[[1]]@covariate), 0)
  expect_equal(lpar[[2]]@distribution, "lognormal")
  expect_equal(lpar[[2]]@mu, 20)
  expect_equal(lpar[[2]]@omega, 1)
  expect_equal(length(lpar[[2]]@covariate), 2)
  expect_equal(lpar[[2]]@rho, list())
  expect_equal(lpar[[3]]@distribution, "lognormal")
  expect_equal(lpar[[3]]@mu, 10)
  expect_equal(lpar[[3]]@omega, 1)
  expect_equal(length(lpar[[3]]@covariate), 1)
  expect_equal(lpar[[3]]@rho, list(c(0.5,0.5)))
})

test_that("Covariates - same covariates for different parameters", {
  lpar <- list(cl=saemixPar(mu.start=20, covariate=c(age=contCov(), sex=binCov())),
               vd=saemixPar(mu.start=10, covariate=c(wt=contCov(), sex=binCov()), rho.param=c("cl","ka")))
  x<-getCovariateModel(lpar)
  expect_equal(length(lpar[[1]]@covariate), 2)
  expect_equal(length(lpar[[2]]@covariate), 2)
})


test_that("Covariates - double covariate definition for one parameter", {
  lpar <- list(cl=saemixPar(mu.start=20, covariate=c(age=contCov(), sex=binCov(name="gender"), gender=binCov(name="gender"))),
               vd=saemixPar(mu.start=10, covariate=c(wt=contCov(), gender=binCov(name="gender")), rho.param=c("cl","ka")))
  x<-removeDuplicateCovDef(lpar)
  expect_equal(length(lpar[[1]]@covariate), 3)
  expect_equal(length(x[[1]]@covariate), 2)
  expect_equal(length(lpar[[2]]@covariate), 2)
})


test_that("Covariates - same covariates for different parameters but with different names", {
  lpar <- list(cl=saemixPar(mu.start=20, covariate=c(age=contCov(), sex=binCov(name="gender"))),
               vd=saemixPar(mu.start=10, covariate=c(wt=contCov(), gender=binCov(name="gender")), rho.param=c("cl","ka")))
  lpar<-removeDuplicateCovDef(lpar)
  x<-getCovariateModel(lpar)
  expect_equal(length(lpar[[1]]@covariate), 2)
  expect_equal(length(x$covariates), 3)
  expect_equal(names(x$covariates), c("age","sex","wt")) # first name found for gender, here in the definition of cl
  expect_equal(sum(x$covariate.model), 4)
  expect_equal(dim(x$covariate.model)[1], 3)
  expect_equal(sum(x$beta.start), 0)
  expect_equal(sum(x$covariate.model.fix), 0)
})

test_that("Covariates - different covariate definitions but same names for different parameters", {
  lpar <- list(cl=saemixPar(mu.start=20, covariate=c(wt=contCov(), sex=binCov())),
               vd=saemixPar(mu.start=10, covariate=c(wt=contCov(transform=log, centering=median), sex=binCov()), rho.param=c("cl","ka")))
  lpar<-removeDuplicateCovDef(lpar)
  x<-getCovariateModel(lpar)
  expect_equal(length(lpar[[1]]@covariate), 2)
  expect_equal(length(lpar[[2]]@covariate), 2)
  expect_equal(sum(x$covariate.model), 4)
  expect_equal(dim(x$covariate.model)[1], 3)
})

test_that("Covariates - different covariate definitions but same names for different parameters, transform and centering given as character strings", {
  lpar <- list(cl=saemixPar(mu.start=20, covariate=c(wt=contCov(), sex=binCov())),
               vd=saemixPar(mu.start=10, covariate=c(wt=contCov(transform="log", centering="median"), sex=binCov()), rho.param=c("cl","ka")))
  lpar<-removeDuplicateCovDef(lpar)
  x<-getCovariateModel(lpar)
  expect_equal(length(lpar[[1]]@covariate), 2)
  expect_equal(length(lpar[[2]]@covariate), 2)
  expect_equal(sum(x$covariate.model), 4)
  expect_equal(dim(x$covariate.model)[1], 3)
})



