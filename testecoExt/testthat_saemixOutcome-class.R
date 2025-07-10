context("Creating SaemixOutcome objects")

test_that("Continuous outcome - implicit", {
  x<-new(Class="SaemixContinuousOutcome")
  expect_equal(x@name.outcome,"y")
  expect_equal(x@type.outcome,"continuous")
  expect_equal(x@distribution,"normal")
  expect_equal(x@error.model,"constant")
  expect_equal(x@error.npar,1)
  expect_equal(x@error.function(1,1),1)
  print(x)
  showall(x)
})


test_that("Continuous outcome - explicit", {
  x<-new(Class="SaemixContinuousOutcome", name.outcome="concentration", error.model="combined1")
  expect_equal(x@name.outcome,"concentration")
  expect_equal(x@type.outcome,"continuous")
  expect_equal(x@distribution,"normal")
  expect_equal(x@error.model,"combined1")
  expect_equal(x@error.npar,1)
  expect_equal(x@error.function(0,c(1,0.5)),1)
  expect_equal(x@error.function(1,c(1,0.5)),1.5)
  expect_equal(x@error.function(2,c(1,0.5)),2)
  print(x)
  showall(x)
})


test_that("Categorical outcome - implicit", {
  x<-new(Class="SaemixDiscreteOutcome")
  expect_equal(x@name.outcome,"y")
  expect_equal(x@type.outcome,"categorical")
  expect_equal(x@distribution,"binomial")
  print(x)
  showall(x)
})

test_that("Categorical outcome - explicit", {
  x<-new(Class="SaemixDiscreteOutcome", name.outcome="bleeding", type.outcome="event", distribution="constantHazard")
  expect_equal(x@name.outcome,"bleeding")
  expect_equal(x@type.outcome,"event")
  expect_equal(x@distribution,"constantHazard")
  print(x)
  showall(x)
})

test_that("List of outcomes", {
  y1<-new(Class="SaemixContinuousOutcome", name.outcome="concentration", error.model="combined1")
  y2<-new(Class="SaemixDiscreteOutcome", name.outcome="hospitalisation", type.outcome="event", distribution="Weibull")
  y3<-new(Class="SaemixDiscreteOutcome", name.outcome="pain", type.outcome="categorical", distribution="binomial")
  y<-list(y1, y2, y3)
  expect_equal(length(y),3)
  print(y)
})

############ Move to specific testthat
context("Creating SaemixVarLevel objects (representing a level of variability)")

test_that("Variability levels", {
  x1<-new(Class="SaemixVarLevel")
  expect_equal(x1@name.level,"iiv")
  expect_equal(x1@variable,"id")
  x2<-new(Class="SaemixVarLevel", name.level="iov", variable="occ")
  expect_equal(x2@name.level,"iov")
  expect_equal(x2@variable,"occ")
})

context("Creating SaemixParameter objects ")
test_that("Parameeters", {
  ka<-new(Class="SaemixParameter", name.par="ka")
  vd<-new(Class="SaemixParameter", name.par="vd", initial=10, covariate="wt")
})


