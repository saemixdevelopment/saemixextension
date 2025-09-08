context("Creating SaemixErrorModel objects via Class constructor")

test_that("Error model - implicit", {
  x<-new(Class="SaemixErrorModel")
  expect_equal(x@name,"constant")
  expect_equal(x@npar,1)
  expect_equal(x@model(1,1),1)
  print(x)
})

test_that("Error model - explicit", {
  x<-new(Class="SaemixErrorModel", name="constant")
  expect_equal(x@name,"constant")
  expect_equal(x@npar,1)
  expect_equal(x@model(1,1),1)
  print(x)
})

test_that("Error model - explicit", {
  x<-new(Class="SaemixErrorModel", name="combined1")
  expect_equal(x@name,"combined1")
  expect_equal(x@npar,2)
  expect_equal(x@model(1,c(1,0.2)),1.2)
  expect_equal(x@model(1,x@param),2)
  print(x)
})

