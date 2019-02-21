context("Testing summary for a likelihood model \n")

test_that("Testing summary for SaemixData object", {
  saemix.fit<-tte.fit
  xlist<-summary(saemix.fit@data)
  expect_length(names(xlist),5)
  expect_equal(xlist$N,10)
  expect_equal(xlist$nobs$ntot,43)
  expect_equal(max(xlist$nobs$nind),18)
  expect_equal(min(xlist$nobs$nind),2)
  expect_equal(sum(xlist$y),23)
})

test_that("Testing summary for SaemixModel object", {
  saemix.fit<-tte.fit
  xlist<-summary(saemix.fit@model)
  expect_length(names(xlist),4)
  expect_equal(xlist$model$modeltype,"likelihood")
  expect_equal(xlist$model$error.model,"constant")
  expect_equal(matrix(c(xlist$covariance.model),dimnames=NULL,ncol=2), diag(2))
})

test_that("Testing summary for SaemixRes object", {
  saemix.fit<-tte.fit
  xlist<-summary(saemix.fit@results)
  expect_equal(xlist$modeltype,"likelihood")
  expect_length(names(xlist),8)
  expect_equal(xlist$sigma,numeric(0))
})

test_that("Testing summary for SaemixObject object", {
  saemix.fit<-tte.fit
  xlist<-summary(saemix.fit)
  expect_equal(xlist$modeltype,"likelihood")
  expect_length(names(xlist),9)
  expect_equal(xlist$sigma,numeric(0))
})

