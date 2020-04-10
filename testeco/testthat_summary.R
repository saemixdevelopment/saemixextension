context("Testing summary for a structural model \n")

test_that("Testing summary for SaemixData object", {
  saemix.fit<-theo.fit1
  xlist<-summary(saemix.fit@data)
  expect_length(names(xlist),7)
  expect_equal(xlist$N,12)
  expect_equal(xlist$nobs$ntot,120)
  expect_equal(max(xlist$nobs$nind),10)
  expect_equal(min(xlist$nobs$nind),10)
})

test_that("Testing summary for SaemixModel object", {
  saemix.fit<-theo.fit1
  xlist<-summary(saemix.fit@model)
  expect_length(names(xlist),4)
  expect_equal(xlist$model$modeltype,"structural")
  expect_equal(xlist$model$error.model,"constant")
  expect_equal(matrix(c(xlist$covariance.model),dimnames=NULL,ncol=3), diag(3))
})

test_that("Testing summary for SaemixRes object", {
  saemix.fit<-theo.fit1
  expect_warning(xlist<-summary(saemix.fit@results))
  expect_equal(xlist$modeltype,"structural")
  expect_length(names(xlist),8)
  expect_lt(abs(xlist$sigma-0.7433444),0.1)
  expect_lt(abs(xlist$logLik[1,2]-xlist$logLik[2,2]),2)
})

test_that("Testing summary for SaemixObject object", {
  saemix.fit<-theo.fit1
  expect_warning(xlist<-summary(saemix.fit))
  expect_equal(xlist$modeltype,"structural")
  expect_length(names(xlist),9)
  expect_lt(abs(xlist$sigma-0.7433444),0.1)
  expect_lt(abs(xlist$logLik[1,2]-xlist$logLik[2,2]),2)
})
