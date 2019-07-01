context("Creating a SaemixSIR object through saemixSIR")

test_that("Creating an empty SaemixSIR object with class", {
  x<-new(Class="SaemixSIR")
  expect_is(x,"SaemixSIR")
})

test_that("Creation of SaemixSIR without SaemixObject", {
  x<-saemixSIR()
  expect_is(x,"character")
  expect_equal(x,"Creation of SaemixSIR object has failed")
})

test_that("Successfull creation of SaemixSIR object with incorrect parameters", {
  x <- new(Class='SaemixSIR',SaemixObject = sim.pd20.saemix, M=2, m=30, cov.mat = matrix(ncol=2, nrow=1), optionll='lkj', est.mu=c(1,2))
  expect_is(x, "SaemixSIR")
  expect_equal(x@optionll, 'linearisation')
  expect_equal(x@M, 5000)
  expect_equal(x@m, 1000)
  cm <- solve(sim.pd20.saemix@results@fim)
  expect_equal(x@cov.mat, cm)
  em <- estpar.vector(sim.pd20.saemix)
  expect_equal(x@est.mu, em)
})



