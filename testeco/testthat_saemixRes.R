context("Testing creation of SaemixRes objects\n")

test_that("Creating an empty SaemixRes object with class", {
  xres<-new(Class="SaemixRes")
  expect_is(xres,"SaemixRes")
  expect_equal(xres@status,"empty")
})

test_that("Attempting to extract residuals from an empty SaemixRes object", {
  xres<-new(Class="SaemixRes")
  expect_equal(length(resid.SaemixRes(xres)),0)
  expect_equal(length(fitted.SaemixRes(xres)),0)
  expect_null(vcov(xres))
})



