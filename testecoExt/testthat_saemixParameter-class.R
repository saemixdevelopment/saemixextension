# Not sure I'm going this way...

context("Creating SaemixParameter objects ")

test_that("Parameters", {
  ka<-new(Class="SaemixParameter", name.par="ka")
  vd<-new(Class="SaemixParameter", name.par="vd", initial=10, covariate="wt")
})
