context("Default control options\n")

test_that("Default option list", {
  #  expect_error(x<-saemixModel());
  x<-saemixControl();
  expect_is(x,"list")
  expect_equal(x$map,TRUE)
  expect_equal(x$fim,TRUE)
  expect_equal(x$ll.is,TRUE)
  expect_equal(x$ll.gq,FALSE)
  expect_equal(x$nb.simpred,100)
})


test_that("Option list - changes", {
  #  expect_error(x<-saemixModel());
  x<-saemixControl(nb.chains=5,nbiter.saemix = c(500,300), ipar.lmcmc = 100)
  expect_is(x,"list")
  expect_equal(x$nb.chains,5)
  expect_equal(x$nbiter.saemix,c(500,300))
  expect_equal(x$ipar.lmcmc,100)
})
