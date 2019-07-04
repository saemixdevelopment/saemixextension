context('Check for correct computation of bic.covariate quantities')

test_that("Testing values of nbeta.random and nbeta.fixed", {
  saemix.fit<-theo.fit1
  nbeta.random1<-saemix.fit@results@nbeta.random
  nbeta.fixed1<-saemix.fit@results@nbeta.fixed
  expect_equal(nbeta.random1,4)
  expect_equal(nbeta.fixed1,0)
  saemix.fit<-theo.fit2
  nbeta.random2<-saemix.fit2@results@nbeta.random
  nbeta.fixed2<-saemix.fit2@results@nbeta.fixed
  expect_equal(nbeta.random2,5)
  expect_equal(nbeta.fixed2,0)
  nbeta.random3<-saemix.fit3@results@nbeta.random
  nbeta.fixed3<-saemix.fit3@results@nbeta.fixed
  expect_equal(nbeta.random3,3)
  expect_equal(nbeta.fixed3,1)
  nbeta.random4<-saemix.fit4@results@nbeta.random
  nbeta.fixed4<-saemix.fit4@results@nbeta.fixed
  expect_equal(nbeta.random4,2)
  expect_equal(nbeta.fixed4,2)
})



test_that("Testing for existence of bic.covariate attributes in SaemixRes objects", {
  saemix.fit<-theo.fit1
  bic.cov.lin1<-saemix.fit@results@bic.covariate.lin
  bic.cov.is1<-saemix.fit@results@bic.covariate.is
  bic.cov.gq1<-saemix.fit@results@bic.covariate.gq
  expect_equal(length(bic.cov.lin1),1)
  expect_equal(length(bic.cov.is1),1)
  expect_equal(length(bic.cov.gq1),0)
  saemix.fit2<-theo.fit2
  bic.cov.lin2<-saemix.fit2@results@bic.covariate.lin
  bic.cov.is2<-saemix.fit2@results@bic.covariate.is
  bic.cov.gq2<-saemix.fit2@results@bic.covariate.gq
  expect_equal(length(bic.cov.lin2),1)
  expect_equal(length(bic.cov.is2),1)
  expect_equal(length(bic.cov.gq2),0)
})
