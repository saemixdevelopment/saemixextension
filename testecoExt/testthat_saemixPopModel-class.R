# covmodel<-matrix(c(1,0,0,0,1,1),nrow=2,byrow=T, dimnames=list(c("sex","wt"), c("ka","vd","cl")))
# parameters<-param4

context("Creating SaemixPopModel objects ")

test_that("Create fixed effect structures from parameters - one parameter", {
  param1<-list(ka=saemixParam())
  xmod <- extractFixedEffectModel(param1)
  expect_equal(xmod$iiv@param, 0)
})


test_that("Create fixed effect structures from parameters - one parameter", {
  param1<-list(ka=saemixParam(mu.init=c(2), sd.init=c(0.5), covariate=c("age")))
  xmod <- extractFixedEffectModel(param1)
  expect_equal(xmod$iiv@param, c(log(param1$ka@mu.init),0))
})

test_that("Create fixed effect structures from parameters - no covariates, one level of variability", {
  param2<-list(ka=saemixParam(mu.init=c(1), sd.init=c(0.8), varlevel=c("iiv","iov")),
               vd=saemixParam(sd.init=0.7),  
               cl=saemixParam(name="cl",mu.init=2, varlevel=c("iiv","iov"), sd.init=c(0.6,0.3), corr = list(iiv=c("ka","vd"),iov=c("vd")))
               )
  xmod <- extractFixedEffectModel(param2)
  print(xmod)
  expect_equal(xmod$iiv@param, c(log(param2$ka@mu.init[1]),log(param2$vd@mu.init),log(param2$cl@mu.init[1])))
  expect_equal(colnames(xmod$iiv@phi.model),c("ka","vd","cl"))
})

test_that("Create fixed effect structures from parameters", {
  param4<-list(ka=saemixParam(mu.init=c(1,3), sd.init=c(0.8,0.5), varlevel=c("iiv","iov")),vd=saemixParam(sd.init=0.7, covariate="wt", covariate.init=c(1), covariate.estim=c("fixed")),  
               cl=saemixParam(name="cl",mu.init=2, varlevel=c("iiv","iov"), sd.init=c(0.6,0.3), corr = list(iiv=c("ka","vd"),iov=c("vd")), covariate=c("wt","sex","age"), covariate.init=c(0.75,0,0), covariate.estim=c("fixed","estimated","estimated"), corr.init=list(iiv=c(-0.5,0.7), iov=0.7), covariate.varlevel=c("iiv","iiv","iov")))
  xmod <- extractFixedEffectModel(param4)
  print(xmod)
  expect_equal(xmod$iiv@param, c(log(param4$ka@mu.init[1]),log(param4$vd@mu.init),1,log(param4$cl@mu.init[1]),0.75,0))
  expect_equal(colnames(xmod$iiv@phi.model),c("ka","vd","cl"))
  expect_equal(rownames(xmod$iiv@phi.model),c("pop","wt","sex"))
  expect_equal(rownames(xmod$iov@phi.model),c("pop","age"))
})

context("Creating SaemixPopModelHat objects ")

test_that("Create fixed effect structures with results from an existing SaemixPopModel object", {
  param4<-list(ka=saemixParam(mu.init=c(1,3), sd.init=c(0.8,0.5), varlevel=c("iiv","iov")),vd=saemixParam(sd.init=0.7, covariate="wt", covariate.init=c(1), covariate.estim=c("fixed")),  
               cl=saemixParam(name="cl",mu.init=2, varlevel=c("iiv","iov"), sd.init=c(0.6,0.3), corr = list(iiv=c("ka","vd"),iov=c("vd")), covariate=c("wt","sex","age"), covariate.init=c(0.75,0,0), covariate.estim=c("fixed","estimated","estimated"), corr.init=list(iiv=c(-0.5,0.7), iov=0.7), covariate.varlevel=c("iiv","iiv","iov")))
  xmod <- extractFixedEffectModel(param4)
  xmod.hat<-new(Class="SaemixPopModelHat", xmod$iiv)
  print(xmod.hat)
  expect_equal(xmod.hat@conf.int$Initial.Value, c(log(param4$ka@mu.init[1]),log(param4$vd@mu.init[1]),1,log(param4$cl@mu.init[1]),0.75,0))
  expect_equal(sum(sum(xmod.hat@conf.int$Status=="fixed")),2)
})

