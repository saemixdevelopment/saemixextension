context("Creating SaemixCovariate objects")

test_that("Creating covariates", {
  cov1<-saemixCov()
  cov2<-saemixCov(name="cov2")
  age<-saemixCov("age", unit="yr")
  comorb1<-saemixCov(name="diabetes", type="binary")
  expect_equal(cov1@name, "")
  expect_equal(cov2@name, "cov2")
  expect_equal(age@type, "continuous")
  expect_equal(comorb1@type, "binary")
})

context("Creating SaemixCovariate objects and applying transformation")

test_that("Transforming continuous covariate ", {
  weight<-rnorm(10, mean=60, sd=10)
  cov.wt<-new(Class="SaemixContinuousCovariate", name="Weight")
  cov.wt2<-new(Class="SaemixContinuousCovariate", name="Weight", transform.function=log, centering.function=mean)
  cov.wt3<-new(Class="SaemixContinuousCovariate", name="Weight", transform.function=log, centering.function=median)
  cov.wt4<-new(Class="SaemixContinuousCovariate", name="Weight", transform.function=log, centering.value=60)
  expect_equal(unique(transformCovariate(cov.wt, weight)-weight), 0)
  expect_equal(unique(transformCovariate(cov.wt2, weight) - log(weight/mean(weight))), 0)
  expect_equal(unique(transformCovariate(cov.wt3, weight) - log(weight/median(weight))), 0)
  expect_equal(unique(transformCovariate(cov.wt4, weight) - log(weight/60)), 0)
})

test_that("Transforming binary covariate - from 0/1", {
  gender<-sample(c("F","M"), 20, replace=TRUE)
  gender<-as.integer(gender=="M")
  cov.gender<-new(Class="SaemixDiscreteCovariate", name="gender")
  gender.trans<-transformCovariate(cov.gender, gender)
  expect_equal(sum(gender.trans), sum(gender))
})


test_that("Transforming binary covariate - from character", {
  gender<-sample(c("F","M"), 20, replace=TRUE)
  cov.gender<-new(Class="SaemixDiscreteCovariate", name="gender", reference="F")
  gender.trans<-transformCovariate(cov.gender, gender)
  expect_equal(sum(gender.trans), sum(gender=="M"))
})

test_that("Transforming binary covariate - from factor", {
  gender<-sample(c("F","M"), 10, replace=TRUE)
  genderF<-factor(gender, labels=c("F","M"))
  cov.gender<-new(Class="SaemixDiscreteCovariate", name="gender", reference="F")
  gender.trans<-transformCovariate(cov.gender, genderF)
  expect_equal(sum(gender.trans), sum(genderF=="M"))
})


test_that("Transforming categorical covariate - character", {
  asthma<-sample(c("none","mild","severe"), 20, replace=TRUE)
  asthmaF<-factor(asthma, labels=c("none","mild","severe"))
  cov<-new(Class="SaemixDiscreteCovariate", name="asthma", type="categorical",reference="none")
  asthmaDum<-transformCovariate(cov, asthmaF)
  expect_equal(sum(asthmaDum[,1]), sum(asthmaF=="mild"))
  expect_equal(sum(asthmaDum[,2]), sum(asthmaF=="severe"))
})

test_that("Transforming categorical covariate to binary by regrouping 2 categories to one", {
  asthma<-sample(c("none","mild","severe"), 20, replace=TRUE)
  asthmaF<-factor(asthma, labels=c("none","mild","severe"))
  covBin<-new(Class="SaemixDiscreteCovariate", name="asthma", type="categorical",reference="none",groups=list(no=c("none","mild"), yes="severe"))
  asthmaBin<-transformCovariate(covBin, asthmaF)
  expect_equal(sum(asthmaBin), sum(asthmaF=="severe"))
})


test_that("Transforming categorical covariate to binary by regrouping 2 categories to one", {
  asthma<-sample(c("none","mild","severe"), 20, replace=TRUE)
  covBin<-new(Class="SaemixDiscreteCovariate", name="asthma", type="categorical",reference="none",groups=list(no=c("none","mild"), yes="severe"))
  asthmaBin<-transformCovariate(covBin, asthma)
  expect_equal(sum(asthmaBin), sum(asthma=="severe"))
})

test_that("Transforming categorical covariate - from integer scores", {
  score<-sample(c(1:5), 20, replace=TRUE)
  cov<-new(Class="SaemixDiscreteCovariate", name="score", type="categorical",groups=list(c(1,2),c(3,4),c(5)))
  score3<-transformCovariate(cov, score)
  expect_equal(sum(score3[,1]), sum(score==3)+sum(score==4))
  expect_equal(sum(score3[,2]), sum(score==5))
})
