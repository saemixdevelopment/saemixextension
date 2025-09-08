context("Creating SaemixOutcome objects via Class constructor")

test_that("Continuous outcome - implicit", {
  x<-new(Class="SaemixContinuousOutcome")
  expect_equal(x@name.outcome,"y")
  expect_equal(x@type.outcome,"continuous")
  expect_equal(x@distribution,"normal")
  expect_equal(x@error.model@name,"constant")
  expect_equal(x@error.model@npar,1)
  expect_equal(x@error.model@model(1,1),1)
  expect_equal(x@density,dnorm)
  expect_equal(do.call(x@density,c(list(0), x@density.param)),dnorm(0))
  print(x)
  showall(x)
})


test_that("Continuous outcome - explicit", {
  x<-new(Class="SaemixContinuousOutcome", name.outcome="concentration", error.model="combined1")
  expect_equal(x@name.outcome,"concentration")
  expect_equal(x@type.outcome,"continuous")
  expect_equal(x@distribution,"normal")
  expect_equal(x@error.model@name,"combined1")
  expect_equal(x@error.model@npar,2)
  expect_equal(x@error.model@model(0,c(1,0.5)),1)
  expect_equal(x@error.model@model(1,c(1,0.5)),1.5)
  expect_equal(x@error.model@model(2,c(1,0.5)),2)
  expect_equal(do.call(x@density,c(list(0), x@density.param)),dnorm(0))
  print(x)
  showall(x)
})


test_that("Categorical outcome - implicit", {
  x<-new(Class="SaemixDiscreteOutcome")
  expect_equal(x@name.outcome,"y")
  expect_equal(x@type.outcome,"categorical")
  expect_equal(x@distribution,"binomial")
  expect_equal(x@density,dbinom)
  expect_equal(do.call(x@density,c(list(0), x@density.param)),0.5)
  print(x)
  showall(x)
})

test_that("Categorical outcome - explicit", {
  x<-new(Class="SaemixDiscreteOutcome", name.outcome="bleeding", type.outcome="event", distribution="exponential")
  expect_equal(x@name.outcome,"bleeding")
  expect_equal(x@type.outcome,"event")
  expect_equal(x@distribution,"exponential")
  expect_equal(do.call(x@density,c(list(0), x@density.param)),dexp(0))
  print(x)
  showall(x)
})

test_that("List of outcomes", {
  y1<-new(Class="SaemixContinuousOutcome", name.outcome="concentration", error.model="combined1")
  y2<-new(Class="SaemixDiscreteOutcome", name.outcome="hospitalisation", type.outcome="event", distribution="Weibull")
  y3<-new(Class="SaemixDiscreteOutcome", name.outcome="pain", type.outcome="categorical", distribution="binomial")
  expect_equal(y2@distribution,"weibull")
  expect_equal(do.call(y2@density,c(list(0), y2@density.param)),dweibull(0,shape=1))
  y<-list(y1, y2, y3)
  expect_equal(length(y),3)
  print(y)
})

context("Creating SaemixOutcome objects via saemixOutcome() function")

test_that("Continuous outcome - implicit", {
  x<-saemixOutcome()
  expect_equal(x@name.outcome,"y")
  expect_equal(x@type.outcome,"continuous")
  expect_equal(x@distribution,"normal")
  expect_equal(x@error.model@name,"constant")
  expect_equal(x@error.model@npar,1)
  expect_equal(x@error.model@model(1,1),1)
  expect_equal(x@density,dnorm)
  expect_equal(do.call(x@density,c(list(0), x@density.param)),dnorm(0))
  print(x)
  showall(x)
})

test_that("Continuous outcome - explicit", {
  x<-saemixOutcome(name.outcome="concentration", error.model="combined1", unit="mg/L")
  expect_equal(x@name.outcome,"concentration")
  expect_equal(x@type.outcome,"continuous")
  expect_equal(x@distribution,"normal")
  expect_equal(x@unit,"mg/L")
  expect_equal(x@error.model@name,"combined1")
  expect_equal(x@error.model@npar,2)
  expect_equal(x@error.model@model(0,c(1,0.5)),1)
  expect_equal(x@error.model@model(1,c(1,0.5)),1.5)
  expect_equal(x@error.model@model(2,c(1,0.5)),2)
  expect_equal(do.call(x@density,c(list(0), x@density.param)),dnorm(0))
  print(x)
  showall(x)
})


test_that("Categorical outcome - implicit", {
  x<-saemixOutcome(type.outcome="categorical")
  expect_equal(x@name.outcome,"y")
  expect_equal(x@type.outcome,"categorical")
  expect_equal(x@distribution,"binomial")
  expect_equal(x@density,dbinom)
  expect_equal(do.call(x@density,c(list(0), x@density.param)),0.5)
  print(x)
  showall(x)
})

test_that("Categorical outcome - explicit", {
  x<-saemixOutcome(name.outcome="bleeding", type.outcome="event", distribution="exponential")
  expect_equal(x@name.outcome,"bleeding")
  expect_equal(x@type.outcome,"event")
  expect_equal(x@distribution,"exponential")
  expect_equal(do.call(x@density,c(list(0), x@density.param)),dexp(0))
  print(x)
  showall(x)
})


test_that("List of outcomes", {
  y1<-saemixOutcome(name.outcome="concentration", error.model="combined1")
  y2<-saemixOutcome(name.outcome="hospitalisation", type.outcome="event", distribution="Weibull")
  y3<-saemixOutcome(name.outcome="pain", type.outcome="categorical", distribution="binomial")
  expect_equal(y2@distribution,"weibull")
  expect_equal(do.call(y2@density,c(list(0), y2@density.param)),dweibull(0,shape=1))
  y<-list(y1, y2, y3)
  names(y)<-c(y1@name.outcome, y2@name.outcome, y3@name.outcome)
  expect_equal(length(y),3)
  print(y)
})

test_that("List of outcomes", {
  y<-list(concentration=saemixOutcome(error.model="combined1", unit="mg/L"), 
        hospitalisation=saemixOutcome( type.outcome="event", distribution="Weibull"), 
        pain=saemixOutcome( type.outcome="categorical", distribution="binomial"))
  for(i in 1:length(y)) y[[i]]@name.outcome <- names(y)[i]
  expect_equal(y[[3]]@distribution,"binomial")
  print(y)
})

