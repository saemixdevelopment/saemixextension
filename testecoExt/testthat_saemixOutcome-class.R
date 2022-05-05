context("Creating typed SaemixOutcome objects")

test_that("Continuous outcome - implicit", {
  x<-new(Class="SaemixContinuousOutcome")
  expect_equal(x@name,"")
  expect_equal(x@type,"continuous")
  expect_equal(x@distribution,"normal")
  expect_equal(x@error.model,"constant")
  expect_equal(x@error.npar,1)
  expect_equal(x@error.function(1,1),1)
  print(x)
  showall(x)
})

test_that("Continuous outcome - explicit", {
  x<-new(Class="SaemixContinuousOutcome", name="concentration", error.model="combined1")
  expect_equal(x@name,"concentration")
  expect_equal(x@type,"continuous")
  expect_equal(x@distribution,"normal")
  expect_equal(x@error.model,"combined1")
  expect_equal(x@error.npar,2)
  expect_equal(x@error.function(0,c(1,0.5)),1)
  expect_equal(x@error.function(1,c(1,0.5)),1.5)
  expect_equal(x@error.function(2,c(1,0.5)),2)
  print(x)
  showall(x)
})


test_that("Continuous outcome - explicit, fixing a.concentration and c.concentration", {
  x<-new(Class="SaemixContinuousOutcome", name="concentration", error.model="power", error.fix=c(T,F,T))
  expect_equal(x@name,"concentration")
  expect_equal(x@type,"continuous")
  expect_equal(x@distribution,"normal")
  expect_equal(x@error.model,"power")
  expect_equal(x@error.npar,3)
  expect_equal(x@error.function(0,c(1,0.5,2)),1)
  expect_equal(x@error.function(1,c(1,0.5,2)),sqrt(1.25))
  expect_equal(x@error.function(2,c(1,0.5,2)),sqrt(2))
  expect_equal(x@error.fix, c(1,0,1))
  print(x)
  showall(x)
})


test_that("Categorical outcome - implicit binary", {
  x<-try(new(Class="SaemixDiscreteOutcome"))
  expect_type(x, "character")
  x<-new(Class="SaemixDiscreteOutcome", type="binary")
  expect_equal(x@name,"")
  expect_equal(x@type,"binary")
  expect_equal(x@distribution,"bernouilli")
  showall(x)
})


test_that("Categorical outcome - implicit categorical", {
  x<-new(Class="SaemixDiscreteOutcome", levels=c("none","moderate","severe"))
  expect_equal(x@name,"")
  expect_equal(x@type,"categorical")
  expect_equal(x@distribution,"categorical")
  expect_equal(as.character(x@levels[2]),"moderate")
  expect_equal(length(x@levels),3)
  print(x)
})

test_that("Categorical outcome - explicit categories", {
  x<-new(Class="SaemixDiscreteOutcome", distribution="categorical", levels=c("none","moderate","severe"))
  expect_equal(x@name,"")
  expect_equal(x@type,"categorical")
  expect_equal(x@distribution,"categorical")
  expect_equal(as.character(x@levels[2]),"moderate")
  expect_equal(length(x@levels),3)
  print(x)
  x1<-new(Class="SaemixDiscreteOutcome", type="categorical", levels=c("none","moderate","severe"))
  expect_equal(x1,x)
})


test_that("Categorical outcome - explicit categories", {
  x<-new(Class="SaemixDiscreteOutcome", type="categorical", levels=c("none","moderate","severe"))
  expect_equal(x@name,"")
  expect_equal(x@type,"categorical")
  expect_equal(x@distribution,"categorical")
  expect_equal(as.character(x@levels[2]),"moderate")
  expect_equal(length(x@levels),3)
  print(x)
})


test_that("Event outcome - explicit", {
  x<-new(Class="SaemixEventOutcome", name="bleeding", distribution="constantHazard")
  expect_equal(x@name,"bleeding")
  expect_equal(x@type,"event")
  expect_equal(x@distribution,"constantHazard")
  print(x)
  showall(x)
})

test_that("List of outcomes", {
  y1<-new(Class="SaemixContinuousOutcome", name="concentration", error.model="combined1")
  y2<-new(Class="SaemixEventOutcome", name="hospitalisation", distribution="Weibull")
  y3<-new(Class="SaemixDiscreteOutcome", name="pain", type="categorical", levels=c("none","mild","moderate","severe"))
  y<-list(y1, y2, y3)
  expect_equal(length(y),3)
  expect_is(y1,"SaemixContinuousOutcome")
  expect_is(y2,"SaemixEventOutcome")
  expect_is(y3,"SaemixDiscreteOutcome")
  print(y)
})




