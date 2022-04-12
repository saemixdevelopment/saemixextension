context("Creating unspecified SaemixOutcome objects")

test_that("Outcome for data", {
  x <- saemixOutcome("concentration", unit="mg/L")
  expect_equal(x@name,"concentration")
  expect_equal(x@type,"continuous")
  print(x)
  pain <- saemixOutcome("pain score",type="categorical", unit="(-)")
  resp <- saemixOutcome("response",type="binary")
  expect_equal(pain@name,"pain score")
  expect_equal(resp@name,"response")
})

context("Creating typed SaemixOutcome objects")

test_that("Continuous outcome - implicit", {
  x<-new(Class="SaemixContinuousOutcome")
  expect_equal(x@name,"y")
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
  expect_equal(x@error.fix, c(TRUE,FALSE,TRUE))
  print(x)
  showall(x)
})


test_that("Categorical outcome - implicit binary", {
  x<-try(new(Class="SaemixDiscreteOutcome"))
  expect_type(x, "character")
  x<-new(Class="SaemixDiscreteOutcome", type="binary")
  expect_equal(x@name,"y")
  expect_equal(x@type,"binary")
  expect_equal(x@distribution,"bernouilli")
  showall(x)
})


test_that("Categorical outcome - implicit categorical", {
  x<-new(Class="SaemixDiscreteOutcome", levels=c("none","moderate","severe"))
  expect_equal(x@name,"y")
  expect_equal(x@type,"categorical")
  expect_equal(x@distribution,"categorical")
  expect_equal(as.character(x@levels[2]),"moderate")
  expect_equal(length(x@levels),3)
  print(x)
})

test_that("Categorical outcome - explicit categories", {
  x<-new(Class="SaemixDiscreteOutcome", distribution="categorical", levels=c("none","moderate","severe"))
  expect_equal(x@name,"y")
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
  expect_equal(x@name,"y")
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

####################
context("Testing functions to create outcomes lists")

test_that("Continuous outcome defined through continuousOutcome - predefined models", {
  x1<-continuousOutcome()
  x2<-continuousOutcome(model="combined1")
  x3<-continuousOutcome(model="combined2")
  x4<-continuousOutcome(model="proportional", start=c(0.2))
  x4<-continuousOutcome(model="proportional", start=c(0.2))
  x5<-continuousOutcome(start=c(4))
  x6<-continuousOutcome(start=c(1,4)) # too many starting values
  x7<-continuousOutcome(model="combined1", start=c(2)) # not enough starting values
  x8<-continuousOutcome(model="power", start=c(0.5, 1, 2))
  expect_equal(x1$error.function(c(0,1), x1$start), c(1,1))
  expect_equal(x1$error.function(c(0,1), c(0.5)), c(0.5,0.5))
  expect_equal(x1$error.function(c(0,1), c(0.5,0.5)), c(0.5,0.5))
  expect_equal(x6$start, c(1))
  expect_equal(x2$error.function(c(0,1), x2$start), c(1,1.5))
  expect_equal(x2$error.function(c(0,1), c(0.5,0.5)), c(0.5,1))
  expect_equal(x3$error.function(c(0,1), c(1,1)), c(1,sqrt(2)))
  expect_equal(x4$error.function(c(0,1), x4$start), c(0,0.2))
  expect_equal(x5$error.function(c(0,1), x5$start), c(4,4))
  expect_equal(x7$start, c(2, 0.5))
  expect_equal(x8$error.function(c(0,0.5), x8$start), c(0.5, sqrt(0.5)))
})

test_that("Continuous outcome defined through continuousOutcome - user models", {
  userError<-function(f,ab) {
    g<-cutoff(sqrt((ab[1]+ab[2]*f)^ab[3]))
    return(g)
  }
  x1<-continuousOutcome(model="user")
  x2<-continuousOutcome(model="user", error.function=proportionalErrorModel, start=c(0.2))
  x3<-continuousOutcome(model="user", error.function=combined1ErrorModel)
  x4<-continuousOutcome(model="user", error.function=userError, start=c(1,1,2))
  expect_equal(x1$error.function, constantErrorModel)
  expect_equal(x1$error.function(c(0,1), x1$start), c(1,1))
  expect_equal(x1$error.npar, 1)
  expect_equal(x2$error.function, proportionalErrorModel)
  expect_equal(x2$start, c(0.2))
  expect_equal(x2$error.npar, 1)
  expect_equal(x2$error.function(c(0,1), x2$start), c(0,0.2))
  expect_equal(x3$error.function(c(0,1), x3$start), c(1,1.5))
  expect_equal(x3$error.function, combined1ErrorModel)
  expect_equal(x3$start, c(1,0.5))
  expect_equal(x3$error.npar, 2)
  expect_equal(x4$error.function(c(0,1), x4$start), c(1,2))
  expect_equal(length(x4$start), 3)
  expect_equal(x4$error.npar, 3)
})

test_that("Discrete outcome defined through discreteOutcome", {
  x1<-discreteOutcome()
  x2<-discreteOutcome(type="event")
  x3<-discreteOutcome(type="categorical", levels=c(1:5))
  expect_equal(x1$type, "binary")
  expect_equal(x1$levels, c(0,1))
  expect_equal(x2$type, "event")
  expect_equal(x2$maxEvents, 1)
  expect_equal(x3$type, "categorical")
  expect_equal(length(x3$levels),5)
})

test_that("Continuous outcome object defined through continuousOutcome - predefined models", {
  x1<-continuousOutcome()
  x<-new(Class="SaemixContinuousOutcome", error.model=x1$error.model, error.npar=x1$error.npar, error.function=x1$error.function, error.parameters=x1$start, error.fix=x1$error.fix)
  expect_equal(x@name,"y")
  expect_equal(x@type,"continuous")
  expect_equal(x@distribution,"normal")
  expect_equal(x@error.model,"constant")
  expect_equal(x@error.npar,1)
  expect_equal(x@error.function(1,1),1)
  print(x)
  showall(x)
})

test_that("SaemixContinuousOutcome object defined through continuousOutcome - user models", {
  userError<-function(f,ab) {
    g<-cutoff(sqrt((ab[1]+ab[2]*f)^ab[3]))
    return(g)
  }
  x1<-continuousOutcome("user", error.function=userError, start=c(1,1,2))
  x4<-new(Class="SaemixContinuousOutcome", error.model=x1$error.model, error.npar=x1$error.npar, error.function=x1$error.function, error.parameters=x1$start, error.fix=x1$error.fix)
  expect_equal(x4@error.function(c(0,1), x4@error.parameters), c(1,2))
  expect_equal(length(x4@error.parameters), 3)
  expect_equal(x4@error.npar, 3)
})


test_that("SaemixDiscreteOutcome object defined through discreteOutcome ", {
  x3<-discreteOutcome(type="categorical", levels=c(1:5))
  x<-new(Class="SaemixDiscreteOutcome",  type=x3$type, levels=x3$levels)
  expect_equal(x@type,"categorical")
  expect_equal(x@distribution,"categorical")
  expect(length(x@levels),5)
  x1<-new(Class="SaemixDiscreteOutcome",  distribution=x3$distribution, levels=x3$levels)
  expect_equal(x,x1)
})

test_that("SaemixEventOutcome object defined through discreteOutcome ", {
  x3<-discreteOutcome(type="event", distribution="Weibull")
  x<-new(Class="SaemixEventOutcome", distribution=x3$distribution)
  expect_equal(x@type,"event")
  expect_equal(x@distribution,"Weibull")
  expect(x@maxEvents, 1)
})

test_that("Using createSaemixOutcome function", {
  x<-createSaemixOutcome(continuousOutcome(model="combined1"))
  expect_equal(x@name,"y")
  expect_equal(x@type,"continuous")
  expect_equal(x@distribution,"normal")
  expect_equal(x@error.model,"combined1")
  expect_equal(x@error.npar,2)
  expect_equal(x@error.function(0,c(1,1)),1)
})

