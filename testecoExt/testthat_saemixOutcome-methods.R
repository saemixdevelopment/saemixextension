context("Creating unspecified SaemixOutcome objects")

test_that("Outcome for data", {
  x <- saemixDataOutcome("concentration", unit="mg/L")
  expect_equal(x@name,"concentration")
  expect_equal(x@type,"continuous")
  print(x)
  pain <- saemixDataOutcome("pain score",type="categorical", unit="(-)")
  resp <- saemixDataOutcome("response",type="binary")
  expect_equal(pain@name,"pain score")
  expect_equal(resp@name,"response")
})


context("Testing continuousOutcome() and discreteOutcome() functions")

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
  expect_equal(x@name,"")
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

context("Testing createSaemixOutcome() function")

test_that("Using createSaemixOutcome function", {
  x<-createSaemixOutcome(continuousOutcome(model="combined1"), name="y")
  expect_equal(x@name,"y")
  expect_equal(x@type,"continuous")
  expect_equal(x@distribution,"normal")
  expect_equal(x@error.model,"combined1")
  expect_equal(x@error.npar,2)
  expect_equal(x@error.function(0,c(1,1)),1)
  expect_equal(x@error.nameparameters,c("a.y","b.y"))
})

context("Creating full SaemixOutcome objects for use in SaemixModel using saemixOutcome with arguments")

test_that("Outcome for model with defaults and units", {
  x <- saemixOutcome("concentration", unit="mg/L")
  expect_equal(x@name,"concentration")
  expect_equal(x@type,"continuous")
  expect_equal(x@error.model,"constant")
  expect_equal(x@error.npar,1)
  expect_equal(x@error.parameters,1)
  pain <- saemixOutcome("pain score",type="categorical", unit="(-)")
  resp <- saemixOutcome("response",type="binary")
  expect_equal(pain@name,"pain score")
  expect_equal(resp@name,"response")
})

test_that("Outcome for model with defaults and units", {
  x <- saemixOutcome("concentration", unit="mg/L", model="combined1", start=c(0.5, 0.3), fix=c(1,0))
  expect_equal(x@name,"concentration")
  expect_equal(x@type,"continuous")
  expect_equal(x@error.model,"combined1")
  expect_equal(x@error.npar,2)
  expect_equal(x@error.parameters,c(0.5,0.3))
  pain <- saemixOutcome("pain score",type="categorical", unit="(-)")
  resp <- saemixOutcome("response",type="binary")
  expect_equal(pain@name,"pain score")
  expect_equal(resp@name,"response")
})


context("Creating outcomes lists and converting them with convert2Outcome")

test_that("3 outcomes defined by saemixOutcome()", {
  lout<-list(conc=saemixOutcome(unit="mg/L", model="combined1", start=c(0.5, 0.3), fix=c(1,0)), 
             effect=saemixOutcome(unit="%", model="proportional"),
             pain=saemixOutcome(type="categorical", levels=1:3)) 
  x<-convertArg2Outcome(lout)
  expect_equal(x[[1]]@name,"conc")
  expect_equal(x[[1]]@type,"continuous")
  expect_equal(x[[1]]@error.model,"combined1")
  expect_equal(x[[1]]@error.npar,2)
  expect_equal(x[[1]]@error.parameters,c(0.5,0.3))  
  expect_equal(x[[2]]@name,"effect")
  expect_equal(x[[2]]@type,"continuous")
  expect_equal(x[[2]]@error.model,"proportional")
  expect_equal(x[[3]]@name,"pain")
  expect_equal(x[[3]]@type,"categorical")
  expect_equal(length(x[[3]]@levels),3)
})


test_that("3 outcomes defined by continousOutcome() or discreteOutcome()", {
  lout2<-list(conc=continuousOutcome(unit="mg/L",model="combined1", start=c(0.5, 0.3), fix=c(1,0)), 
              effect=continuousOutcome(unit="%", model="proportional"),
              pain=discreteOutcome(levels=1:3),
              discreteOutcome(type="event"))
  x<-convertArg2Outcome(lout2)
  expect_equal(x[[1]]@name,"conc")
  expect_equal(x[[1]]@type,"continuous")
  expect_equal(x[[1]]@error.model,"combined1")
  expect_equal(x[[1]]@error.npar,2)
  expect_equal(x[[1]]@error.parameters,c(0.5,0.3))  
  expect_equal(x[[2]]@name,"effect")
  expect_equal(x[[2]]@type,"continuous")
  expect_equal(x[[2]]@error.model,"proportional")
  expect_equal(x[[2]]@error.npar,1)
  expect_equal(x[[2]]@error.parameters,c(1))  
  expect_equal(x[[3]]@name,"pain")
  expect_equal(x[[3]]@type,"categorical")
  expect_equal(length(x[[3]]@levels),3)
  expect_equal(x[[4]]@name,"out4")
  expect_equal(x[[4]]@type,"event")
})


test_that("3 outcomes defined by the parameter names+types", {
  lout3<-c(conc="continuous", effect="continuous", pain="categorical")
  x<-convertArg2Outcome(lout3)
  expect_equal(x[[1]]@name,"conc")
  expect_equal(x[[1]]@type,"continuous")
  expect_equal(x[[1]]@error.model,"constant")
  expect_equal(x[[1]]@error.npar,1)
  expect_equal(x[[1]]@error.parameters,1)  
  expect_equal(x[[2]]@name,"effect")
  expect_equal(x[[2]]@type,"continuous")
  expect_equal(x[[2]]@error.model,"constant")
  expect_equal(x[[3]]@name,"pain")
  expect_equal(x[[3]]@type,"categorical")
  expect_equal(length(x[[3]]@levels),0)
})

test_that("3 outcomes defined by the names alone (defaults to continuous outcome)", {
  lout4<-c("conc","effect","pain") # pain will be considered as continuous
  x<-convertArg2Outcome(lout4)
  expect_equal(x[[1]]@name,"conc")
  expect_equal(x[[1]]@type,"continuous")
  expect_equal(x[[1]]@error.model,"constant")
  expect_equal(x[[1]]@error.npar,1)
  expect_equal(x[[1]]@error.parameters,1)  
  expect_equal(x[[2]]@name,"effect")
  expect_equal(x[[2]]@type,"continuous")
  expect_equal(x[[2]]@error.model,"constant")
  expect_equal(x[[3]]@name,"pain")
  expect_equal(x[[3]]@type,"continuous")
  expect_equal(x[[3]]@error.model,"constant")
  expect_equal(x[[3]]@error.npar,1)
})

test_that("3 outcomes defined by the name=type or type alone", {
  lout5<-c(conc="continuous","continuous",pain="categorical")
  x<-convertArg2Outcome(lout5)
  expect_equal(x[[1]]@name,"conc")
  expect_equal(x[[1]]@type,"continuous")
  expect_equal(x[[1]]@error.model,"constant")
  expect_equal(x[[1]]@error.npar,1)
  expect_equal(x[[1]]@error.parameters,1)  
  expect_equal(x[[2]]@name,"out2")
  expect_equal(x[[2]]@type,"continuous")
  expect_equal(x[[2]]@error.model,"constant")
  expect_equal(x[[2]]@error.npar,1)
  expect_equal(x[[2]]@error.parameters,1)  
  expect_equal(x[[3]]@name,"pain")
  expect_equal(x[[3]]@type,"categorical")
  expect_equal(length(x[[3]]@levels),0)
})

test_that("Mixing types", {
  lout6<-c(conc=saemixOutcome(unit="mg/L",model="combined1", start=c(0.5, 0.3), fix=c(1,0)),"continuous",pain="categorical")
  x<-convertArg2Outcome(lout6)
  expect_equal(x[[1]]@name,"conc")
  expect_equal(x[[1]]@type,"continuous")
  expect_equal(x[[1]]@error.model,"combined1")
  expect_equal(x[[1]]@error.npar,2)
  expect_equal(x[[1]]@error.parameters,c(0.5,0.3))  
  expect_equal(x[[2]]@name,"out2")
  expect_equal(x[[2]]@type,"continuous")
  expect_equal(x[[2]]@error.model,"constant")
  expect_equal(x[[2]]@error.npar,1)
  expect_equal(x[[2]]@error.parameters,1)  
  expect_equal(x[[3]]@name,"pain")
  expect_equal(x[[3]]@type,"categorical")
  expect_equal(length(x[[3]]@levels),0)
})
