context("Creating SaemixModel objects ")

test_that("Create model from a list of parameters - basic", {
  param1<-list(ka=saemixParam(mu.init=c(1,3), sd.init=c(0.5)),
               vd=saemixParam(sd.init=0.7),  
               cl=saemixParam(name="cl",mu.init=2))
  model1cpt<-function(psi,id,xidep) { 
    dose<-xidep[,1]
    tim<-xidep[,2]  
    ka<-psi[id,1]
    V<-psi[id,2]
    CL<-psi[id,3]
    k<-CL/V
    ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
    return(ypred)
  }
  saemix.model <- new(Class="SaemixModel", parameters=param1, model=model1cpt)
  expect_equal(saemix.model@noutcome,1)
  expect_equal(saemix.model@nphi,3)
  print(saemix.model)
})


test_that("Create model from a list of parameters", {
  param4<-list(ka=saemixParam(mu.init=c(1,3), sd.init=c(0.8,0.5), varlevel=c("iiv","iov")),vd=saemixParam(sd.init=0.7, covariate="wt", covariate.init=c(1), covariate.estim=c("fixed")),  
               cl=saemixParam(name="cl",mu.init=2, varlevel=c("iiv","iov"), sd.init=c(0.6,0.3), corr = list(iiv=c("ka","vd"),iov=c("vd")), covariate=c("wt","sex","age"), covariate.init=c(0.75,0,0), covariate.estim=c("fixed","estimated","estimated"), corr.init=list(iiv=c(-0.5,0.7), iov=0.7), covariate.varlevel=c("iiv","iiv","iov")))
  model1cpt<-function(psi,id,xidep) { 
    dose<-xidep[,1]
    tim<-xidep[,2]  
    ka<-psi[id,1]
    V<-psi[id,2]
    CL<-psi[id,3]
    k<-CL/V
    ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
    return(ypred)
  }
  
  saemix.model <- new(Class="SaemixModel", parameters=param4, model=model1cpt)
  expect_equal(saemix.model@noutcome,1)
  expect_equal(saemix.model@nphi,3)
  print(saemix.model)
})

# ToDo: debug print and show
## decide on outcome structure:
### separate count and categorical ?
### 
## tests on outcome
### computing error for continuous models

test_that("Create model from a list of parameters, responses from a list of outcome", {
  param <- list(ka=saemixParam(mu.init=c(1,3), sd.init=c(0.8,0.5), varlevel=c("iiv","iov")),
                vd=saemixParam(sd.init=0.7, covariate="wt", covariate.init=c(1), covariate.estim=c("fixed")),  
                cl=saemixParam(name="cl",mu.init=2, varlevel=c("iiv","iov"), sd.init=c(0.6,0.3), corr = list(iiv=c("ka","vd"),iov=c("vd")), covariate=c("wt","sex","age"), covariate.init=c(0.75,0,0), covariate.estim=c("fixed","estimated","estimated"), corr.init=list(iiv=c(-0.5,0.7), iov=0.7), covariate.varlevel=c("iiv","iiv","iov")))
  y<-list(concentration=saemixOutcome(error.model="combined1", unit="mg/L"), 
        hospitalisation=saemixOutcome( type.outcome="event", distribution="Weibull"), 
        pain=saemixOutcome( type.outcome="categorical", distribution="binomial"))
  model1cpt<-function(psi,id,xidep) { 
    dose<-xidep[,1]
    tim<-xidep[,2]  
    ka<-psi[id,1]
    V<-psi[id,2]
    CL<-psi[id,3]
    k<-CL/V
    ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
    return(ypred)
  }
  saemix.model <- new(Class="SaemixModel", parameters=param, model=model1cpt, outcome=y)
  expect_equal(saemix.model@noutcome,3)
  for(i in 1:length(y)) y[[i]]@name.outcome <- names(y)[i]
  expect_equal(saemix.model@outcome,y)
  expect_equal(saemix.model@nphi,3)
  expect_equal(saemix.model@outcome[[1]]@unit,"mg/L")
  print(saemix.model)
})

test_that("Create model using saemixModel() with default parameters given as a vector", {
  model1cpt<-function(psi,id,xidep) { 
    dose<-xidep[,1]
    tim<-xidep[,2]  
    ka<-psi[id,1]
    V<-psi[id,2]
    CL<-psi[id,3]
    k<-CL/V
    ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
    return(ypred)
  }
  saemix.model <-saemixModel(model=model1cpt, parameters=c("ka","V","CL"))
  expect_equal(saemix.model@noutcome,1)
  expect_equal(saemix.model@outcome[[1]]@name.outcome,"y")
  expect_equal(saemix.model@outcome[[1]]@error.model@name,"constant")
  expect_equal(saemix.model@nphi,3)
  expect_equal(saemix.model@outcome[[1]]@unit,"")
  print(saemix.model)
  expect_equal(saemix.model@log,"Parameter given as names only\n")
})


test_that("Create model using saemixModel(), outcome not specified (defaults to 1 continuous outcome)", {
  param <- list(ka=saemixParam(mu.init=c(1,3), sd.init=c(0.8,0.5), varlevel=c("iiv","iov")),
                vd=saemixParam(sd.init=0.7, covariate="wt", covariate.init=c(1), covariate.estim=c("fixed")),  
                cl=saemixParam(name="cl",mu.init=2, varlevel=c("iiv","iov"), sd.init=c(0.6,0.3), corr = list(iiv=c("ka","vd"),iov=c("vd")), covariate=c("wt","sex","age"), covariate.init=c(0.75,0,0), covariate.estim=c("fixed","estimated","estimated"), corr.init=list(iiv=c(-0.5,0.7), iov=0.7), covariate.varlevel=c("iiv","iiv","iov")))
  model1cpt<-function(psi,id,xidep) { 
    dose<-xidep[,1]
    tim<-xidep[,2]  
    ka<-psi[id,1]
    V<-psi[id,2]
    CL<-psi[id,3]
    k<-CL/V
    ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
    return(ypred)
  }
  saemix.model <-saemixModel(model=model1cpt, parameters=param)
  expect_equal(saemix.model@noutcome,1)
  expect_equal(saemix.model@outcome[[1]]@name.outcome,"y")
  expect_equal(saemix.model@outcome[[1]]@error.model@name,"constant")
  expect_equal(saemix.model@nphi,3)
  expect_equal(saemix.model@outcome[[1]]@unit,"")
  print(saemix.model)
})


test_that("Create model using saemixModel() with multiple responses", {
  param <- list(ka=saemixParam(mu.init=c(1,3), sd.init=c(0.8,0.5), varlevel=c("iiv","iov")),
                vd=saemixParam(sd.init=0.7, covariate="wt", covariate.init=c(1), covariate.estim=c("fixed")),  
                cl=saemixParam(name="cl",mu.init=2, varlevel=c("iiv","iov"), sd.init=c(0.6,0.3), corr = list(iiv=c("ka","vd"),iov=c("vd")), covariate=c("wt","sex","age"), covariate.init=c(0.75,0,0), covariate.estim=c("fixed","estimated","estimated"), corr.init=list(iiv=c(-0.5,0.7), iov=0.7), covariate.varlevel=c("iiv","iiv","iov")))
  y<-list(concentration=saemixOutcome(error.model="combined1", unit="mg/L"), 
          hospitalisation=saemixOutcome( type.outcome="event", distribution="Weibull"), 
          pain=saemixOutcome( type.outcome="categorical", distribution="binomial"))
  model1cpt<-function(psi,id,xidep) { 
    dose<-xidep[,1]
    tim<-xidep[,2]  
    ka<-psi[id,1]
    V<-psi[id,2]
    CL<-psi[id,3]
    k<-CL/V
    ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
    return(ypred)
  }
  saemix.model <-saemixModel(model=model1cpt, parameters=param, outcome=y)
  expect_equal(saemix.model@noutcome,3)
  for(i in 1:length(y)) y[[i]]@name.outcome <- names(y)[i]
  expect_equal(saemix.model@outcome,y)
  expect_equal(saemix.model@nphi,3)
  expect_equal(saemix.model@outcome[[1]]@unit,"mg/L")
  print(saemix.model)
})

