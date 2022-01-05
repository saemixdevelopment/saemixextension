context("Plots for SaemixModel")

test_that("Plot model for default psi, no additional predictors", {
  expmod<-function(psi,id,xidep) { 
    x<-xidep[,1]  
    a<-psi[id,1]
    k<-psi[id,2]
    ypred<-a*(1-exp(-k*x))
    return(ypred)
  }
  x1<-saemixModel(model=expmod,description="Exponential model", 
                 psi0=matrix(c(10, 1), ncol=2,byrow=TRUE, dimnames=list(NULL, c("a","k"))))
  plot(x1, range=c(0,10))
  plot(x1, range=c(0,10), verbose=TRUE)
  plot(x1, range=c(0,10), predictors=10, verbose=TRUE) # Additional predictor ignored here, plot unchanged
}

test_that("Plot model for psi defined by user", {
  expmod<-function(psi,id,xidep) { 
    x<-xidep[,1]  
    a<-psi[id,1]
    k<-psi[id,2]
    ypred<-a*(1-exp(-k*x))
    return(ypred)
  }
  x1<-saemixModel(model=expmod,description="Exponential model", 
                  psi0=matrix(c(10, 1), ncol=2,byrow=TRUE, dimnames=list(NULL, c("a","k"))))
  plot(x1, range=c(0,10), psi=c(5, 0.5), verbose=TRUE)
}

test_that("Plot model for default psi, given additional predictors", {
  model1cpt<-function(psi,id,xidep) { 
    tim<-xidep[,1]  
    dose<-xidep[,2]
    ka<-psi[id,1]
    V<-psi[id,2]
    CL<-psi[id,3]
    k<-CL/V
    ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
    return(ypred)
  }
  x2<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", 
                 psi0=matrix(c(1.5,30,1), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))))
  # without the predictors, no plot
  plot(x2, verbose=TRUE) # Prints a bunch of messages pointing to possible errors
  # with the predictors, plots the predictions
  plot(x2, range=c(0,24), predictors=300, verbose=TRUE) #
})

test_that("Plot model for default psi, across the wrong predictor (order of tim and dose in model)", {
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
  x2<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", 
                  psi0=matrix(c(1.5,30,1), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))))
  plot(x2, range=c(0,24), predictors=300, verbose=TRUE) # should give a straight line as the range is over dose, not tim, ie C(t=300) over a dose range of 0-24
})

context("Plots for SaemixModel+SaemixData created separately")

test_that("Plot model, for data in an SaemixData object, for default psi and selected psi", {
  model1cpt<-function(psi,id,xidep) { 
    tim<-xidep[,1]  
    dose<-xidep[,2]
    ka<-psi[id,1]
    V<-psi[id,2]
    CL<-psi[id,3]
    k<-CL/V
    ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
    return(ypred)
  }
  x3<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", 
                  psi0=matrix(c(1.5,30,1), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))))
  x<-saemixData(name.data=file.path(datDir,"theo.saemix.tab"),header=T,na=".", name.group=c("Id"),name.predictors=c("Time","Dose"),name.response=c("Concentration"),units=list(x="hr",y="mg/L"), name.X="Time",verbose=F)
  plot(x3, x)
  plot(x3, x, ilist=1:9)
  plot(x3, x, psi=c(1.5, 30, 2), ilist=1:9, xlim=c(0,20))
})

context("Plots for SaemixModel+SaemixData extracted from a fitted object")

test_that("Plot model, for data in an SaemixData object, for default psi and selected psi", {
  model1cpt<-function(psi,id,xidep) { 
    tim<-xidep[,1]  
    dose<-xidep[,2]
    ka<-psi[id,1]
    V<-psi[id,2]
    CL<-psi[id,3]
    k<-CL/V
    ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
    return(ypred)
  }
  x3<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", 
                  psi0=matrix(c(1.5,30,1), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))))
  x<-saemixData(name.data=file.path(datDir,"theo.saemix.tab"),header=T,na=".", name.group=c("Id"),name.predictors=c("Time","Dose"),name.response=c("Concentration"),units=list(x="hr",y="mg/L"), name.X="Time",verbose=F)
  opt<-list(seed=123456, save=FALSE, displayProgress=FALSE)
  yfit<-saemix(x3, x, opt)
  plot(yfit@model, yfit@data) # initial values
  plot(yfit@model, yfit@data, psi=yfit@results@fixed.psi) # population fits
  plot(yfit@model, yfit@data, psi=yfit@results@cond.mean.psi) # individual fits
})

test_that("Plot model, for data in an SaemixData object, for default psi and selected psi, when dose is the first predictor (works because predictors are matched here)", {
  model1cpt<-function(psi,id,xidep) { 
    tim<-xidep[,2]  
    dose<-xidep[,1]
    ka<-psi[id,1]
    V<-psi[id,2]
    CL<-psi[id,3]
    k<-CL/V
    ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
    return(ypred)
  }
  x3<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", 
                  psi0=matrix(c(1.5,30,1), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))))
  x<-saemixData(name.data=file.path(datDir,"theo.saemix.tab"),header=T,na=".", name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),units=list(x="hr",y="mg/L"), name.X="Time",verbose=F)
  opt<-list(seed=123456, save=FALSE, displayProgress=FALSE)
  yfit<-saemix(x3, x, opt)
  plot(yfit@model, yfit@data) # initial values
  plot(yfit@model, yfit@data, psi=yfit@results@fixed.psi) # population fits
  plot(yfit@model, yfit@data, psi=yfit@results@cond.mean.psi) # individual fits
})

