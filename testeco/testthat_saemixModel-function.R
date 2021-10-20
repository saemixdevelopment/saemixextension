context("Auxiliary functions for SaemixModel objects - covariance model\n")

test_that("Covariance models - test matrix containing elements other than 0/1", {
  xmat<-matrix(c(1,0,0,2),ncol=2)
  expect_false(validate.covariance.model(xmat))
})

test_that("Covariance models - not a square matrix", {
  xmat<-matrix(c(1,0,0,2,3,4),ncol=2)
  expect_false(validate.covariance.model(xmat))
})

test_that("Covariance models - non symmetrical", {
  xmat<-matrix(c(1,0,1,1),ncol=2)
  expect_false(validate.covariance.model(xmat))
})

test_that("Covariance models - test covariance model within an object", {
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
  x<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))))
  expect_true(validate.covariance.model(x@covariance.model))
})

context("Predict function for SaemixModel objects\n")

test_that("Predictions given an SaemixModel object, using the psi0 slot", {
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
  x<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))))
  
  xidep<-data.frame(dose=100, tim=seq(0,24,2))
  xtim<-xidep[,2]
  ypred1<-predict.saemixmodel(x, xidep)
  ycomp1<-xidep[1,1]*1/(20*(1-0.5/20))*(exp(-0.5/20*xtim)-exp(-1*xtim))
  expect_null(rownames(ypred1$param))
  expect_equal(colnames(ypred1$param)[-c(1)], colnames(x@psi0))
  expect_identical(ypred1$predictions$pred,ycomp1)
})

test_that("Predictions given an SaemixModel object for a different set of parameters", {
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
  x<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))))
  
  xidep<-data.frame(dose=100, tim=seq(0,24,2))
  xtim<-xidep[,2]
  
  ypred2<-predict.saemixmodel(x, xidep, psi=c(2, 25, 0.5))
  ycomp2<-xidep[1,1]*2/(25*(2-0.5/25))*(exp(-0.5/25*xtim)-exp(-2*xtim))
  expect_identical(ypred2$predictions$pred,ycomp2)
})

test_that("Wrong number of parameters", {
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
  x<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))))
  
  xidep<-data.frame(dose=100, tim=seq(0,24,2))
  xtim<-xidep[,2]
  
  ypred2<-predict.saemixmodel(x, xidep, psi=c(2, 25))
  expect_null(ypred2)
})


test_that("Predictions given an SaemixModel object for several subjects", {
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
  x<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))))
  
  xidep<-data.frame(dose=100, tim=seq(0,24,2))
  xtim<-xidep[,2]
  id<-rep(1:3, each=length(xidep[,2]))
  xidep<-do.call(rbind,rep(list(xidep),3))
  psi1<-do.call(rbind,list(c(2, 25, 0.5), c(1,20,0.5),c(1.5, 20, 1)))
  
  ypred2<-predict.saemixmodel(x, xidep, psi=psi1, id=id)
  ycomp2<-c(xidep[1,1]*2/(25*(2-0.5/25))*(exp(-0.5/25*xtim)-exp(-2*xtim)),
            xidep[1,1]*1/(20*(1-0.5/20))*(exp(-0.5/20*xtim)-exp(-1*xtim)),
            xidep[1,1]*1.5/(20*(1.5-1/20))*(exp(-1/20*xtim)-exp(-1.5*xtim)))
  expect_identical(ypred2$predictions$pred,ycomp2)
})

# object<-x
# predictors<-xidep
# psi<-psi1
