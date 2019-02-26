context("Testing creation of SaemixModel for structural models\n")

test_that("Creating an empty SaemixModel object with class", {
  x<-new(Class="SaemixModel")
  expect_is(x,"SaemixModel")
  expect_null(body(x@model))
  expect_is(x@psi0,"matrix")
  expect_equal(dim(x@psi0),c(0,0))
  #  expect_error(new(Class="SaemixModel"));    # did not test for particular error message (not sure if good practice), but could be changed if desired
})

test_that("Errors in using saemixModel - no model", {
#  expect_error(x<-saemixModel());
  x<-saemixModel();
  expect_is(x,"character")
  expect_equal(x,"Creation of SaemixModel failed")
})

test_that("Errors in using saemixModel - no initial values for parameters", {
  tte.model<-function(psi,id,xidep) {
    T<-xidep[,1]
    N <- nrow(psi)
    Nj <- length(T)
    censoringtime = 6
    lambda <- psi[id,1]
    beta <- psi[id,2]
    init <- which(T==0)
    cens <- which(T==censoringtime)
    ind <- setdiff(1:Nj, append(init,cens))
    hazard <- (beta/lambda)*(T/lambda)^(beta-1)
    H <- (T/lambda)^beta
    logpdf <- rep(0,Nj)
    logpdf[cens] <- -H[cens] + H[cens-1]
    logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind])
    return(logpdf)
  }
  x<-saemixModel(model=tte.model);    
  expect_is(x,"character")
  expect_equal(x,"Creation of SaemixModel failed")
})

test_that("Minimal SaemixModel object", {
  tte.model<-function(psi,id,xidep) {
    T<-xidep[,1]
    N <- nrow(psi)
    Nj <- length(T)
    censoringtime = 6
    lambda <- psi[id,1]
    beta <- psi[id,2]
    init <- which(T==0)
    cens <- which(T==censoringtime)
    ind <- setdiff(1:Nj, append(init,cens))
    hazard <- (beta/lambda)*(T/lambda)^(beta-1)
    H <- (T/lambda)^beta
    logpdf <- rep(0,Nj)
    logpdf[cens] <- -H[cens] + H[cens-1]
    logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind])
    return(logpdf)
  }
  x<-saemixModel(model=tte.model, psi0=matrix(c(1.,20,0.5), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))), modeltype="structural")
  expect_is(x, "SaemixModel") # tests for particular class
  expect_equal(x@nb.parameters,3)
})

test_that("Successful creation of a SaemixModel object of type structural", {
  tte.model<-function(psi,id,xidep) {
    T<-xidep[,1]
    N <- nrow(psi)
    Nj <- length(T)
    censoringtime = 6
    lambda <- psi[id,1]
    beta <- psi[id,2]
    init <- which(T==0)
    cens <- which(T==censoringtime)
    ind <- setdiff(1:Nj, append(init,cens))
    hazard <- (beta/lambda)*(T/lambda)^(beta-1)
    H <- (T/lambda)^beta
    logpdf <- rep(0,Nj)
    logpdf[cens] <- -H[cens] + H[cens-1]
    logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind])
    return(logpdf)
  }
  x<-saemixModel(model=tte.model,description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))))
  expect_is(x, "SaemixModel") # tests for particular class
  expect_equal(x@nb.parameters,3)
})

test_that("Successful creation of a SaemixModel object with 2 responses", {
  tte.model<-function(psi,id,xidep) {
    T<-xidep[,1]
    N <- nrow(psi)
    Nj <- length(T)
    censoringtime = 6
    lambda <- psi[id,1]
    beta <- psi[id,2]
    init <- which(T==0)
    cens <- which(T==censoringtime)
    ind <- setdiff(1:Nj, append(init,cens))
    hazard <- (beta/lambda)*(T/lambda)^(beta-1)
    H <- (T/lambda)^beta
    logpdf <- rep(0,Nj)
    logpdf[cens] <- -H[cens] + H[cens-1]
    logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind])
    return(logpdf)
  }
  x<-saemixModel(model=tte.model,description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))), name.response=c("PK","PD"),name.sigma=c("add","prop"), error.model=c("combined","proportional"))
  expect_is(x, "SaemixModel") # tests for particular class
  expect_equal(x@nb.parameters,3)
  expect_equal(length(x@indx.res),3)
})

test_that("Successful creation of a SaemixModel object with all arguments", {
  tte.model<-function(psi,id,xidep) {
    T<-xidep[,1]
    N <- nrow(psi)
    Nj <- length(T)
    censoringtime = 6
    lambda <- psi[id,1]
    beta <- psi[id,2]
    init <- which(T==0)
    cens <- which(T==censoringtime)
    ind <- setdiff(1:Nj, append(init,cens))
    hazard <- (beta/lambda)*(T/lambda)^(beta-1)
    H <- (T/lambda)^beta
    logpdf <- rep(0,Nj)
    logpdf[cens] <- -H[cens] + H[cens-1]
    logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind])
    return(logpdf)
  }
  x<-saemixModel(model=tte.model,description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))), transform.par=c(1,1,1), covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE), fixed.estim=c(1,1,1), covariance.model=matrix(c(1,0,0,0,1,1,0,1,1),ncol=3,byrow=TRUE), omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), error.model="combined",error.init=c(1,0.5))
  expect_is(x, "SaemixModel") # tests for particular class
  expect_equal(x@nb.parameters,3)
  expect_identical(colnames(x@omega.init),colnames(x@psi0))
  expect_equal(dim(x@omega.init),c(3,3))
  expect_identical(x@error.init,c(1,0.5))
  expect_identical(sum(x@covariate.model),1)
  expect_identical(sum(x@betaest.model),4)
  expect_identical(sum(x@transform.par),3)
  expect_identical(sum(x@fixed.estim),3)
})
