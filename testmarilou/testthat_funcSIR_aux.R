context('Testing auxiliary functions of SIR')

################################## DONNÉES ET MODÈLES À TESTER ############################################
####### X1 - invalid covariance model ######
xdata1<-saemixData(name.data=simPDEmax20, 
                   name.group=c("id"),name.predictors=c("dose"),
                   name.response=c("ysim"), name.X="dose") 
xmodel1<-saemixModel(model=modelPD,
                     description="PD Emax model",
                     psi0=matrix(c(20,200,10),ncol=3,byrow=TRUE, 
                                 dimnames=list(NULL, c("E0","Emax","ED50"))), fixed.estim = c(1,1,1),transform.par=c(1,1,1), 
                     error.model = "proportional", covariance.model = matrix(c(1,0,0,0,0,1,0,1,1),nrow=3, ncol=3, byrow=T))
###invalid covariance model
x1<-saemix(xmodel1,xdata1)

######### X2 - full covariance model ##########
xmodel2<-saemixModel(model=modelPD,
                     description="PD Emax model",
                     psi0=matrix(c(20,200,10),ncol=3,byrow=TRUE, 
                                 dimnames=list(NULL, c("E0","Emax","ED50"))), fixed.estim = c(1,1,1),transform.par=c(1,1,1), 
                     error.model = "proportional", covariance.model = matrix(c(1,1,1,1,1,1,1,1,1),nrow=3, ncol=3, byrow=T))
x2<-saemix(xmodel2,xdata1)


########### X3 - invalid covariance model 2 ###########
xmodel3<-saemixModel(model=modelPD,
                     description="PD Emax model",
                     psi0=matrix(c(20,200,10),ncol=3,byrow=TRUE, 
                                 dimnames=list(NULL, c("E0","Emax","ED50"))), fixed.estim = c(1,1,1),transform.par=c(1,1,1), 
                     error.model = "proportional", covariance.model = matrix(c(1,1,0,0,0,0,0,0,1),nrow=3, ncol=3, byrow=T))
xmodel4<-saemixModel(model=modelPD,
                     description="PD Emax model",
                     psi0=matrix(c(20,200,10),ncol=3,byrow=TRUE, 
                                 dimnames=list(NULL, c("E0","Emax","ED50"))), fixed.estim = c(1,1,1),transform.par=c(1,1,1), 
                     error.model = "proportional", covariance.model = matrix(c(1,0,0,1,0,0,0,0,1),nrow=3, ncol=3, byrow=T))

xmodel5<-saemixModel(model=modelPD,
                     description="PD Emax model",
                     psi0=matrix(c(20,200,10),ncol=3,byrow=TRUE, 
                                 dimnames=list(NULL, c("E0","Emax","ED50"))), fixed.estim = c(1,1,1),transform.par=c(1,1,1), 
                     error.model = "proportional", covariance.model = matrix(c(1,0,0,0,0,0,0,0,1),nrow=3, ncol=3, byrow=T))

###invalid covariance models
x3<-saemix(xmodel3,xdata1)
x4<-saemix(xmodel4,xdata1)
x5 <- saemix(xmodel5,xdata1)
#########################################################################################################################



######### indx.covomega #########
test_that("Non-valid covariance model in SaemixObject do not give wrong indx.covomega", {
  indx <- indx.covomega(x1)
  expect_equal(indx, NULL)
})

test_that("indx.covomega gives the right indx with full covariance model", {
  indx <- indx.covomega(x2)
  expect_equivalent(indx, matrix(c(2,1,3,1,3,2), ncol=2,byrow=T))
})


######### lllin.saemix #########
test_that("llin.saemix computes correctly log-likelihood by linearisation", {
  x <- test1
  ll <- lllin.saemix(x)[1,1]
  expect_equal(ll,x@results@ll.lin)
})

test_that("Replacing est.mu with replacePopPar does not change de ll.lin",{
  x <- test1
  est.mu <- estpar.vector(x)
  llx <- lllin.saemix(x)
  y <- replacePopPar.saemixObject(x,est.mu)
  lly <- lllin.saemix(y)
  expect_equal(llx,lly)
})




