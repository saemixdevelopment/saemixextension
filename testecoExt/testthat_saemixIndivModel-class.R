context("Creating SaemixIndivModel objects ")

test_that("Create individual model from parameters", {
  param4<-list(ka=saemixParam(mu.init=c(1,3), sd.init=c(0.8,0.5), varlevel=c("iiv","iov")),vd=saemixParam(sd.init=0.7, covariate="wt", covariate.init=c(1), covariate.estim=c("fixed")),  
               cl=saemixParam(name="cl",mu.init=2, varlevel=c("iiv","iov"), sd.init=c(0.6,0.3), corr = list(iiv=c("ka","vd"),iov=c("vd")), covariate=c("wt","sex","age"), covariate.init=c(0.75,0,0), covariate.estim=c("fixed","estimated","estimated"), corr.init=list(iiv=c(-0.5,0.7), iov=0.7), covariate.varlevel=c("iiv","iiv","iov")))
  indivmodel <- new(Class="SaemixIndivModel", param4)
  print(indivmodel)
  expect_equal(indivmodel@covariate, param4$cl@covariate)
  expect_equal(indivmodel@varlevel, param4$ka@varlevel)
})

context("Generating parameters using individual model")

# Generating mu+beta.cov
test_that("Generating  phibar (population value of phi adjusted for covariates) and psibar for one parameter, with covariates", {
  param1<-list(ka=saemixParam(mu.init=c(2), sd.init=c(0.5), covariate=c("age"),covariate.init=c(0.2)))
  indivmodel <- new(Class="SaemixIndivModel", param1)
  
  nsuj<-10
  cdesign <- matrix(c(rep(1,nsuj), log(seq(from=50,to=(50+2*(nsuj-1)), by=2)/60)), ncol=2)
  colnames(cdesign)<-c("pop","lage")
  
  phipop <- cdesign %*% indivmodel@popmodel[[1]]@phi
  psipop<-exp(phipop)
  expect_equal(min(psipop), param1$ka@mu.init*exp(param1$ka@covariate.init*log(50/60)))
  expect_equal(max(psipop), param1$ka@mu.init*exp(param1$ka@covariate.init*log((50+2*(nsuj-1))/60)))
})

test_that("Generating individual phi for one parameter, with covariates", {
  param1<-list(ka=saemixParam(mu.init=c(2), sd.init=c(0.5), covariate=c("age"),covariate.init=c(0.2)))
  indivmodel <- new(Class="SaemixIndivModel", param1)
  
  nsuj<-10
  cdesign <- matrix(c(rep(1,nsuj), log(seq(from=50,to=(50+2*(nsuj-1)), by=2)/60)), ncol=2)
  colnames(cdesign)<-c("pop","lage")
  
  phipop <- cdesign %*% indivmodel@popmodel[[1]]@phi

  # pas sûre de ce qu'on essaie de faire ici...  simplifier ?
  omega.eta<-indivmodel@varmodel[[1]]@subomega-mydiag(mydiag(indivmodel@varmodel[[1]]@subomega))+mydiag(cutoff(mydiag(indivmodel@varmodel[[1]]@subomega),.Machine$double.eps))
  chol.omega<-try(chol(omega.eta))
  somega<-solve(omega.eta)

  eta<-matrix(rnorm(nsuj*indivmodel@nphi),ncol=indivmodel@nphi) %*% chol.omega
  phipop[,indivmodel@varmodel[[1]]@idcol.eta] <-phipop[,indivmodel@varmodel[[1]]@idcol.eta]+eta
  expect_equal(dim(phipop),c(nsuj,1))
  print(indivmodel@transform[[1]](phipop))
})

context("Testing computational functions")

test_that("Converting from phi to psi and back using transphi", {
  param2<-list(ka=saemixParam(mu.init=c(2), sd.init=c(0.5), covariate=c("lage"),covariate.init=c(0.2)),
               vd=saemixParam(mu.init=c(20), sd.init=0.7, covariate="lwt", covariate.init=c(1), covariate.estim=c("fixed")))
  indivmodel <- new(Class="SaemixIndivModel", param2)
  
  cdesign <- matrix(c(rep(1,nsuj), log(seq(from=50,to=(50+2*(nsuj-1)), by=2)/60),
                      log(seq(from=90, to=(90-4*(nsuj-1)), by=-4)/70)), ncol=3)
  colnames(cdesign)<-c("pop","lage","lwt")
  phipop <- cdesign %*% indivmodel@popmodel[[1]]@phi
  psipop<- matrix(c(indivmodel@transform[[1]](phipop[,1]),indivmodel@transform[[2]](phipop[,2])), ncol=2)
  
  expect_equivalent(transphi(phipop, indivmodel@transform), psipop)
  expect_equivalent(transphi(psipop, indivmodel@invtransform), phipop)
})
