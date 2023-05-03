#########################################################################################
# Jackknife saemix

#' Jackknife for saemix fits
#' 
#' This function performs leave-one-out analyses for a saemixObject run, removing one subject at a time and computing parameter estimates.
#' 
#' @param saemixObject an object returned by the \code{\link{saemix}} function
#' @param nmax maximum number of subjects in the leave-one-out approach (when nmax is smaller than the number of subjects in the dataset, a random sample of size nmax is used to compute the jackknife distribution), defaults to NULL (all the subjects, can last a long time with a large number of subjects)
#' @param compute.likelihood whether to estimate the likelihood for each jackknife run. Defaults to FALSE, in which case only the population parameters are estimated (caution: estimating the likelihood with IS is very time-consuming to do repeatedly)
#' @param saemix.options list of options to run the saemix algorithm. Defaults to the options in the object, but suppressing the estimation of individual parameters (map=FALSE) and likelihood (ll.is=FALSE), and the options to print out intermediate results and graphs (displayProgress=FALSE, save.graphs=FALSE, save=FALSE)
#' 
#' @details 
#'   
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>
#' 
#' @references 
#' 
#' @examples 
#' 
#' # Bootstrap for the theophylline data 
#' data(theo.saemix)
#' 
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
#'   name.group=c("Id"),name.predictors=c("Dose","Time"),
#'   name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#'   units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
#' 
#' model1cpt<-function(psi,id,xidep) { 
#' 	  dose<-xidep[,1]
#' 	  tim<-xidep[,2]  
#' 	  ka<-psi[id,1]
#' 	  V<-psi[id,2]
#' 	  CL<-psi[id,3]
#' 	  k<-CL/V
#' 	  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
#' 	  return(ypred)
#' }
#' 
#' saemix.model<-saemixModel(model=model1cpt,
#'   description="One-compartment model with first-order absorption", 
#'   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,
#'   dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),
#'   covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
#'   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#'   omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")
#' 
#' saemix.options<-list(algorithm=c(1,0,0),seed=632545,save=FALSE,save.graphs=FALSE, 
#' displayProgress=FALSE)
#' 
#' # Not run (strict time constraints for CRAN)
#' \donttest{
#' saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' jack.distribution <- saemix.jackknife(saemix.fit)
#' print(jack.distribution)
#' }
#' 
#' 
#' @export saemix.jackknife


saemix.jackknife<-function(saemixObject, nmax=NULL, saemix.options=NULL, compute.likelihood=FALSE) {
  # Default options for jackknife: run estimation step only
  if(is.null(saemix.options)) {
    saemix.options<-saemixObject["options"]
    saemix.options$directory<-"current"
    saemix.options$fix.seed<-FALSE
    saemix.options$map<-FALSE   # Only parameter estimates are required for bootstrap
    saemix.options$fim<-FALSE
    saemix.options$displayProgress<-FALSE 
    saemix.options$save.graphs<-FALSE
    saemix.options$save<-FALSE
    saemix.options$print<-FALSE
  } 
  if(!compute.likelihood) {
    saemix.options$ll.is<-FALSE
    saemix.options$ll.gq<-FALSE
  } else {
    saemix.options$fim<-TRUE
    if(!is.null(saemix.options$ll.is) && saemix.options$ll.is)
      saemix.options$map<-TRUE # need individual estimates to compute the FIM
  }
  verbose <- saemix.options$warnings
  if(saemixObject@model@modeltype=="structural") idx.eps<-saemixObject@model@indx.res else idx.eps<-integer(0)
  idx.iiv<-saemixObject@model@indx.omega
  idx.rho<-which(saemixObject@model@covariance.model[lower.tri(saemixObject@model@covariance.model)]==1)
  bootstrap.distribution<-failed.runs<-data.frame()
  nelements <- length(saemixObject@results@fixed.effects)+length(idx.iiv)+length(idx.rho)+length(idx.eps)
  # Starting point: estimates from the fit 
  model.jack<-saemixObject["model"]
  model.jack@psi0 <- model.jack["betaest.model"]
  model.jack@psi0[model.jack["betaest.model"]==1]<-saemixObject@results@fixed.effects
  
  saemix.data <- saemixObject@data
  if(!is.null(nmax) && nmax<saemix.data@N)
    samplesuj<-sort(sample(1:saemix.data@N, nmax)) else samplesuj<-1:saemix.data@N
  zesuj<-unique(saemix.data@data$index)

  jackknife.distribution<-NULL
  for(isuj in samplesuj) {
    sdata <- subset(saemix.data, index!=zesuj[isuj])
    yfit <- try(saemix(model.jack, sdata, control=saemix.options))
    if(is(yfit,"try-error")) {
      l1<-c(isuj,rep(NA,nelements))
      failed.runs <- rbind(failed.runs, c(iboot, yfit))
    } else {
      res<-yfit@results
      l1<-c(isuj,res@fixed.effects, diag(res@omega)[idx.iiv])
      if(length(idx.rho)>0) l1<-c(l1,res@omega[lower.tri(res@omega)][idx.rho])
      if(length(idx.eps)>0) l1<-c(l1, res@respar[idx.eps])
      if(length(res@ll.lin)>0) l1<-c(l1, res@ll.lin)
      if(length(res@ll.is)>0) l1<-c(l1, res@ll.is)
      if(length(res@ll.gq)>0) l1<-c(l1, res@ll.gq)
    }
    jackknife.distribution<-rbind(jackknife.distribution,l1) 
  }
  # Names
  nampar<-colnames(saemixObject@model@covariance.model)
  namcol<-c(saemixObject@results@name.fixed, saemixObject@results@name.random)
  if(length(idx.rho)>0) {
    for(i in 1:(length(nampar)-1)) {
      for(j in (i+1):length(nampar)) {
        if(saemixObject@model@covariance.model[i,j]==1) {
          namcol<-c(namcol,paste("cov.",nampar[i],nampar[j],sep=""))
        }
      }
    }
  }
  if(length(idx.eps)>0) namcol<-c(namcol,saemixObject@model@name.sigma[idx.eps])
  if(length(res@ll.lin)>0) namcol<-c(namcol,"LL.lin")
  if(length(res@ll.is)>0) namcol<-c(namcol,"LL.is")
  if(length(res@ll.gq)>0) namcol<-c(namcol,"LL.gq")
  colnames(jackknife.distribution)<-c("Replicate",namcol)
  if(verbose && dim(failed.runs)[1]>0) {
    cat(dim(failed.runs)[1],"failed:\n")
    print(head(failed.runs))
  }
  return(jackknife.distribution)
}

#########################################################################################
# BCa corrected percentile for bootstrap

#' Compute BCa corrected percentiles for bootstrap
#' 
#' This function computes the BCa corrected percentile.
#' 
#' @param saemixObject an object returned by the \code{\link{saemix}} function
#' @param bootstrap.distribution the bootstrap distribution from which BCa percentiles are to be computed (obtained by a call to the boostrap function \code{\link{bootstrap.saemix}})
#' @param jackknife.distribution the jackknife distribution to compute the acceleration factor (if NULL, the BC percentile with an acceleration of 1 will be reported instead). The jackknife distribution can be obtained by a call to the function \code{\link{jackknife.saemix}}
#' @param compute.likelihood whether to compute BCa percentiles also for likelihood (if the initial fit or either of the two distributions doesn't have these then only percentiles for the population parameters are obtained) 
#' 
#' @details Formula
#'   
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>
#' 
#' @references 
#' 
#' Efron 
#' 
#' @examples 


bca.percentile <- function(saemixObject, bootstrap.distribution, jackknife.distribution=NULL, alpha=0.05, compute.likelihood=FALSE) {
  # check there are 
  if(compute.likelihood) {
    idx.lin <- length(grep("ll.lin", colnames(bootstrap.distribution)))*length(grep("ll.lin", colnames(jackknife.distribution)))*length(saemixObject@results@ll.lin)
    if(idx.lin==1) dist.lin<-TRUE else dist.lin<-FALSE
    idx.is <- length(grep("ll.is", colnames(bootstrap.distribution)))*length(grep("ll.is", colnames(jackknife.distribution)))*length(saemixObject@results@ll.is)
    if(idx.is==1) dist.is<-TRUE else dist.is<-FALSE
    idx.gq <- length(grep("ll.gq", colnames(bootstrap.distribution)))*length(grep("ll.gq", colnames(jackknife.distribution)))*length(saemixObject@results@ll.gq)
    if(idx.gq==1) dist.gq<-TRUE else dist.gq<-FALSE
    if(idx.lin+idx.gq+idx.is==0) compute.likelihood<-FALSE # no matching information
  } 
  if(!compute.likelihood) {
    dist.lin <- FALSE
    dist.is <- FALSE
    dist.gq <- FALSE
  }
  if(alpha<=0 | alpha>=1) alpha <- 0.05
  u <- c(alpha/2, 1-alpha/2) 
  zu <- qnorm(u)
  # Extracting original parameters from the model
  if(saemixObject@model@modeltype=="structural") idx.eps<-saemixObject@model@indx.res else idx.eps<-integer(0)
  idx.iiv<-saemixObject@model@indx.omega
  idx.rho<-which(saemixObject@model@covariance.model[lower.tri(saemixObject@model@covariance.model)]==1)
  res<-saemixObject@results
  origpar<-c(res@fixed.effects, diag(res@omega)[idx.iiv])
  if(length(idx.rho)>0) origpar<-c(origpar,res@omega[lower.tri(res@omega)][idx.rho])
  if(length(idx.eps)>0) origpar<-c(origpar, res@respar[idx.eps])
  nelements <- length(saemixObject@results@fixed.effects)+length(idx.iiv)+length(idx.rho)+length(idx.eps)
  
  # Computing Ii = (\hat(\theta) - \hat(\theta_{-i}))*(N-1)
  if(is.null(jackknife.distribution)) {
    accel <- rep(1,nelements+idx.is+idx.lin+idx.gq+1)
  } else {
    nsuj <- dim(jackknife.distribution)[1]
    for(icol in 1:nelements) {
      jackknife.distribution[,(icol+1)] <- origpar[icol]-jackknife.distribution[,(icol+1)]
    }
    jackknife.distribution<-jackknife.distribution*(nsuj-1)
    if(dist.lin) jackknife.distribution[,"ll.lin"] <- res@ll.lin-jackknife.distribution[,"ll.lin"]
    if(dist.is) jackknife.distribution[,"ll.is"] <- res@ll.is-jackknife.distribution[,"ll.is"]
    if(dist.gq) jackknife.distribution[,"ll.gq"] <- res@ll.gq-jackknife.distribution[,"ll.gq"]
    
    # Computing acceleration for all columns
    accel <- colSums(jackknife.distribution**3)/(colSums(jackknife.distribution**2)**(3/2))/6
  }
  
  # correcting the bootstrap distribution
  boot.bca <- NULL
  for(icol in 1:nelements) {
    z0 <- qnorm(mean(bootstrap.distribution[,(icol+1)] <= origpar[icol]))
    uadj <- pnorm(z0 + (z0+zu)/(1-accel[icol+1]*(z0+zu))) 
    boot.bca<-cbind(boot.bca,quantile(bootstrap.distribution[,(icol+1)], uadj))
  }
  colnames(boot.bca)<-colnames(bootstrap.distribution)[2:(nelements+1)]
  if(dist.lin) {
    z0 <- qnorm(mean(bootstrap.distribution[,"ll.lin"] <= res@ll.lin))
    uadj <- pnorm(z0 + (z0+zu)/(1-accel["ll.lin"]*(z0+zu))) 
    boot.bca<-cbind(boot.bca,quantile(bootstrap.distribution[,"ll.lin"], uadj))
    colnames(boot.bca)[dim(boot.bca)[2]]<-"ll.lin"
  }
  if(dist.is) {
    z0 <- qnorm(mean(bootstrap.distribution[,"ll.is"] <= res@ll.is))
    uadj <- pnorm(z0 + (z0+zu)/(1-accel["ll.is"]*(z0+zu))) 
    boot.bca<-cbind(boot.bca,quantile(bootstrap.distribution[,"ll.is"], uadj))
    colnames(boot.bca)[dim(boot.bca)[2]]<-"ll.is"
  }
  if(dist.gq) {
    z0 <- qnorm(mean(bootstrap.distribution[,"ll.gq"] <= res@ll.gq))
    uadj <- pnorm(z0 + (z0+zu)/(1-accel["ll.gq"]*(z0+zu))) 
    boot.bca<-cbind(boot.bca,quantile(bootstrap.distribution[,"ll.gq"], uadj))
    colnames(boot.bca)[dim(boot.bca)[2]]<-"ll.gq"
  }
  rownames(boot.bca)<-paste0(u*100,c("%"))
  
  return(boot.bca)
}

