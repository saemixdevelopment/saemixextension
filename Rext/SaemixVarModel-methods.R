#' @include aaa_generics.R
#' @include SaemixParameter.R
#' @include SaemixParameter-methods.R
NULL

####################################################################################
####  Auxiliary function

#' Matrix diagonal
#' 
#' Extract or replace the diagonal of a matrix, or construct a diagonal matrix (replace diag function from R-base)
#' 
#' @param x	a matrix, vector or 1D array, or missing.
#' @param nrow Optional number of rows for the result when x is not a matrix. 
#' @param ncol Optional number of columns for the result when x is not a matrix. 
#' 
#' @return If x is a matrix then diag(x) returns the diagonal of x. The resulting vector will have names if the matrix x has matching column and rownames.
#' @seealso \code{diag}
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @keywords models
#' @examples
#' 
#' mydiag(1)
#' mydiag(c(1,2))
#' 
#' @export mydiag

# Redefining diag function, too many problems with the R version
mydiag <- function (x = 1, nrow, ncol) {
  if (is.matrix(x)) {
    if (nargs() > 1L) 
      stop("'nrow' or 'ncol' cannot be specified when 'x' is a matrix")
    if ((m <- min(dim(x))) == 0L) 
      return(vector(typeof(x), 0L))
    y <- c(x)[1L + 0L:(m - 1L) * (dim(x)[1L] + 1L)]
    nms <- dimnames(x)
    if (is.list(nms) && !any(sapply(nms, is.null)) && identical((nm <- nms[[1L]][seq_len(m)]), 
                                                                nms[[2L]][seq_len(m)])) 
      names(y) <- nm
    return(y)
  }
  if (is.array(x) && length(dim(x)) != 1L) 
    stop("'x' is an array, but not 1D.")
  if (missing(x)) 
    n <- nrow
  else n <- length(x)
  if (!missing(nrow)) 
    n <- nrow
  if (missing(ncol)) 
    ncol <- n
  p <- ncol
  y <- array(0, c(n, p))
  if ((m <- min(n, p)) > 0L) 
    y[1L + 0L:(m - 1L) * (n + 1L)] <- x
  y
}

############ VALIDITY OF COVARIANCE MODEL 


#' Validate the structure of the covariance model
#' 
#' Check that a matrix corresponds to a structure defining a covariance model for a non-linear mixed effect model.
#' Such a matrix should be composed of only 0s and 1s, with at least one element set to 1, and should be square and symmetrical.
#' 1s on the diagonal indicate that the corresponding parameter has interindividual variability and that its variance will be estimated.
#' 1s as off-diagonal elements indicate that a covariance between the two corresponding parameters will be estimated.
#' 
#' @param x	a matrix
#' @param verbose	a boolean indicating whether warnings should be output if x is not a valid covariance model
#' 
#' @return a boolean, TRUE if x is an acceptable structure and FALSE if not. Messages will be output to describe why x isn't a valid covariance model if the argument verbose is TRUE.
#' @seealso \code{SaemixModel}
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Belhal Karimi
#' @keywords models
#' @examples
#' 
#' covarmodel<-diag(c(1,1,0))
#' validate.covariance.model(covarmodel) # should return TRUE
#' 
#' @export validate.covariance.model
#' 
validate.covariance.model <- function(x, verbose=TRUE){
  #non-square matrix
  if(dim(x)[1]!=dim(x)[2]) {
    if(verbose) message("Error initialising SaemixModel object:\n   The covariance model needs to be a square matrix, please check dimensions.\n")
    return(FALSE)
  }
  # only 0s
  s <- sum(abs(x))
  if(s==0) {
    if(verbose) message("At least one parameter should have IIV in the model, the covariance model may not be only 0s.")
    #   return(FALSE)
  }
  
  #values other than 1 or 0
  s <- sum(x[x!=1 & x!=0])
  if (s>0){
    if(verbose) message("Error initialising SaemixModel object:\n  Invalid covariance model, only 0 or 1 values accepted, please change covariance model.\n")
    return(FALSE)
  }
  
  #asymmetrical
  if (!all(t(x)==x)){
    if(verbose) message("Error initialising SaemixModel object:\n  The matrix defining the covariance model is not symmetrical, please change covariance model.\n")
    return(FALSE)
  }
  
  #values other than 0 when diagonal number is 0
  for (i in 1:nrow(x)){
    for(j in 1:ncol(x)){
      if (x[i,j]!=0){
        if(x[i,i]==0 |x[j,j]==0){
          if(verbose) message("Error initialising SaemixModel object:\n  The covariance model is invalid, please change covariance model.\n")
          return(FALSE)
        }
      }
    }
  }
  return(TRUE)
}

################################################################################

#' Function to create variability levels when given a list of parameters
#' 
#' Constructs the statistical model when given a list of parameters. 
#' This function is used by XXXTODO and is not intended to be called by the user.
#' 
#' @name extractVarModel
#' @aliases validateSaemixVar addSaemixIndices
#' 
#' @param parameters a list of parameters of class SaemixParam
#' @param varModel if given, the list of variability levels to extract. If omitted, the variability levels will be derived from the variability levels contained in the parameters.
#' 
#' @details Starting from saemix 4.0, parameters are specified through their distribution.
#' 
#' TODO: list of parameters
#' 
#' TODO: briefly describe other helper functions
#' 
#' @return An SaemixVarModel object containing the levels of variability in varModel (if given) or in the parameters (see \code{\link{SaemixVarModel}}).
#' 
#' @references E Comets, A Lavenu, M Lavielle M (2017). Parameter estimation in nonlinear mixed effect models using saemix,
#' an R implementation of the SAEM algorithm. Journal of Statistical Software, 80(3):1-41.
#' 
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}
#' @seealso \code{\link{SaemixData}},\code{\link{SaemixModel}}, \code{\link{saemixControl}},\code{\link{saemix}}
#' 
#' @examples
#' # Statistical model including two levels of variability
#' saemix.parameters <- list(ka=saemixParam(sd.init=c(0.8,0.5), varlevel=c("iiv","iov")),
#'      vd=saemixParam(sd.init=0.7),  
#'      cl=saemixParam(name="cl",mu.init=2, varlevel=c("iiv","iov"), sd.init=c(0.6,0.3), corr = list(iiv=c("ka","vd"),iov=c("vd")), covariate=c("wt","sex","age"), covariate.init=c(0.75,0,0), covariate.estim=c("fixed","estimated","estimated"), corr.init=list(iiv=c(-0.5,0.7), iov=0.7)))
#' saemix.statmodel <- extractVarModel(saemix.parameters)
#' print(saemix.statmodel)
#' 
#' @export 


extractVarModel <- function(parameters, varlevel=c()) {
  xcheck<-checkParameters(parameters=parameters, varlevel=varlevel)
  if(inherits(xcheck,"character")) {
    return("Extraction failed: parameters should be a list of SaemixParameter objects\n")
  }
  parameters <- xcheck$parameters
  varlevel.order<-xcheck$varlevel.order
  logmsg<-xcheck$logmsg
  
  nampar <- names(parameters)
  npar<-length(nampar)
  list.varModel<-vector(length=length(varlevel.order),mode="list")
  names(list.varModel)<-varlevel.order
  for(ilev in 1:length(varlevel.order)) {
    xlev <- varlevel.order[ilev] # name of the variability level 
    omega.model<-omega<-matrix(0, ncol=npar, nrow=npar, dimnames=list(nampar, nampar))
    omega.estim<-matrix("estimated", ncol=npar, nrow=npar, dimnames=list(nampar, nampar))
    for(ipar in 1:npar) {
      # Is there an omega for this parameter at variability level xlev
      if(length(parameters[[ipar]]@varlevel)>0) idx <- grep(xlev, parameters[[ipar]]@varlevel)[1] else idx<-integer(0)
      if(length(idx)>0 && !is.na(idx)) {
        omega.model[ipar,ipar]<-1
        omega[ipar,ipar]<-(parameters[[ipar]]@sd.init[idx])**2
        omega.estim[ipar,ipar]<-parameters[[ipar]]@sd.estim[idx]
        # Are there correlations associated
        if(length(parameters[[ipar]]@corr)>0) idx.corr<-grep(xlev, names(parameters[[ipar]]@corr))[1] else idx.corr<-integer(0)
        if(length(idx.corr)>0 && !is.na(idx.corr)) {
          corpar <- parameters[[ipar]]@corr[[idx.corr]]
          for(j in 1:length(corpar)) {
            jpar<-grep(corpar[j],rownames(omega.model))
            idx.jpar <- grep(xlev, parameters[[jpar]]@varlevel)
            if(length(idx.jpar)>0 && !is.na(idx.jpar)) { # only if the other parameter has variability, associate a covariance and its initial value
              omega.model[ipar,jpar]<-omega.model[jpar,ipar]<-1
              omega[ipar,jpar]<-omega[jpar,ipar]<-((parameters[[ipar]]@corr.init[[idx.corr]][j])**2)*(parameters[[ipar]]@sd.init[idx])*(parameters[[jpar]]@sd.init[idx.jpar])**2
              omega.estim[ipar,jpar]<-omega.estim[jpar,ipar]<-parameters[[ipar]]@corr.estim[[idx.corr]][j]
            }
          }
        } 
      }
    }
    vlev <-  new(Class="SaemixVarModel", omega.model=omega.model, name.level=varlevel.order[ilev])
    vlev@omega <- omega
    vlev@nphi <- nrow(omega)
    vlev@log <- logmsg
    vlev@omega.estim <- omega.estim
    # Sanitise matrices to check conformity
    ## include warning messages
    vlev <- validateSaemixVar(vlev)
    # Add indices
    vlev <- addSaemixIndices(vlev)
    list.varModel[[ilev]] <- vlev
    
  }
  
  return(list.varModel)
}

validateSaemixVar <- function(xlev, verbose=FALSE) {
  logmsg<-""
  if(length(xlev@omega.model)>0) {
    npar<-nrow(xlev@omega.model)
    nampar<-colnames(xlev@omega)
    # remove correlations associated with a parameter without variability
    idx<-which(mydiag(xlev@omega.model)==0)
    if(length(idx)>0) {
      msg<-paste0("Removing correlations (if any) for parameters without variability at level ",xlev@name.level,"\n")
      logmsg<-paste0(logmsg,msg)
      if(verbose) cat(msg)
      xlev@omega.model[,idx]<-xlev@omega.model[idx,]<-0
      xlev@omega[,idx]<-xlev@omega[idx,]<-0
    }
    # add correlation terms if missing (eg a correlated with b correlated with c but missing correlation between a and c)
    for(ipar in 1:npar) {
      setcor<-which(xlev@omega.model[,ipar]==1)
      setcor<-setcor[-c(ipar)]
      if(length(setcor)>1) {
        for(i in 1:length(setcor)) {
          for(j in 2:length(setcor)) {
            if(xlev@omega.model[i,j]==0) {
              msg<-paste0("Added missing correlation between ",nampar[i]," and ",nampar[j],"\n")
              logmsg<-paste0(logmsg,msg)
              if(verbose) cat(msg)
            }
            xlev@omega.model[i,j]<-xlev@omega.model[j,i]<-1
          }
        }
      }
    }
  }
  xlev@log<-paste0(xlev@log,logmsg)
  return(xlev)
}

## idcol.eta: parameters with variability
## idcol.eta.fix: parameters with fixed variability
## idvec.var: variance parameters in param
## idvec.cov: covariance parameters 

## idmat.estim: which var and cov parameters are estimated (from lower triangular matrix omega)
## subomega: initialised to the submatrix of omega
## ToDo: check initial values of omega as >0


setGeneric(name="addSaemixIndices",
           def=function(object,verbose=FALSE) standardGeneric("addSaemixIndices")
)

setMethod("addSaemixIndices",
          signature="SaemixVarModel",
          function(object, verbose=FALSE) {
            object@idcol.eta<-which(mydiag(object@omega.model)==1)
            object@idcol.eta.fix<-which(mydiag(object@omega.model)==1 & mydiag(object@omega.estim=="fixed"))
            object@subomega <- object@omega[object@idcol.eta, object@idcol.eta,drop=FALSE]
            # Parameter names and initial values
            nampar<-colnames(object@omega.model)
            name.rand1<-xpar<-c()
            myidx.var<-myidx.cov<-c()
            ipar<-0
            for(iom in 1:object@nphi) {
              for(jom in iom:object@nphi) {
                if(object@omega.model[iom, jom]==1) {
                  ipar<-ipar+1
                  if(iom==jom) {
                    myidx.var<-c(myidx.var,ipar)
                    xpar<-c(xpar,object@omega[iom,jom])
                    name.rand1<-c(name.rand1,paste("Var",nampar[iom],sep="."))
                  }
                  else {
                    myidx.cov<-c(myidx.cov,ipar)
                    xpar<-c(xpar,object@omega[iom,jom])
                    name.rand1<-c(name.rand1,paste("Cov",nampar[iom],nampar[jom],sep="."))
                  }
                }
              }
            }
            if(length(myidx.var)>0) {
              object@param.names <- name.rand1
              object@param <- xpar
              xvec<-object@omega.estim[lower.tri(object@omega.estim, diag=TRUE)]
              object@idvec.var <-myidx.var
              object@idmat.estim<-c(myidx.var, myidx.cov)
            }
            if(length(myidx.cov)>0) object@idvec.cov <-myidx.cov
            return(object)
          }
)

