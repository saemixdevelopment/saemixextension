
#' Class "SaemixVarLevel"
#' 
#' An object of the SaemixVarLevel class, representing the variability for a given level of random effects
#' It contains the model itself (the structure of the variance-covariance matrix), which elements are
#' fixed (not estimated) and the initial values for the variance-covariance matrix parameters.
#' 
#' @name SaemixVarLevel-class 
#' @docType class
#' @aliases SaemixVarLevel SaemixVarLevel-class 
#' @aliases print,SaemixVarLevel showall,SaemixVarLevel show,SaemixVarLevel
#' 
#' @section Objects from the Class: 
#' An object of the SaemixVarLevel class contains the following slots:
#' @slot name.level Object of class \code{"character"}: name associated with the variability level. If not given, will use the name of the variable the variability is associated to (eg: id, occ, ...)
#' @slot variable name of the variable (in the dataset) to which the variability is associated (eg: id, occ, ...)
#' @slot omega a variance-covariance matrix giving the variances and covariances of the random effects associated with the current variability level. The matrix should be block positive definite (ie variances may be missing from the model, see omega.model, but the blocks containing elements present in the model should be positive definite)
#' @slot omega.model a matrix of 0/1 giving the structure of the variance-covariance matrix (1 indicates elements present in the model, 0 denotes elements not in the model). Defaults to the identity matrix.
#' @slot omega.model.fix a matrix of 0/1 indicating elements of the variance-covariance matrix which are fixed (not estimated). Defaults to a matrix of 0's.
#' @slot omega.tri (internal) lower triangular matrix in vector form
#' @slot omega.names (internal) names of the variance-covariance parameters, in the same order as in omega.tri
#' @slot idxmat.triomega (internal) index of estimated elements in omega.tri
#' @slot idxmat.triomega.fix (internal) index of estimated elements in omega.tri that have been fixed by the user
#' @slot idxmat.triomega.var (internal) index of elements in omega.tri corresponding to variance components
#' @slot idxmat.triomega.covar (internal) index of elements in omega.tri corresponding to covariance components
#' @slot index.eta (internal) index of parameters with IIV (1's in omega.model)
#' @section Methods:
#'   \describe{
#'     \item{[<-}{\code{signature(x = "SaemixVarLevel")}: replace elements of object}
#'     \item{[}{\code{signature(x = "SaemixVarLevel")}: access elements of object}
#'     \item{initialize}{\code{signature(.Object = "SaemixVarLevel")}: internal function to initialise object, not to be used}
#'     \item{print}{\code{signature(x = "SaemixVarLevel")}: prints details about the object (more extensive than show)}
#'     \item{showall}{\code{signature(object = "SaemixVarLevel")}: shows all the elements in the object}
#'     \item{show}{\code{signature(object = "SaemixVarLevel")}: prints details about the object}
#' 	 }
#' 
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}
#' 
#' @seealso \code{\link{saemixData}} \code{\link{SaemixModel}} \code{\link{saemixControl}} \code{\link{saemix}}
#' @examples
#' showClass("SaemixVarLevel")
#' 
#' @keywords classes
#' @exportClass SaemixVarLevel


# Variability level class - generic
setClass(Class = "SaemixVarLevel",
         representation=representation(
           name.level = "character", # name of variability level
           variable = "character", # which variable (in the dataset) is the variability associated to
           omega = "matrix", # variance-covariance matrix (values)
           chol.omega = "matrix", # Cholesky decomposition of omega
           omega.model = "matrix", # structure of the variance-covariance matrix (1=present in the model, 0=absent)
           omega.model.fix = "matrix", # fixed elements of the variance-covariance matrix (1=fixed, 0=estimated)
           omega.tri = "numeric", # lower triangular matrix of omega, in vector form
           omega.names = "character", # names of the variance-covariance parameters in omega.tri
           # indices in omega (columns/lines)
           index.omega.novar = "integer", # index of parameters without IIV ## i0.omega2
           index.omega.var = "integer", # index of parameters with IIV ## i1.omega2
           # --- Indices in matrices
           # indices in omega
           idxmat.omega = "integer",  # index of 1's in omega.model ## indest.omega           
           # indices in triangular matrices
           idxmat.triomega = "integer", # index of estimated elements in omega.tri (=1 in omega.model)
           idxmat.triomega.fix = "integer", # index of fixed elements in omega.tri (=1 in omega.model.fix)
           idxmat.triomega.var = "integer", # index of variances in omega.tri
           idxmat.triomega.covar = "integer", # index of covariances in omega.tri
           index.eta = "integer" #  index of parameters with estimated variances
         ),
         validity=function(object){
           # Check all sizes are commensurate & check symmetry using validate.covariance.matrix TODO
           validate.covariance.model(object@omega.model)
           return(TRUE)
         }
)

setMethod( 
  f="initialize",
  signature="SaemixVarLevel",
  definition=function(.Object, size=1, name.level="", variable="", omega=NULL, omega.model=NULL, omega.model.fix=NULL){
    if(name.level=="") {
      if(variable=="") {
        name.level<-"iiv"
        variable<-"id"
      } else {
        name.level<-paste0("var.",variable)
        if(tolower(variable)=="id") name.level<-"iiv"
        if(tolower(variable)=="occ") name.level<-"iov"
      }
    }
    .Object@name.level <- name.level
    if(variable=="") variable<-"id"
    .Object@variable <- variable
    if(!is.null(omega)) {  # check omega is conform if given
      if(!(is.matrix(omega)) || nrow(omega)!=ncol(omega)) {
        message(paste0("Variance matrix omega not conform, ignoring and setting size to ",size))
        omega<-NULL
      }
    } 
    if(!is.null(omega)) {
      size<-ncol(omega)
      if(!is.null(omega.model) && validate.covariance.model(omega.model)) {
        if(ncol(omega.model)==size) omega <- omega*(omega.model) else omega.model<-1*(omega!=0) # if right size, multiply omega and omega.model, if omega.model isn't the same size, 
      } else omega.model<-1*(omega!=0)
    } else { # omega not given or not conform but omega.model given
      if(!is.null(omega.model) && validate.covariance.model(omega.model)) { # if omega.model is conform, adjust omega
        size<-ncol(omega.model)
        omega<-mydiag(size)
      } else omega.model<-NULL
    }
    if(size<1 & is.null(omega) & is.null(omega.model)) {
      message("Please give at least one of size (size of variance model), omega (value of variance-covariance matrix, a symmetric square matrix) or omega.model (a non-negative symmetric square matrix), returning NULL")
      return()
    }
    
    if(is.null(omega.model.fix)) .Object@omega.model.fix<-diag(x=0, nrow=size, ncol=size) else .Object@omega.model.fix<-omega.model.fix
    if(is.null(omega.model)) .Object@omega.model<-diag(x=1, nrow=size, ncol=size) else .Object@omega.model<-omega.model
    .Object@index.eta<-which(mydiag(.Object@omega.model)>0)
    
    .Object@index.omega.var <- which(mydiag(omega.model)>0) ## i1.omega2
    .Object@index.omega.novar <- which((1-mydiag(omega.model))>0) ## i0.omega2
    .Object@idxmat.omega <-  which(omega.model>0) # index of 1's in omega.model ## indest.omega
    
    if(is.null(omega)) omega<-diag(x=1, nrow=size, ncol=size) 
    # TBD:
    # omega<-omega*.Object@omega.model
    .Object@omega<-omega  
    .Object@omega.tri <- omega[lower.tri(omega, diag=TRUE)] # vector of the lower triangular elements, including the diagonal
    .Object@idxmat.triomega <- which(.Object@omega.model[lower.tri(.Object@omega.model, diag=TRUE)]==1) # elements of omega.tri in the model
    .Object@idxmat.triomega.fix <- which(.Object@omega.model.fix[lower.tri(.Object@omega.model.fix, diag=TRUE)]==1) # elements of omega.tri fixed by the user
    # Elements corresponding to variances
    mat1<-vec2mat(1:length(.Object@omega.tri))
    idx<-diag(mat1)
    .Object@idxmat.triomega.var<-intersect(.Object@idxmat.triomega,as.integer(idx))
    # Elements corresponding to covariances
    idx<-mat1[lower.tri(mat1)]
    .Object@idxmat.triomega.covar<-intersect(.Object@idxmat.triomega,as.integer(idx))
    
    validObject(.Object)
    return(.Object)
  }
)


# Getteur
setMethod(
  f ="[",
  signature = "SaemixVarLevel" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "name.level"={return(x@name.level)},
            "variable"={return(x@variable)},
            "omega"={return(x@omega)},
            "chol.omega"={return(x@chol.omega)},
            "omega.model"={return(x@omega.model)},
            "omega.model.fix"={return(x@omega.model.fix)},
            "omega.tri"={return(x@omega.tri)},
            "omega.names"={return(x@omega.names)},
            "index.omega.novar"={return(x@index.omega.novar)},
            "index.omega.var"={return(x@index.omega.var)},
            "idxmat.omega"={return(x@idxmat.omega)},
            "idxmat.triomega"={return(x@idxmat.triomega)},
            "idxmat.triomega.tri"={return(x@idxmat.triomega.tri)},
            "idxmat.triomega.fix"={return(x@idxmat.triomega.fix)},
            "idxmat.triomega.var"={return(x@idxmat.triomega.var)},
            "idxmat.triomega.covar"={return(x@idxmat.triomega.covar)},
            "index.eta"={return(x@index.eta)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixVarLevel" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "name.level"={x@name.level<-value},
            "variable"={x@variable<-value},
            "omega"={x@omega<-value},
            "chol.omega"={x@chol.omega<-value},
            "omega.model"={x@omega.model<-value},
            "omega.model.fix"={x@omega.model.fix<-value},
            "omega.tri"={x@omega.tri<-value},
            "omega.names"={x@omega.names<-value},
            "index.omega.novar"={x@index.omega.novar<-value},
            "index.omega.var"={x@index.omega.var<-value},
            "idxmat.omega"={x@idxmat.omega<-value},
            "idxmat.triomega"={x@idxmat.triomega<-value},
            "idxmat.triomega.tri"={x@idxmat.triomega<-value},
            "idxmat.triomega.fix"={x@idxmat.triomega.fix<-value},
            "idxmat.triomega.var"={x@idxmat.triomega.var<-value},
            "idxmat.triomega.covar"={x@idxmat.triomega.covar<-value},
            "index.eta"={x@index.eta<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

############################### function to create a SaemixVarLevel 
# Resizing a second matrix based on a first (helper function)
resizeMatrix<-function(mat1, mat2) {
  size<-dim(mat1)[1]
  if(dim(mat2)[1]>size) mat2<-mat2[1:size,1:size] else {
    mat21<-diag(nrow=size, ncol=size)
    mat21[1:dim(mat2)[1],1:dim(mat2)[1]]<-mat2
    mat2<-mat21
  }
  mat2
}

# Manipulating matrices to vector form (by column)
vec2mat <- function(x){
  p <- sqrt(1 + 8 * length(x))/ 2 - 0.5
  m <- matrix(0, p, p)
  m[lower.tri(m, diag=TRUE)] <- x
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}
# Note: faster version using same functions from Rfast if needed (but mostly for large matrices)



# creator for an SaemixVarLevel object
saemixVarModel<-function(size=1, name.level="", variable="", omega=NULL, omega.model=NULL, omega.model.fix=NULL, parameter.names=NULL, verbose=FALSE) {
  # the size of omega or omega.model takes precedence over size
  if(!is.null(omega)) {
    if(size!=dim(omega)[1]) {
      if(verbose & size!=1) message("Changing size to adjust the number of parameters to the size of omega\n")
      size<-dim(omega)[1]
    }
  } else {
    if(!is.null(omega.model)) {
      if(size!=dim(omega.model)[1]) {
        if(verbose & size!=1) message("Changing size to adjust the number of parameters to the size of omega.model\n")
        size<-dim(omega.model)[1]
      }
    }
  }
  if(!is.null(omega.model) & !is.null(omega)) { # Checking consistency between omega and omega.model
    omega.model1<-resizeMatrix(omega,omega.model)
    omega.model1[omega!=0]<-1
    if(!identical(omega.model1, omega.model)) {
      if(verbose) message("Mismatch between omega.model and omega, adjusting omega.model to the non-zero elements in omega")
      omega.model<-omega.model1
    }
  }
  if(is.null(omega.model) & !is.null(omega)) {
    # omega.model defaults to the non-null elements in omega
    omega.model<-omega
    omega.model[omega!=0]<-1
    if(verbose) message("Inferring covariance model from omega")
  }
  # Check the sizes of omega, omega.model and omega.model.fix are the same (if given)
  # Adjusting w/r to omega (if given) or to omega.model (if omega not given and omega.model is given)
  if(!is.null(omega.model) & !is.null(omega)) {
    if(dim(omega.model)[1]!=dim(omega)[1]) {
      if(verbose) message("Adjusting size of omega.model to the size of omega, please check the input\n")
      omega.model<-resizeMatrix(omega, omega.model)
    }
  }
  if(!is.null(omega.model.fix) & !is.null(omega.model)) {
    if(dim(omega.model.fix)[1]!=dim(omega.model)[1]) {
      if(verbose) message("Adjusting size of omega.model.fix to the size of omega.model, please check the input\n")
      omega.model.fix<-resizeMatrix(omega.model, omega.model.fix)
    }
  }
  if(is.null(omega.model) & !is.null(omega) & !is.null(omega.model.fix)) {
    if(dim(omega)[1]!=dim(omega.model.fix)[1]) {
      if(verbose) message("Adjusting size of omega.model.fix to the size of omega, please check the input\n")
      omega.model.fix<-resizeMatrix(omega, omega.model.fix)
    }
  }
  x<-new(Class="SaemixVarLevel", size=size, name.level=name.level, variable=variable, omega=omega, omega.model=omega.model, omega.model.fix=omega.model.fix)
  if(!is.null(parameter.names) & length(parameter.names==size)) 
    x<-saemixVarNames(x, parameter.names)
  return(x)
}

# Attribution of parameter names
saemixVarNames<-function(object, parameter.names) {
  if(missing(parameter.names))
    parameter.names<-paste0("theta",1:dim(object@omega)[1])
  npar<-length(parameter.names)
  varnames<-c()
  for(jnam in 1:npar) {
    for(inam in jnam:npar) {
      if(inam==jnam) varnames<-c(varnames, paste0("omega2.",parameter.names[inam])) else 
        varnames<-c(varnames, paste0("cov.",parameter.names[jnam],".",parameter.names[inam]))
    }
  }
  object@omega.names<-varnames
  return(object)
}



############################### Show/print
setMethod("show","SaemixVarLevel",
          function(object) {
            cat("Variability level:",object@name.level,"(associated with",object@variable,")\n")
            cat("    variance-covariance model\n")
            print(object@omega.model)
          }
)

setMethod("showall","SaemixVarLevel",
          function(object) {
            cat("Variability level:",object@name.level,"(associated with",object@variable,")\n")
            cat("    model\n")
            print(object@omega.model)
            cat("    variance-covariance matrix\n")
            print(object@omega.model*object@omega)
            if(sum(object@omega.model.fix)>0 & length(object@omega.names)==0) {
              cat("    some elements are fixed:\n")
              print(object@omega.model.fix)
            }
            if(length(object@omega.names)>0) {
              cat("    estimated parameters: ")
              if(length(object@idxmat.triomega.fix)>0) cat(object@omega.names[object@idxmat.triomega[!(object@idxmat.triomega %in% c(object@idxmat.triomega.fix))]],"\n") else  cat(object@omega.names[object@idxmat.triomega],"\n")
              if(length(object@idxmat.triomega.fix)>0)
                cat("    fixed parameters: ",object@omega.names[object@idxmat.triomega[(object@idxmat.triomega %in% c(object@idxmat.triomega.fix))]],"\n")
            }
          }
)

setMethod("print","SaemixVarLevel",
          function(x,nlines=10,...) {
            show(x)
          }
)

########################################################  Auxiliary functions
# Redefining diag function, too many problems with the R version

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

########################################################  Validity of covariance model

# Testing the Validity of covariance model

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
          if(verbose) message("Error initialising SaemixModel object:\n  Covariances can only be included between 2 parameters with variances.\n")
          return(FALSE)
        }
      }
    }
  }
  # Check that the matrix has a block structure - doesn't work, fails for simple block :-/
  # indx1<-which(diag(x)==0)
  # if(length(indx1)>0) x1<-x[-c(indx1),-c(indx1)] else x1<-x # removing the lines without variances
  # could maybe work by changing the off-diagonal elements to 0.5 instead of 1...
  # xchol<-try(chol(x1))
  # if(is(xchol, 'try-error')) {
  #   if(verbose) message("Error initialising SaemixModel object:\n  Covariance matrices should be block-diagonal.\n")
  #   return(FALSE)
  # }
  return(TRUE)
}

# TODO: a function to make sure levels of variability are nested and compatible with one another

