
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
#' @slot index.omega (internal) index of estimated elements in omega.tri
#' @slot index.omega.fix (internal) index of estimated elements in omega.tri that have been fixed by the user
#' @slot index.omega.var (internal) index of elements in omega.tri corresponding to variance components
#' @slot index.omega.covar (internal) index of elements in omega.tri corresponding to covariance components
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
           omega.model = "matrix", # structure of the variance-covariance matrix (1=present in the model, 0=absent)
           omega.model.fix = "matrix", # fixed elements of the variance-covariance matrix (1=fixed, 0=estimated)
           omega.tri = "numeric", # lower triangular matrix of omega, in vector form
           omega.names = "character", # names of the variance-covariance parameters in omega.tri
           index.omega = "numeric", # index of estimated elements in omega.tri (=1 in omega.model)
           index.omega.fix = "numeric", # index of fixed elements in omega.tri (=1 in omega.model.fix)
           index.omega.var = "numeric", # index of variances in omega.tri
           index.omega.covar = "numeric" # index of covariances in omega.tri
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
    if(is.null(omega.model.fix)) .Object@omega.model.fix<-diag(x=0, nrow=size, ncol=size) else .Object@omega.model.fix<-omega.model.fix
    if(is.null(omega.model)) .Object@omega.model<-diag(x=1, nrow=size, ncol=size) else .Object@omega.model<-omega.model
    if(is.null(omega)) omega<-diag(x=1, nrow=size, ncol=size) 
    # TBD:
    # omega<-omega*.Object@omega.model
    .Object@omega<-omega
    .Object@omega.tri <- omega[lower.tri(omega, diag=TRUE)] # vector of the lower triangular elements, including the diagonal
    .Object@index.omega <- which(.Object@omega.model[lower.tri(.Object@omega.model, diag=TRUE)]==1) # elements of omega.tri in the model
    .Object@index.omega.fix <- which(.Object@omega.model.fix[lower.tri(.Object@omega.model.fix, diag=TRUE)]==1) # elements of omega.tri fixed by the user
    # Elements corresponding to variances
    mat1<-vec2mat(1:length(.Object@omega.tri))
    idx<-diag(mat1)
    .Object@index.omega.var<-intersect(.Object@index.omega,idx)
    # Elements corresponding to covariances
    idx<-mat1[lower.tri(mat1)]
    .Object@index.omega.covar<-intersect(.Object@index.omega,idx)
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
            "omega.model"={return(x@omega.model)},
            "omega.model.fix"={return(x@omega.model.fix)},
            "omega.tri"={return(x@omega.tri)},
            "omega.names"={return(x@omega.names)},
            "index.omega"={return(x@index.omega)},
            "index.omega.tri"={return(x@index.omega.tri)},
            "index.omega.fix"={return(x@index.omega.fix)},
            "index.omega.var"={return(x@index.omega.var)},
            "index.omega.covar"={return(x@index.omega.covar)},
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
            "omega.model"={x@omega.model<-value},
            "omega.model.fix"={x@omega.model.fix<-value},
            "omega.tri"={x@omega.tri<-value},
            "omega.names"={x@omega.names<-value},
            "index.omega"={x@index.omega<-value},
            "index.omega.tri"={x@index.omega<-value},
            "index.omega.fix"={x@index.omega.fix<-value},
            "index.omega.var"={x@index.omega.var<-value},
            "index.omega.covar"={x@index.omega.covar<-value},
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
              if(length(object@index.omega.fix)>0) cat(object@omega.names[object@index.omega[!(object@index.omega %in% c(object@index.omega.fix))]],"\n") else  cat(object@omega.names[object@index.omega],"\n")
              if(length(object@index.omega.fix)>0)
                cat("    fixed parameters: ",object@omega.names[object@index.omega[(object@index.omega %in% c(object@index.omega.fix))]],"\n")
            }
          }
)

setMethod("print","SaemixVarLevel",
          function(x,nlines=10,...) {
            show(x)
          }
)


