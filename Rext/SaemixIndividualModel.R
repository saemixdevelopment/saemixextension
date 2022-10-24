# More structured definition of an individual model matching covariate model and data for a given variability level
# defining design variables and setting indices to access items
# in the end, should be created by a function associating a variability model, a covariate model (par-cov relationship and cov transformation) and data (find cov, find level of variability)
# Contains most of what used to be Uargs

# SaemixVarLevel
# name.level = "character", # name of variability level
# variable = "character", # which variable (in the dataset) is the variability associated to
# omega = "matrix", # variance-covariance matrix (values)
# omega.model = "matrix", # structure of the variance-covariance matrix (1=present in the model, 0=absent)
# omega.model.fix = "matrix", # fixed elements of the variance-covariance matrix (1=fixed, 0=estimated)
# omega.tri = "numeric", # lower triangular matrix of omega, in vector form
# omega.names = "character", # names of the variance-covariance parameters in omega.tri
# idxmat.triomega = "numeric", # index of estimated elements in omega.tri (=1 in omega.model)
# idxmat.triomega.fix = "numeric", # index of fixed elements in omega.tri (=1 in omega.model.fix)
# idxmat.triomega.var = "numeric", # index of variances in omega.tri
# idxmat.triomega.covar = "numeric", # index of covariances in omega.tri
# index.eta = "numeric" #  index of parameters with estimated variances


# Variability level class - generic
setClass(Class = "SaemixIndivModel",
         contains = "SaemixVarLevel",
         representation=representation(
           # --- names
           name.modpar = "character", # name of parameters in the model (psi, ex: ka, cl)
           name.fixedpar = "character", # name of fixed parameters (mu+beta) for the current level of variability
           name.covariates = "character", # name of covariates in the model for this level of variability
           # --- number of parameters, covariates,...
           N="integer", # number of 'subjects' associated to this variability level (eg: N=nb of subjects for IIV, N=sum(n_{occ,i}) for IOV,...)
           nb.modpar="integer", # number of model parameters (eg ka, cl) = number of mu ## nb.parameters
           nb.etas = "integer", # number of parameters with random effects (non-null elements in diag(omega))
           nb.fixedpar = "integer", # number of parameters mu+beta (mu: at highest level, eg IIV=population parameter, then eg systematic(fixed) occasion effect) ## was nb.betas
           nb.cov = "integer", # number of covariates (=nrow(covariate.model)) ## TODO: transform covariates, possibly adding columns
           nb.mubeta.optim = "integer", # nb of fixed parameters to optimise on
           # --- parameters
           fixedpar="numeric", # value for mu+beta
           # --- parameter and covariate models
           covariate.model="matrix",	# covariate model (covariate design matrix)
           covariate.model.fix="matrix",	# fixed parameters in the covariate model (defaults to all 0's)
           # mu+covariate model (not sure if still needed), maybe remove
           betaest.model="matrix",	# design matrix (row for mu then covariate.model)
           # for IIV: rbind(rep(1, nb.modpar), covariate.model); for other levels, covariate.model or rbind(rep(0, nb.modpar), covariate.model) ?
           # --- design matrices
           Mcovariates="matrix", # Mcovariate: matrix with N rows and columns corresponding to the (transformed) covariate values for each parameter x covariate relationship (including mu for IIV, a column of 1's)
           LCOV="matrix", # LCOV: (nb.fixedpar x nb.modpar) design matrix to pass from the vector of fixed parameters (mu+beta) to the model parameters; used in M-step to compute d1.omega
           MCOV="matrix", # MCOV: values of fixed parameters in matrix form (same structure as LCOV); used to compute mean.phi (population value of phi accounting for individual covariates)
           COV="matrix", # N x nb(mu+beta) containing 1's for mu and covariate values for each subject;  used in M-step to compute mean.phi as MCOV %*% COV (phi = MCOV %*% COV + eta)
           # --- helper matrices
           #    COV1="matrix", # initialised as Uargs$COV[,Uargs$ind.fix1] = Uargs$COV[,index.fixedpar.fix] # maybe remove
           # modified in main to Uargs$COV[,Uargs$ind.fix11] = Uargs$COV[,index.fixedpariiv.estim] when flag.fmin is TRUE and kiter=saemix.options$nbiter.sa (ie end of simulated annealing and minimisation) => removed, make the change in M-step
           COV2="matrix", # used in M-step to compute comega; square matrix of size n(mu+beta)
           # block structure consisting of N (for cells muxmu), sum_i cov_i^(c) (for cells mu x beta) and sum_i cov_i^(c) cov_i(c) (c' can be equal to c, for cells beta x beta)
           # maybe remove
           dstatCOV="matrix", # used in M-step to compute betas[Uargs$ind.fix1]= elements of (mu+beta) fixed to their starting value
           # defined as COV[,index.fixedpar.fix,drop=FALSE]%*%MCOV[index.fixedpar.fix,,drop=FALSE]
           # TBD if we keep, could use the full form, also not sure about the indices :-/
           MCOV0="matrix", # used in compute.Uy as args$MCOV0[args$j0.covariate]<-b0
           # phi0<-args$COV0 %*% args$MCOV0
           # defined in initialisation as MCOV[index.fixedparnoiiv.estim,index.omega.novar,drop=FALSE]
           COV0="matrix", # used in compute.Uy, COV[,index.fixedparnoiiv.estim,drop=FALSE]
           # # phi0<-args$COV0 %*% args$MCOV0 # same line as MCOV0
           # --- indices in vectors
           # indices in fixedpar
           index.mu = "integer", # index of mu's in fixedpar ## indx.betaI 
           index.beta = "integer", # index of beta's in fixedpar ## indx.betaC
           index.fixedpar.fix = "integer",   # fixed fixed effects (fixed elements of mu+beta) ## ind.fix0
           index.fixedpar.estim= "integer",  # estimated fixed effects (estimated elements of mu+beta) ## ind.fix1
           index.fixedpariiv.estim= "integer",  # Index of beta effects on parameters with IIV ## ind.fix11
           index.fixedparnoiiv.estim= "integer",  # Index of beta effects on parameters without IIV ## ind.fix10
           # --- indices in matrices
            # indices in betaest.model
           idxmat.fixedpar = "integer", # index of fixedpar in betaest.model ## ind.covariates
           idxmat.mu = "integer", # index of mu's in betaest.model 
           idxmat.beta = "integer", # index of beta's in betaest.model 
           # indices in LCOV/MCOV
           idxmat.mcov.par = "integer", #  index of elements in LCOV present in the model ##  j.covariate
           idxmat.mcov.fixedpar.optim ="integer" #  index in the matrix MCOV (or LCOV) of fixed parameters (mu+beta) on parameters without IIV that are fixed to their starting value, to optimise ## j0.covariate
         ),
         validity=function(object){
           # Check all sizes are commensurate & check symmetry using validate.covariance.matrix TODO
           validate.covariance.model(object@omega.model)
           return(TRUE)
         }
)


# Getteur
setMethod(
  f ="[",
  signature = "SaemixIndivModel" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            # From SaemixVarLevel
            "name.level"={return(x@name.level)},
            "variable"={return(x@variable)},
            "omega"={return(x@omega)},
            "omega.model"={return(x@omega.model)},
            "omega.model.fix"={return(x@omega.model.fix)},
            "chol.omega"={return(x@chol.omega)},
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
            # Elements added in the child class
            "name.modpar"={return(x@name.modpar)},
            "name.fixedpar"={return(x@name.fixedpar)},
            "name.covariates"={return(x@name.covariates)},
            "nb.modpar"={return(x@nb.modpar)},
            "N"={return(x@N)},
            "ncov"={return(x@ncov)},
            "nb.etas"={return(x@nb.etas)},
            "nb.fixedpar"={return(x@nb.fixedpar)},
            "fixedpar"={return(x@fixedpar)},
            "covariate.model"={return(x@covariate.model)},
            "covariate.model.fix"={return(x@covariate.model.fix)},
            "betaest.model"={return(x@betaest.model)},
            "Mcovariates"={return(x@Mcovariates)},
            "LCOV"={return(x@LCOV)},
            "MCOV"={return(x@MCOV)},
            "COV"={return(x@COV)},
            "COV2"={return(x@COV2)},
            "dstatCOV"={return(x@dstatCOV)},
            "MCOV0"={return(x@MCOV0)},
            "COV0"={return(x@COV0)},
            "index.mu"={return(x@index.mu)},
            "index.beta"={return(x@index.beta)},
            "index.fixedpar.fix"={return(x@index.fixedpar.fix)},
            "index.fixedpar.estim"={return(x@index.fixedpar.estim)},
            "index.fixedpariiv.estim"={return(x@index.fixedpariiv.estim)},
            "index.fixedparnoiiv.estim"={return(x@index.fixedparnoiiv.estim)},
            "idxmat.fixedpar"={return(x@idxmat.fixedpar)},
            "idxmat.mu"={return(x@idxmat.mu)},
            "idxmat.beta"={return(x@idxmat.beta)},
            "idxmat.mcov.par"={return(x@idxmat.mcov.par)},
            "idxmat.mcov.fixedpar.optim"={return(x@idxmat.mcov.fixedpar.optim)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixIndivModel" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "name.level"={x@name.level<-value},
            "variable"={x@variable<-value},
            "omega"={x@omega<-value},
            "omega.model"={x@omega.model<-value},
            "omega.model.fix"={x@omega.model.fix<-value},
            "chol.omega"={x@chol.omega<-value},
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
            # Elements added in the child class
            "nb.modpar"={x@nb.modpar<-value},
            "ncov"={x@ncov<-value},
            "N"={x@N<-value},
            "nb.etas"={x@nb.etas<-value},
            "nb.fixedpar"={x@nb.fixedpar<-value},
            "name.modpar"={x@name.modpar<-value},
            "name.fixedpar"={x@name.fixedpar<-value},
            "name.covariates"={x@name.covariates<-value},
            "fixedpar"={x@fixedpar<-value},
            "covariate.model"={x@covariate.model<-value},
            "covariate.model.fix"={x@covariate.model.fix<-value},
            "betaest.model"={x@betaest.model<-value},
            "Mcovariates"={x@Mcovariates<-value},
            "LCOV"={x@LCOV<-value},
            "MCOV"={x@MCOV<-value},
            "COV"={x@COV<-value},
            "COV2"={x@COV2<-value},
            "dstatCOV"={x@dstatCOV<-value},
            "MCOV0"={x@MCOV0<-value},
            "COV0"={x@COV0<-value},
            "index.mu"={x@index.mu<-value},
            "index.beta"={x@index.beta<-value},
            "index.fixedpar.fix"={x@index.fixedpar.fix<-value},
            "index.fixedpar.estim"={x@index.fixedpar.estim<-value},
            "index.fixedpariiv.estim"={x@index.fixedpariiv.estim<-value},
            "index.fixedparnoiiv.estim"={x@index.fixedparnoiiv.estim<-value},
            "idxmat.fixedpar"={x@idxmat.fixedpar<-value},
            "idxmat.mu"={x@idxmat.mu<-value},
            "idxmat.beta"={x@idxmat.beta<-value},
            "idxmat.mcov.par"={x@idxmat.mcov.par<-value},
            "idxmat.mcov.fixedpar.optim"={x@idxmat.mcov.fixedpar.optim<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

# 
# setMethod( 
#   f="initialize",
#   signature="SaemixIndivModel",
#   definition=function(.Object, size=1, name.level="", variable="", omega=NULL, omega.model=NULL, omega.model.fix=NULL){
#     .Object <- callNextMethod(.Object, size=size, name.level=name.level, variable=variable, omega=omega, omega.model=omega.model, omega.model.fix=omega.model.fix)
#     # .Object@omega.tri <- omega[lower.tri(omega, diag=TRUE)] # vector of the lower triangular elements, including the diagonal
#     # .Object@idxmat.triomega <- which(.Object@omega.model[lower.tri(.Object@omega.model, diag=TRUE)]==1) # elements of omega.tri in the model
#     # .Object@idxmat.triomega.fix <- which(.Object@omega.model.fix[lower.tri(.Object@omega.model.fix, diag=TRUE)]==1) # elements of omega.tri fixed by the user
#     # # Elements corresponding to variances
#     # mat1<-vec2mat(1:length(.Object@omega.tri))
#     # idx<-diag(mat1)
#     # .Object@idxmat.triomega.var<-intersect(.Object@idxmat.triomega,idx)
#     # # Elements corresponding to covariances
#     # idx<-mat1[lower.tri(mat1)]
#     # .Object@idxmat.triomega.covar<-intersect(.Object@idxmat.triomega,idx)
#     # 
#     validObject(.Object)
#     return(.Object)
#   }
# )

setMethod( 
  f="initialize",
  signature="SaemixIndivModel",
  definition=function(.Object, var.model, covariate.model=matrix(nrow=0, ncol=0), covariate.model.fix=matrix(nrow=0, ncol=0)){
    # Parameters and covariance model
    .Object@name.level <- var.model@name.level
    .Object@variable <- var.model@variable
    .Object@omega <- var.model@omega
    .Object@omega.model <- var.model@omega.model
    .Object@omega.model.fix <- var.model@omega.model.fix
    .Object@nb.modpar <- ncol(var.model@omega.model)
    if(!is.null(colnames(var.model@omega))) .Object@name.modpar <- colnames(var.model@omega)
    if(length(.Object@name.modpar)==0 & !is.null(colnames(var.model@omega.model))) .Object@name.modpar <- colnames(var.model@omega.model)
    if(length(var.model@omega.names)>0) .Object@omega.names<-var.model@omega.names
    # Indices filled in from variability model
    .Object@idxmat.triomega <- var.model@idxmat.triomega
    .Object@idxmat.triomega.var <- var.model@idxmat.triomega.var
    .Object@idxmat.triomega.covar <- var.model@idxmat.triomega.covar
    .Object@idxmat.triomega.fix <- var.model@idxmat.triomega.fix
    .Object@idxmat.omega <- var.model@idxmat.omega
    .Object@omega <- var.model@omega
    .Object@chol.omega <- var.model@chol.omega
    .Object@index.omega.var <-var.model@index.omega.var
    .Object@index.omega.novar <-var.model@index.omega.novar
    # Covariate model
    if(length(covariate.model)>0) {
      if(ncol(covariate.model)!=.Object@nb.modpar) {
        message("Mismatch between covariate.model and var.model, ignoring covariate.model\n")
        ncov<-0
        } else {
        .Object@covariate.model <- covariate.model
        if(!is.null(rownames(covariate.model))) .Object@name.covariates <- rownames(covariate.model)
        ncov<-nrow(covariate.model)
        if(length(covariate.model.fix)==0) covariate.model.fix <- covariate.model*0
        if(ncol(covariate.model.fix)!=.Object@nb.modpar || nrow(covariate.model.fix)!=nrow(covariate.model)) {
          message("Mismatch between covariate.model and covariate.model.fix\n") 
          covariate.model.fix <- covariate.model*0
        } 
        .Object@covariate.model.fix <- covariate.model.fix
        }
    } else ncov<-0
    .Object@nb.cov <- as.integer(ncov)
    
    validObject(.Object)
    return(.Object)
  }
)


############################### Show/print
setMethod("show","SaemixIndivModel",
          function(object) {
            cat("Variability level:",object@name.level,"(associated with",object@variable,")\n")
            cat("    variance-covariance model\n")
            print(object@omega.model)
            if(length(object@covariate.model)>0)
              print(object@covariate.model)
          }
)

setMethod("showall","SaemixIndivModel",
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

setMethod("print","SaemixIndivModel",
          function(x,nlines=10,...) {
            show(x)
          }
)
