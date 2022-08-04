###########################################################
# Class SaemixData
###########################################################
setClass(Class="SaemixData",
  representation=representation(
    name.data="character",	# name of dataset
    messages="logical",		# whether to print messages when creating the object
    outcome="list", # data outcomes (a vector or list of elements with type SaemixOutcome including name, type and units; defaults to y1, y2,... with units="" if not given)
    covariates="list", # covariates in the dataset (a list of elements with type SaemixCovariate including name, type and units)
    # names
    name.group="character",	# name of column with ID
    name.predictors="character",# name of column(s) with predictors 
    name.response="character",	# name of column with response
    name.covariates="character",# name of column(s) with covariates [name.covariates=names(object@covariates)]
    name.X="character",		# name of predictor used on X axis for graphs
    name.mdv="character", # name of column indicating a missing response
    name.cens="character", # name of column indicating a censored response
    name.occ="character", # name of column with the occasion
    name.ytype="character", # name of column with type of response (1,2,...)
    #    trans.cov="list",	# a list of transformations applied to the covariates [REMOVE, now associated to the covariate model in the model object]
    units="character",		# units for predictors (a vector)
    # data
    data="data.frame",		# the data (data frame with columns name.group (subject id), index (id renamed to 1:N), name.predictors (predictors), name.response (possibly transformed during fit), name.covariates), mdv (missing data), cens (censored data, 1=censored & value in column response is the LOQ, ytype (type of the response), occ (occasion)); binary covariates are modified to 0/1
    ocov="data.frame",		# original covariates from the dataset, before transformation
    yorig="numeric",		# vector of responses in original dataset
    # indices for the variability levels
    var.match="matrix",  # matrix to match the different levels of variability to the corresponding lines in the dataset
    # for ivarlevel, defined as rep(1:Nunit[ivarlevel], times=nvar.obs[[ivarlevel]])
    # maybe remove
    ind.gen="logical",	# vector of booleans (same size as name.covariates); TRUE=genetic covariate, FALSE=non-genetic covariates # MAYBE remove (only needed when creating the object)
    # summary statistics
    N="integer",		# number of subjects
    nvarlevel="integer",		# number of variability levels
    Nunit="integer", # vector with the number of units for each variability level Nunits[1]=N, Nunits[2]=N*sum_i(nocc_i),...
    nrep.unit="list", # list with for each variability level, the number of lines in the dataset corresponding to this variability level
    # nrep.unit[[1]] = nind.obs
    id.unit="matrix", # matrix with nvarlevel columns
    ntot.obs="numeric",		# total number of observations (=dim(tab)[1])
    nind.obs="numeric"		# total number of observations for each subject [maybe replace with nrep.unit[[1]]]
    # TODO: change all nind.obs to nrep.unit[[1]] in the algorithm
  ),
  validity=function(object){
    #    cat ("--- Checking SaemixData object ---\n")
    if (length(object@name.data)==0) {
      if(object@messages) message("[ SaemixData : validation ] Please provide a name for the data (dataset or datafile on disk).")
      return("No dataset provided")
    }
    return(TRUE)
  }
)


setMethod(
  f="initialize",
  signature="SaemixData",
  definition= function (.Object, name.data, outcome, name.group="id", name.predictors="time", name.response="y",
                        name.covariates=c(), name.mdv="mdv", name.cens="cens",name.occ="occ",name.ytype="ytype",
                        units=c(), verbose=TRUE){
    #    cat ("--- initialising SaemixData Object --- \n")
    if(missing(name.data)) stop("Please provide a data.frame containing the data.")
    if(missing(outcome)) outcome<-list(saemixDataOutcome(name=name.response)) else {
      # check
      if(!is(outcome,"list")) outcome<-list(outcome)
      for(i in 1:length(outcome)) 
        if(!is(outcome,"SaemixOutcome")) stop("If given, outcome must be a list of SaemixDataOutcome objects.")
    }
    .Object@name.data<-name.data
    .Object@outcome<-outcome
    .Object@name.group<-name.group
    .Object@name.predictors<-name.predictors
    .Object@name.response<-name.response
    .Object@name.mdv<-name.mdv
    .Object@name.cens<-name.cens
    .Object@name.ytype<-name.ytype
    .Object@name.occ<-name.occ
    if(length(units)==0) units<-rep("-",length(name.predictors))
    if(is(units,"list")) units<-units$x
    .Object@units <- units
    .Object@messages<-verbose # remove ?
    validObject(.Object)
    return (.Object )
  }
)


###########################################################
# Classes SaemixSimData and SaemixRepData
###########################################################

#' @rdname SaemixData-class
#' @exportClass SaemixSimData

setClass(
  Class="SaemixRepData", # Saemix data, replicated for different chains
  representation=representation(
    N="integer",		# number of subjects
    NM="integer",		# number of subjects, replicated
    dataM="data.frame",		# replicated data with columns IdM, xM, yM
    #    IdM="numeric",		# subject id
    #    XM="data.frame",		# matrix of predictors
    #    yM="numeric",		# vector of responses 
    nrep="integer"		# number of replicates
  ),
  validity=function(object){
    #    cat ("--- Checking SaemixData object ---\n")
    return(TRUE)
  }
)

#' @rdname SaemixData-class
#' @exportClass SaemixRepData

setClass(
  Class="SaemixSimData", # Saemix predicted and simulated data
  representation=representation(
    N="integer",		# number of subjects
    name.group="character", # name of column with ID element
    name.response="character",	# name of column with response
    name.predictors="character",# name of column(s) with predictors 
    name.X="character",		# name of predictor used on X axis for graphs
    units="character",		# units for predictors (a vector)
#    units="list",		# units (list with components for x, y, and cov)
    data="data.frame",		# ECO TODO: do we need to keep it here ?
    nsim="numeric",		# number of simulations
    datasim="data.frame",	# simulated data with columns idsim (id in replications), irep (replication number), ypred (simulated predictions, without error), ysim (simulated data, with error)
    sim.psi="data.frame"	# simulated parameters
  ),
  validity=function(object){
    #    cat ("--- Checking saemixSimData object ---\n")
    return(TRUE)
  }
)

# Initialize method for saemixRepData and saemixSimData
#' @rdname initialize-methods
#' 
#' @param nb.chains number of chains used in the algorithm
#' 
#' @exportMethod initialize

setMethod(
  f="initialize",
  signature="SaemixRepData",
  definition= function (.Object,data=NULL,nb.chains=1){
    #    cat ("--- initialising SaemixData Object --- \n")
    if(is.null(data)) {
      .Object@N<-.Object@NM<-integer(0)
      .Object@dataM<-data.frame()
    } else {
      N<-data@N
      .Object@N<-N
      .Object@nrep<-as.integer(nb.chains)
      .Object@NM<-as.integer(N*nb.chains)
      IdM<-kronecker(c(0:(nb.chains-1)),rep(N,data@ntot.obs))+rep(data@data[,"index"], nb.chains)
      yM<-rep(data@data[,data@name.response],nb.chains)
      XM<-do.call(rbind,rep(list(data@data[,c(data@name.predictors,data@name.cens,data@name.mdv,data@name.ytype),drop=FALSE]), nb.chains))
      .Object@dataM<-data.frame(IdM=c(IdM),XM,yM=yM)
    }
    # Object validation
    #    validObject(.Object)
    return (.Object )
  }
)

#' @rdname initialize-methods
#' 
#' @param datasim dataframe containing the simulated data
#' 
#' @exportMethod initialize

setMethod(
  f="initialize",
  signature="SaemixSimData",
  definition= function (.Object,data=NULL,datasim=NULL) {
    #    cat ("--- initialising SaemixData Object --- \n")
    if(!is.null(data)) {
      .Object@N<-data@N
      .Object@name.group<-data@name.group
      .Object@name.response<-data@name.response
      .Object@name.predictors<-data@name.predictors
      .Object@name.X<-data@name.X
      .Object@units<-data@units
      .Object@data<-data@data
    }
    if(is.null(data) || is.null(datasim) || dim(datasim)[1]==0) {
      .Object@datasim<-data.frame()
      .Object@nsim<-0
    } else {
      .Object@datasim<-datasim
      .Object@nsim<-dim(datasim)[1]/dim(data@data)[1]
    }
    # Object validation
    #    validObject(.Object)
    return (.Object )
  }
)

###########################################################
# Set/Get methods for SaemixData class and derived classes
###########################################################

#' Get/set methods for SaemixData object
#' 
#' Access slots of a SaemixData object using the object\["slot"\] format
#' 
#' @name [
#' @aliases [<-,SaemixData-method [,SaemixData-method
#' 
#' @param x object
#' @param i element to be replaced
#' @param j element to replace with
#' @param value value to replace with
#' @param drop whether to drop unused dimensions
#' @docType methods
#' @keywords methods
#' @exportMethod [
#' @exportMethod [<-
#' @exportPattern "^[[:alpha:]]+"
#' @rdname extract-methods


# Getteur
setMethod(
  f ="[",
  signature = "SaemixData" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "name.data"={return(x@name.data)},
            "messages"={return(x@messages)},
            "outcome"={return(x@outcome)},
            "covariates"={return(x@covariates)},
            "name.group"={return(x@name.group)},
            "name.predictors"={return(x@name.predictors)},
            "name.response"={return(x@name.response)},
            "name.covariates"={return(x@name.covariates)},
            "name.X"={return(x@name.X)},
            "name.mdv"={return(x@name.mdv)},
            "name.cens"={return(x@name.cens)},
            "name.occ"={return(x@name.occ)},
            "name.ytype"={return(x@name.ytype)},
            "units"={return(x@units)},
            "data"={return(x@data)},
            "ocov"={return(x@ocov)},
            "var.match"={return(x@var.match)},
            "N"={return(x@N)},
            "yorig"={return(x@yorig)},
            "ind.gen"={return(x@ind.gen)},
            "ntot.obs"={return(x@ntot.obs)},
            "nind.obs"={return(x@nind.obs)},
            "nvarlevel"={return(x@nvarlevel)},
            "Nunit"={return(x@Nunit)},
            "nrep.unit"={return(x@nrep.unit)},
            "id.unit"={return(x@id.unit)},
            stop("No such attribute\n")
    )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixData" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "name.data"={x@name.data<-value},
            "messages"={x@messages<-value},
            "outcome"={x@outcome<-value},
            "covariates"={x@covariates<-value},
            "name.group"={x@name.group<-value},
            "name.predictors"={x@name.predictors<-value},
            "name.response"={x@name.response<-value},
            "name.covariates"={x@name.covariates<-value},
            "name.X"={x@name.X<-value},
            "name.mdv"={x@name.mdv<-value},
            "name.cens"={x@name.cens<-value},
            "name.occ"={x@name.occ<-value},
            "name.ytype"={x@name.ytype<-value},
            "trans.cov"={x@trans.cov<-value},
            "units"={x@units<-value},
            "data"={x@data<-value},
            "ocov"={x@ocov<-value},
            "var.match"={x@var.match<-value},
            "N"={x@N<-value},
            "ind.gen"={x@ind.gen<-value},
            "yorig"={x@yorig<-value},
            "ntot.obs"={x@ntot.obs<-value},
            "nind.obs"={x@nind.obs<-value},
            "nvarlevel"={x@nvarlevel<-value},
            "Nunit"={x@Nunit<-value},
            "nrep.unit"={x@nrep.unit<-value},
            "id.unit"={x@id.unit<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

# For saemixRepData

#' extract parts of SaemixRepData
#'
#' @name [
#' @aliases [,SaemixRepData-method
#' @docType methods
#' @rdname extract-methods
#' 
setMethod(
  f ="[",
  signature = "SaemixRepData" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "N"={return(x@N)},
            "NM"={return(x@NM)},
            "dataM"={return(x@dataM)},
            "nrep"={return(x@nrep)},
            stop("No such attribute\n")
    )
  }
)

#' replace names of SaemixRepData
#'
#' @name [
#' @aliases [<-,SaemixRepData-method
#' @docType methods
#' @rdname extract-methods

setReplaceMethod(
  f ="[",
  signature = "SaemixRepData" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "N"={x@N<-value},
            "NM"={x@NM<-value},
            "dataM"={x@dataM<-value},
            "nrep"={x@nrep<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

# For saemixSimData

#' extract parts of SaemixSimData
#'
#' @name [
#' @aliases [,SaemixSimData-method
#' @docType methods
#' @rdname extract-methods
#' 

setMethod(
  f ="[",
  signature = "SaemixSimData" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "N"={return(x@N)},
            "name.group"={return(x@name.group)},
            "name.response"={return(x@name.response)},
            "name.predictors"={return(x@name.predictors)},
            "name.X"={return(x@name.X)},
            "units"={return(x@units)},
            "data"={return(x@data)},
            "nsim"={return(x@nsim)},
            "sim.psi"={return(x@sim.psi)},
            "datasim"={return(x@datasim)},
            stop("No such attribute\n")
    )
  }
)

#' replace names of SaemixSimData
#'
#' @name [
#' @aliases [<-,SaemixSimData-method
#' @docType methods
#' @rdname extract-methods

setReplaceMethod(
  f ="[",
  signature = "SaemixSimData" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "N"={x@N<-value},
            "name.group"={x@name.group<-value},
            "name.response"={x@name.response<-value},
            "name.predictors"={x@name.predictors<-value},
            "name.X"={x@name.X<-value},
            "units"={x@units<-value},
            "data"={x@data<-value},
            "sim.psi"={x@sim.psi<-value},
            "datasim"={x@datasim<-value},
            "nsim"={x@nsim<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

####################################################################################
####			SaemixData class - method to print/show data		####
####################################################################################


#' @rdname print-methods
#' 
#' @param x an object of type SaemixData, SaemixModel, SaemixRes or SaemixObject
#' @param nlines maximum number of lines of data to print (defaults to 10)
#' @param ... additional arguments passed on the print function
#' 
#' @exportMethod print

setMethod("print","SaemixData",
          function(x,nlines=10,...) {
            digits<-2;nsmall<-2
            cat("Object of class SaemixData\n")
            cat("    longitudinal data for use with the SAEM algorithm\n")
            cat("Dataset",x@name.data,"\n")
            st1<-paste(x@name.response," ~ ",paste(x@name.predictors,collapse=" + ")," | ", x@name.group,sep="")
            cat("    Structured data:",st1,"\n")
            cat(paste0("Outcome",ifelse(length(x@outcome)>1,"s:",":")), paste(names(x@outcome), collapse=", "),"\n")
            if(length(x@name.predictors)>1) {
              cat("    X variable for graphs:",x@name.X,paste("(",x@units[1],")",sep=""),"\n")
            } else  cat("    Predictor:",x@name.X,paste("(",x@units[1],")",sep=""),"\n")
            # Eco change to adjust to SaemixCovariate object (maybe call specific method)
            # ncov<-length(x@name.covariates)
            # if(ncov>0) {
            #   cat("    covariates:",paste(paste(x@name.covariates," (",x@units$covariates,")",sep=""),collapse=", "),"\n")
            #   if(length(x@ocov)>0) {
            #     for(icov in 1:ncov) {
            #       if(is.factor(x@ocov[,icov]) | length(unique(x@ocov[,icov]))==2) cat("      reference class for covariate",x@name.covariates[icov],": ",levels(as.factor(x@ocov[,icov]))[1],"\n")
            #     }
            #   }
            # }
            if(length(x@data)>0) {
              if(nlines==0) return()
              cat("Dataset characteristics:\n")
              cat("    number of subjects:    ",x@N,"\n")
              if(x@N>0) {
                cat("    number of observations:",x@ntot.obs,"\n")
                cat("    average/min/max nb obs:",format(mean(x@nind.obs),digits=digits, nsmall=nsmall), " / ", min(x@nind.obs)," / ",max(x@nind.obs),"\n")
                #    if(length(x@data)>0) print(x@data)
              }
              xdat<-x@data
              if(length(x@ocov)>0) xdat[,x@name.covariates]<-x@ocov
              if(nlines==(-1)) {
                cat("Data:\n")
                print(xdat)
              } else {
                cat("First",nlines,"lines of data:\n")
                nrowShow <- min (nlines , nrow(xdat ))
                print(xdat[1:nrowShow,-c(1)])
              }
            } else cat("No data.\n")
          }
)

#' @rdname show-methods
#' 
#' @param object an object of type SaemixData, SaemixModel, SaemixRes or SaemixObject
#' 
#' @exportMethod show

setMethod("show","SaemixData",
          function(object) {
            cat("Object of class SaemixData\n")
            cat("    longitudinal data for use with the SAEM algorithm\n")
            cat("Dataset",object@name.data,"\n")
            if(length(object@outcome)==0) {
              cat("    empty object\n")
              return(invisible())
            }
            st1<-paste(object@name.response," ~ ",paste(object@name.predictors,collapse=" + ")," | ", object@name.group,sep="")
            cat("    Structured data:",st1,"\n")
            cat(paste0("Outcome",ifelse(length(object@outcome)>1,"s:",":")), paste(names(object@outcome), collapse=", "),"\n")
            if(length(object@name.predictors)>1) {
              cat("    X variable for graphs:",object@name.X, paste("(",object@units[1],")",sep=""),"\n")
            }
            ncov<-length(object@name.covariates)
            # Eco change to adjust to SaemixCovariate object (maybe call specific method)
            # if(ncov>0) {
            #   cat("    covariates:",paste(paste(object@name.covariates," (",object@units$covariates,")",sep=""),collapse=", "),"\n")
            #   if(length(object@ocov)>0) {
            #     for(icov in 1:ncov) {
            #       if(is.factor(object@ocov[,icov])) cat("      reference class for covariate",object@name.covariates[icov],": ",levels(object@ocov[,icov])[1],"\n")
            #     }
            #     if(length(object@data)>0) object@data[,object@name.covariates]<-object@ocov
            #   }
            # }
            if(length(object@data)>0) {
              if(object@N>0) cat(object@ntot.obs,"    observations in",object@N,"subjects\n")
              cat("First lines of data:\n")
              nrowShow <- min (10 , nrow(object@data ))
              print(object@data[1:nrowShow,-c(1)])
            } else cat("No data.\n")
          }
)

#' @rdname showall-methods
#' @aliases showall
#' @exportMethod showall

# Could be print, with only head of data
setMethod("showall",signature="SaemixData",
          function(object) {
            digits<-2;nsmall<-2
            cat("Object of class SaemixData\n")
            cat("    longitudinal data for use with the SAEM algorithm\n")
            cat("Dataset",object@name.data,"\n")
            cat("    header:",object@header,"\n")
            cat("    sep:",object@sep,"\n")
            cat("    na:",object@na,"\n")
            st1<-paste(object@name.response," ~ ",paste(object@name.predictors,collapse=" + ")," | ", object@name.group,sep="")
            cat("    Structured data:",st1,"\n")
            cat("    subject identifier:    ",object@name.group,"\n")
            cat("    predictors:       ",object@name.predictors,"\n")
            # Eco: adjust to SaemixOutcome
#            cat("    response:         ",object@name.response,paste("(",object@units$y,")",sep=""),"\n")
            cat("    X variable for graphs:",object@name.X,paste("(",object@units[1],")",sep=""),"\n")
            # Eco change to adjust to SaemixCovariate object (maybe call specific method)
            # ncov<-length(object@name.covariates)
            # if(ncov>0) {
            #   cat("    covariates:",paste(paste(object@name.covariates," (",object@units$covariates,")",sep=""),collapse=", "),"\n")
            #   if(length(object@ocov)>0) {
            #     for(icov in 1:ncov) {
            #       if(is.factor(object@ocov[,icov])) cat("      reference class for covariate",object@name.covariates[icov],": ",levels(object@ocov[,icov])[1],"\n")
            #     }
            #     if(length(object@data)>0) object@data[,object@name.covariates]<-object@ocov
            #   }
            # }
            cat("Dataset characteristics:\n")
            cat("    number of subjects:    ",object@N,"\n")
            if(object@N>0) {
              cat("    number of observations:",object@ntot.obs,"\n")
              cat("    average/min/max nb obs:",format(mean(object@nind.obs),digits=digits, nsmall=nsmall), " / ", min(object@nind.obs)," / ",max(object@nind.obs),"\n")
              #    if(length(object@data)>0) print(object@data)
            }
            if(length(object@data)>0) {
              cat("First lines of data:\n")
              nrowShow <- min (10 , nrow(object@data ))
              ncolShow <- min (10 , ncol(object@data))
              print(object@data[1:nrowShow,-c(1)])
            } else cat("No data.\n")
          }
)

# SaemixRepData
#' @rdname show-methods
#' @exportMethod show

setMethod("show","SaemixRepData",
          function(object) {
            cat("Object of class saemixRepData\n")
            if(length(object@N)>0) {
              cat("    replicated data used in the SAEM algorithm\n")
              cat("    number of subjects in initial dataset",object@N,"\n")
              cat("    number of replications",object@nrep,"\n")
              cat("    number of subjects in replicated dataset",object@NM,"\n")
            } else cat("Empty object \n")
          } 
)

# SaemixSimData
#' @rdname show-methods
#' @exportMethod show

setMethod("show","SaemixSimData",
          function(object) {
            cat("Object of class SaemixSimData\n")
            cat("    data simulated according to a non-linear mixed effect model\n")
            if(length(object@N)>0) {
              cat("Characteristics of original data\n")
              cat("    number of subjects:",object@N,"\n")
              cat("    summary of response:\n")
              print(summary(object@data[,object@name.response]))
              cat("Characteristics of simulated data\n")
              if(dim(object@datasim)[1]>0) {
                cat("    number of simulated datasets:",object@nsim,"\n")
                cat("    summary of simulated response\n")
                print(summary(object@datasim$ysim))
              } else cat("    no simulations performed yet\n")
            }}
)

