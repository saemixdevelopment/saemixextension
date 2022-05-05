########################################################################
# Defines covariate transformation models
## parent class: SaemixCovariateTransform
## child classes: covmodelCont2Cont, covmodelCont2Cat, covmodelCont2Cat

setClass(Class = "SaemixCovariateTransform",
         representation=representation(
           name = "character", # name of transformed covariate
           type = "character", # name of transformed covariate
           name.orig = "character", # name of original covariate
           type.orig = "character" # type of original covariate
         ),
         validity=function(object){
           return(TRUE)
         }
)

setClass(Class = "covmodelCont2Cont",
         contains="SaemixCovariateTransform",
         representation=representation(
           transform.function = "function", # Transformation applied to the covariate (a function)
           centering.function = "function", # if given, the function is applied (eg: median, mean)
           centering.value = "numeric" # if given and centering.function not given, the (fixed) value used to center
         ),
         validity=function(object){
           return(TRUE)
         }
)

setClass(Class = "covmodelCont2Cat",
         contains="SaemixCovariateTransform",
         representation=representation(
           name.cat = "character", # Category names after transformation
           reference = "character", # reference category (if "", defaults to the first category)
           ncat = "numeric", # if given, the number of categories to dichotomise to (otherwise determined by breaks)
           breaks = "numeric" # if given, the breaks used to create the categories
         ),
         validity=function(object){
           return(TRUE)
         }
)

setClass(Class = "covmodelCat2Cat",
         contains="SaemixCovariateTransform",
         representation=representation(
           name.cat = "character", # Category names after transformation
           reference = "character", # reference category (if "", defaults to the first category)
           ncat = "numeric", # number of categories
           groups = "list" # the groups
         ),
         validity=function(object){
           return(TRUE)
         }
)

########################################################################
# Initialize

setMethod( 
  f="initialize",
  signature="SaemixCovariateTransform",
  definition=function(.Object, name, type="continuous", name.orig=""){
    .Object@name <- name
    .Object@type <- type
    if(name.orig=="") .Object@name.orig <- name else .Object@name.orig <- name.orig
    validObject(.Object)
    return(.Object)
  }
)

setMethod( 
  f="initialize",
  signature="covmodelCont2Cont",
  definition=function(.Object, name, name.orig="", transform.function, centering.function, centering.value, verbose=FALSE){
    .Object <- callNextMethod(.Object, name=name, name.orig=name.orig)
    .Object@type<-"continuous"
    .Object@type.orig<-"continuous"
    if(!missing(transform.function) && is(transform.function,"function")) .Object@transform.function <- transform.function else .Object@transform.function<-function(x) x
    if(!missing(centering.value)) {
      .Object@centering.value<-centering.value
    } else {
      if(!missing(centering.function) && is(centering.function,"function")) .Object@centering.function<-centering.function else .Object@centering.function<-function(x) 1
    }
    validObject(.Object)
    return(.Object)
  }
)


setMethod( 
  f="initialize",
  signature="covmodelCont2Cat",
  definition=function(.Object, name, breaks, name.orig="", ncat=2, reference=character(), name.cat=character(), verbose=FALSE){
    .Object <- callNextMethod(.Object, name=name, name.orig=name.orig)
    .Object@type<-"categorical"
    .Object@type.orig<-"continuous"
    if(!missing(breaks)) {
      ncat<-length(breaks)+1
      .Object@breaks<-breaks
      if(length(name.cat)>0) {
        if(length(name.cat)!=ncat) {
          if(verbose) message("Size mismatch between name.cat and the number of categories, ignoring")
          name.cat<-character()
        }
      }
      if(length(name.cat)==0) {
        x1<-c("min",breaks,"max")
        name.cat<-paste(x1[-length(x1)],x1[-c(1)],sep="-")
      }
    } else { # only nb of categories given
      if(length(name.cat)!=ncat) {
        if(verbose) message("Size mismatch between name.cat and the number of categories, ignoring")
#        name.cat<-paste0("G",1:ncat)
        name.cat<-character()
      }
    }
    .Object@name.cat<-name.cat
    .Object@ncat<-ncat
    if(length(reference)==0 && length(.Object@name.cat)>0)
      reference<-.Object@name.cat[1]
    if(length(.Object@name.cat)>0 && is.na(match(reference, .Object@name.cat))) { # reference name not found in the categories
      if(verbose) message(paste("Category",reference,"not found in", .Object@name.cat,": changing reference to category",.Object@name.cat[1]))
      reference<-.Object@name.cat[1]
    }
    .Object@reference<-reference #min(reference, .Object@ncat)
    validObject(.Object)
    return(.Object)
  }
)

setMethod( 
  f="initialize",
  signature="covmodelCat2Cat",
  definition=function(.Object, name, groups, name.orig="", reference=character(), name.cat=character(), verbose=FALSE){
    .Object <- callNextMethod(.Object, name=name, name.orig=name.orig)
    .Object@type<-"categorical"
    .Object@type.orig<-"categorical"
    if(missing(groups) & length(name.cat)==0) { # Initialise to all empty, we'll use the categorical covariate as is
      if(verbose) message("Initialising to empty number of categories")
      .Object@groups<-list()
      .Object@name.cat<-character()
      ncat<-numeric(0)
    }
    if(!missing(groups)) {
      ncat<-length(groups)
      .Object@groups<-groups
      if(length(name.cat)>0 && length(groups)==length(name.cat))
        .Object@name.cat <- name.cat else {
          name1<-1:length(groups)
          for(i in 1:length(name1)) name1[i]<-paste(unlist(groups[[i]]),collapse="-")
          .Object@name.cat <- name1
        }
    } else {
      if(length(name.cat)>0) {
        ncat<-length(name.cat)
        .Object@name.cat <- name.cat
      } # else .Object@name.cat <- paste0("G",1:ncat)
    }
    .Object@ncat<-ncat
    if(length(reference)==0 && length(.Object@name.cat)>0)
      reference<-.Object@name.cat[1]
    if(length(.Object@name.cat)>0 && is.na(match(reference, .Object@name.cat))) { # reference name not found in the categories
      if(verbose) message(paste("Category",reference,"not found in", .Object@name.cat,": changing reference to category",.Object@name.cat[1]))
      reference<-.Object@name.cat[1]
    }
    .Object@reference<-reference #min(reference, .Object@ncat)
    validObject(.Object)
    return(.Object)
  }
)


########################################################################
# Getteur
setMethod(
  f ="[",
  signature = "SaemixCovariateTransform" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "name"={return(x@name)},
            "type"={return(x@type)},
            "name.orig"={return(x@name.orig)},
            "type.orig"={return(x@type.orig)},
            stop("No such attribute\n")
    )
  }
)
setMethod(
  f ="[",
  signature = "covmodelCont2Cont" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "transform.function"={return(x@transform.function)},
            "centering.function"={return(x@centering.function)},
            "centering.value"={return(x@centering.value)},
            stop("No such attribute\n")
    )
  }
)
setMethod(
  f ="[",
  signature = "covmodelCont2Cat" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "name.cat"={return(x@name.cat)},
            "reference"={return(x@reference)},
            "ncat"={return(x@ncat)},
            "breaks"={return(x@breaks)},
            stop("No such attribute\n")
    )
  }
)
setMethod(
  f ="[",
  signature = "covmodelCat2Cat" ,
  definition = function (x,i,j,drop ){
    switch (EXPR=i,
            "name.cat"={return(x@name.cat)},
            "reference"={return(x@reference)},
            "ncat"={return(x@ncat)},
            "groups"={return(x@groups)},
            stop("No such attribute\n")
    )
  }
)


# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixCovariateTransform" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "name"={x@name<-value},
            "type"={x@type<-value},
            "name.orig"={x@name.orig<-value},
            "type.orig"={x@type.orig<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)
setReplaceMethod(
  f ="[",
  signature = "covmodelCont2Cont" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "transform.function"={x@transform.function<-value},
            "centering.function"={x@centering.function<-value},
            "centering.value"={x@centering.value<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

setReplaceMethod(
  f ="[",
  signature = "covmodelCont2Cat" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "name.cat"={x@name.cat<-value},
            "reference"={x@reference<-value},
            "ncat"={x@ncat<-value},
            "breaks"={x@breaks<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

setReplaceMethod(
  f ="[",
  signature = "covmodelCat2Cat" ,
  definition = function (x,i,j,value){
    switch (EXPR=i,
            "name.cat"={x@name.cat<-value},
            "reference"={x@reference<-value},
            "ncat"={x@ncat<-value},
            "groups"={x@groups<-value},
            stop("No such attribute\n")
    )
    validObject(x)
    return(x)
  }
)

########################################################################
# Print and show methods

setMethod("show","covmodelCont2Cont",
          function(object) {
            cat(object@type,"covariate")
            if(object@name!="") cat(": ",object@name)
            cat("\n")
            if(!identical(object@transform.function, function(x) x)) {
              cat("    transformation applied\n")
            }
            if(!identical(object@centering.function, function(x) 1)) {
              cat("    centering function applied\n")
            } else {
              if(length(object@centering.value)>0) cat("    centering on value:",object@centering.value,"\n")
            }
          }
)
setMethod("showall","covmodelCont2Cont",
          function(object) {
            cat(object@type,"covariate")
            if(object@name!="") cat(": ",object@name)
            cat("\n")
            if(!identical(object@transform.function, function(x) x)) {
              cat("    transformation:\n")
              print(object@transform.function)
            }
            if(!identical(object@centering.function, function(x) 1)) {
              cat("    centering function:\n")
              print(object@centering.function)
            } else {
              if(length(object@centering.value)>0) cat("    centering on value:",object@centering.value,"\n")
            }
          }
)

setMethod("print","covmodelCont2Cont",
          function(x,nlines=10,...) {
            show(x)
          }
)

# Categorical covariate transformed from a continuous covariate
setMethod("show","covmodelCont2Cat",
          function(object) {
            cat(object@type,"covariate")
            if(object@name!="") {
              cat(": ",object@name)
              if(object@name.orig!="") cat(" transformed from ",object@name.orig)
            }
            cat("\n")
            if(length(object@ncat)>0) 
              cat("      in",object@ncat,"categories\n")
            if(length(object@breaks)>0) {
              cat("      breaks used to cut covariate: ",object@breaks,"\n")
            }
          }
)

setMethod("showall","covmodelCont2Cat",
          function(object) {
            show(object)
          }
)

setMethod("print","covmodelCont2Cat",
          function(x,nlines=10,...) {
            show(x)
          }
)

# Categorical covariate transformed from a categorical covariate
setMethod("show","covmodelCat2Cat",
          function(object) {
            cat(object@type,"covariate")
            if(object@name!="") {
              cat(": ",object@name)
              if(object@name.orig!="") cat(" transformed from ",object@name.orig)
            }
            cat("\n")
            if(length(object@ncat)>0) 
              cat("      in",object@ncat,"categories\n")
            if(length(object@groups)>0)
              for(i in 1:length(object@groups)) cat("Group",object@name.cat[i],":",object@groups[[i]],"\n")
          }
)

setMethod("showall","covmodelCat2Cat",
          function(object) {
            show(object)
          }
)

setMethod("print","covmodelCat2Cat",
          function(x,nlines=10,...) {
            show(x)
          }
)

########################################################################

