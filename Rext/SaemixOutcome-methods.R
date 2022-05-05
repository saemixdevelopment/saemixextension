################################################################
# functions to create continuous and discrete outcomes

#' Define SaemixOutcome objects for use in data or model objects
#' 
#' These functions provide easy constructors for single or multiple outcomes used in 
#' the definition of data objects (parent class, with only name, type and unit) or 
#' model objects (children classes, with more detailed information such as the residual 
#' error model for continuous outcomes).
#' 
#' @name saemixOutcome
#' @aliases saemixDataOutcome 
#' @aliases discreteOutcome continuousOutcome
#' @examples
#' # Defining a typed outcome for SaemixData objects
#' conc <- saemixDataOutcome("concentration", unit="mg/L")
#' conc
#' pain <- saemixDataOutcome("pain score",type="categorical", unit="(-)")
#' pain
#' resp <- saemixDataOutcome("response",type="binary")
#' resp
#' 
#' # Defining a typed outcome for SaemixModel objects
#' conc <- saemixOutcome("concentration", unit="mg/L")
#' conc
#' pain <- saemixOutcome("pain score",type="categorical", unit="(-)")
#' pain
#' 

# For use when defining a data object (parent)

saemixDataOutcome<-function(name="", type="continuous", unit="") {
  new(Class = "SaemixOutcome", name=name, type=type, unit=unit)
}

# For use when defining a model object (child)

saemixOutcome<-function(name="",type="continuous", unit="", distribution="", model="constant", start=c(1), error.function=constantErrorModel, npar=1, fix=c(0), ncat=NULL, levels=NULL, maxEvents=NULL, intervalCensored=FALSE, intervalLength=NULL) {
  if(length(fix)>0 & is(fix,"logical")) fix<-as.integer(fix!=0)
  if(type=="continuous") saemix.outcome<-new(Class="SaemixContinuousOutcome", name=name, error.model=model, error.npar=npar, error.function=error.function, error.parameters=start, error.fix=fix, unit=unit)
  if(type %in% c("binary","categorical","count")) saemix.outcome<-new(Class="SaemixDiscreteOutcome", name=name, levels=levels, distribution=distribution, type=type, unit=unit)
  if(type=="event") saemix.outcome<-new(Class="SaemixEventOutcome", name=name, distribution=distribution, maxEvents=maxEvents, intervalCensored=intervalCensored, intervalLength=intervalLength, unit=unit)
  # if(name!="" && type=="continuous") {
  #   saemix.outcome@name<-name
  #   rootpar<-strsplit(saemix.outcome@error.nameparameters,".",fixed=T)
  #   vec<-c()
  #   for(i in 1:length(rootpar)) vec<-rootpar[[i]][1]
  #   saemix.outcome@error.nameparameters <- paste(vec,name,sep=".")
  # }
  return(saemix.outcome)
  
}


continuousOutcome<-function(model="constant", start=c(1), error.function=constantErrorModel, npar=1, fix=c(0), unit="") {
  if(model=="user") {
    ierr<-0
    if(identical(error.function, constantErrorModel)) {
      model<-"constant"
      ierr<-1
    }
    if(identical(error.function, proportionalErrorModel)) {
      model<-"proportional"
      ierr<-1
    }
    if(identical(error.function, combined1ErrorModel)) {
      model<-"combined1"
      ierr<-1
    }
    if(identical(error.function, combined2ErrorModel)) {
      model<-"combined2"
      ierr<-1
    }
    if(identical(error.function, powerErrorModel)) {
      model<-"power"
      ierr<-1
    }
    if(ierr==1) 
      message("User error model set to ",model," for consistency with the error.function given\n")
  }
  outlist<-list(type="continuous", error.model=model, start=start, error.function=error.function, error.npar=npar, error.fix=fix)
  if(outlist$error.model=="constant" | outlist$error.model=="exponential") {
    outlist$error.function<-constantErrorModel
    outlist$error.npar<-1
    outlist$error.fix<-outlist$error.fix[1]
    outlist$start<-outlist$start[1]
    if(length(outlist$start)==0) outlist$start<-c(1)
    if(outlist$start<=0) {
      outlist$start<-1
      message("Starting value for constant error model must be positive, changed to 1\n")
    }
  }
  if(outlist$error.model=="proportional") {
    outlist$error.npar<-1
    outlist$error.function<-proportionalErrorModel
    outlist$start<-outlist$start[1]
    outlist$error.fix<-outlist$error.fix[1]
    if(length(outlist$start)==0) outlist$start<-c(0.5)
    if(outlist$start<=0) {
      outlist$start<-0.5
      message("Starting value for proportional error model must be positive, changed to 0.5\n")
    }
  }
  if(outlist$error.model=="combined1") {
    outlist$error.function<-combined1ErrorModel
    outlist$error.npar<-2
    if(length(outlist$start)>2) outlist$start<-abs(outlist$start[1:2])
    if(length(outlist$start)==1) outlist$start<-c(abs(outlist$start),0.5)
    if(length(outlist$start)==0) outlist$start<-c(1,0.5)
    if(length(outlist$error.fix)<2) outlist$error.fix<-rep(outlist$error.fix, length.out=2) else outlist$error.fix<-outlist$error.fix[1:2]
    x<-try(outlist$error.function(c(1), outlist$start))
    if(is.na(x)) message("Check starting value for error function, returns NA for x=1\n")
  }
  if(outlist$error.model=="combined2") {
    outlist$error.function<-combined2ErrorModel
    outlist$error.npar<-2
    if(length(outlist$start)>2) outlist$start<-abs(outlist$start[1:2])
    if(length(outlist$start)==1) outlist$start<-c(abs(outlist$start),0.5)
    if(length(outlist$start)==0) outlist$start<-c(1,0.5)
    if(length(outlist$error.fix)<2) outlist$error.fix<-rep(outlist$error.fix, length.out=2) else outlist$error.fix<-outlist$error.fix[1:2]
    x<-try(outlist$error.function(c(1), outlist$start))
    if(is.na(x)) message("Check starting value for error function, returns NA for x=1\n")
  }
  if(outlist$error.model=="power") {
    outlist$error.function<-powerErrorModel
    outlist$error.npar<-3
    if(length(outlist$start)>3) outlist$start<-outlist$start[1:3]
    if(length(outlist$start)==0) outlist$start<-c(1,0.5,2)
    if(length(outlist$start)==1) outlist$start<-c(outlist$start,0.5,2)
    if(length(outlist$start)==2) outlist$start<-c(outlist$start,2)
    outlist$start[1:2]<-abs(outlist$start[1:2]) # a and b must be positive
    if(length(outlist$error.fix)<3) outlist$error.fix<-rep(outlist$error.fix, length.out=3) else outlist$error.fix<-outlist$error.fix[1:3]
    x<-try(outlist$error.function(c(1), outlist$start))
    if(is.na(x)) message("Check starting value for error function, returns NA for x=1\n")
  }
  if(outlist$error.model=="user") {
    if(length(start)!=npar) outlist$error.npar<-length(start)
    x<-try(error.function(c(1), start))
    if(length(outlist$error.fix)<outlist$error.npar) outlist$error.fix<-rep(outlist$error.fix, length.out=outlist$error.npar) else outlist$error.fix<-outlist$error.fix[1:outlist$error.npar]
    if(is.na(x)) {
      message("Check user-defined error function, returns NA for x=1\n")
      return()
    }
  }
  outlist$unit<-unit
  return(outlist)
}

discreteOutcome<-function(type="binary", distribution="", ncat=NULL, levels=NULL, maxEvents=NULL, intervalCensored=FALSE, intervalLength=NULL, unit="") {
  if(is.null(levels) & !is.null(ncat)) levels<-c(1:ncat)
  if(length(levels)>2 & type=="binary") type<-distribution<-"categorical"
  outlist<-list(type=type, distribution=distribution, levels=levels, maxEvents=maxEvents, intervalCensored=intervalCensored, intervalLength=intervalLength)
  if(outlist$distribution=="") {
    if(type=="categorical") outlist$distribution<-"categorical"
    if(type=="binary") outlist$distribution<-"bernouilli"
  }
  if(is.null(levels) & type=="binary")
    outlist$levels<-c(0,1)
  if(outlist$type=="event") {
    if(is.null(outlist$maxEvents)) outlist$maxEvents<-1
  }
  if(is.null(outlist$intervalLength)) outlist$intervalLength<-0
  outlist$unit<-unit
  return(outlist)
}


################################################
## Function to transform outcomes to a list of SaemixOutcome objects

#' @aliases createSaemixOutcome convertArg2Outcome
#' @example 
#' # Defining a list of outcomes for SaemixModel objects
#' # first outcome: only name is given, defaults to continuous outcome with constant error model
#' # second outcome: defined as a continuous outcome with options
#' # third outcome: categorical, defined using the name="type" version
#' lout <- convertArg2Outcome(list("concentration",effect=continuousOutcome(model="combined1", start=c(2,0.5)), pain="categorical"))
#' lout
#' @export createSaemixOutcome

############ Create outcomes from a list

createSaemixOutcome<-function(response, name="") {
  if(response$type=="continuous") saemix.outcome<-new(Class="SaemixContinuousOutcome", name=name, error.model=response$error.model, error.npar=outcome$error.npar, error.function=response$error.function, error.parameters=response$start, error.fix=response$error.fix, unit=response$unit)
  if(response$type %in% c("binary","categorical","count")) saemix.outcome<-new(Class="SaemixDiscreteOutcome", name=name, levels=response$levels, distribution=response$distribution, type=response$type, unit=response$unit)
  if(response$type=="event") saemix.outcome<-new(Class="SaemixEventOutcome", name=name, distribution=response$distribution, maxEvents=response$maxEvents, intervalCensored=response$intervalCensored, intervalLength=response$intervalLength, unit=response$unit)
  return(saemix.outcome)
}

############ Convert a possibly mixed list of outcomes into SaemixOutcome of different types
# outcome can be a vector or a list
# the elements of outcome can be:
## SaemixOutcome objects
## a character string giving the type of the outcome (one of "continuous","binary","categorical","event")
### example: "continuous" 
## a character string giving the name of the outcome (not one of "continuous","binary","categorical","event")
### example: "concentration" 
## a named character string outcomeName="type" where "type" is one of of "continuous","binary","categorical","event" 
### example: concentration="continuous" 

#' @export convertArg2Outcome

convertArg2Outcome <- function(outcome, verbose=FALSE) {
  if(is(outcome,"numeric")) { # only a number of responses => assume they are all continuous
    smx.out<-vector(outcome[1], mode="list")
    for(i in 1:outcome[1]) smx.out[[i]]<-createSaemixOutcome(continuousOutcome(), name=paste0("y",i))
    names(smx.out)<-paste0("y",1:outcome[1])
    return(smx.out)
  }
  saemix.outcome<-vector(mode="list", length=length(outcome))
  nameout<-paste0("out",1:length(outcome))
  if(!is.null(names(outcome))) names(saemix.outcome)<-names(outcome) else names(saemix.outcome)<-nameout
  names(saemix.outcome)[names(saemix.outcome)==""]<-nameout[names(saemix.outcome)==""]
  if(!is(outcome,"list")) { # parameters are given as names only, names+type, or type only
    for(i in 1:length(outcome)) {
      if(outcome[i] %in% c("continuous","binary","count","categorical","event")) {
        saemix.outcome[[i]]<-saemixOutcome(type=outcome[i])
        if(is.null(names(saemix.outcome)[i]) || names(saemix.outcome)[i]=="") names(saemix.outcome)[i]<-paste0("out",i)
      } else {
        saemix.outcome[[i]]<-saemixOutcome(type="continuous")
        names(saemix.outcome)[i]<-outcome[i]
      }
    }
  } else {
    for(i in 1:length(outcome)) {
      if(is(outcome[[i]],"SaemixOutcome")) {
        saemix.outcome[[i]]<-outcome[[i]]
      } else {
        if(is(outcome[[i]],"list")) {
          x1<-try(createSaemixOutcome(outcome[[i]], name=names(saemix.outcome)[i]))
          if(is(x1, "try-error")) {
            if(verbose) message("Outcome number",i,"of invalid type, assuming a continuous outcome")
            saemix.outcome[[i]]<-saemixOutcome(name=names(saemix.outcome)[i])
          } else saemix.outcome[[i]]<-x1
        } else {
          if(outcome[i] %in% c("continuous","binary","count","categorical","event")) {
            saemix.outcome[[i]]<-saemixOutcome(type=outcome[i])
            if(is.null(names(saemix.outcome)[i]) || names(saemix.outcome)[i]=="") names(saemix.outcome)[i]<-paste0("out",i)
          } else {
            saemix.outcome[[i]]<-saemixOutcome(type="continuous")
            names(saemix.outcome)[i]<-outcome[i]
          }
        }
      }
    }
  }
  for(i in 1:length(outcome)) {
    if(saemix.outcome[[i]]@name=="") saemix.outcome[[i]]@name<-names(saemix.outcome)[i]
    if(saemix.outcome[[i]]@type=="continuous") {
      idx<-which(saemix.outcome[[i]]@error.nameparameters %in% c("a.","b.","c."))
      if(length(idx)>0) saemix.outcome[[i]]@error.nameparameters[idx]<-paste0(saemix.outcome[[i]]@error.nameparameters[idx], names(saemix.outcome)[i])
    }
  }
  return(saemix.outcome)
}

