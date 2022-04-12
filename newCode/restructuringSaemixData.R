saemixDir<-"/home/eco/work/saemix/saemixextension"
source(file.path(saemixDir,"Rext","SaemixData.R"))

theo<-read.table(file.path(saemixDir,"data","theo.saemix.tab"), header=TRUE)
xdat<-new(Class="SaemixData", name.data="theo")

xdat<-new(Class="SaemixData", name.data="theo", name.predictors=c(3,2), name.response=4, name.covariates=c(5,6))
xdat<-new(Class="SaemixData", name.data="theo", name.predictors=c(3,2), name.response=4, covariates=c(wt=saemixCov(name="Weight",unit="kg"), sex=saemixCov(name="Sex",type="binary")))

mytest<-function(objet) {
  x1<-try(objet)
  if(is(x1,"try-error")) return(NULL) else return(x1)
}

saemixData<-function(name.data, name.group="", name.predictors=c(), name.response="",  name.X="", outcome=c(), covariates=c(), name.covariates=c(), name.genetic.covariates=c(), name.mdv="", name.cens="", name.occ="", name.ytype="", units=c(), verbose=TRUE, header=TRUE, sep="\t", na=c("NA","."), automatic=TRUE) {
  # name.data can be a data.frame or the path to a file on disk
  if(missing(name.data)) {
    if(verbose) cat("Error in saemixData: please provide the name of the datafile or dataframe (between quotes)\n")
    return("Creation of SaemixData object failed")
  }
  x1<-try(name.data)
  if(is(x1,"try-error")) return(invisible())
  if(is.data.frame(name.data)) {
    data_from_name.data <- TRUE
    dat <- name.data
    name.data<-deparse(substitute(name.data))
  } else {
    data_from_name.data <- FALSE
    dat <- try(read.table(object@name.data,header=header, sep=sep, na.strings=na))
    if(inherits(dat,"try-error")) stop("The file ",object@name.data," does not exist. Please check the name and path.\n")      
    if(object@messages) {
      cat("These are the first lines of the dataset as read into R. Please check the format of the data is appropriate, if not, modify the na and/or sep items and retry:\n")
      print(head(dat))
    }
  }
  if(is.null(colnames(dat))) colnames(dat)<-paste0("x",1:dim(dat)[2])
  # Creating an SaemixData object
  xdat<-new(Class="SaemixData", name.data=name.data)
  xdat@data<-dat
  if(name.group=="") {
    if(!automatic) {
      if(verbose) cat("   Missing ID identifier, assuming the ID is in column 1 of the dataset.\n")
      name.group<-colnames(dat)[1]
    } else {
      vnames<-validate.names(name.group,colnames(dat),recognisednames=c("id","subject","sujet","group","groupe"),verbose = verbose, automatic=automatic)
      if(length(vnames)==0) {
        if(verbose) message("Please provide a valid name for the ID column.\n")
        return("Creation of SaemixData object failed")
      }
      xdat@name.group<-vnames
    }
  }
  # # ECO TODO: reconnaissance automatique (avant affectation a la valeur 2) ?
  # if(missing(name.predictors)) {
  #   if(automatic) {
  #     name.predictors<-"2"      
  #     if(verbose) cat("   Missing predictors identifier, assuming there is one predictor in column 2 of the dataset.\n")
  #   } else name.predictors<-character()
  # }
  # if(missing(name.response)) {
  #   if(automatic) {
  #     if(verbose) cat("   Missing response identifier, assuming the response is in column 3 of the dataset.\n")
  #     name.response<-"3"
  #   } else name.response<-character()
  # }
  # if(missing(name.covariates)) name.covariates<-character()
  # if(missing(name.mdv)) name.mdv<-character()
  # if(missing(name.cens)) name.cens<-character()
  # if(missing(name.occ)) name.occ<-character()
  # if(missing(name.ytype)) name.ytype<-character()
  # if(missing(name.X)) name.X<-character()
  # .Object@name.group<-name.group
  # .Object@name.predictors<-name.predictors
  # .Object@name.response<-name.response
  # .Object@name.covariates<-name.covariates
  # .Object@name.mdv<-name.mdv
  # .Object@name.cens<-name.cens
  # .Object@name.occ<-name.occ
  # .Object@name.ytype<-name.ytype
  # .Object@outcome<-character()
  # if(missing(units)) units<-rep("-",length(name.predictors))
  # .Object@units<-units
  # .Object@name.X<-name.X
  
  return(xdat)
}

# testing it works
x1<-saemixData(theo)

# testing it works inside a loop
mytest2<-function(dat) {
  for(i in 1:2) {
    name.data<-dat
    x1<-saemixData(name.data)
  }
  return(x1)
}
x1<-mytest2(theo)

#################

setMethod(
  f="initialize",
  signature="SaemixData",
  definition= function (.Object, name.data, outcome, name.group=c("1"), name.predictors=c("2"), name.response=c("3"), name.covariates=c(), name.X, units, name.mdv, name.cens, name.occ, name.ytype, verbose=TRUE){
    #    cat ("--- initialising SaemixData Object --- \n")
    if(missing(name.data)) stop("Please provide a name for the data (dataset or datafile on disk).")
    .Object@name.data<-name.data
    .Object@messages<-verbose
    return (.Object )
  }
)
    xdat<-new(Class="SaemixData", name.data="theo")
    
    

    # Object validation
    #    validObject(.Object)
    return (.Object )
  }
)


## Create a longitudinal data structure from a file or a dataframe
## Helper function not intended to be called by the user
## @exportMethod read

#' Create a longitudinal data structure from a file or a dataframe
#'  
#' Helper function not intended to be called by the user
#'  
#' @param object an SaemixData object
#' @param dat the name of a dataframe in the R environment, defaults to NULL; if NULL, the function will
#' attempt to read the file defined by the slot name.data.
#'  
#' @rdname read-methods
#' @aliases readSaemix,SaemixData readSaemix,SaemixData-method
#'  
#' @exportMethod readSaemix

setMethod("readSaemix",
          signature="SaemixData",
          function(object, dat = NULL) {
            ow <- options("warn")
            options("warn"=-1)
            # ce test devrait aller dans la definition de la classe
            if(class(object@name.data)!="character") {
              if(object@messages) message("Please provide the name of the data (data.frame or path to file on disk) as a character string.\n")
              return("Creation of SaemixData object failed")
            }
            if(is.null(dat)) {
              if(object@messages) cat("Reading data from file",object@name.data,"\n")
              header<-object@header
              if(is.null(header)) header<-TRUE
              sep<-object@sep
              if(is.null(sep)) sep<-""
              na.strings<-object@na
              if(is.null(na.strings)) na.strings<-"NA"
              dat<-try(read.table(object@name.data,header=header,sep=sep,na.strings=na.strings))
              if(inherits(dat,"try-error")) stop("The file ",object@name.data," does not exist. Please check the name and path.\n")      
              if(object@messages) {
                cat("These are the first lines of the dataset as read into R. Please check the format of the data is appropriate, if not, modify the na and/or sep items and retry:\n")
                print(head(dat))
              }
            }
            if(dim(dat)[2]<3) {
              if(object@messages) message("The dataset does not contain enough data. The non-linear mixed effect model requires at least 3 columns, with subject ID, predictor (at least one) and response. \nPlease check the field separator, currently given as:", paste("sep=\"",object@sep,"\"",sep=""),"\n")
              return("Creation of SaemixData object failed")
            }
            # Automatic recognition of columns 
            #    ID (one of id, subject, sujet, group, groupe regardless of case)
            #    response (one of Y, conc, concentration, resp, response, y, dv regardless of case)
            #    predictors (time and/or dose, x, regardless of case)
            # ECO TODO: improve automatic recognition ?
            # check that we have at least a column id, response and X
            
            vnames<-validate.names(object@name.group,colnames(dat),recognisednames=c("id","subject","sujet","group","groupe"),verbose = object@messages, automatic=object@automatic)
            if(length(vnames)==0) {
              if(object@messages) message("Please provide a valid name for the ID column.\n")
              return("Creation of SaemixData object failed")
            }
            object@name.group<-vnames
            
            vnames<-validate.names(object@name.predictors,colnames(dat),recognisednames=c("time","temps","tps","tim","x","dose"),verbose = object@messages, automatic=object@automatic)
            if(length(vnames)==0) {
              if(object@messages) message("Please provide a valid name for the predictor(s).\n")
              return("Creation of SaemixData object failed")
            }
            object@name.predictors<-vnames
            
            vnames<-validate.names(object@name.response,colnames(dat),recognisednames=c("response","resp","conc","concentration","y","dv"),verbose = object@messages, automatic=object@automatic)
            if(length(vnames)==0) {
              if(object@messages) message("Please provide a valid name for the response.\n")
              return("Creation of SaemixData object failed")
            }
            object@name.response<-vnames[1]
            if(length(vnames)>1 & object@messages) cat("Using the response",object@name.response,"as dependent variable.\n")
            
            if(length(object@name.covariates)>0) {
              if(object@name.covariates[1]!="") {
                i1<-suppressWarnings(as.integer(object@name.covariates[!is.na(suppressWarnings(as.integer(object@name.covariates)))]))
                object@name.covariates[!is.na(suppressWarnings(as.integer(object@name.covariates)))]<- colnames(dat)[i1]
              }
              idx<-object@name.covariates[!(object@name.covariates %in% colnames(dat))]
              if(length(idx)>0) {
                if(object@messages) cat("Covariates",object@name.covariates[idx],"not found.\n") 
                object@units$covariates<-object@units$covariates[object@name.covariates %in% colnames(dat)]
                object@name.covariates<-object@name.covariates[object@name.covariates %in% colnames(dat)]
              }
            }
            if(nchar(object@name.group)*length(object@name.predictors)*nchar(object@name.response)<=0) {
              stop("Please check the structure of the data file and provide information concerning which columns specify the group structure (ID), the predictors (eg dose, time) and the response (eg Y, conc). See documentation for automatic recognition of column names for these elements.\n")
            }
            if(nchar(object@name.X)==0)
              object@name.X<-object@name.predictors[1]
            if(!is.na(suppressWarnings(as.integer(object@name.X)))) {
              if(dim(dat)[2]<suppressWarnings(as.integer(object@name.X))) {
                if(object@messages) cat("Attribute name.X",object@name.X,"does not correspond to a valid column in the dataset, setting the X axis for graphs to",object@name.predictors[1],".\n")
                object@name.X<-object@name.predictors[1]
              } else object@name.X<-colnames(dat)[suppressWarnings(as.integer(object@name.X))]
            } 
            if(match(object@name.X,object@name.predictors,nomatch=0)==0) {
              if(object@messages) cat("Attribute name.X",object@name.X,"does not correspond to a valid column in the dataset, setting the X axis for graphs to",object@name.predictors[1],".\n")
              object@name.X<-object@name.predictors[1]
            }
            if(nchar(object@name.mdv)==0) mdv<-rep(0,dim(dat)[1]) else {mdv<-dat[,object@name.mdv]}
            mdv[is.na(dat[,object@name.response])]<-1
            #    if(sum(mdv)>0) 
            object@name.mdv<-"mdv"
            if(nchar(object@name.cens)==0) cens<-rep(0,dim(dat)[1]) else cens<-dat[,object@name.cens]
            object@name.cens <-"cens"
            if(nchar(object@name.occ)==0) occ<-rep(1,dim(dat)[1]) else occ<-dat[,object@name.occ]
            object@name.occ<-"occ"
            if(nchar(object@name.ytype)==0) ytype<-rep(1,dim(dat)[1]) else ytype<-dat[,object@name.ytype]
            object@name.ytype<-"ytype"
            all.names<-c(object@name.group,object@name.predictors, object@name.response,object@name.covariates)
            
            dat<-dat[,all.names,drop=FALSE]
            dat<-cbind(dat,mdv=mdv,cens=cens,occ=occ,ytype=ytype)
            
            if(!is(dat,"data.frame")) dat<-as.data.frame(dat)
            # Saving covariates in the original format in ocov, transforming binary covariates in dat to factors
            object@ocov<-dat[,object@name.covariates,drop=FALSE]
            for(icov in object@name.covariates) {
              if(length(unique(dat[,icov]))==2) dat[,icov]<-suppressWarnings(as.integer(factor(dat[,icov])))-1
            }
            # Removing missing values in predictor columns
            # dat<-dat[!is.na(dat[,object@name.response]),]
            for(i in object@name.predictors) {
              if(sum(is.na(dat[,i]))>0) {
                if(object@messages) cat("Removing missing values for predictor",i,"\n")
                if(!is.null(dim(object@ocov)[1])) object@ocov<-object@ocov[!is.na(dat[,i]),,drop=FALSE]
                dat<-dat[!is.na(dat[,i]),]
              }
            }
            # Removing subjects with only MDV in responses
            idx<-c();inull<-c()
            for(isuj in unique(dat[,object@name.group])) {
              if(sum(1-dat$mdv[dat[,object@name.group]==isuj])==0) {
                inull<-c(inull,isuj)
                idx<-c(idx,which(dat[,object@name.group]==isuj))
              }
            }
            #  if(object@messages) print(idx)
            if(length(inull)>0) {
              if(object@messages) cat("Some subjects have no observations, removing them:",inull,"\n")
              dat<-dat[-idx,]
              if(!is.null(dim(object@ocov)[1])) object@ocov<-object@ocov[-idx,,drop=FALSE]
            }
            if(!is.null(dim(object@ocov)[1])) object@ocov <- object@ocov[order(dat[,object@name.group], dat[,object@name.X]),,drop=FALSE]
            
            # ECO TODO: missing data in covariates kept for the moment, only excluded depending on the model
            #    for(i in object@name.covariates) dat<-dat[!is.na(dat[,i]),]
            object@ntot.obs<-dim(dat)[1] # total number of observations
            dat <- dat[order(dat[,object@name.group], dat[,object@name.X]),]
            id<-dat[,object@name.group]
            object@N<-length(unique(id))
            nind.obs<-tapply(id,id,length) # individual numbers of observations (1xN)
            nind.obs<-nind.obs[match(unique(id),names(nind.obs))]
            object@nind.obs<-c(nind.obs)
            dat<-cbind(index=rep(1:object@N,times=nind.obs),dat)
            object@data<-dat
            
            options(ow) # reset
            validObject(object)
            return(object)
          }
)
