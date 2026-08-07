#' @include aaa_generics.R
#' @include SaemixData.R
NULL
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
            st1<-paste(x@name.response," ~ ",paste(x@name.predictors,collapse=" + "),sep="")
            for(i in length(x@varlevel):1) st1<-paste(st1," | ",x@varlevel[i])
            cat("    Structured data:",st1,"\n")
            if(length(x@name.predictors)>1) {            
              cat("    X variable for graphs:",x@name.X,paste("(",x@units$x,")",sep=""),"\n")
            } else  cat("    Predictor:",x@name.X,paste("(",x@units$x,")",sep=""),"\n")
            ncov<-length(x@name.covariates)
            if(ncov>0) {
              cat("    covariates:",paste(paste(x@name.covariates," (",x@units$covariates,")",sep=""),collapse=", "),"\n")
              # if(length(x@ocov)>0) {
              #   for(icov in 1:ncov) {
              #     if(is.factor(x@ocov[,icov]) | length(unique(x@ocov[,icov]))==2) cat("      reference class for covariate",x@name.covariates[icov],": ",levels(as.factor(x@ocov[,icov]))[1],"\n")
              #   }
              # }
            }
            if(FALSE) {
              cat("    Group column:",x@name.group,"\n")
              cat("    Predictors:",x@name.predictors,"\n")
              cat("    X variable for graphs:",x@name.X,paste("(",x@units$x,")",sep=""),"\n")
              cat("    Response column:",x@name.response, paste("(",x@units$y,")",sep=""),"\n")
              cat("    Covariates:",x@name.covariates,"\n")
            }
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
              #              if(length(x@ocov)>0) try(xdat[,x@name.covariates]<-x@ocov)
              if(nlines==(-1)) {
                cat("Data:\n")
                print(xdat)
              } else {
                cat("First",nlines,"lines of data:\n")
                nrowShow <- min (nlines , nrow(xdat ))
                print(xdat[1:nrowShow,-c(1)])
              }
            } else cat("No data.\n")
            if(x@messages) cat(x@log)
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
            st1<-paste(object@name.response," ~ ",paste(object@name.predictors,collapse=" + "),sep="")
            for(i in length(object@varlevel):1) st1<-paste(st1," | ",object@varlevel[i])
            cat("    Structured data:",st1,"\n")
            if(length(object@name.predictors)>1) {
              cat("    X variable for graphs:",object@name.X, paste("(",object@units$x,")",sep=""),"\n")
            }
            ncov<-length(object@name.covariates)
            if(ncov>0) {
              cat("    covariates:",paste(paste(object@name.covariates," (",object@units$covariates,")",sep=""),collapse=", "),"\n")
              # if(length(object@ocov)>0) {
              #   for(icov in 1:ncov) {
              #     if(is.factor(object@ocov[,icov])) cat("      reference class for covariate",object@name.covariates[icov],": ",levels(object@ocov[,icov])[1],"\n")
              #   }
              #   if(length(object@data)>0) try(object@data[,object@name.covariates]<-object@ocov)
              # }
            }
            if(length(object@data)>0) {
              if(object@N>0) cat(object@ntot.obs,"    observations in",object@N,"subjects\n")
              cat("First lines of data:\n")
              nrowShow <- min (10 , nrow(object@data ))
              print(object@data[1:nrowShow,-c(1)])
            } else cat("No data.\n")
            if(object@messages) cat(object@log)
            
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
            st1<-paste(object@name.response," ~ ",paste(object@name.predictors,collapse=" + "),sep="")
            for(i in length(object@varlevel):1) st1<-paste(st1," | ",object@varlevel[i])
            cat("    Structured data:",st1,"\n")
            cat("    subject identifier:    ",object@name.group,"\n")
            cat("    predictors:       ",object@name.predictors,"\n")
            cat("    response:         ",object@name.response,paste("(",object@units$y,")",sep=""),"\n")
            cat("    X variable for graphs:",object@name.X,paste("(",object@units$x,")",sep=""),"\n")
            ncov<-length(object@name.covariates)
            if(ncov>0) {
              cat("    covariates:",paste(paste(object@name.covariates," (",object@units$covariates,")",sep=""),collapse=", "),"\n")
              # if(length(object@ocov)>0) {
              #   for(icov in 1:ncov) {
              #     if(is.factor(object@ocov[,icov])) cat("      reference class for covariate",object@name.covariates[icov],": ",levels(object@ocov[,icov])[1],"\n")
              #   }
              #   if(length(object@data)>0) try(object@data[,object@name.covariates]<-object@ocov)
              # }
            }
            cat("Dataset characteristics:\n")
            cat("    number of subjects:    ",object@N,"\n")
            if(object@N>0) {
              cat("    number of observations:",object@ntot.obs,"\n")
              cat("    average/min/max nb obs:",format(mean(object@nind.obs),digits=digits, nsmall=nsmall), " / ", min(object@nind.obs)," / ",max(object@nind.obs),"\n")
              if(length(object@varlevel) > 1 && length(object@var.N) >1) {
                nocc <- object@var.N[[2]]
                cat("    variability levels:    ",paste(object@varlevel, collapse=" > "),"\n")
                cat("    occasions per subject: ",format(mean(nocc),digits=digits, nsmall=nsmall), " / ", min(nocc)," / ",max(nocc)," (mean/min/max)\n")
                cat("    total subject x occasion pairs:",sum(nocc),"\n")
              }
              #    if(length(object@data)>0) print(object@data)
            }
            if(length(object@data)>0) {
              cat("First lines of data:\n")
              nrowShow <- min (10 , nrow(object@data ))
              ncolShow <- min (10 , ncol(object@data))
              print(object@data[1:nrowShow,-c(1)])
            } else cat("No data.\n")
            cat("Warning messages tracked during the creation of the object:\n")
            cat(object@log)
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

####################################################################################
####				Summary method for SaemixData			####
####################################################################################

#' Summarising longitudinal data
#' 
#' summary method for class SaemixData
#' @name summary-methods
#' @aliases summary summary,SaemixData summary,SaemixData-method
#' 
#' @param object an object of class SaemixData
#' @param print a boolean controlling whether to print the output or return it silently
#' @param ... additional arguments (ignored)
#'     
#' @return a list with a number of elements extracted from the dataset
#' \describe{
#' \item{N}{ number of subjects}
#' \item{nobs}{ the total number of observations} 
#' \item{nind.obs}{a vector giving the number of observations for each subject}
#' \item{id}{subject ID; x: predictors; y: response, and, if present in the data, covariates: the covariates (as many lines as observations) and ind.covariates: the individual covariates (one}
#' }
#' @exportMethod summary

setMethod("summary","SaemixData",
          function(object, print=TRUE, ...) {
            digits<-2;nsmall<-2
            if(print) {
              cat("Object of class SaemixData\n")
              cat("    longitudinal data for use with the SAEM algorithm\n")
              cat("Dataset",object@name.data,"\n")
              st1<-paste(object@name.response," ~ ",paste(object@name.predictors,collapse=" + "),sep="")
              for(i in length(object@varlevel):1) st1<-paste(st1," | ",object@varlevel[i])
              cat("    Structured data:",st1,"\n")
              if(length(object@name.predictors)>1) cat("    X variable for graphs:",object@name.X,paste("(",object@units$x,")",sep=""),"\n")
              if(length(object@name.covariates)>0) {
                cat("    covariates:",paste(paste(object@name.covariates," (",object@units$covariates,")",sep=""),collapse=", "),"\n")
              }
              cat("Dataset characteristics:\n")
              cat("    number of subjects:    ",object@N,"\n")
              if(object@N>0) {
                cat("    number of observations:",object@ntot.obs,"\n")
                cat("    average/min/max nb obs:",format(mean(object@nind.obs),digits=digits, nsmall=nsmall), " / ", min(object@nind.obs)," / ",max(object@nind.obs),"\n")
                if(length(object@varlevel) > 1 && length(object@var.N) >1) {
                  nocc <- object@var.N[[2]]
                  cat("    variability levels:    ",paste(object@varlevel, collapse=" > "),"\n")
#                  cat("    occasions per subject: ",format(mean(nocc),digits=digits, nsmall=nsmall), " / ", min(nocc)," / ",max(nocc)," (mean/min/max)\n")
#                  cat("    total subject x occasion pairs:",sum(nocc),"\n")
                }
              }
            }
            res<-list(N=object@N,nobs=list(ntot=object@ntot.obs,nind=object@nind.obs), id=object@data[,object@name.group],x=object@data[,object@name.predictors,drop=FALSE], y=object@data[,object@name.response])
            if(length(object@varlevel) > 1 && length(object@var.N) >1) {
              res$varlevel <- object@varlevel
              res$nocc <- object@var.N[[2]]
              res$ntot.occ <- sum(object@var.N[[2]])
            }
            if(length(object@name.covariates)>0) {
              res$covariates<-object@data[,object@name.covariates]
              ucov<-object@data[,c(object@name.group,object@name.covariates)]
              ucov<-ucov[match(unique(object@data$index),object@data$index),]
              res$ind.covariates<-ucov
            }
            invisible(res)
          }
)
