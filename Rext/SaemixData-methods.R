#' @include aaa_generics.R
#' @include SaemixData.R
#' @include SaemixData-printmethods.R
#' @include SaemixData-covariatemethods.R
NULL

####################################################################################
####	      		SaemixData class - method to read data		                    	####
####################################################################################

#' Name validation (## )Helper function not intended to be called by the user)
#' 
#' Helper function, checks if the names given by the user match to the names in the dataset. If not, automatic recognition is attempted when automatic=TRUE.
#'
#' @name validate.names
#' 
#' @param usernames vector of strings
#' @param datanames vector of strings
#' @param recognisednames vector of strings, values for automatic recognition
#' @param verbose logical, whether to print warning messages
#' @param automatic a boolean indicating whether to attempt automatic name recognition when some colum names are missing or wrong (defaults to TRUE)
#' @return a vector with valid names
#' @examples 
#' # TODO
#' @keywords methods
#' @export 
NULL

validate.names<-function(usernames,datanames,recognisednames=c(),verbose=TRUE, automatic=TRUE) {
  valnames<-usernames[usernames!=""]
  if(length(valnames)>0) {
    remcol<-c() # keep track of columns to remove
    # Detect names given as column numbers
    convnames<-suppressWarnings(as.integer(usernames))
    icol.int<-which(!is.na(convnames))
    if(length(icol.int)>0) {
      namcol.int<-suppressWarnings(as.integer(usernames[icol.int]))
      ioutrange<-namcol.int[namcol.int>length(datanames) | namcol.int<0]
      if(length(ioutrange)>0) {
        if(verbose) cat("Column number(s)",ioutrange,"do(es) not exist in the dataset, please check\n")
        remcol<-icol.int[namcol.int %in% ioutrange]
        icol.int<-icol.int[!namcol.int %in% ioutrange]
        namcol.int<-namcol.int[!namcol.int %in% ioutrange]
      }
      valnames[icol.int]<-datanames[namcol.int]
    }
    # Validate name columns given as strings
    icol.str<-which(is.na(convnames))
    if(length(icol.str)>0) {
      namcol.str<-usernames[icol.str]
      inotfound<-which(is.na(match(namcol.str,datanames)))
      if(length(inotfound)>0) {
        remcol<-c(remcol,icol.str[inotfound])
        if(verbose & namcol.str[inotfound]!="") cat("Column name(s)",namcol.str[inotfound],"do(es) not exist in the dataset, please check\n")
      }
    }
    if(length(remcol)>0) {
      if(verbose) cat("Remove columns",remcol,"(",usernames[remcol],")","\n")
      valnames<-valnames[-remcol]
    }
  }

  # Automatic column detection (only if no name given or none recognised)
  if(length(valnames)==0) {
    if(length(recognisednames)==0) {
      if(verbose) cat("Please check input\n")
      return(valnames)
    } 
    if(automatic) {
      if(verbose) cat("No valid name given, attempting automatic recognition\n")
      valnames<-datanames[match(recognisednames,tolower(datanames),nomatch=0)]
      if(length(valnames)==0) {
        if(verbose) cat("Automatic recognition failed, please check input\n")} else {
          if(verbose) cat("Automatic recognition of columns",valnames,"successful \n")
        }
    } else return(NULL)
  }
  # Remove columns
  return(valnames)
}

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
            if(!inherits(object@name.data,"character")) {
              if(object@messages) message("Please provide the name of the data (data.frame or path to file on disk) as a character string.\n")
              return("Creation of SaemixData object failed")
            }
            if(!object@automatic) {
              if(length(object@name.group)==0 | length(object@name.predictors)==0 | length(object@name.response)==0) {
                if(object@messages) message("Please provide the name of the relevant columns in the dataset (minimally: grouping level, predictors and response) or set automatic to TRUE to attempt automatic column recognition\n")
                return("Creation of SaemixData object failed")
              }
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
            
            varlevel<-object@varlevel
            object@name.group<-as.character(varlevel[1])
            if(length(object@name.group)>0)
              if(!is.na(as.numeric(object@name.group))) object@name.group <- colnames(dat)[as.numeric(object@name.group)]
            vnames<-validate.names(object@name.group,colnames(dat),recognisednames=c("id","subject","sujet","group","groupe"),verbose = object@messages, automatic=object@automatic)
            if(length(vnames)==0) {
              if(object@messages) message("Please provide a valid name for the ID column or set automatic to TRUE to attempt automatic column recognition.\n")
              return("Creation of SaemixData object failed")
            }
            if(!identical(object@name.group,vnames)) {
              msg<-paste0("   Name of the grouping factor changed from ",object@name.group," to ",vnames,"\n")
              object@log<-paste0(object@log,msg)
              
            }
            object@name.group<-vnames
            varlevel<-c(vnames)

            if(length(object@name.predictors)>0) {
              for(i1 in 1:length(object@name.predictors)) {
              if(!is.na(as.numeric(object@name.predictors[i1]))) object@name.predictors[i1] <- colnames(dat)[as.numeric(object@name.predictors[i1])]
            }}
            vnames<-validate.names(object@name.predictors,colnames(dat),recognisednames=c("time","temps","tps","tim","x","dose"),verbose = object@messages, automatic=object@automatic)
            if(length(vnames)==0) {
              if(object@messages) message("Please provide a valid name for the predictor(s) or set automatic to TRUE to attempt automatic column recognition.\n")
              return("Creation of SaemixData object failed")
            }
            if(!identical(object@name.predictors,vnames)) {
              msg<-paste0("   Name of the predictors changed from ",object@name.predictors," to ",vnames,"\n")
              object@log<-paste0(object@log,msg)
              
            }
            object@name.predictors<-vnames
            
            if(length(object@name.response)>0)
              if(!is.na(as.numeric(object@name.response))) object@name.response <- colnames(dat)[as.numeric(object@name.response)]
            vnames<-validate.names(object@name.response,colnames(dat),recognisednames=c("response","resp","conc","concentration","y","dv"),verbose = object@messages, automatic=object@automatic)
            if(length(vnames)==0) {
              if(object@messages) message("Please provide a valid name for the response.\n")
              return("Creation of SaemixData object failed")
            }
            if(!identical(object@name.response,vnames)) {
              msg<-paste0("   Name of the response changed from ",object@name.response," to ",vnames,"\n")
              object@log<-paste0(object@log,msg)
            }
            object@name.response<-vnames[1]
            if(length(vnames)>1) {
              if(object@messages) cat("Using the response",object@name.response,"as dependent variable.\n")
            } 
            
            if(length(object@name.covariates)>0) {
              if(object@name.covariates[1]!="") {
                i1<-suppressWarnings(as.integer(object@name.covariates[!is.na(suppressWarnings(as.integer(object@name.covariates)))]))
                object@name.covariates[!is.na(suppressWarnings(as.integer(object@name.covariates)))]<- colnames(dat)[i1]
              }
              idx<-object@name.covariates[!(object@name.covariates %in% colnames(dat))]
              if(length(idx)>0) {
                if(object@messages) cat("Covariates",object@name.covariates[idx],"not found.\n") 
                msg<-paste0("  Covariates ",object@name.covariates[idx],"not found\n")
                object@log<-paste0(object@log,msg)
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
            if(length(object@name.occ)==0) {
              occ<-rep(1,dim(dat)[1])
              if(!("occ" %in% colnames(dat))) dat<-cbind(dat,occ=occ)
              } else {
              occ<-dat[,object@name.occ]
              varlevel[2]<-object@name.occ
            }
            object@varlevel <- varlevel
#            object@name.occ<-"occ" # changed... now second level of variability may have another name
            if(nchar(object@name.ytype)==0) ytype<-rep(1,dim(dat)[1]) else ytype<-dat[,object@name.ytype]
            object@name.ytype<-"ytype"
            all.names<-c(object@varlevel,object@name.predictors, object@name.response,object@name.covariates)
            
            dat<-dat[,all.names,drop=FALSE]
            dat<-cbind(dat,mdv=mdv,cens=cens,ytype=ytype)
            
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
                msg<-paste0("   Removing missing values for predictor ",i,"\n")
                object@log<-paste0(object@log,msg)
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
              msg<-paste0("   Removing subjects with no observations, id= ",inull,"\n")
              object@log<-paste0(object@log,msg)
              dat<-dat[-idx,]
              if(!is.null(dim(object@ocov)[1])) object@ocov<-object@ocov[-idx,,drop=FALSE]
            }
            # ECO TODO: missing data in covariates kept for the moment, only excluded depending on the model
            #    for(i in object@name.covariates) dat<-dat[!is.na(dat[,i]),]
            object@ntot.obs<-dim(dat)[1] # total number of observations
            # Sorting in the order of levels of variability then name.X (=usually time)
            datord<-list(dat[,object@varlevel[[1]]])
            if(length(varlevel)>1){
              for(ilev in 2:length(varlevel)) datord[[ilev]]<-dat[,object@varlevel[[ilev]]]
            }
            datord[[(length(varlevel)+1)]]<-dat[,object@name.X]
            dat <- dat[do.call("order",datord),]
            if(!is.null(dim(object@ocov)[1])) object@ocov <- object@ocov[do.call("order",datord),,drop=FALSE]
            id<-dat[,object@name.group]
            object@N<-length(unique(id))
            nind.obs<-tapply(id,id,length) # individual numbers of observations (1xN)
            nind.obs<-nind.obs[match(unique(id),names(nind.obs))]
            object@nind.obs<-c(nind.obs)
            dat<-cbind(index=rep(1:object@N,times=nind.obs),dat)
            object@varlevel <- varlevel
            index<-list(dat$index)
            nb.obs<-list(nind.obs)
            Nunit<-list(object@N) # or list(rep(1,object@N)) ?
            if(length(varlevel)>1) {
              ilab<-dat[,object@varlevel[1]]
              for(ilev in 2:length(varlevel)) {
                ilab<-cbind(ilab,dat[,object@varlevel[ilev]])
                ilab1<- ilab[!duplicated(ilab),]
                Nunit[[ilev]] <- tapply(ilab1[,1],ilab1[,1],length)
                idx<-c(which(!duplicated(ilab)),dim(ilab)[1]+1)
                x1<-diff(idx) # individual numbers of observations per unit at that varlevel
                print(head(x1),20)
                #  x1<-tapply(ilab,ilab,length) # weird sort issue !!
                index[[ilev]]<-rep(1:length(x1),times=x1)
                nb.obs[[ilev]]<-x1
              }
              ilab<-paste0(ilab,"-",dat[,object@name.X])
              dat$index <- index[[ilev]] # if more than 1 level of variability, index needs to be reset to the last variability level
            } 
            object@data<-dat
            names(index)<-varlevel
            names(nb.obs)<-varlevel
            object@var.index<-index
            object@var.nobs<-nb.obs
            object@var.N<-Nunit
            options(ow) # reset
            validObject(object)
            return(object)
          }
)

####################################################################################
####		Creating an object of SaemixData class - User-level function	####
####################################################################################

#' Function to create an SaemixData object
#' 
#' This function creates an SaemixData object. The only mandatory argument is
#' the name of the dataset. If the dataset has a header (or named columns), the
#' program will attempt to detect which column correspond to ID, predictor(s)
#' and response. Warning messages will be printed during the object creation
#' and should be read for details.
#' 
#' Multiple responses can be specified through the ytype argument (1 for the first response, 2 for the second, etc...)
#' and each can be given a name and a type (continuous, categorical or event)
#' 
#' This function is the user-friendly constructor for the SaemixData object
#' class. The read.saemixData is a helper function, used to read the dataset,
#' and is not intended to be called directly.
#' 
#' @name saemixData
#' 
#' @param name.data name of the dataset (can be a character string giving the name of a file on disk or of a dataset in the R session, or the name of a dataset
#' @param header whether the dataset/file contains a header. Defaults to TRUE
#' @param sep the field separator character. Defaults to any number of blank spaces ("")
#' @param na a character vector of the strings which are to be interpreted as NA values. Defaults to c(NA)
#' @param name.group name (or number) of the column containing the subject id
#' @param name.predictors name (or number) of the column(s) containing the predictors (the algorithm requires at least one predictor x)
#' @param name.response name (or number) of the column containing the response variable y modelled by predictor(s) x
#' @param outcome a vector of the form c(y1="continuous", y2="categorical", y3="event",...) specifying the name (y1, y2,..) and type (one of "continuous","categorical" (for count, binary, categorical) or "event" (survival-type data)) of each outcome. If not given, saemix will set the response names to y1, y2, etc... and the type to "continuous" unless the corresponding response only has 2 values, in which case that response will be set to "categorical"
#' @param name.covariates name (or number) of the column(s) containing the covariates, if present (otherwise missing)
#' @param name.genetic.covariates name (or number) of the column(s) containing the covariates, if present (otherwise missing)
#' @param name.mdv name of the column containing the indicator for missing variable
#' @param name.cens name of the column containing the indicator for censoring
#' @param name.occ name of the column containing the occasion
#' @param name.ytype name of the column containing the index of the response
#' @param name.X name of the column containing the regression variable to be used on the X axis in the plots (defaults to the first predictor)
#' @param units list with up to three elements, x, y and optionally covariates, containing the units for the X and Y variables respectively, as well as the units for the different covariates (defaults to empty)
#' @param verbose a boolean indicating whether messages should be printed out during the creation of the object
#' @param automatic a boolean indicating whether to attempt automatic name recognition when some colum names are missing or wrong (defaults to TRUE)
#' @details This function is the user-friendly constructor for the SaemixData object class. The read is a helper function, used to read the dataset, and is not intended to be called directly.
#' @return An SaemixData object (see \code{\link{saemixData}}).
#' @references E Comets, A Lavenu, M Lavielle M (2017). Parameter estimation in nonlinear mixed effect models using saemix,
#' an R implementation of the SAEM algorithm. Journal of Statistical Software, 80(3):1-41.
#' 
#' E Kuhn, M Lavielle (2005). Maximum likelihood estimation in nonlinear mixed effects models. 
#' Computational Statistics and Data Analysis, 49(4):1020-1038.
#' 
#' E Comets, A Lavenu, M Lavielle (2011). SAEMIX, an R version of the SAEM algorithm. 20th meeting of the 
#' Population Approach Group in Europe, Athens, Greece, Abstr 2173.
#' 
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}, Audrey Lavenu, Marc Lavielle.
#' @seealso \code{\link{SaemixData}},\code{\link{SaemixModel}}, \code{\link{saemixControl}},\code{\link{saemix}}
#' @examples
#' 
#' data(theo.saemix)
#' # Legacy
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
#'   name.group=c("Id"),name.predictors=c("Dose","Time"),
#'   name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#'   units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
#'   
#' # New version
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
#'   varlevel="Id",name.predictors=c("Dose","Time"),
#'   name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#'   units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
#' 
#' print(saemix.data)
#' 
#' plot(saemix.data)
#' @export 

saemixData<-function(name.data,header,sep,na,name.group, name.predictors, name.response,  name.X, outcome=character(), name.covariates=c(), name.genetic.covariates=c(), name.mdv="", name.cens="", name.occ=character(0), varlevel=c(), covariate.varlevel=c(), name.ytype="", units=list(x="",y="",covariates=c()), verbose=TRUE, automatic=TRUE) {
  # setting proper types for the SaemixData class
  if(missing(name.data)) {
    if(verbose) cat("Error in saemixData: please provide the name of the datafile or dataframe (between quotes)\n")
    return("Creation of SaemixData object failed")
  }
  if(is.data.frame(name.data)) {
    data_from_name.data <- TRUE
    dat <- name.data
    name.data<-deparse(substitute(name.data))
  } else {
    data_from_name.data <- FALSE
  }
  if(length(varlevel)==0) {
    is.legacy<-TRUE
    if(missing(name.group)) name.group<-"" else {
      name.group<-as.character(name.group)
      varlevel<-c(name.group)
    }
    if(missing(name.occ)) name.occ<-character(0) else  {
      name.occ<-as.character(name.occ)
      varlevel<-c(varlevel,name.occ)
    }
  } else {
    is.legacy<-FALSE
    name.group <- as.character(varlevel[1])
  }
  if(missing(header)) header<-TRUE
  if(missing(sep)) sep<-""
  if(missing(na)) na<-"NA" else {na<-as.character(na);na[is.na(na)]<-"NA"}
  if(missing(name.predictors)) name.predictors<-"" else name.predictors<-as.character(name.predictors)
  if(missing(name.response)) name.response<-"" else  name.response<-as.character(name.response)
  if(missing(name.mdv)) name.mdv<-"" else  name.mdv<-as.character(name.mdv)
  if(missing(name.cens)) name.cens<-"" else  name.cens<-as.character(name.cens)
  if(missing(name.ytype)) name.ytype<-"" else  name.ytype<-as.character(name.ytype)
  if(missing(name.X)) name.X<-"" else name.X<-as.character(name.X)
  if(!is.character(outcome)) outcome<-as.character(outcome)
  name.covariates<-c(as.character(name.covariates),as.character(name.genetic.covariates))
  x<-new(Class="SaemixData",name.data=name.data,header=header,sep=sep,na=na, name.predictors=name.predictors,name.X=name.X, name.response=name.response,name.covariates=name.covariates, units=units, name.mdv=name.mdv, name.cens=name.cens, name.ytype=name.ytype, varlevel=varlevel, covariate.varlevel=covariate.varlevel, verbose=verbose, automatic=automatic)
  #  showall(x)
  if(data_from_name.data) {
    x1<-readSaemix(x, dat)
  } else
    x1<-readSaemix(x)
  if(is(x1,"SaemixData")) {
    if(x@name.X=="") {
      x@name.X<-x@name.predictors[1]
      if(length(x@name.predictors)>1) x@log<-paste0(x@log, "   Assuming ",x@name.X," is the predictor on the X-axis\n")
    }
    # Genetic covariates
    igen<-rep(FALSE,length(name.covariates))
    igen[match(name.genetic.covariates,name.covariates)]<-TRUE
    x1@ind.gen<-igen
    # Outcomes
    if(identical(outcome,character())) { # Outcome not given
      if(length(unique(x1@data$ytype[x1@data$ytype>0]))==1) { # Only one response
        if(length(unique(x1@data[,x1@name.response]))>2) {
          x1@outcome<-c(y1="continuous")
          names(x1@outcome)<-x1@name.response
          if(verbose) cat("Assuming response is continuous\n")
        } else {
          x1@outcome<-c(y1="categorical")
          names(x1@outcome)<-x1@name.response
          if(verbose) cat("Only two modalities in the response, assuming the response is binary\n")
        }
      } else { # More than one response
        nresp<-length(unique(x1@data$ytype[x1@data$ytype>0])) # potentially we want to pass doses as ytype=0
        vec<-c()
        for(iresp in 1:nresp)
          if(length(unique(x1@data[x1@data$ytype==iresp,x1@name.response]))>2) vec<-c(vec,"continuous") else vec<-c(vec,"categorical")
        names(vec)<-paste("y",1:nresp,sep="")
        x1@outcome<-vec
      }
    } else { # Use the outcome vector given to specify name and type of outcomes
      if(is.null(names(outcome))) { # if only type is given, set outcome names to y1, y2,...
        if(length(outcome)==1) names(vec)<-x1@name.response else names(outcome)<-paste("y",1:length(outcome),sep="")
      }
      i1<-which(names(outcome)=="")
      if(length(i1)>0) names(outcome)[i1]<-paste("y",i1,sep="") # complete missing outcome names
      nresp<-length(unique(x1@data$ytype[x1@data$ytype>0]))
      if(length(outcome)>nresp) {
        if(verbose) cat(paste("Size of argument outcome longer than the number of responses in the data, using only the first", nresp,"first elements\n"))
        outcome<-outcome[1:nresp]
      }
      if(length(outcome)<nresp) {
        for(iresp in (length(outcome)+1):nresp) {
          if(length(unique(x1@data[x1@data$ytype==iresp,name.response]))>2) outcome<-c(outcome,"continuous") else outcome<-c(outcome,"categorical")
          names(outcome)[iresp]<-paste("y",iresp,sep="")
        }
      }
      x1@outcome<-outcome
    }
    # Extract covariate matrices corresponding to variability levels defined in the data
    # ToDo: automated sorting of variability levels (eg: if no variability at occ level, set to iiv)
    x1 <- extractVarlevelCovariates(x1,verbose=verbose)
    if(verbose) cat("\n\nThe following SaemixData object was successfully created:\n\n")
  }
  if(verbose) print(x1,nlines=0)
  return(x1)
}
