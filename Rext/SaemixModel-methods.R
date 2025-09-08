## ToDo
### adapt summary, plots, predictions to the new Class SaemixModel

####################################################################################
####				Summary method for SaemixModel			####
####################################################################################

#' @rdname summary-methods
#' @exportMethod summary

setMethod("summary","SaemixModel",
          function(object, print=TRUE, ...) {
            if(print) {
              cat("Nonlinear mixed-effects model\n")
              cat("  Model function")
              if(length(object@description)>0 && nchar(object@description)>0) {
                cat(": ",object@description,"\n")}
              else {
                cat("\n")
                print(object@model)
              }
              fix1<-ifelse(object@fixed.estim==1,""," [fixed]")
              cat("    ",object@nb.parameters,"parameters:", paste(object@name.modpar,fix1,sep=""),"\n")
              cat("     error model:",object@error.model,"\n")
              if(dim(object@covariate.model)[1]>0) {
                cat("     covariate model:\n")
                print(object@covariate.model) 
              } else cat("No covariate\n")
            }
            distrib<-c("normal","log-normal","probit","logit")
            tab.par<-data.frame(Parameter=object@name.modpar, Distribution=distrib[object@transform.par+1], Estimated=ifelse(as.numeric(object@betaest.model[1,])==1,"estimated","fixed"), Initial.value=object@psi0[1,])
            tab.res<-data.frame(parameters=object@name.sigma,Initial.value=object@error.init)   
            res<-list(model=list(modeltype=object@modeltype,model.function=object@model, error.model=object@error.model),parameters=list(fixed=tab.par, residual.error=tab.res),covariance.model=object@covariance.model, covariate.model=object@covariate.model)
            invisible(res)
          }
)

####################################################################################
####			SaemixModel class - method to plot			####
####################################################################################

#' Plot model predictions using an SaemixModel object
#' 
#' This function will plot predictions obtained from an SaemixModel object over a given range of X. Additional predictors may be passed on to the function using the predictors argument.
#' 
## #' @name plot-SaemixModel
#' 
#' @param x an SaemixData object or an SaemixSimData object
#' @param y unused, present for compatibility with base plot function
#' @param range range of X over which the model is to be plotted. Important note: the *first* predictor will be used for the X-axis, the other
#' predictors when present need to be passed sequentially in the predictors argument, in the order in which they appear in the model
#' Less important note: please use explicitely range=XXX where XXX is of the form c(a,b) to pass the plotting range on the X-axis)
#' @param psi parameters of the model 
#' @param predictors additional predictors needed to pass on to the model
#' @param ... additional arguments to be passed on to plot (titles, legends, ...). Use verbose=TRUE to print some messages 
#' concerning the characteristics of the plot
#' 
#' @aliases plot,SaemixModel-methods 
#' @aliases plot,SaemixModel
#' @aliases plot-SaemixModel
#' @keywords methods
### #' @docType methods
#' @exportMethod plot
#' @rdname plot-SaemixModel
#' 
#'@examples
#' # Note that we have written the PK model so that time is the first predictor (xidep[,1]) 
#' # and dose the second
#' model1cpt<-function(psi,id,xidep) { 
#'      tim<-xidep[,1]  
#'      dose<-xidep[,2]
#'      ka<-psi[id,1]
#'      V<-psi[id,2]
#'      CL<-psi[id,3]
#'      k<-CL/V
#'      ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
#'      return(ypred)
#'      }
#'  x<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", 
#'                psi0=matrix(c(1.5,30,1), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))))
#'  # Plot the model over 0-24h, using the parameters given in psi0 and a dose of 300
#'  plot(x, range=c(0,24), predictors=300, verbose=TRUE)
#'  # Plot the model over 0-24h, using another set of parameters and a dose of 350
#'  plot(x, range=c(0,24), psi=c(1.5,20,2), predictors=350, verbose=TRUE)


# Plot simulations from the model
# ECO TODO: test for graphical parameters and set them properly
# ECO TODO: adjust to multiple responses

setMethod("plot","SaemixModel",
          function(x, y , range=c(0,1), psi=NULL, predictors=NULL, ...) {
            # If verbose=TRUE, print messages
            args1<-match.call(expand.dots=TRUE)
            list.args <- list(...)
            i1<-match("verbose",names(args1))
            if(is.na(i1)) verbose<-FALSE else verbose<-eval(args1[[i1]])
            # Set psi by default to the starting parameters given in the model (if not given as arguments)
            if(is.null(psi)) psi<-x@psi0[1,,drop=FALSE]
            if(is.null(dim(psi)[1])) psi<-matrix(psi,nrow=1) else psi<-psi[1,,drop=FALSE]
            npred<-length(x@name.predictors)
            if(npred==0 & is.null(predictors)) npred<-1 else {
              if(npred==0 & !missing(predictors)) {
                npred<-1+length(predictors)
              } else {
                if(npred>1 & (missing(predictors) || length(predictors)<(npred-1))) {
                  if(verbose) message("Please provide the value of the predictors other than X\n")
                  return("Missing predictors")
                }
              }
            }
            if(length(x@name.response)>1) {
              if(verbose) message("Currently the plot can only be obtained for single-response models.\n")
              return()
            }
            if(length(x@name.X)>0 & length(x@name.predictors)>0 && x@name.X != x@name.predictors[1]){
              if(verbose) message("Warning: X predictor supposed to be on the first axis, exiting without plot\n")
              return()
            }
            npts<-100
            #    id<-rep(1,npts+1)
            psi<-matrix(rep(psi, npts+1), byrow=T, nrow=(npts+1))
            id<-matrix(rep(1,npts+1), ncol=1)
            xval<-range[1]+(range[2]-range[1])*c(0:100)/100
            if(npred==1) {
              xdep<-matrix(xval,ncol=1)
            } else {
              xdep<-cbind(xval,matrix(rep(predictors[1:(npred-1)],(npts+1)), byrow=T,nrow=(npts+1)))
              if(length(x@name.X)>0) {
                colnames(xdep)<-c(x@name.X,x@name.predictors[x@name.predictors!=x@name.X])
                xdep<-xdep[,match(x@name.predictors,colnames(xdep))]
              } else colnames(xdep)<-paste("Predictor",1:npred)
            }
            ypred<-try(x@model(psi,id,xdep))
            if(!is.numeric(ypred) & verbose) {
              message("Problem when attempting to obtain predictions from the model.\n")
              message("Usage: plot(x, range=c(0,1), psi, predictors) \n")
              message("Possible solutions can be:\n")
              message("   1. provide suitable values for X (option range=c(<lower bound>, <upper bound>))\n")
              message("   2. provide values for additional predictors (option predictors=c(<value for predictor 1>, <value for predictor 2>, ...)).\n")
              message("   3. check values for the model parameters (defaults to component psi0[1,] of the model).\n")
              message("   4. the predictor used the X-axis is assumed to be in the first column; please check your model is written in a compatible way.\n")
            } else {
              if(length(x@name.X)==0 | length(x@name.predictors)==0) {
                if(verbose) message("Warning: X predictor supposed to be on the first axis\n")}
              if(verbose) message("Plot characteristics:\n")
              if(npred>1) {
                for(j in 1:dim(xdep)[2]) {
                  if(length(x@name.X)==0) {
                    if(j>1) message("   predictor:",colnames(xdep)[j],"=",xdep[1,j],"\n")
                  } else {
                    if(colnames(xdep)[j]!=x@name.X) message("    predictor:",colnames(xdep)[j],"=",xdep[1,j],"\n")
                  }
                }}
              if(verbose) message("   range for X-axis:",min(xval),"-",max(xval),"\n")
              if(verbose) message("   parameters used: ", paste(x@name.modpar,"=",psi[1,],collapse=", "),"\n")
              plot(xval,ypred,type="l",xlab=ifelse(length(x@name.X)==0, "X",x@name.X),ylab=ifelse(length(x@name.response)==0, "Response",x@name.response),...)
            }
          }
)

####################################################################################
####			SaemixModel class - method to obtain predictions given a set of predictors and parameters			####
####################################################################################

#' Predictions for a new dataset
#' 
#' @param object an SaemixModel object
#' @param psi a vector or a dataframe giving the parameters for which predictions are to be computed (defaults to empty). 
#' The number of columns in psi (or the number of elements of psi, if psi is given as a vector) should match the number of
#' parameters in the model, otherwise an error message will be shown and the function will return empty.
#' If psi is NA, the predictions are computed for the population parameters in the model (first line of the psi0 slot).
#' Covariates are not taken into account in the prediction. 
#' If psi is a dataframe, each line will be used for a separate 'subject' in the predictors dataframe, as 
#' indicated by the id argument; if id is not given, only the first line of psi will be used. 
#' @param predictors a dataframe with the predictors for the model (must correspond to the predictors used by the model function), or an SaemixData object (the predictors will then be extracted from the object).
#' @param id a vector of indices of length equal to the number of lines in predictors, matching each line of predictors to the 
#' corresponding line in psi, ie the parameters for this predictors (defaults to empty). If id is given, the unique values in id must be equal
#' to the number of lines in psi, otherwise id will be set to 1. If id is given and its values do not take the consecutive values 1:N, the
#' indices will be matched to 1:N to follow the lines in psi.
#' 
#' @param \dots unused argument, for consistency with the generic
#' 
#' @details The function uses the model slot of the SaemixModel object to obtain predictions, using the predictors object. The
#' user is responsible for giving all the predictors needed by the model function.
#' if psi is not given, the predictions will be computed for the population parameters (first line of the psi0 slot) of the object.
#' 
#' @details The predictions correspond to the structure of the model; for models defined in terms of their likelihood, the predictions 
#' are the log-pdf of the model (see documentation for details).
#' 
#' @details Warning: this function is currently under development and the output may change in future versions of the package 
#' to conform to the usual predict functions.
#' 
#' @return a list with two components
#' \describe{
#' \item{param}{a dataframe with the estimated parameters}
#' \item{predictions}{a dataframe with the population predictions}
#' }
#' 
#' @examples 
#' data(theo.saemix)
#' xpred<-theo.saemix[,c("Dose","Time")]
#' 
#' model1cpt<-function(psi,id,xidep) { 
#' 	  dose<-xidep[,1]
#' 	  tim<-xidep[,2]  
#' 	  ka<-psi[id,1]
#' 	  V<-psi[id,2]
#' 	  CL<-psi[id,3]
#' 	  k<-CL/V
#' 	  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
#' 	  return(ypred)
#' }
#' 
#' saemix.model<-saemixModel(model=model1cpt,modeltype="structural",
#'   description="One-compartment model with first-order absorption", 
#'   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,
#'   dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),
#'   covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
#'   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#'   omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")
#' 
#' head(predict(saemix.model, xpred)$predictions)
#' head(predict(saemix.model, xpred, psi=c(2, 40, 0.5))$predictions)
#' indpsi<-data.frame(ka=2, V=seq(25,47,2), CL=seq(2.5,4.7, 0.2))
#' head(predict(saemix.model, xpred, psi=indpsi)$predictions)
#' 
#' @importFrom stats predict
#' @method predict SaemixModel
#' @export 
#' 

predict.SaemixModel<-function(object, predictors, psi=c(), id=c(), ...) {
  if(!is(predictors,"SaemixData")  & !is(predictors,"data.frame")) {
    message("The predictors argument should be either a dataframe or an SaemixData object to extract the predictors from.\n")
    return()
  }
  if(is(predictors,"SaemixData")) {
    id<-predictors@data[,predictors@name.group]
    predictors <- predictors@data[,predictors@name.predictors,drop=FALSE]
  }
  xidep<-predictors
  if(length(id)==0 || length(id)!=dim(predictors)[1]) 
    id<-rep(1,dim(xidep)[1]) 
  idkeep<-id
  if(max(id)>length(unique(id))) { # indexes need to go from 1 to N
    id1<-1:length(unique(id))
    id2<-unique(id)
    id<-id1[match(id,id2)]
  }
  if(length(psi)==0) psi<-object["psi0"][1,,drop=FALSE]
  if(is.null(dim(psi))) psi<-as.data.frame(t(psi)) # psi given as a vector
  if(dim(psi)[2] != object@nb.parameters) {
    message(paste0("psi must have a number of columns equal to the number of parameters in the model (",object@nb.parameters,")\n")
    )
    return()
  }
  if(dim(psi)[1]==1 & length(unique(id))>1)
    psi<-do.call(rbind,rep(list(psi),length(unique(id))))
  colnames(psi)<-colnames(object["psi0"])
  rownames(psi)<-NULL
  ypred<-object["model"](psi, id, xidep)
  return(list(param=cbind(id=1:dim(psi)[1],psi), predictions=data.frame(id=idkeep, xidep, pred=unname(ypred))))
}

#' Check initial fixed effects for an SaemixModel object applied to an SaemixData object
#' 
#' @param model an SaemixModel object
#' @param data an SaemixData object (the predictors will then be extracted from the object using the name.predictors slot of the object)
#' @param psi a vector or a dataframe giving the parameters for which predictions are to be computed (defaults to empty). 
#' The number of columns in psi (or the number of elements of psi, if psi is given as a vector) should match the number of
#' parameters in the model, otherwise an error message will be shown and the function will return empty.
#' If psi is NA, the predictions are computed for the population parameters in the model (first line of the psi0 slot).
#' Covariates are not taken into account in the prediction. 
#' If psi is a dataframe, each line will be used for a separate 'subject' in the predictors dataframe, as 
#' indicated by the id argument; if id is not given, only the first line of psi will be used. 
#' @param id the vector of subjects for which individual plots will be obtained. If empty, the first 12 subjects in the dataset will be used (subject id's are taken from the name.group slot in the data object). If id is given, individual plots will be shown for the matching subjects in the dataset (eg if id=c(1:6), the first 6 subjects in the dataframe will be used for the plots, retrieving their ID from the data object)
#' 
#' @param \dots unused argument, for consistency with the generic
#' 
#' @details The function uses the model slot of the SaemixModel object to obtain predictions, using the predictors object. The
#' user is responsible for giving all the predictors needed by the model function.
#' if psi is not given, the predictions will be computed for the population parameters (first line of the psi0 slot) of the object.
#' 
#' @details The predictions correspond to the structure of the model. For models defined as a structural model, 
#' individual plots for the subjects in id overlaying the predictions for the parameters psi and the individual data
#' are shown, and the predictions correspond to f(t_ij, psi).
#' For models defined in terms of their likelihood, the predictions returned correspond to the log-likelihood.
#' No individual graphs are currently available for discrete data models.
#' 
#' @details Warning: this function is currently under development and the output may change in future versions of the package 
#' 
# #' @seealso \code{\link[predict.SaemixModel]}, \code{\link[plotDiscreteData]},  \code{\link[ggplot]}
#' 
#' @return the predictions corresponding to the values for each observation in the predictors of either the model f or log-likelihood. 
#' For Gaussian data models, the function also plots the data overlayed with the model predictions for each subject in id 
#' (where id is the index in the N subjects).
#' 
#' @examples 
#' data(theo.saemix)
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
#'   name.group=c("Id"),name.predictors=c("Dose","Time"),
#'   name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#'   units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
#' 
#' model1cpt<-function(psi,id,xidep) { 
#' 	  dose<-xidep[,1]
#' 	  tim<-xidep[,2]  
#' 	  ka<-psi[id,1]
#' 	  V<-psi[id,2]
#' 	  CL<-psi[id,3]
#' 	  k<-CL/V
#' 	  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
#' 	  return(ypred)
#' }
#' 
#' saemix.model<-saemixModel(model=model1cpt,modeltype="structural",
#'   description="One-compartment model with first-order absorption", 
#'   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,
#'   dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),
#'   covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
#'   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#'   omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")
#' 
#' checkInitialFixedEffects(saemix.model, saemix.data, id=c(1:6))
#' checkInitialFixedEffects(saemix.model, saemix.data, id=c(1:6), psi=c(0.5, 30, 2)) # better fit
#' 
#' @export 

checkInitialFixedEffects<-function(model, data, psi=c(), id=c(), ...) {
  addarg <- list(...)
  verbose<-FALSE
  if("verbose" %in% names(addarg)) verbose <- as.logical(addarg["verbose"])
  if(is.na(verbose)) verbose<-FALSE
  if(!is(data,"SaemixData") ) {
    if(verbose) message("The data argument should be an SaemixData object to extract the predictors from.\n")
    return()
  }
  if(length(psi)==0) psi<-model["psi0"][1,,drop=FALSE]
  if(is.null(dim(psi))) psi<-as.data.frame(t(psi)) # psi given as a vector
  if(dim(psi)[2] != model@nb.parameters) {
    message(paste0("psi must have a number of columns equal to the number of parameters in the model (",model@nb.parameters,")\n")
    )
    return()
  }
  idall<-data@data[,"index"]
  xidep <- data@data[,data@name.predictors,drop=FALSE]
  #  if(dim(psi)[1]==1 & length(unique(idall))>1)
  psi<-do.call(rbind,rep(list(psi),length(unique(idall))))
  colnames(psi)<-colnames(model["psi0"])
  rownames(psi)<-NULL
  # Model predictions
  ypred<-model["model"](psi, idall, xidep)
  
  # For the plot, select subjects corresponding to number id, or use the first 12 subjects
  idplot <- intersect(idall, id)
  if(length(idplot)==0) idplot<-1:12
  idkeep <- which(data@data$index %in% idplot) # retrieve data for these subjects
  
  if(model@modeltype=="structural") {
    # Individual graphs
    obspl<-data.frame(id=data@data[idkeep,data@name.group], x=data@data[idkeep,data@name.X], y=data@data[idkeep,data@name.response])
    predpl<-data.frame(id=data@data[idkeep,data@name.group], x=data@data[idkeep,data@name.X], y=ypred[idkeep])
    myplot <- ggplot(data=obspl, aes(x=.data$x, y=.data$y, group=.data$id)) + geom_point() + geom_line(data=predpl) + 
      xlab(paste0(data@name.X," (",data@units$x,")")) + ylab(paste0(data@name.response," (",data@units$y,")")) + 
      theme_bw() + facet_wrap(.~id, nrow=3, ncol=4)
    print(myplot)
  } else {
    if(verbose) message("Individual plots are only available for models dealing with continuous outcomes.\n")
  }
  invisible(ypred)
}

####################################################################################
####			SaemixModel & SaemixData class - method to plot	predictions from a model for the data in a dataset		####
####################################################################################

#' Plot model predictions for a new dataset. If the dataset is large, only the first 20 subjects (id's) will be shown.
#' 
#' @param x an SaemixModel object
#' @param y an SaemixData object
#' @param ... additional arguments. Passing psi=X where X is a vector or a dataframe will allow
#' changing the parameters for which predictions are to be computed (defaults to the population parameters
#' defined by the psi element of x) (see details)
#' 
#' @details The function uses the model slot of the SaemixModel object to obtain predictions, using the dataset contained in the 
#' SaemixData object. The user is responsible for making sure data and model match.
#' If psi is not given, the predictions will be computed for the population parameters (first line of the psi0 slot) of the object.
#' If psi is given, the number of columns in psi (or the number of elements of psi, if psi is given as a vector) should match 
#' the number of parameters in the model, otherwise an error message will be shown and the function will return empty.
#' If psi is a dataframe, each line will be used for a separate subject of the smx.data object. Elements of psi will be recycled 
#' if psi has less lines than the number of subjects in the dataset.
#' 
#' @details Currently this function only works for models defined as 'structural'.
#' 
#' @return a ggplot object
#' 
#' @aliases plot.SaemixModel
#' 
#' @examples 
#' data(theo.saemix)
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA,
#'    name.group=c("Id"),name.predictors=c("Dose","Time"),
#'    name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#'    units=list(x="hr",y="mg/L", covariates=c("kg","-")), name.X="Time")
#' 
#' model1cpt<-function(psi,id,xidep) { 
#' 	  dose<-xidep[,1]
#' 	  tim<-xidep[,2]  
#' 	  ka<-psi[id,1]
#' 	  V<-psi[id,2]
#' 	  CL<-psi[id,3]
#' 	  k<-CL/V
#' 	  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
#' 	  return(ypred)
#' }
#' 
#' saemix.model<-saemixModel(model=model1cpt,modeltype="structural",
#'   description="One-compartment model with first-order absorption", 
#'   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,
#'   dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),
#'   covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
#'   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#'   omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")
#' 
#' plot(saemix.model, saemix.data)
#' plot(saemix.model, saemix.data, psi=c(2, 40, 3))
#' indpsi<-data.frame(ka=2, V=seq(25,47,2), CL=seq(2.5,4.7, 0.2))
#' plot(saemix.model, saemix.data, psi=indpsi)
#' 
#' @exportMethod plot
#' 

# Plot the data, either as points or as lines grouped by x@name.group
setMethod("plot",c("SaemixModel","SaemixData"),
          function(x, y, ...) {
            if(x@modeltype[1]!="structural") {
              message("Currently plots of the model are only available for continuous response models\n")
              return()
            }
            args1<-match.call(expand.dots=TRUE)
            list.args <- list(...)
            i1<-match("verbose",names(args1))
            if(!is.na(i1)) verbose<-FALSE else verbose<-eval(args1[[i1]])
            i1<-match("psi",names(args1))
            if(!is.na(i1)) psi<-eval(args1[[i1]]) else psi<-x["psi0"][1,,drop=FALSE]
            if(is.null(dim(psi))) psi<-as.data.frame(t(psi)) # psi given as a vector
            if(dim(psi)[2] != x@nb.parameters) {
              message(paste0("psi must have a number of columns equal to the number of parameters in the model (",x@nb.parameters,")\n"))
              return()
            }
            if(dim(psi)[1]==1 || dim(psi)[1]<y@N)
              psi<-do.call(rbind,rep(list(psi),length.out=y@N))
            i1<-match("ilist",names(args1))
            if(!is.na(i1)) ilist<-eval(args1[[i1]]) else ilist<-NA
            
            nvalues<-100
            xt<-seq(min(y@data[,y@name.X]), max(y@data[,y@name.X]), length.out=nvalues)
            xidep<-data.frame(x=xt)
            colnames(xidep)<-y@name.X
            id<-y@data[,y@name.group]
            if(length(y@name.predictors)>1) {
              otherpred<-y@name.predictors[y@name.predictors != y@name.X]
              x1<-y@data[match(unique(id), id), otherpred, drop=FALSE]
              dat1<-NULL
              for(i in 1:length(unique(id)))
                dat1<-rbind(dat1, 
                            do.call(rbind,rep(list(x1[i,,drop=FALSE]), nvalues)))
              xidep<-cbind(xidep, dat1)
              colnames(xidep[2:dim(xidep)[2]])<-otherpred
              xidep<-xidep[,y["name.predictors"]] # Sort the predictors back in the correct order...
            }
            id<-rep(1:length(unique(id)), each=nvalues)
            ypred<-predict(x, predictors=xidep, psi=psi, id=id)
            gpred<-cbind(id=id,xidep,y=ypred$predictions$pred)
            colnames(gpred)[colnames(gpred)==y@name.X]<-"x"
            
            gdat<-y@data
            colnames(gdat)[colnames(gdat)==y@name.X]<-"x"
            colnames(gdat)[colnames(gdat)==y@name.response]<-"y"
            colnames(gdat)[colnames(gdat)==y@name.group]<-"id"
            zesuj<-unique(gdat$id)
            if(!is.na(ilist)[1]) {
              if(is(ilist,"numeric")) ilist<-zesuj[intersect(ilist, 1:length(zesuj))] else ilist<-intersect(ilist, zesuj)
            } else ilist<-zesuj[1:min(20, length(zesuj))]
            gdat1<-gdat[gdat$id %in% ilist,]
            gpred1<-gpred[gpred$id %in% ilist,]
            if(length(unique(gdat$id))>20) {
              nrow<-4
              ncol<-5
            } else {
              nrow<-NULL
              ncol<-NULL
            }
            
            g1<-ggplot(data=gdat1, aes(x=.data$x, y=.data$y, group=.data$id)) + geom_point() + geom_line(data=gpred1,aes(x=.data$x, y=.data$y)) +
              facet_wrap(.~id, nrow=nrow, ncol=ncol) + 
              labs(x=y@name.X, y=y@name.response) + theme_bw()
            return(g1)
          } 
)


