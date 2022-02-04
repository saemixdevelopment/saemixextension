###########  Choose a model for covariates and random effects simultaneously by BIC in a stepwise algorithm  ###########

#' Stepwise procedure for joint selection of covariates and random effects
#' 
#' Joint selection of covariates and random effects in a nonlinear mixed effects model by a stepwise-type
#' algorithm based on two different versions of BIC for covariate selection and random effects selection
#' respectively. Selection is made among the covariates as such specified in the SaemixData object.
#' Only uncorrelated random effects structures are considered.
#' 
#' @param saemixObject An object returned by the \code{\link{saemix}} function
#' @param trace If TRUE, a table summarizing the steps of the algorithm is printed. Default "TRUE"  
#' @param direction The mode of stepwise search, can be one of "both", "backward", or "forward", 
#' with a default of "forward". 
#' @param covariate.init A matrix specifying the initial covariate structure to be considered in the algorithm.
#' covariate.init is only used if the direction argument is "both". 
#' @return An object of the SaemixObject class storing the covariate model and the covariance structure of 
#' random effects of the final model.
#' @author Maud Delattre
#' @references M Delattre, M Lavielle, MA Poursat (2014) A note on BIC in mixed effects models. 
#' Electronic Journal of Statistics 8(1) p. 456-475
#' M Delattre, MA Poursat (2017) BIC strategies for model choice in a population approach. 
#' (arXiv:1612.02405)
#' @keywords selection stepwise covariate
#' @examples 
#' data(theo.saemix)
#' 
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA,
#'   name.group=c("Id"),name.predictors=c("Dose","Time"),
#'   name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#'   units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
#' 
#' # Definition of models to be compared
#' model1cpt<-function(psi,id,xidep) { 
#'    dose<-xidep[,1]
#'    tim<-xidep[,2]  
#'    ka<-psi[id,1]
#'    V<-psi[id,2]
#'    CL<-psi[id,3]
#'    k<-CL/V
#'    ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
#'    return(ypred)
#' }
#' 
#' saemix.model1<-saemixModel(model=model1cpt,modeltype="structural", 
#'      description="One-compartment model with first-order absorption",
#'      psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, 
#'      dimnames=list(NULL, c("ka","V","CL"))), transform.par=c(1,1,1),
#'      covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE))
#'                           
#' saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
#' saemix.fit1<-saemix(saemix.model1,saemix.data,saemix.options)     
#' 
#' \dontrun{
#' res.forward <- step.saemix(saemix.fit1, direction = "forward")
#' res.backward <- step.saemix(saemix.fit1, direction = "backward")
#' covariate.init <- matrix(c(1,0,0,0,1,0),ncol=3,nrow=2)
#' res.stepwise <- step.saemix(saemix.fit1, direction="both")
#' }
#' @export step.saemix


step.saemix <- function(saemixObject, trace = TRUE, direction = "forward", covariate.init = NULL){
  
  saemix.data <- saemixObject["data"]
  
  if (class(saemixObject)=="SaemixObject"){
  
  if (length(saemix.data@name.covariates)==0){
    stop("The saemixData object should contain covariates.")
  } else{
    
    if (direction %in% c("forward","backward","both")){
      if (direction=="forward"){
        res <- forward.procedure(saemixObject,trace)
      }
      if (direction=="backward"){
        warning("Backward procedures are not recommended if initial complete model is too complex.")
        res <- backward.procedure(saemixObject,trace)
      }
      if (direction=="both"){
        res <- stepwise.procedure(saemixObject,covariate.init,trace)
      }
    } else (
      stop('Invalid argument direction.')
    )
  }
  } else{stop("Invalid class for argument saemixObject.")}
  
  return(res)
}





