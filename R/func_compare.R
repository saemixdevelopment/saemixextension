###########################  Model Comparison with BIC or AIC 	#############################

#' Model comparison with information criteria (AIC, BIC).
#'
#' A specific penalty is used for BIC (BIC.cov) when the compared models have in common the
#' structural model and the covariance structure for the random effects (see Delattre et al., 2014).
#'
#' Note that the comparison between two or more models will only be valid if they are
#' fitted to the same dataset.
#'
#' @param \dots Two or more objects returned by the \code{\link{saemix}} function
#' @param method The method used for computing the likelihood : "is" (Importance Sampling),
#' "lin" (Linearisation) or "gq" (Gaussian quadrature). The default is to use importance sampling "is".
#' If the requested likelihood is not available in all model objects, the method stops with an error.
#' @return A matrix of information criteria is returned, with at least two columns containing respectively
#' AIC and BIC values for each of the compared models. When the models have in common the structural model
#' and the covariance structure for the random effects, the matrix includes an additional column with BIC.cov
#' values that are more appropriate when the comparison only concerns the covariates.
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Maud Delattre
#' @references Delattre, M., Lavielle, M. and Poursat, M.A. (2014) A note on BIC in mixed effects models.
#' Electronic Journal of Statistics 8(1) p. 456-475
#' @keywords model comparison, AIC, BIC
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
#' # Model with one covariate
#' saemix.model1<-saemixModel(model=model1cpt,modeltype="structural",
#'   description="One-compartment model, clearance dependent on weight",
#'   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),
#'   transform.par=c(1,1,1),covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE))
#' # Model with two covariates
#' saemix.model2<-saemixModel(model=model1cpt,modeltype="structural",
#'    description="One-compartment model, clearance dependent on weight and volume dependent on sex",
#'    psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),
#'    transform.par=c(1,1,1),covariate.model=matrix(c(0,0,1,0,1,0),ncol=3,byrow=TRUE))
#' # Model with three covariates
#' saemix.model3<-saemixModel(model=model1cpt,modeltype="structural",
#'   description="One-compartment model, clearance and absorption dependent on weight, volume dependent on sex",
#'   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),
#'   transform.par=c(1,1,1),covariate.model=matrix(c(1,0,1,0,1,0),ncol=3,byrow=TRUE))
#'
#' # Running the main algorithm to estimate the population parameters
#' saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
#' saemix.fit1<-saemix(saemix.model1,saemix.data,saemix.options)
#' saemix.fit2<-saemix(saemix.model2,saemix.data,saemix.options)
#' saemix.fit3<-saemix(saemix.model3,saemix.data,saemix.options)
#'
#' # Model comparison
#'
#' compare.saemix(saemix.fit1, saemix.fit2, saemix.fit3)
#' compare.saemix(saemix.fit1, saemix.fit2, saemix.fit3, method = "lin")
#' # We would need to run llgq.saemix first to make this work
#' #compare.saemix(saemix.fit1, saemix.fit2, saemix.fit3, method = "gq")
#'
#' @export compare.saemix

compare.saemix<-function(..., method = c("is", "lin", "gq")) {

  mod.list <- list(...)
  nb.mod <- length(mod.list)

  if (nb.mod < 2){
    stop("'compare.saemix' requires at least two models.")
  }

  list.class <- sapply(mod.list, class)
  if (!all(list.class=="SaemixObject")) {
    stop("All inputs should have class 'SaemixObject'.")
  }

  method <- match.arg(method)

  # This returns TRUE if the requested likelihood is available
  ll_available <- function(x, method) {
    length(slot(x@results, paste0("ll.", method))) > 0
  }

  for (k in 1:nb.mod){
    x <- mod.list[[k]]
    if (!ll_available( x, method)) {
      stop("The likelihood calculated by 'll", method, ".saemix' is not available for model ", k)
    }
  }

  aic_slot <- paste0("aic.", method)
  bic_slot <- paste0("bic.", method)
  bic_cov_slot <- paste0("bic.covariate.", method)

  info <- matrix(0,nb.mod,2)

  rownames(info) <- 1:nb.mod
  colnames(info) <- c("AIC","BIC")

  str.model <- list()
  cov.model <- list()

  # Comparison of datasets : if models are estimated on different datasets, comparison does
  # not make sense.

  same.data <- rep(NA,(nb.mod-1))

  for (k in 2:nb.mod){
    same.data[k-1] <- identical(mod.list[[1]]@data,mod.list[[k]]@data)
  }

  if ("FALSE" %in% same.data){
    stop('Compared models should be fitted on the same data.')
  }

  # Comparison of the statistical models : BIC.covariate does only make sense if model type,
  # structural model (and residual model if appropriate), and covariance structure for the random
  # effects are the same

  for (k  in 1:nb.mod){
    info[k,] <- c(
      slot(mod.list[[k]]@results, aic_slot),
      slot(mod.list[[k]]@results, bic_slot))
  }

  same.model.type <- rep(NA,(nb.mod-1))
  same.str.model <- rep(NA,(nb.mod-1))

  for (k in 2:nb.mod){
    same.model.type[k-1] <- identical(mod.list[[1]]@model@modeltype,mod.list[[k]]@model@modeltype)
    same.str.model[k-1] <- identical(mod.list[[1]]@model@model,mod.list[[k]]@model@model)
  }

  if (!(FALSE %in% same.model.type)){

    if (!(FALSE %in% same.str.model)){

      same.cov.model <- rep(NA,(nb.mod-1))

      for (k in 2:nb.mod){
        same.cov.model[k-1] <- identical(mod.list[[1]]@model@covariance.model,mod.list[[k]]@model@covariance.model)
      }

      if (mod.list[[1]]@model@modeltype == "structural"){

        same.res.model <- rep(NA,(nb.mod-1))

        for (k in 2:nb.mod){
          same.res.model[k-1] <- identical(mod.list[[1]]@model@error.model,mod.list[[k]]@model@error.model)
        }

        if ((!(FALSE %in% same.cov.model)) && (!(FALSE %in% same.res.model))){
          bic.cov <- rep(NA,nb.mod)
          for (k in 1:nb.mod){
            bic.cov[k] <- slot(mod.list[[k]]@results, bic_cov_slot)
          }
          info <- cbind(info, bic.cov)
          colnames(info) <- c("AIC","BIC","BIC.cov")

        }

      } else {
        if (!(FALSE %in% same.cov.model)){
          bic.cov <- rep(NA,nb.mod)
          for (k in 1:nb.mod){
            bic.cov[k] <- slot(mod.list[[k]]@results, bic_cov_slot)
          }
          info <- cbind(info, bic.cov)
          colnames(info) <- c("AIC","BIC","BIC.cov")
        }
      }
    }
  }

  message("Likelihoods calculated by ",
    switch(method,
      is = "importance sampling",
      lin = "linearisation",
      gq = "Gaussian quadrature"))
  return(info)
}
