###########################  Model Comparison with BIC or AIC 	#############################

#' Compares two or more models with information criteria.
#' 
#' A specific penalty can be used in BIC when the compared models have in common the 
#' structural model and the covariance structure for the random effects (BIC.cov). 
#' 
#' @param mod.list A list of two objects returned by the \code{\link{saemix}} function
#' @param method The method used for computing the likelihood : "is" (Importance Sampling), 
#' "ll" (Linearisation) or "gq" (Gaussian quadrature). Default "is"
#' @return A matrix of information criteria is returned
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{SaemixObject}},\code{\link{saemix}}
#' @references COMPLETER
#' @keywords models
#' @examples
#' @export compare.saemix


compare.saemix<-function(mod.list,...) {
  
  if (!is.list(mod.list)){
    stop("'mod.list' must be a list.")
  } else {
    if (length(mod.list)<=1){
      stop("'compare.saemix' requires at least two models.") 
    } else {
      list.class <- sapply(mod.list, class)
      if (!all(list.class=="SaemixObject")) {
        stop("All inputs should have class 'SaemixObject'.")
      }
    }
  }
  
  nb.mod <- length(mod.list)
  
  
  args1<-match.call(expand.dots=TRUE)
  i1<-match("method",names(args1))
  
  if(!is.na(i1)) {
    str1<-as.character(args1[[i1]])
    if(str1 %in% c("is","lin","gq")) {
      method<-str1
    } else method<-"is"
  } else method<-"is"
  
  
  for (k in 1:nb.mod){
    if(method=="is" & length(mod.list[[k]]@results@ll.is)==0) {
      namObj<-deparse(substitute(mod.list[[k]]))
      mod.list[[k]]<-llis.saemix(mod.list[[k]])
      assign(namObj,mod.list[[k]],envir=parent.frame())
    }
  }
  
  for (k in 1:nb.mod){
    if(method=="gq" & length(mod.list[[k]]@results@ll.gq)==0) {
      namObj<-deparse(substitute(mod.list[[k]]))
      mod.list[[k]]<-llgq.saemix(mod.list[[k]])
      assign(namObj,mod.list[[k]],envir=parent.frame())
    }    
  }
  
  
  for (k in 1:nb.mod){
    if(method=="lin" & length(mod.list[[k]]@results@ll.lin)==0) {
      namObj<-deparse(substitute(mod.list[[k]]))
      mod.list[[k]]<-fim.saemix(mod.list[[k]])
      assign(namObj,mod.list[[k]],envir=parent.frame())
    }
  }
  
  
  if (method=="is"){
    
    info <- matrix(0,nb.mod,2)
    rownames(info) <- as.character(seq(1,nb.mod))
    colnames(info) <- c("AIC","BIC")
    
    str.model <- list()
    cov.model <- list()
    
    for (k in 1:nb.mod){
      info[k,] <- c(mod.list[[k]]@results@aic.lin,mod.list[[k]]@results@bic.lin)
      str.model[[k]] <- deparse(mod.list[[k]]@model@model)
      cov.model[[k]] <- deparse(mod.list[[k]]@model@covariance.model)
    }
    
    same.str.model <- rep(NA,(nb.mod-1))
    same.cov.model <- rep(NA,(nb.mod-1))
    
    for (k in 2:nb.mod){
      same.str.model[k] <- all(str.model[[1]]==str.model[[k]])
      same.cov.model[k] <- all(cov.model[[1]]==cov.model[[k]])
    }
    
    # Si même modèle structurel et même structure d'effets aléatoires
    if ((!("FALSE" %in% same.str.model)) && (!("FALSE" %in% same.cov.model))){
      bic.cov <- rep(NA,nb.mod)
      for (k in 1:nb.mod){
        bic.cov[k] <- mod.list[[k]]@results@bic.covariate.is
      }
      info <- cbind(info, bic.cov)
      colnames(info) <- c("AIC","BIC","BIC.cov")
    }
    cat("Likelihoods computed by importance sampling \n")
  }
  
  if (method=="lin"){
    
    info <- matrix(0,nb.mod,2)
    rownames(info) <- as.character(seq(1,nb.mod))
    colnames(info) <- c("AIC","BIC")
    
    str.model <- list()
    cov.model <- list()
    
    for (k in 1:nb.mod){
      info[k,] <- c(mod.list[[k]]@results@aic.lin,mod.list[[k]]@results@bic.lin)
      str.model[[k]] <- deparse(mod.list[[k]]@model@model)
      cov.model[[k]] <- deparse(mod.list[[k]]@model@covariance.model)
    }
    
    same.str.model <- rep(NA,(nb.mod-1))
    same.cov.model <- rep(NA,(nb.mod-1))
    
    for (k in 2:nb.mod){
      same.str.model[k] <- all(str.model[[1]]==str.model[[k]])
      same.cov.model[k] <- all(cov.model[[1]]==cov.model[[k]])
    }
    
    # Si même modèle structurel et même structure d'effets aléatoires
    if ((!("FALSE" %in% same.str.model)) && (!("FALSE" %in% same.cov.model))){
      bic.cov <- rep(NA,nb.mod)
      for (k in 1:nb.mod){
        bic.cov[k] <- mod.list[[k]]@results@bic.covariate.lin
      }
      info <- cbind(info, bic.cov)
      colnames(info) <- c("AIC","BIC","BIC.cov")
    }
    cat("Likelihoods computed by linearisation \n")
  }
  
  
  if (method=="gq"){
    
    info <- matrix(0,nb.mod,2)
    rownames(info) <- as.character(seq(1,nb.mod))
    colnames(info) <- c("AIC","BIC")
    
    str.model <- list()
    cov.model <- list()
    
    for (k in 1:nb.mod){
      info[k,] <- c(mod.list[[k]]@results@aic.lin,mod.list[[k]]@results@bic.lin)
      str.model[[k]] <- deparse(mod.list[[k]]@model@model)
      cov.model[[k]] <- deparse(mod.list[[k]]@model@covariance.model)
    }
    
    same.str.model <- rep(NA,(nb.mod-1))
    same.cov.model <- rep(NA,(nb.mod-1))
    
    for (k in 2:nb.mod){
      same.str.model[k] <- all(str.model[[1]]==str.model[[k]])
      same.cov.model[k] <- all(cov.model[[1]]==cov.model[[k]])
    }
    
    # Si même modèle structurel et même structure d'effets aléatoires
    if ((!("FALSE" %in% same.str.model)) && (!("FALSE" %in% same.cov.model))){
      bic.cov <- rep(NA,nb.mod)
      for (k in 1:nb.mod){
        bic.cov[k] <- mod.list[[k]]@results@bic.covariate.gq
      }
      info <- cbind(info, bic.cov)
      colnames(info) <- c("AIC","BIC","BIC.cov")
    }
    
    cat("Likelihoods computed by Gaussian quadrature \n")
  }
  
  return(info)
}