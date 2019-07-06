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
    
    # Comparison of datasets : if models are estimated on different datasets, comparison does 
    # not make sense. 
    
    same.data <- rep(NA,(nb.mod-1))
    
    for (k in 2:nb.mod){
      same.data[k-1] <- identical(mod.list[[1]]@data,mod.list[[k]]@data)
    }
    
    
    if ("FALSE" %in% same.data){
      stop('Compared models should be fitted on the same data.')
    } else{
      # Comparison of the statistical models : BIC.covariate does only make sense if model type,
      # structural model (and residual model if appropriate), and covariance structure for the random 
      # effects are the same
      
      for (k  in 1:nb.mod){
        info[k,] <- c(mod.list[[k]]@results@aic.is,mod.list[[k]]@results@bic.is)
      }
      
      same.model.type <- rep(NA,(nb.mod-1))
      same.str.model <- rep(NA,(nb.mod-1))
      
      for (k in 2:nb.mod){
        same.model.type[k-1] <- identical(mod.list[[1]]@model@modeltype,mod.list[[k]]@model@modeltype)
        same.str.model[k-1] <- identical(mod.list[[1]]@model@model,mod.list[[k]]@model@model)
      }
      
      if (!("FALSE" %in% same.model.type)){
        
        if (!("FALSE" %in% same.str.model)){
          
          same.cov.model <- rep(NA,(nb.mod-1))
          
          for (k in 2:nb.mod){
            same.cov.model[k-1] <- identical(mod.list[[1]]@model@covariance.model,mod.list[[k]]@model@covariance.model)
          }
          
          if (mod.list[[1]]@model@modeltype == "structural"){
            
            same.res.model <- rep(NA,(nb.mod-1))
            
            
            
            for (k in 2:nb.mod){
              same.res.model[k-1] <- identical(mod.list[[1]]@model@error.model,mod.list[[k]]@model@error.model)
            }
            
            if ((!("FALSE" %in% same.cov.model)) && (!("FALSE" %in% same.res.model))){
              bic.cov <- rep(NA,nb.mod)
              for (k in 1:nb.mod){
                bic.cov[k] <- mod.list[[k]]@results@bic.covariate.is
              }
              info <- cbind(info, bic.cov)
              colnames(info) <- c("AIC","BIC","BIC.cov")
              
            }
            
          } else {
            if (!("FALSE" %in% same.cov.model)){
              bic.cov <- rep(NA,nb.mod)
              for (k in 1:nb.mod){
                bic.cov[k] <- mod.list[[k]]@results@bic.covariate.is
              }
              info <- cbind(info, bic.cov)
              colnames(info) <- c("AIC","BIC","BIC.cov")
            }
          }
        }
        
      }
      cat("Likelihoods computed by importance sampling \n")
    }
  }
  
  if (method=="lin"){
    
    info <- matrix(0,nb.mod,2)
    rownames(info) <- as.character(seq(1,nb.mod))
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
    } else{
      # Comparison of the statistical models : BIC.covariate does only make sense if model type,
      # structural model (and) residual model if appropriate), and covariance structure for the random 
      # effects are the same
      
      for (k  in 1:nb.mod){
        info[k,] <- c(mod.list[[k]]@results@aic.lin,mod.list[[k]]@results@bic.lin)
      }
      
      same.model.type <- rep(NA,(nb.mod-1))
      same.str.model <- rep(NA,(nb.mod-1))
      
      for (k in 2:nb.mod){
        same.model.type[k-1] <- identical(mod.list[[1]]@model@modeltype,mod.list[[k]]@model@modeltype)
        same.str.model[k-1] <- identical(mod.list[[1]]@model@model,mod.list[[k]]@model@model)
      }
      
      if (!("FALSE" %in% same.model.type)){
        
        if (!("FALSE" %in% same.str.model)){
          
          same.cov.model <- rep(NA,(nb.mod-1))
          
          for (k in 2:nb.mod){
            same.cov.model[k-1] <- identical(mod.list[[1]]@model@covariance.model,mod.list[[k]]@model@covariance.model)
          }
          
          if (mod.list[[1]]@model@modeltype == "structural"){
            
            same.res.model <- rep(NA,(nb.mod-1))
            
            
            
            for (k in 2:nb.mod){
              same.res.model[k-1] <- identical(mod.list[[1]]@model@error.model,mod.list[[k]]@model@error.model)
            }
            
            if ((!("FALSE" %in% same.cov.model)) && (!("FALSE" %in% same.res.model))){
              bic.cov <- rep(NA,nb.mod)
              for (k in 1:nb.mod){
                bic.cov[k] <- mod.list[[k]]@results@bic.covariate.lin
              }
              info <- cbind(info, bic.cov)
              colnames(info) <- c("AIC","BIC","BIC.cov")
              
            }
            
          } else {
            warning('Linearisation is not appropriate for computing likelihoods in discrete models.')
            if (!("FALSE" %in% same.cov.model)){
              bic.cov <- rep(NA,nb.mod)
              for (k in 1:nb.mod){
                bic.cov[k] <- mod.list[[k]]@results@bic.covariate.lin
              }
              info <- cbind(info, bic.cov)
              colnames(info) <- c("AIC","BIC","BIC.cov")
            }
          }
        }
        
      }
      
    
      
      cat("Likelihoods computed by linearisation \n")
    }
  }
  
  if (method=="gq"){
    
    info <- matrix(0,nb.mod,2)
    rownames(info) <- as.character(seq(1,nb.mod))
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
    } else{
      # Comparison of the statistical models : BIC.covariate does only make sense if model type,
      # structural model (and residual model if appropriate), and covariance structure for the random 
      # effects are the same
      
      for (k  in 1:nb.mod){
        info[k,] <- c(mod.list[[k]]@results@aic.gq,mod.list[[k]]@results@bic.gq)
      }
      
      same.model.type <- rep(NA,(nb.mod-1))
      same.str.model <- rep(NA,(nb.mod-1))
      
      for (k in 2:nb.mod){
        same.model.type[k-1] <- identical(mod.list[[1]]@model@modeltype,mod.list[[k]]@model@modeltype)
        same.str.model[k-1] <- identical(mod.list[[1]]@model@model,mod.list[[k]]@model@model)
      }
      
      if (!("FALSE" %in% same.model.type)){
        
        if (!("FALSE" %in% same.str.model)){
          
          same.cov.model <- rep(NA,(nb.mod-1))
          
          for (k in 2:nb.mod){
            same.cov.model[k-1] <- identical(mod.list[[1]]@model@covariance.model,mod.list[[k]]@model@covariance.model)
          }
          
          if (mod.list[[1]]@model@modeltype == "structural"){
            
            same.res.model <- rep(NA,(nb.mod-1))
            
            
            
            for (k in 2:nb.mod){
              same.res.model[k-1] <- identical(mod.list[[1]]@model@error.model,mod.list[[k]]@model@error.model)
            }
            
            if ((!("FALSE" %in% same.cov.model)) && (!("FALSE" %in% same.res.model))){
              bic.cov <- rep(NA,nb.mod)
              for (k in 1:nb.mod){
                bic.cov[k] <- mod.list[[k]]@results@bic.covariate.gq
              }
              info <- cbind(info, bic.cov)
              colnames(info) <- c("AIC","BIC","BIC.cov")
              
            }
            
          } else {
            if (!("FALSE" %in% same.cov.model)){
              bic.cov <- rep(NA,nb.mod)
              for (k in 1:nb.mod){
                bic.cov[k] <- mod.list[[k]]@results@bic.covariate.gq
              }
              info <- cbind(info, bic.cov)
              colnames(info) <- c("AIC","BIC","BIC.cov")
            }
          }
        }
        
      }
      
      cat("Likelihoods computed by Gaussian quadrature \n")
    }
  }
  return(info)
}