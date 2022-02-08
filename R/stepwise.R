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
#' @param covariate.init A matrix specifying the initial covariate structure to be considered in the algorithm.
#' @return An object of the SaemixObject class storing the covariate model and the covariance structure of 
#' random effects of the final model.
#' @author Maud Delattre
#' @references M Delattre, M Lavielle, MA Poursat (2014) A note on BIC in mixed effects models. 
#' Electronic Journal of Statistics 8(1) p. 456-475
#' M Delattre, MA Poursat (2017) BIC strategies for model choice in a population approach. 
#' (arXiv:1612.02405)
#' @keywords selection stepwise covariate

stepwise.procedure <-
  function(saemixObject,
           covariate.init = NULL,
           trace = TRUE) {
    warning("Likelihood computed by importance sampling.")
    
    method <- "is"
    
    saemix.data <- saemixObject["data"]
    saemix.data["messages"] <- FALSE
    saemix.model <- saemixObject["model"]
    
    
    saemix.options <- saemixObject["options"]
    saemix.options$displayProgress <- FALSE
    saemix.options$map <- FALSE
    
    
    saemix.options$ll.is <- TRUE
    saemix.options$ll.gq <- FALSE
    saemix.options$fim <- FALSE
    
    
    nb.psi <- saemix.model["nb.parameters"]
    name.psi <- saemix.model["name.modpar"]
    psi.init <- saemix.model["psi0"]
    transform.par <- saemix.model["transform.par"]
    fixed.estim <- saemix.model["fixed.estim"]
    
    
    
    nb.covariates <- length(saemix.data["name.covariates"])
    
    if (is.null(covariate.init)) {
      covariate.init <- matrix(0, ncol = nb.psi, nrow = nb.covariates)
      covariate.model <- covariate.init
      colnames(covariate.init) <- name.psi
      rownames(covariate.init) <- saemix.data["name.covariates"]
    }
    else{
      dim.cov.init <- dim(covariate.init)
      if ((dim.cov.init[1] != length(saemix.data@name.covariates)) ||
          (dim.cov.init[2] != nb.psi)) {
        stop("Incorrect initialisation for covariates.")
      } else{
        covariate.model <- covariate.init
        colnames(covariate.init) <- name.psi
        rownames(covariate.init) <- saemix.data["name.covariates"]
      }
    }
    
    
    
    MM <-
      matrix(0, saemix.model@nb.parameters, saemix.model@nb.parameters)
    cov.str <- list()
    
    jj <- 0
    for (nb.re in 1:saemixObject@model@nb.parameters) {
      covv <- combn(seq(1, saemixObject@model@nb.parameters), nb.re)
      for (j in 1:dim(covv)[2]) {
        cov.str[[jj + j]] <- MM
        diag(cov.str[[jj + j]])[covv[, j]] <- 1
      }
      jj <- jj + dim(covv)[2]
    }
    
    
    nb.cov.str <- length(cov.str)
    
    # Initializing the covariate structure of the model
    
    
    if (saemixObject@model@modeltype == "structural") {
      ## 1 Continuous models
      
      error.model <- saemix.model["error.model"]
      
      BIC <-
        rep(0, nb.cov.str) # Values of the appropriate BIC at each step of the procedure (either the one for the
      # covariates or the one for the random effects)
      BIC.tot <-
        rep(0, nb.cov.str) # Values of the BIC for the global models (covariates and random effects) at each
      # step of the procedure
      BIC.step <- list()
      BIC.tot.step <- list()
      model.step <- list()
      best.model <- list()
      best.model[[1]] <- covariate.model
      
      
      ind.step <- 1
      
      for (j in 1:nb.cov.str) {
        model <- saemixModel(
          model = saemix.model@model,
          modeltype = saemix.model@modeltype,
          description = "",
          psi0 = psi.init,
          transform.par = transform.par,
          covariate.model = covariate.init,
          fixed.estim = fixed.estim,
          covariance.model = cov.str[[j]],
          error.model = error.model
        )
        res <- saemix(model, saemix.data, saemix.options)
        BIC[j] <- BIC(res)
        BIC.tot[j] <- BIC(res)
      }
      
      best.cov.str <- cov.str[[which.min(BIC)[1]]]
      
      best.BIC <- BIC.tot[which.min(BIC)[1]]
      
      BIC.step[[ind.step]] <- BIC
      BIC.tot.step[[ind.step]] <- BIC.tot
      model.step[[ind.step]] <- seq(1, nb.cov.str)
      best.model[[ind.step + 1]] <- best.cov.str
      
      cov.model <- list()
      BIC <- rep(0, nb.psi * nb.covariates)
      BIC.tot <- rep(0, nb.psi * nb.covariates)
      k <- 0
      ind.step <- ind.step + 1
      
      for (j in 1:(nb.psi * nb.covariates)) {
        k <- k + 1
        covariate.new.model <- covariate.model
        covariate.new.model[j] <- 1 - covariate.model[j]
        cov.model[[k]] <- covariate.new.model
        model <- saemixModel(
          model = saemix.model@model,
          modeltype = saemix.model@modeltype,
          description = "",
          psi0 = psi.init,
          transform.par = transform.par,
          covariate.model = covariate.new.model,
          fixed.estim = fixed.estim,
          covariance.model = best.cov.str,
          error.model = error.model
        )
        
        res <- saemix(model, saemix.data, saemix.options)
        BIC[j] <- BIC.covariate(res)
        BIC.tot[j] <- BIC(res)
      }
      
      min.BIC <- min(BIC)
      indice.best <- which.min(BIC)[1]
      select.BIC.tot <- BIC.tot[indice.best]
      
      
      if (select.BIC.tot < best.BIC) {
        best.covariate.model <- cov.model[[indice.best]]
      } else{
        best.covariate.model <-
          matrix(c(rep(0, nb.psi * nb.covariates)), ncol = nb.psi, byrow = TRUE)
      }
      
      BIC.step[[ind.step]] <- BIC
      BIC.tot.step[[ind.step]] <- BIC.tot
      model.list <- rep('', length(cov.model))
      best.model[[ind.step + 1]] <- best.covariate.model
      
      
      for (m in 1:length(cov.model)) {
        for (p in 1:length(name.psi)) {
          cov.list <- c('')
          nb.cov.psi <-
            length(saemix.data@name.covariates[which(cov.model[[m]][, p] == 1)])
          if (nb.cov.psi != 0) {
            for (q in 1:nb.cov.psi) {
              cov.list <-
                paste(cov.list, saemix.data@name.covariates[which(cov.model[[m]][, p] ==
                                                                    1)][q], '')
            }
            model.list[m] <-
              paste(model.list[m], name.psi[p], '(', cov.list, ')')
          }
          
          if (model.list[m] == ""){model.list[m]<-"<none>"}
          
        }
      }
      
      model.list <- rep('', length(cov.model))
      
      for (m in 1:length(cov.model)) {
        for (p in 1:length(name.psi)) {
          cov.list <- c('')
          nb.cov.psi <-
            length(saemix.data@name.covariates[which(cov.model[[m]][, p] == 1)])
          if (nb.cov.psi != 0) {
            for (q in 1:nb.cov.psi) {
              cov.list <-
                paste(cov.list, saemix.data@name.covariates[which(cov.model[[m]][, p] ==
                                                                    1)][q], '')
            }
            model.list[m] <-
              paste(model.list[m],
                    name.psi[p],
                    '(',
                    substring(cov.list, 2, nchar(cov.list) - 1),
                    ')',
                    sep = "")
          } 
          
          if (model.list[m] == ""){model.list[m]<-"<none>"}
        }
      }
      model.step[[ind.step]] <- model.list
      
      while (select.BIC.tot <= best.BIC) {
        best.BIC <- select.BIC.tot
        
        BIC <- rep(0, nb.cov.str)
        BIC.tot <- rep(0, nb.cov.str)
        ind.step <- ind.step + 1
        
        for (j in 1:nb.cov.str) {
          model <- saemixModel(
            model = saemix.model@model,
            modeltype = saemix.model@modeltype,
            description = "",
            psi0 = psi.init,
            transform.par = transform.par,
            covariate.model = best.covariate.model,
            fixed.estim = fixed.estim,
            covariance.model = cov.str[[j]],
            error.model = error.model
          )
          
          res <- saemix(model, saemix.data, saemix.options)
          BIC[j] <- BIC(res)
          BIC.tot[j] <- BIC(res)
        }
        
        
        min.BIC <- min(BIC)
        indice.best <- which(BIC == min.BIC)[1]
        best.cov.str <- cov.str[[indice.best]]
        best.BIC <- BIC.tot[indice.best]
        BIC.step[[ind.step]] <- BIC
        BIC.tot.step[[ind.step]] <- BIC.tot
        model.step[[ind.step]] <- seq(1, nb.cov.str)
        best.model[[ind.step + 1]] <- best.cov.str
        
        
        cov.model <- list()
        BIC <- rep(0, nb.psi * nb.covariates)
        BIC.tot <- rep(0, nb.psi * nb.covariates)
        k <- 0
        ind.step <- ind.step + 1
        
        for (j in 1:(nb.psi * nb.covariates)) {
          k <- k + 1
          covariate.new.model <- best.covariate.model
          covariate.new.model[j] <- 1 - covariate.new.model[j]
          cov.model[[j]] <- covariate.new.model
          model <- saemixModel(
            model = saemix.model@model,
            modeltype = saemix.model@modeltype,
            description = "",
            psi0 = psi.init,
            transform.par = transform.par,
            covariate.model = covariate.new.model,
            fixed.estim = fixed.estim,
            covariance.model = best.cov.str,
            error.model = error.model
          )
          
          res <- saemix(model, saemix.data, saemix.options)
          BIC[j] <- BIC.covariate(res)
          BIC.tot[j] <- BIC(res)
        }
        
        BIC.non.null <- BIC[which(BIC != 0)]
        min.BIC <- min(BIC.non.null)
        BIC.tot.non.null <- BIC.tot[which(BIC.tot != 0)]
        select.BIC.tot <- BIC.tot.non.null[which.min(BIC.non.null)]
        
        if (select.BIC.tot <= best.BIC) {
          best.BIC <- select.BIC.tot
          best.covariate.model <- cov.model[[which.min(BIC.non.null)]]
        }
        BIC.step[[ind.step]] <- BIC.non.null
        BIC.tot.step[[ind.step]] <- BIC.tot.non.null
        model.list <- rep('', length(cov.model))
        
        for (m in 1:length(cov.model)) {
          for (p in 1:length(name.psi)) {
            cov.list <- c('')
            nb.cov.psi <-
              length(saemix.data@name.covariates[which(cov.model[[m]][, p] == 1)])
            if (nb.cov.psi != 0){
              if (nb.cov.psi > 1) {
                for (q in 1:(nb.cov.psi - 1)) {
                  cov.list <-
                    paste(cov.list,
                          saemix.data@name.covariates[which(cov.model[[m]][, p] == 1)][q],
                          ',',
                          sep = "")
                }
                cov.list <-
                  paste(cov.list, saemix.data@name.covariates[which(cov.model[[m]][, p] ==
                                                                      1)][nb.cov.psi], sep = "")
              }
              else {
                cov.list <-
                  paste(cov.list, saemix.data@name.covariates[which(cov.model[[m]][, p] ==
                                                                      1)][1], sep = "")
              }
              model.list[m] <-
                paste(model.list[m], name.psi[p], '(', cov.list, ')', sep = "")
            } 
            
            if (model.list[m] == ""){model.list[m]<-"<none>"}
          }
        }
        model.step[[ind.step]] <- model.list
        best.model[[ind.step + 1]] <- best.covariate.model
        
      }
      
      
      
      best.fit <- saemixModel(
        model = saemixObject@model@model,
        modeltype = saemixObject@model@modeltype,
        description = saemixObject@model@description,
        psi0 = saemixObject@model@psi0,
        fixed.estim = saemixObject@model@fixed.estim,
        covariate.model = best.model[[length(best.model)]],
        covariance.model = best.model[[length(best.model) -
                                         1]],
        error.model = saemixObject@model@error.model
      )
      
      nb.steps = length(BIC.step)
      
      # Print main steps of the procedure on the console
      if (trace == TRUE) {
        cat('\n')
        cat('\n')
        cat('------------------\n')
        cat('---- Summary  ----\n')
        cat('------------------\n')
        res.summary <- matrix(NA, nb.steps, 3)
        colnames(res.summary) <- c("Covariates", "R.E.", "BIC")
        
        name.init.cov <- c('')
        for (p in 1:length(name.psi)) {
          nb.cov.psi <-
            length(saemix.data@name.covariates[which(covariate.init[, p] == 1)])
          init.cov <- c('')
          if (nb.cov.psi != 0) {
            for (q in 1:nb.cov.psi) {
              init.cov <-
                paste(init.cov, saemix.data@name.covariates[which(covariate.init[, p] ==
                                                                    1)][q], '')
            }
            name.init.cov <-
              paste(name.init.cov,
                    name.psi[p],
                    '(',
                    substring(init.cov, 2, nchar(init.cov) - 1),
                    ')',
                    sep = "")
          }
        }
        
        if (name.init.cov == "") {
          name.init.cov = "<none>"
        }
        
        re <-
          as.character(name.psi[which(diag(cov.str[[which.min(BIC.step[[1]])]]) ==
                                        1)])
        re <- paste(re, collapse = ",")
        res.summary[1, ] <-
          c(name.init.cov, re, as.character(round(BIC.tot.step[[1]][which.min(BIC.step[[1]])], 2)))
        for (j in 2:nb.steps) {
          if ((j / 2 - floor(j / 2)) != 0) {
            re <-
              as.character(name.psi[which(diag(cov.str[[which.min(BIC.step[[j]])]]) ==
                                            1)])
            re <- paste(re, collapse = ",")
            res.summary[j, ] <-
              c("----", re, as.character(round(BIC.tot.step[[j]][which.min(BIC.step[[j]])], 2)))
          } else{
            res.summary[j, ] <-
              c(model.step[[j]][which.min(BIC.step[[j]])], "----", as.character(round(BIC.tot.step[[j]][which.min(BIC.step[[j]])], 2)))
          }
        }
        
        rownames(res.summary) <- rep("", nb.steps)
        
        print(res.summary, quote = FALSE)
        
        cat('\n')
        cat('\n')
        cat('---------------------\n')
        cat('---- Final Model ----\n')
        cat('---------------------\n')
        cat('Covariate model \n')
        colnames(best.model[[nb.steps + 1]]) <- name.psi
        rownames(best.model[[nb.steps + 1]]) <-
          saemix.data["name.covariates"]
        print(best.model[[nb.steps + 1]])
        cat('Random effects structure \n')
        mat <-  best.model[[nb.steps]]
        colnames(mat) <- rownames(mat) <- name.psi
        print(mat)
        cat('BIC=', (BIC.tot.step[[nb.steps - 1]][which.min(BIC.step[[nb.steps -
                                                                        1]])]))
      }
    }
    else {
      ## 2 Discrete models
      
      
      
      BIC <-
        rep(0, nb.cov.str) # Values of the appropriate BIC at each step of the procedure (either the one for the
      # covariates or the one for the random effects)
      BIC.tot <-
        rep(0, nb.cov.str) # Values of the BIC for the global models (covariates and random effects) at each
      # step of the procedure
      BIC.step <- list()
      BIC.tot.step <- list()
      model.step <- list()
      best.model <- list()
      best.model[[1]] <- covariate.model
      
      
      ind.step <- 1
      
      for (j in 1:nb.cov.str) {
        model <- saemixModel(
          model = saemix.model@model,
          modeltype = saemix.model@modeltype,
          description = "",
          psi0 = psi.init,
          transform.par = transform.par,
          covariate.model = covariate.init,
          fixed.estim = fixed.estim,
          covariance.model = cov.str[[j]]
        )
        res <- saemix(model, saemix.data, saemix.options)
        BIC[j] <- BIC(res)
        BIC.tot[j] <- BIC(res)
      }
      
      best.cov.str <- cov.str[[which.min(BIC)[1]]]
      
      best.BIC <- BIC.tot[which.min(BIC)[1]]
      
      BIC.step[[ind.step]] <- BIC
      BIC.tot.step[[ind.step]] <- BIC.tot
      model.step[[ind.step]] <- seq(1, nb.cov.str)
      best.model[[ind.step + 1]] <- best.cov.str
      
      cov.model <- list()
      BIC <- rep(0, nb.psi * nb.covariates)
      BIC.tot <- rep(0, nb.psi * nb.covariates)
      k <- 0
      ind.step <- ind.step + 1
      
      for (j in 1:(nb.psi * nb.covariates)) {
        k <- k + 1
        covariate.new.model <- covariate.model
        covariate.new.model[j] <- 1 - covariate.model[j]
        cov.model[[k]] <- covariate.new.model
        model <- saemixModel(
          model = saemix.model@model,
          modeltype = saemix.model@modeltype,
          description = "",
          psi0 = psi.init,
          transform.par = transform.par,
          covariate.model = covariate.new.model,
          fixed.estim = fixed.estim,
          covariance.model = best.cov.str
        )
        
        res <- saemix(model, saemix.data, saemix.options)
        BIC[j] <- BIC.covariate(res)
        BIC.tot[j] <- BIC(res)
      }
      
      min.BIC <- min(BIC)
      indice.best <- which.min(BIC)[1]
      select.BIC.tot <- BIC.tot[indice.best]
      
      
      if (select.BIC.tot < best.BIC) {
        best.covariate.model <- cov.model[[indice.best]]
      } else{
        best.covariate.model <-
          matrix(c(rep(0, nb.psi * nb.covariates)), ncol = nb.psi, byrow = TRUE)
      }
      
      BIC.step[[ind.step]] <- BIC
      BIC.tot.step[[ind.step]] <- BIC.tot
      model.list <- rep('', length(cov.model))
      best.model[[ind.step + 1]] <- best.covariate.model
      
      
      for (m in 1:length(cov.model)) {
        for (p in 1:length(name.psi)) {
          cov.list <- c('')
          nb.cov.psi <-
            length(saemix.data@name.covariates[which(cov.model[[m]][, p] == 1)])
          if (nb.cov.psi != 0) {
            for (q in 1:nb.cov.psi) {
              cov.list <-
                paste(cov.list, saemix.data@name.covariates[which(cov.model[[m]][, p] ==
                                                                    1)][q], '')
            }
            model.list[m] <-
              paste(model.list[m], name.psi[p], '(', cov.list, ')')
          } 
          
          if (model.list[m] == ""){model.list[m]<-"<none>"}
          
        }
      }
      
      model.list <- rep('', length(cov.model))
      
      for (m in 1:length(cov.model)) {
        for (p in 1:length(name.psi)) {
          cov.list <- c('')
          nb.cov.psi <-
            length(saemix.data@name.covariates[which(cov.model[[m]][, p] == 1)])
          if (nb.cov.psi != 0) {
            for (q in 1:nb.cov.psi) {
              cov.list <-
                paste(cov.list, saemix.data@name.covariates[which(cov.model[[m]][, p] ==
                                                                    1)][q], '')
            }
            model.list[m] <-
              paste(model.list[m],
                    name.psi[p],
                    '(',
                    substring(cov.list, 2, nchar(cov.list) - 1),
                    ')',
                    sep = "")
          }
        }
      }
      model.step[[ind.step]] <- model.list
      
      while (select.BIC.tot <= best.BIC) {
        best.BIC <- select.BIC.tot
        
        BIC <- rep(0, nb.cov.str)
        BIC.tot <- rep(0, nb.cov.str)
        ind.step <- ind.step + 1
        
        for (j in 1:nb.cov.str) {
          model <- saemixModel(
            model = saemix.model@model,
            modeltype = saemix.model@modeltype,
            description = "",
            psi0 = psi.init,
            transform.par = transform.par,
            covariate.model = best.covariate.model,
            fixed.estim = fixed.estim,
            covariance.model = cov.str[[j]]
          )
          
          res <- saemix(model, saemix.data, saemix.options)
          BIC[j] <- BIC(res)
          BIC.tot[j] <- BIC(res)
        }
        
        
        min.BIC <- min(BIC)
        indice.best <- which(BIC == min.BIC)[1]
        best.cov.str <- cov.str[[indice.best]]
        best.BIC <- BIC.tot[indice.best]
        BIC.step[[ind.step]] <- BIC
        BIC.tot.step[[ind.step]] <- BIC.tot
        model.step[[ind.step]] <- seq(1, nb.cov.str)
        best.model[[ind.step + 1]] <- best.cov.str
        
        
        cov.model <- list()
        BIC <- rep(0, nb.psi * nb.covariates)
        BIC.tot <- rep(0, nb.psi * nb.covariates)
        k <- 0
        ind.step <- ind.step + 1
        
        for (j in 1:(nb.psi * nb.covariates)) {
          k <- k + 1
          covariate.new.model <- best.covariate.model
          covariate.new.model[j] <- 1 - covariate.new.model[j]
          cov.model[[j]] <- covariate.new.model
          model <- saemixModel(
            model = saemix.model@model,
            modeltype = saemix.model@modeltype,
            description = "",
            psi0 = psi.init,
            transform.par = transform.par,
            covariate.model = covariate.new.model,
            fixed.estim = fixed.estim,
            covariance.model = best.cov.str
          )
          
          res <- saemix(model, saemix.data, saemix.options)
          BIC[j] <- BIC.covariate(res)
          BIC.tot[j] <- BIC(res)
        }
        
        BIC.non.null <- BIC[which(BIC != 0)]
        min.BIC <- min(BIC.non.null)
        BIC.tot.non.null <- BIC.tot[which(BIC.tot != 0)]
        select.BIC.tot <- BIC.tot.non.null[which.min(BIC.non.null)]
        
        if (select.BIC.tot <= best.BIC) {
          best.BIC <- select.BIC.tot
          best.covariate.model <- cov.model[[which.min(BIC.non.null)]]
        }
        BIC.step[[ind.step]] <- BIC.non.null
        BIC.tot.step[[ind.step]] <- BIC.tot.non.null
        model.list <- rep('', length(cov.model))
        
        for (m in 1:length(cov.model)) {
          for (p in 1:length(name.psi)) {
            cov.list <- c('')
            nb.cov.psi <-
              length(saemix.data@name.covariates[which(cov.model[[m]][, p] == 1)])
            if (nb.cov.psi != 0){
              if (nb.cov.psi > 1) {
                for (q in 1:(nb.cov.psi - 1)) {
                  cov.list <-
                    paste(cov.list,
                          saemix.data@name.covariates[which(cov.model[[m]][, p] == 1)][q],
                          ',',
                          sep = "")
                }
                cov.list <-
                  paste(cov.list, saemix.data@name.covariates[which(cov.model[[m]][, p] ==
                                                                      1)][nb.cov.psi], sep = "")
              }
              else {
                cov.list <-
                  paste(cov.list, saemix.data@name.covariates[which(cov.model[[m]][, p] ==
                                                                      1)][1], sep = "")
              }
              model.list[m] <-
                paste(model.list[m], name.psi[p], '(', cov.list, ')', sep = "")
            } 
            
            if (model.list[m] == ""){model.list[m]<-"<none>"}
          }
        }
        model.step[[ind.step]] <- model.list
        best.model[[ind.step + 1]] <- best.covariate.model
        
      }
      
      
      
      best.fit <- saemixModel(
        model = saemixObject@model@model,
        modeltype = saemixObject@model@modeltype,
        description = saemixObject@model@description,
        psi0 = saemixObject@model@psi0,
        fixed.estim = saemixObject@model@fixed.estim,
        covariate.model = best.model[[length(best.model)]],
        covariance.model = best.model[[length(best.model) -
                                         1]]
      )
      
      nb.steps = length(BIC.step)
      
      # Print main steps of the procedure on the console
      if (trace == TRUE) {
        cat('\n')
        cat('\n')
        cat('------------------\n')
        cat('---- Summary  ----\n')
        cat('------------------\n')
        res.summary <- matrix(NA, nb.steps, 3)
        colnames(res.summary) <- c("Covariates", "R.E.", "BIC")
        
        name.init.cov <- c('')
        for (p in 1:length(name.psi)) {
          nb.cov.psi <-
            length(saemix.data@name.covariates[which(covariate.init[, p] == 1)])
          init.cov <- c('')
          if (nb.cov.psi != 0) {
            for (q in 1:nb.cov.psi) {
              init.cov <-
                paste(init.cov, saemix.data@name.covariates[which(covariate.init[, p] ==
                                                                    1)][q], '')
            }
            name.init.cov <-
              paste(name.init.cov,
                    name.psi[p],
                    '(',
                    substring(init.cov, 2, nchar(init.cov) - 1),
                    ')',
                    sep = "")
          }
        }
        
        if (name.init.cov == "") {
          name.init.cov = "<none>"
        }
        
        re <-
          as.character(name.psi[which(diag(cov.str[[which.min(BIC.step[[1]])]]) ==
                                        1)])
        re <- paste(re, collapse = ",")
        res.summary[1, ] <-
          c(name.init.cov, re, as.character(round(BIC.tot.step[[1]][which.min(BIC.step[[1]])], 2)))
        for (j in 2:nb.steps) {
          if ((j / 2 - floor(j / 2)) != 0) {
            re <-
              as.character(name.psi[which(diag(cov.str[[which.min(BIC.step[[j]])]]) ==
                                            1)])
            re <- paste(re, collapse = ",")
            res.summary[j, ] <-
              c("----", re, as.character(round(BIC.tot.step[[j]][which.min(BIC.step[[j]])], 2)))
          } else{
            res.summary[j, ] <-
              c(model.step[[j]][which.min(BIC.step[[j]])], "----", as.character(round(BIC.tot.step[[j]][which.min(BIC.step[[j]])], 2)))
          }
        }
        
        rownames(res.summary) <- rep("", nb.steps)
        
        print(res.summary, quote = FALSE)
        
        cat('\n')
        cat('\n')
        cat('---------------------\n')
        cat('---- Final Model ----\n')
        cat('---------------------\n')
        cat('Covariate model \n')
        colnames(best.model[[nb.steps + 1]]) <- name.psi
        rownames(best.model[[nb.steps + 1]]) <-
          saemix.data["name.covariates"]
        print(best.model[[nb.steps + 1]])
        cat('Random effects structure \n')
        mat <-  best.model[[nb.steps]]
        colnames(mat) <- rownames(mat) <- name.psi
        print(mat)
        cat('BIC=', (BIC.tot.step[[nb.steps - 1]][which.min(BIC.step[[nb.steps -
                                                                        1]])]))
      }
    }
    
    return(best.fit)
    
  }