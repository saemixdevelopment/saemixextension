Model <- R6Class("Model", 
  portable = FALSE,
  private = list(   
    omega_sqrt = 0,  # cholesky decomposition of omega matrix
    omega_derivatives = list(), # list of omega derivative matricies 
    par_count = 0L,   # number of parameters for the log-likelihood function
    eta_zero = 0,       # zero random effect vector with size re_count
    ll_max = 0,       # largest likelihood value encountered
    min_ll_allowed = 0, # smallest likelihood value allowed (smaller will be truncated)
    fe_index = 0L,    # index vector for the mu vector
    re_index = 0L,    # index vector for the b vector
    par_index = 0L,    # index vector for the full parameter vector,
    adaptive = NULL,
    gq_settings = NULL,
    mc_settings = NULL,
    h = NULL,
    adaptation_mode = NULL,
    lhs_sampling = NULL,
    settings_y_integration = NULL,
    
    inverse_integration_function = NULL,
    
    fim_std_sampling = function(design){
      matrix(integrate_mc(function(eta.vec){
        apply(eta.vec,1,function(eta){
          sim <-simulation_eta_function(design, eta)
          observed_fim(sim, design, eta)
        })
      }, re_count, settings=settings_y_integration), nrow=theta_count)
    },
    
    fim_inverse_sampling = function(design){
      y <- inverse_simulation_eta_function(design, eta_zero, urand=NULL)
      length_y <- length(y)
      matrix(inverse_integration_function(function(urand_matrix){
        apply(urand_matrix, 1, function(urand){
          eta <- qnorm(urand[seq_len(re_count)]) 
          y <- inverse_simulation_eta_function(design, eta, urand=urand[-seq_len(re_count)])
          observed_fim(y, design, eta)
        })
      }, re_count+length_y, settings=settings_y_integration), nrow=theta_count)
    },
    
    calc_fim = NULL
    ),
  public = list(
    mu = 0,         # parameter values
    omega = 0,        # variance of IIV random effects
    theta = 0,   # vector of parameters values
    theta_names = "", # names for the model parameters 
    fe_count = 0L,   # number of fixed effect parameters
    re_count = 0L,   # number of random effect parameters
    theta_count = 0L,
    parameter_function = NULL,
    log_likelihood_function = NULL,
    simulation_function = NULL,
    inverse_simulation_function = NULL, 
    initialize = function(parameter_function, log_likelihood_function, 
                          simulation_function, 
                          inverse_simulation_function=function(...) stop("Inverse simulation function needs to be implemented for this sampling method."),
                          mu, omega=diag(0, nrow=NROW(mu)), settings=defaults.agq()){
      # Creates a new model object  
      # 
      # Args:
      #  parameter_function: function with signature f(mu, b) returning a named list list(p1=,p2=,...)
      #                       with the parameter values given mu and b 
      #  log_likelihood_function: function with signature f(y, design, p1, p2, ...) returning a 
      #                           row vector of length y  with log_likelihood for each y_j
      #  simulation_function: function with signature f(design, p1, p2, ...) returning a row vector of 
      #                       simulated data
      #  mu: row vector of fixed effect values 
      #  omega: covariance matrix for the random effects (automatically transformed to a matrix if 
      #         given as scalar)
      #  debug: if true, the model object exposes all private functions publically otherwise only the ones 
      #         used by external routines
      #
      # Returns:
      #  A model object (implemented as list) 
      parameter_function <<- parameter_function
      log_likelihood_function <<- log_likelihood_function
      simulation_function <<- simulation_function
      inverse_simulation_function <<- inverse_simulation_function
      # convert omega to matrix if a scalar was supplied
      if(!is.matrix(omega)&NROW(omega==1)) omega <- diag(omega, 1)
      
      mu <<- mu
      omega <<- omega
      # create list of matrix derivatives of the omega matrix
      # list of unique entries
      unique_entries <- which(lower.tri(omega, T) & omega!=0.0, arr.ind = T)
      # reorder to have diagnoal elements first
      unique_entries <- unique_entries[order(unique_entries[,1]!=unique_entries[,2], unique_entries[,1], unique_entries[,2]),,drop=F]

      zero_matrix <- omega
      zero_matrix[] <- 0
      omega_derivatives <<- mapply(function(row, col) {zero_matrix[row, col] <- 1
                                                       zero_matrix[col, row] <- 1
                                                       zero_matrix}, unique_entries[,1], unique_entries[,2], SIMPLIFY=F, USE.NAMES=F)
      # vector of parameters valuesv
      theta <<- c(mu, omega[unique_entries])
            
      # number of fixed effect parameters
      fe_count <<- NROW(mu)
      
      # number of random effect parameters
      re_count <<- NROW(omega)
      
      # number of total parameters
      theta_count <<- NROW(theta)

      # number of parameters for the log-likelihood function (#arguments minus y and design)
      par_count <<- length(formals(log_likelihood_function)) - 2 
      
      # cholesky decomposition of omega
      if(re_count>0){
        omega_sqrt <<- t(chol(omega))
      }else{
        omega_sqrt <<- diag(nrow=0)
      }
      
      # zero random effect vector with size re_count
      eta_zero <<- rep(0, re_count)
      
      ll_max <<- 0
      
      min_ll_allowed <<- log(.Machine$double.xmin)/2 
      
      
      # vector 1:fe_count
      fe_index <<- seq_len(fe_count)
      
      # vector 1:re_count
      re_index <<- seq_len(re_count)
      
      # vector 1:par_count
      par_index <<- seq_len(par_count)
      
      # names for the fixed effects
      names_mu <- names(mu)
      if(is.null(names_mu)){
        fe_names  <- paste0("mu", fe_index) 
      }else{
        fe_names  <- c(paste0("mu", fe_index)[names_mu==""], names_mu[names_mu!=""]) 
      }
      
      theta_names <<- c(fe_names, paste0("omega",unique_entries[,1], ".", unique_entries[,2]))
      
      update_settings(settings)
    },
    
    
    # log likelihood function of eta
    log_likelihood_eta_function = function(y, design, eta){
      do.call(log_likelihood_function, c(list(y=y, design=design),
                                         parameter_function(mu, omega_sqrt %*% eta)))    
    },
    
    # simulation_function of eta
    simulation_eta_function = function(design, eta){
      do.call(simulation_function, c(list(design=design),
                                     parameter_function(mu, omega_sqrt %*% eta)))
    },
    
    inverse_simulation_eta_function = function(design, eta, urand){
      do.call(inverse_simulation_function, c(list(design=design, urand=urand),
                                             parameter_function(mu, omega_sqrt %*% eta)))
    },
    
    dpardmu = function(b){
      # Calculates the Jacbian of the parameter function wrt mu 
      # 
      # Args:
      #  b:  parameter vector of random effects to take the derivative at
      #
      # Returns:
      #  Jacobian (n.parXn.mu) wrt mu 
      par0 <- unlist(parameter_function(mu, b))
      vapply(seq_along(mu), function(i) (unlist(parameter_function(mu+h*(fe_index==i), b))-par0), 
             rep(0,NROW(par0)))/h
    },
    
    dpardb = function(b){
      # Calculates the Jacbian of the parameter function wrt b 
      # 
      # Args:
      #  b:  parameter vector of random effects to take the derivative at
      #
      # Returns:
      #  Jacobian (n.parXn.b) wrt b 
      par0 <- unlist(parameter_function(mu, b))
      vapply(seq_along(b), 
             function(i) (unlist(parameter_function(mu, b+h*(re_index==i)))-par0), 
             rep(0,NROW(par0)))/h
    },
    
    
    gradient_dll_dpar = function(y, design, ...){
      # Calculates the gradient of the log likelihood function wrt to its parameters
      # 
      # Args:
      #  y:       data
      #  design:  design
      #  ...:     parameters to the log-likelihood function to calculate derivative for
      #
      # Returns:
      #  Matrix (n.parxn.obs) of derivatives
      
      params <- list(...)
      ll0 <- do.call(log_likelihood_function, c(list(y=y, design=design), params))
      vapply(par_index, function(i) {
        params[[i]] <- params[[i]]+h
        (do.call(log_likelihood_function, c(list(y=y, design=design), params))-ll0)/h
      }, rep(0, NROW(y)))
    },
    
    jacobian_dbdens_domega = function(eta){
      inv_sqrt_omega <- solve(omega_sqrt)
      domega <- omega
      domega[] <- 0
      sapply(omega_derivatives, function(domega) {
        M <- inv_sqrt_omega%*%domega%*%t(inv_sqrt_omega)
        M[upper.tri(M)] <- 0
        diag(M) <- 0.5*diag(M)
        d_sqrt_omega_d_omega <- omega_sqrt%*%M
        (-tr(d_sqrt_omega_d_omega%*%inv_sqrt_omega)+eta%*%inv_sqrt_omega%*%d_sqrt_omega_d_omega%*%eta)
      })
    },
    
    dchol_doemgau = function(eta){
      inv_sqrt_omega <- solve(omega_sqrt)
      domega <- omega
      domega[] <- 0
      sapply(re_index, function(i) {
        domega[i,i] <- 1
        M <- inv_sqrt_omega%*%domega%*%t(inv_sqrt_omega)
        M[upper.tri(M)] <- 0
        diag(M) <- 0.5*diag(M)
        d_sqrt_omega_d_omega <- omega_sqrt%*%M
        (d_sqrt_omega_d_omega%*%eta)
      })
    }, 
    
    gradient_joint_likelihood = function(y, design, eta){
      # Calculates the gradient of the joint likelihood function wrt to mu and b
      # 
      # Args:
      #  y:       data
      #  design:  design
      #  b:  parameter vector of random effects to take the derivative at
      #  mu: parameter vector of fixed effects to take the derivative at
      #
      # Returns:
      #  Vector (n.mu+n.b) of derivatives
      b <- omega_sqrt %*% eta
      dLLdP <- do.call(gradient_dll_dpar, c(list(y, design), parameter_function(mu,b)))
      dPdMu <- dpardmu(b)
      #dPdb <- dpardb(b)
      #dCholdOmegau <- dchol_doemgau(eta)
      dLLdMu <- .colSums(dLLdP, NROW(y), par_count) %*% dPdMu
      #dLLdb <- .colSums(dLLdP, NROW(y), fe_count) %*% dPdb 
      #dLLdOmega <- dLLdb%*%dCholdOmegau
      dLLdOmega <- jacobian_dbdens_domega(eta)
      # likelihood scaling
      ll <- sum(log_likelihood_eta_function(y, design, eta))-ll_max
      like <- ifelse(ll>min_ll_allowed,exp(ll),0)
      return(c(dLLdMu, dLLdOmega)*like)
    },
    
    d2_log_likelihood_deta2 = function(y, design, eta){
      hessian(function(x) sum(log_likelihood_eta_function(y, design, x)))
    },
    
    simulation_eta_function_vec = function(design, eta){
      sim <- apply(eta, 1, FUN=simulation_eta_function, design=design)
      max_length <- max(unlist(lapply(sim, NROW))) 
      filled_sim <- lapply(sim, function(x) c(x, rep(0, max_length-NROW(x))))
      matrix(unlist(filled_sim), nrow=NROW(eta), byrow=T)
    },
    
    log_likelihood_eta_function_vec = function(y, design, eta) {
      .colSums(apply(eta,1,FUN=log_likelihood_eta_function, y=y, design=design), NROW(y), NROW(eta))
    },
    
    likelihood_eta_function_vec = function(y, design, eta) {
      ll <- log_likelihood_eta_function_vec(y, design, eta)-ll_max  
      ifelse(ll>min_ll_allowed, exp(ll), 0)
    },
    
    marginal_likelihood = function(y, design, eta) {
      if(adaptive){
        integrate_gq(likelihood_eta_function_vec,re_count, center=eta, 
                     settings = gq_settings, y=y, design=design)[1]
      }else{
        integrate_gq(likelihood_eta_function_vec,re_count,
                     settings = gq_settings,  y=y, design=design)[1]      
      }
    },
    
    
    gradient_joint_likelihood_vec = function(y, design, eta) {
      apply(eta, 1, FUN=function(eta.col) gradient_joint_likelihood(y, design, eta.col))
    },
    
    observed_fim= function(y, design, eta_sim) {
      # prevent numerical underflow by shifting all likelihood values upwards i.e. normalizing the joint likelihood to 1 at the simulated eta value
      if(use_scaled_likelihood){
        ll_max <<- sum(log_likelihood_eta_function(y, design, eta_sim))  
      }
      if(adaptive){
        if(adaptation_mode == 'mode' | adaptation_mode == 'center'){
          f <- function(eta) -sum(log_likelihood_eta_function(y, design, eta)) - log(prod(dnorm(eta))/det(omega_sqrt))
          res <- optim(eta_sim, f, method = "BFGS", hessian = T)
          center <- res$par
          hess <- -hessian(function(x) sum(log_likelihood_eta_function(y, design, x)), center)+
            diag(re_count)
          scale <- solve(hess)
          #scale <- solve(res$hessian)
          if(adaptation_mode == 'center') scale <- diag(re_count)
        }else{
          center <- eta_sim
          hess <- -hessian(function(x) sum(log_likelihood_eta_function(y, design, x)), center)+
            diag(re_count)
          scale <- solve(hess)
        }  
        if(use_scaled_likelihood){
          ll_max <<- sum(log_likelihood_eta_function(y, design, center))  
        }
        
      }else{
        center <- eta_zero
        scale <- diag(re_count)
      }
      gradient <- integrate_gq(gradient_joint_likelihood_vec, re_count, 
                               center=center, 
                               scale=scale, 
                               settings = gq_settings, y=y, design=design)
      ml <- integrate_gq(likelihood_eta_function_vec,re_count, 
                         center=center, 
                         scale=scale, 
                         settings = gq_settings, y=y, design=design)[1]
      if(is.na(ml) || ml==0 || is.nan(gradient)) return(matrix(0, nrow = theta_count, ncol = theta_count))
      if(any(is.nan(gradient %*% t(gradient)/ml^2))) browser()   
      gradient %*% t(gradient)/ml^2
    },
 

    fim = function(design){
      fim <- calc_fim(design)
      rownames(fim) <- theta_names
      colnames(fim)<- theta_names
      return(fim)
    },
    
    update_settings = function(new_settings=defaults.agq()){
      # option for adaptive quadrature
      adaptive <<- new_settings$adaptive
      # step length for numerical derivatives
      h <<- new_settings$derivative_step
      
      gq_settings  <<- filter_settings(new_settings, "gq")
      
      mc_settings <<- filter_settings(new_settings, "mc")
      
      adaptation_mode  <<- new_settings$adaptation_mode
      
      lhs_sampling <<- new_settings$lhs_sampling
      
      settings_y_integration <<- filter_settings(new_settings, "y_integration")
      
      switch(settings_y_integration$method,
             mc = calc_fim <<- fim_std_sampling,
             lhs = {
               calc_fim <<- fim_inverse_sampling
               inverse_integration_function <<- integrate_mc_lhs
             },
             qrmc = {
               calc_fim <<- fim_inverse_sampling
               inverse_integration_function <<- integrate_qrmc
             },
             stop("Unknown data sampling method '", y_sampling$method,"'"))
      
      # option to switch likelihood scaling on or off to prevent numerical underflows
      use_scaled_likelihood <<- new_settings$scaled_likelihood
      
      if(!is.null(new_settings$seed)) set.seed(new_settings$seed)
    }
    
))


calc_fim <- function(model, design, settings=NULL){
  if(!is.null(settings)) model$update_settings(settings)
  return(model$fim(design))
}

calc_rse <- function(model, fim){
  sqrt(diag(solve(fim)))/model$theta*100  
}

simulate_data <- function(model, design, subjects=1000, long_format=F){
  y <- model$simulation_eta_function_vec(design, eta=matrix(rnorm(subjects*model$re_count), nrow=subjects))
  if(!long_format) {
    return(y)
  }else{
    do.call(rbind, lapply(seq_len(subjects), function(i) transform(subject=i, design, y=y[i,])))
  }
}

plot_model <- function(model, design, subjects=100){
  library(ggplot2)
  data <- simulate_data(model, design, subjects, long_format = T)
  data$y <- factor(data$y)
  facet_formula <- sprintf("~%s", paste(colnames(design), sep="+"))
  ggplot(data, aes(x=y))+geom_histogram()+facet_wrap(as.formula(facet_formula))
}


