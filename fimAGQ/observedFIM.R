
observed_fim= function(y, design, eta_sim, options=list(use_scaled_likelihood=TRUE, adaptive=TRUE, adaptation_mode="mode")) {
  # prevent numerical underflow by shifting all likelihood values upwards i.e. normalizing the joint likelihood to 1 at the simulated eta value
  if(!is.null(options$use_scaled_likelihood)) options$use_scaled_likelihood <-TRUE 
  if(!is.null(options$adaptive)) options$adaptive <-TRUE
  if(!is.null(options$adaption_mode)) options$adaption_mode <-"mode" 
  
  if(options$use_scaled_likelihood){
    ll_max <<- sum(log_likelihood_eta_function(y, design, eta_sim))  
  }
  if(options$adaptive){
    if(options$adaptation_mode == 'mode' | options$adaptation_mode == 'center'){
      f <- function(eta) -sum(log_likelihood_eta_function(y, design, eta)) - log(prod(dnorm(eta))/det(omega_sqrt))
      res <- optim(eta_sim, f, method = "BFGS", hessian = T)
      center <- res$par
      hess <- -hessian(function(x) sum(log_likelihood_eta_function(y, design, x)), center)+
        diag(re_count)
      scale <- solve(hess)
      #scale <- solve(res$hessian)
      if(options$adaptation_mode == 'center') scale <- diag(re_count)
    }else{
      center <- eta_sim
      hess <- -hessian(function(x) sum(log_likelihood_eta_function(y, design, x)), center)+
        diag(re_count)
      scale <- solve(hess)
    }  
    if(options$use_scaled_likelihood){
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
}
