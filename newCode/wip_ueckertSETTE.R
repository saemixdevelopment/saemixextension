###### Exact FIM by AGQ (code by Sebastian Ueckert)

For non-Gaussian models, the exact FIM should be computed, and two approaches have been proposed using either numerical integration by a combination of MC and adaptive Gaussian quadrature (MC/AGQ, Ueckert et al 2017) or stochastic integration by MCMC (Rivière et al. 2017).

Both these approaches are computationally intensive. 

**Eco** pas sûre que la méthode pour calculer la FIM exacte dans un protocole de population soit utilisable pour calculer la FIM empirique d'un jeu de données... 

**Eco** no example for TTE models in Sebastian's code so need to define the different functions by hand...

Here we use code provided by Sebastian Ueckert implementing the MC/AGQ approach, as the MCMC requires the installation of rStan. In this approach, the information matrix (FIM) over the population is first decomposed the sum of the individual FIM:
  $$
  FIM(\Psi, \Xi) = \sum_{i=1}^{N} FIM(\Psi, \xi_i)
$$
  where $\xi_i$ denotes the individual design in subject $i$. Assuming $Q$ different elementary designs, the FIM can also be summed over the different designs weighted by the number of subjects $N_q$ in design $q$ as:
  $$
  FIM(\Psi, \Xi) = \sum_{q=1}^{Q} N_q FIM(\Psi, \xi_q)
$$
  
  In the following, we first load the functions needed to compute the exact FIM. We then define a model object with the following components:
  
  - *parameter_function*: a function returning the list of parameters as the combination of fixed and random effects
- *log_likelihood_function*: using the parameters, computes the log-likelihood for all y in the dataset
- *simulation_function*: using the parameters, computes the log-likelihood and produces a random sample from the corresponding distribution
- *inverse_simulation_function*: supposed to be the quantile function but not quite sure :-/ (here, returns the category in which is urand)
- *mu*: the fixed parameters
- *omega*: the variance-covariance matrix

For *mu* and *omega*, we use the results from the saemix fit.

**TODO Alexandra** ?
  
  ```{r TTEmodelFIMSEbyAGQ}
# Code Sebastian
source(file.path(dirAGQ,"default_settings.R"))
source(file.path(dirAGQ,"helper_functions.R"))
source(file.path(dirAGQ,"integration.R"))
source(file.path(dirAGQ,"model.R"))

saemix.fit <- tte.fit

# TODO - adapt to TTE model ???
# model <- Model$new(
#   parameter_function = function(mu, b) list(alp1=mu[1]+b[1], alp2=mu[2], alp3=mu[3], alp4=mu[4], beta=mu[5] + b[2]),
#   log_likelihood_function = function(y, design, alp1, alp2, alp3, alp4, beta) {
#     logit1<-alp1 + beta*design$time
#     logit2<-logit1+alp2
#     logit3<-logit2+alp3
#     logit4<-logit3+alp4
#     pge1<-exp(logit1)/(1+exp(logit1))
#     pge2<-exp(logit2)/(1+exp(logit2))
#     pge3<-exp(logit3)/(1+exp(logit3))
#     pge4<-exp(logit4)/(1+exp(logit4))
#     pobs = (y==1)*pge1+(y==2)*(pge2 - pge1)+(y==3)*(pge3 - pge2)+(y==4)*(pge4 - pge3)+(y==5)*(1 - pge4)
#     log(pobs)
#   },
#   simulation_function = function(design, alp1, alp2, alp3, alp4, beta) {
#     logit1<-alp1 + beta*design$time
#     logit2<-logit1+alp2
#     logit3<-logit2+alp3
#     logit4<-logit3+alp4
#     pge1<-exp(logit1)/(1+exp(logit1))
#     pge2<-exp(logit2)/(1+exp(logit2))
#     pge3<-exp(logit3)/(1+exp(logit3))
#     pge4<-exp(logit4)/(1+exp(logit4))
#     x<-runif(length(time))
#     ysim<-1+as.integer(x>pge1)+as.integer(x>pge2)+as.integer(x>pge3)+as.integer(x>pge4)
#   },
#   inverse_simulation_function = function(design, urand,alp1, alp2, alp3, alp4, beta) {
#     if(is.null(urand)) return(seq_along(design$time))
#     logit1<-alp1 + beta*design$time
#     logit2<-logit1+alp2
#     logit3<-logit2+alp3
#     logit4<-logit3+alp4
#     pge1<-exp(logit1)/(1+exp(logit1))
#     pge2<-exp(logit2)/(1+exp(logit2))
#     pge3<-exp(logit3)/(1+exp(logit3))
#     pge4<-exp(logit4)/(1+exp(logit4))
#     1+as.integer(urand>pge1)+as.integer(urand>pge2)+as.integer(urand>pge3)+as.integer(urand>pge4)
#   },
#   mu = saemix.fit@results@fixed.effects,
#   omega = saemix.fit@results@omega[c(1,5),c(1,5)]
#   )
#   log_likelihood_function = function(y, design, lambda, beta) {
#     hazard <- (beta/lambda)*(design$time/lambda)^(beta-1) # ln(H')
#     H <- (design$time/lambda)^beta # ln(H)
#     logpdf<-rep(0,length(design$time))
#     logpdf[]
# 
#   },
# 
# weibulltte.model<-function(psi,id,xidep) {
#   T<-xidep[,1]
#   y<-xidep[,2] # events (1=event, 0=no event)
#   cens<-which(xidep[,3]==1) # censoring times (subject specific)
#   init <- which(T==0)
#   lambda <- psi[id,1] # Parameters of the Weibull model
#   beta <- psi[id,2]
#   Nj <- length(T)
#   
#   ind <- setdiff(1:Nj, append(init,cens)) # indices of events
#   hazard <- (beta/lambda)*(T/lambda)^(beta-1) # ln(H')
#   H <- (T/lambda)^beta # ln(H)
#   logpdf <- rep(0,Nj) # ln(l(T=0))=0
#   logpdf[cens] <- -H[cens] + H[cens-1] # ln(l(T=censoring time))
#   logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind]) # ln(l(T=event time))
#   return(logpdf)
# }
# 
# simulateWeibullTTE <- function(psi,id,xidep) {
#   T<-xidep[,1]
#   y<-xidep[,2] # events (1=event, 0=no event)
#   cens<-which(xidep[,3]==1) # censoring times (subject specific)
#   init <- which(T==0)
#   lambda <- psi[,1] # Parameters of the Weibull model
#   beta <- psi[,2]
#   Nj <- length(T)
#   ind <- setdiff(1:Nj, append(init,cens)) # indices of events
#   tevent<-T
#   Vj<-runif(dim(psi)[1])
#   tsim<-lambda*(-log(Vj))^(1/beta) # nsuj events
#   tevent[T>0]<-tsim
#   tevent[tevent[cens]>T[cens]] <- T[tevent[cens]>T[cens]]
#   return(tevent)
# }


# TBContinued
```