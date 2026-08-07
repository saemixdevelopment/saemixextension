This is the development github for the saemix package for R (https://cran.r-project.org/web/packages/saemix/index.html).

Known issues:
- There is currently a bug when using kernel 4 in saemix (see open issues).
  - in R/main_estep.R:247, the fourth MCMC kernel computes the candidate log-likelihood from phiM instead of phiMc, effectively not using the predictions in the computation of the difference in log-likelihood
  - kernel 4 is off in the defaults (nbiter.mcmc = c(2,2,2,0)), so this only affects users who turn it on.
