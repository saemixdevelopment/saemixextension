

# Design matrices for the run
setClass(Class = "SaemixIndivStructure",
         contains = "SaemixIndivModel",
         representation=representation(
           # --- number of parameters, covariates,...
           N="integer", # number of 'subjects' associated to this variability level (eg: N=nb of subjects for IIV, N=sum(n_{occ,i}) for IOV,...)
           nb.cov = "integer", # number of covariates (=nrow(covariate.model)) ## TODO: transform covariates, possibly adding columns
           # --- design matrices
           Mcovariates="matrix", # Mcovariate: matrix with N rows and columns corresponding to the (transformed) covariate values for each parameter x covariate relationship (including mu for IIV, a column of 1's)
           COV="matrix", # N x nb(mu+beta) containing 1's for mu and covariate values for each subject;  used in M-step to compute mean.phi as MCOV %*% COV (phi = MCOV %*% COV + eta)
           # --- helper matrices
           #    COV1="matrix", # initialised as Uargs$COV[,Uargs$ind.fix1] = Uargs$COV[,index.fixedpar.fix] # maybe remove
           # modified in main to Uargs$COV[,Uargs$ind.fix11] = Uargs$COV[,index.fixedpariiv.estim] when flag.fmin is TRUE and kiter=saemix.options$nbiter.sa (ie end of simulated annealing and minimisation) => removed, make the change in M-step
           dstatCOV="matrix", # used in M-step to compute betas[Uargs$ind.fix1]= elements of (mu+beta) fixed to their starting value
           # defined as COV[,index.fixedpar.fix,drop=FALSE]%*%MCOV[index.fixedpar.fix,,drop=FALSE]
           # TBD if we keep, could use the full form, also not sure about the indices :-/
           MCOV0="matrix", # used in compute.Uy as args$MCOV0[args$j0.covariate]<-b0
           # phi0<-args$COV0 %*% args$MCOV0
           # defined in initialisation as MCOV[index.fixedparnoiiv.estim,index.omega.novar,drop=FALSE]
           COV0="matrix" # used in compute.Uy, COV[,index.fixedparnoiiv.estim,drop=FALSE]
           # # phi0<-args$COV0 %*% args$MCOV0 # same line as MCOV0
         ),
         validity=function(object){
           return(TRUE)
         }
)
