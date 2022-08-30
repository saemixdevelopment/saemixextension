####################################################################################
####		Defining class containing data, model and options		####
####################################################################################

#' @include aaa_generics.R
NULL

#' Class "SaemixObject"
#' 
#' An object of the SaemixObject class, storing the input to saemix, and the results obtained by a call
#' to the SAEM algorithm
#' 
#' Details of the algorithm can be found in the pdf file accompanying the package.
#' 
#' @name SaemixObject-class
#' @docType class
#' @aliases SaemixObject-class SaemixObject [<-,SaemixObject-method
#' print,SaemixObject predict,SaemixObject showall,SaemixObject show,SaemixObject summary,SaemixObject
#' @section Objects from the Class:
#' An object of the SaemixObject class is created after a call to \code{\link{saemix}} and contain the following slots:
#'   \describe{
#'     \item{\code{data}:}{Object of class \code{"SaemixData"}: saemix dataset, created by a call to \code{saemixData}}
#'     \item{\code{model}:}{Object of class \code{"SaemixModel"}: saemix model, created by a call to \code{saemixModel}}
#'     \item{\code{results}:}{Object of class \code{"SaemixData"}: saemix dataset, created by a call to \code{saemixData}}
#'     \item{\code{rep.data}:}{Object of class \code{"SaemixRepData"}: (internal) replicated saemix dataset, used the execution of the algorithm}
#'     \item{\code{sim.data}:}{Object of class \code{"SaemixSimData"}: simulated saemix dataset}
#'     \item{\code{options}:}{Object of class \code{"list"}: list of settings for the algorithm}
#'     \item{\code{prefs}:}{Object of class \code{"list"}: list of graphical options for the graphs}
#'   }
#' @section Methods:
#'   \describe{
#'     \item{[<-}{\code{signature(x = "SaemixObject")}: replace elements of object}
#'     \item{[}{\code{signature(x = "SaemixObject")}: access elements of object}
#'     \item{initialize}{\code{signature(.Object = "SaemixObject")}: internal function to initialise object, not to be used}
#'     \item{plot}{\code{signature(x = "SaemixObject")}: plot the data}
#'     \item{print}{\code{signature(x = "SaemixObject")}: prints details about the object (more extensive than show)}
#'     \item{showall}{\code{signature(object = "SaemixObject")}: shows all the elements in the object}
#'     \item{show}{\code{signature(object = "SaemixObject")}: prints details about the object}
#'     \item{summary}{\code{signature(object = "SaemixObject")}: summary of the object. Returns a list with a number of elements extracted from the object.}
#' 	 }
#' @references E Comets, A Lavenu, M Lavielle M (2017). Parameter estimation in nonlinear mixed effect models using saemix,
#' an R implementation of the SAEM algorithm. Journal of Statistical Software, 80(3):1-41.
#' 
#' E Kuhn, M Lavielle (2005). Maximum likelihood estimation in nonlinear mixed effects models. 
#' Computational Statistics and Data Analysis, 49(4):1020-1038.
#' 
#' E Comets, A Lavenu, M Lavielle (2011). SAEMIX, an R version of the SAEM algorithm. 20th meeting of the 
#' Population Approach Group in Europe, Athens, Greece, Abstr 2173.
#' 
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}
#' @author Audrey Lavenu
#' @author Marc Lavielle.
#' @seealso \code{\link{SaemixData}} \code{\link{SaemixModel}} \code{\link{saemixControl}} \code{\link{saemix}}
#' \code{\link{plot.saemix}},
#' @keywords classes
#' @include SaemixData.R
#' @include SaemixModel.R
#' @include SaemixRes.R
#' @exportClass SaemixObject
#' @examples
#' 
#' showClass("SaemixObject")
#' 

###############################
# definition
setClass(Class="SaemixObject",
  representation=representation(
    data="SaemixData",		# Data
    model="SaemixModel",	# Model
    results="SaemixRes",	# Fit results
    rep.data="SaemixRepData",	# Data replicates during algorithm (nb chains)
    sim.data="SaemixSimData", 	# Simulated data
    options="list",		# Options and parameters for algorithm
    prefs="list"		# Options for graphs
  ),
  validity=function(object){
#    cat ("--- Checking SaemixObject object ---\n")
    validObject(object@data)
    validObject(object@model)
    return(TRUE)
  }
)

###############################
# Initialize
#' @rdname initialize-methods
#' 
#' @param .Object a SaemixObject, SaemixRes, SaemixData or SaemixModel object to initialise
#' @param data an SaemixData object
## #' @param model an SaemixModel object
#' @param options a list of options passed to the algorithm
#' 
#' @exportMethod initialize

setMethod(
  f="initialize",
  signature="SaemixObject",
  definition= function (.Object, data, model, options=list()){
#    cat ("--- initialising SaemixObject --- \n")
    ind.exp<-which(model["error.model"]=='exponential')
    if(length(ind.exp)>0) {
      y<-yobs<-data["data"][,data["name.response"]]
      for(ityp in ind.exp) y[data["data"][,data["name.ytype"]]==ityp]<-log(cutoff(yobs[data["data"][,data["name.ytype"]]==ityp]))
      data["yorig"]<-yobs
      data["data"][,data["name.response"]]<-y
    }
    .Object@data<-data
# Adjusting number of covariates
    if(dim(model@covariate.model)[1]>length(data@name.covariates)) {
      cat("The number of covariates in model (",dim(model@covariate.model)[1],") is larger than the number of covariates in the dataset (",length(data@name.covariates),"), keeping only the first",length(data@name.covariates),":",data@name.covariates,".\n")
      model@covariate.model<-model@covariate.model[1:length(data@name.covariates),]
    }
    if(dim(model@covariate.model)[1]<length(data@name.covariates) & dim(model@covariate.model)[1]>0) {
    	cat("The number of covariates in model (",dim(model@covariate.model)[1],") is smaller than the number of covariates in the dataset (",length(data@name.covariates),"), assuming no covariate-parameter relationship for the remaining covariates; please check covariates:",data@name.covariates,".\n")
    	l1<-rep(0,dim(model@covariate.model)[2])
    	n1<-length(data@name.covariates)-dim(model@covariate.model)[1]
    	model@covariate.model<-rbind(model@covariate.model, matrix(rep(l1,n1),nrow=n1))
    }
# setting the names of the fixed effects
    if(dim(model@covariate.model)[1]>0) {
      nam.with.cov<-rep("",length(model@covariate.model))
      row1<-matrix(rep(model@name.modpar,length(data@name.covariates)), ncol=length(model@name.modpar),byrow=TRUE)
      col1<-matrix(rep(data@name.covariates,length(model@name.modpar)), ncol=length(model@name.modpar))
      idcov<-which(model@covariate.model==1)
      nam.with.cov[idcov]<-paste("beta_",col1[idcov],"(",row1[idcov],")",sep="")
      nam1<-rbind(model@name.modpar,matrix(nam.with.cov, ncol=length(model@name.modpar)))
      nam1<-c(nam1)
      model@name.fixed<-nam1[nam1!=""]
    } else model@name.fixed<-model@name.modpar
    i1.omega2<-model@indx.omega
    model@name.random<-paste("omega2",model@name.modpar[model@indx.omega], sep=".")
    .Object@model<-model
    .Object@model@betaest.model<-matrix(c(rep(1,.Object@model@nb.parameters), c(t(.Object@model@covariate.model))),ncol=.Object@model@nb.parameters,byrow=TRUE)
# Aligning the name of the response for the model object to the same as in the data object
    if(model@name.response=="") model@name.response<-data@name.response
    
# Covariates
    .Object@model@name.cov<-.Object@data@name.covariates
    if(length(.Object@model@name.cov)>0 & sum(.Object@model@covariate.model)>0) {
      try(rownames(.Object@model@covariate.model)<-.Object@model@name.cov)
      try(rownames(.Object@model@betaest.model)[2:(1+ length(.Object@model@name.cov))]<-.Object@model@name.cov)
    }
    ucov <- rownames(.Object@model@covariate.model)[ rowSums(.Object@model@covariate.model)>0]
    if(length(ucov)>0) {
      for(icov in ucov) {
      	xdat<-subset(.Object@data@data,is.na(get(icov)))
      	if(dim(xdat)[1]>0) {
      		imis<-unique(xdat[,.Object@data@name.group])
      		cat("Missing values for covariate",as.character(icov),"for which a parameter-covariate relationship is estimated: removing subject(s)",imis,"from the dataset.\n")
      	}
      	.Object@data<-subset(.Object@data,!is.na(get(icov)))
      }
    }
# Initialising options
    opt<-saemixControl()
    if(length(options)>0) {
      for(i in names(options)) opt[i]<-options[i]
      while(length(options["nbiter.mcmc"][[1]])<4) options["nbiter.mcmc"][[1]]<-c(options["nbiter.mcmc"][[1]],0) # nb of kernels now 4, complete if shorter
      if(!opt$fix.seed) {
        rm(.Random.seed)
        runif(1)
        opt$seed<-.Random.seed[5]
      }
      if(is.null(options$nb.chains) & data@N>0) opt$nb.chains<-ceiling(50/data@N)
      if(data@N>0 && data@N<50 & opt$nb.chains<ceiling(50/data@N)) {
        cat("The number of subjects is small, increasing the number of chains to", ceiling(50/data@N),"to improve convergence\n")
        opt$nb.chains<-ceiling(50/data@N)
      }
      # Options that would have been changed in saemixControl and need to be set here from the elements of the options list
      if(is.null(options$nbiter.sa)) opt$nbiter.sa<-opt$nbiter.saemix[1]/2
      if(opt$nbiter.sa>opt$nbiter.saemix[1]) {
        if(opt$warnings) message("The number of iterations for the simulated annealing should be lower or equal to the number of iterations in the first stage of the algorithm, setting it to K1=nbiter.saemix[1]. We advise setting it to nbiter.saemix[1]/2.\n")
        opt$nbiter.sa<-opt$nbiter.saemix[1]
      }
      if(!is.null(options$fix.seed) && !(options$fix.seed)) {
        rm(.Random.seed)
        runif(1)
        opt$seed<-.Random.seed[5]
      }
    } else{
      opt$nb.chains<-opt$nb.chains
    }
    if(opt$ipar.lmcmc<2) {
      opt$ipar.lmcmc<-2
      cat("Value of L_MCMC too small, setting it to 2 (computation of the conditional means and variances of the individual parameters)\n")
    }

    # if(is.null(opt$nbiter.sa)){
    #   opt$nbiter.sa<-round(opt$nbiter.saemix[1]/2)
    # } else{
    #   opt$nbiter.sa<-max(1,opt$nbiter.sa)
    # }

    opt$nbiter.tot<-sum(opt$nbiter.saemix)
    .Object@options<-opt
# Options for plots
    .Object@prefs<-saemix.plot.setoptions(.Object)
# Object validation
    validObject(.Object)
    return (.Object )
  }
)

###########################	Default options		#############################

#' List of options for running the algorithm SAEM
#' 
#' List containing the variables relative to the optimisation algorithm. All
#' these elements are optional and will be set to default values when running
#' the algorithm if they are not specified by the user.
#' 
#' All the variables are optional and will be set to their default value when
#' running \code{\link{saemix}}.
#' 
#' The function \code{\link{saemix}} returns an object with an element options
#' containing the options used for the algorithm, with defaults set for
#' elements which have not been specified by the user.
#' 
#' These elements are used in subsequent functions and are not meant to be used
#' directly.
#' 
#' @param map a boolean specifying whether to estimate the individual parameters (MAP estimates). Defaults to TRUE
#' @param fim a boolean specifying whether to estimate the Fisher Information Matrix and derive the estimation errors 
#' for the parameters. Defaults to TRUE. The linearised approximation to the log-likelihood is also computed in the process
#' @param ll.is a boolean specifying whether to estimate the log-likelihood by importance sampling. Defaults to TRUE
#' @param ll.gq a boolean specifying whether to estimate the log-likelihood by Gaussian quadrature. Defaults to FALSE
#' @param nbiter.saemix nb of iterations in each step (a vector containing 2
#' elements, nbiter.saemix\[1\] for the exploration phase of the algorithm (K1) and nbiter.saemix\[2\]
#' for the smoothing phase (K2))
#' @param nb.chains nb of chains to be run in parallel in the MCMC algorithm
#' (defaults to 1)
#' @param nbiter.burn nb of iterations for burning
#' @param nbiter.map nb of iterations of the MAP kernel (4th kernel) to run at the beginning 
#' of the estimation process (defaults to nbiter.saemix\[1\]/10 if nbiter.mcmc\[4\] is more than 0) 
#' (EXPERIMENTAL, see Karimi et al. 2019 for details)
#' @param nbiter.mcmc nb of iterations in each kernel during the MCMC step
#' @param nbiter.sa nb of iterations subject to simulated annealing (defaults to nbiter.saemix\[1\]/2, 
#' will be cut down to K1=nbiter.saemix\[1\] if greater than that value). We recommend to stop 
#' simulated annealing before the end of the exploration phase (nbiter.saemix\[1\]).
#' @param proba.mcmc probability of acceptance
#' @param stepsize.rw stepsize for kernels q2 and q3 (defaults to 0.4)
#' @param rw.init initial variance parameters for kernels (defaults to 0.5)
#' @param alpha.sa parameter controlling cooling in the Simulated Annealing
#' algorithm (defaults to 0.97)
#' @param fix.seed TRUE (default) to use a fixed seed for the random number
#' generator. When FALSE, the random number generator is initialised using a
#' new seed, created from the current time.  Hence, different sessions started
#' at (sufficiently) different times will give different simulation results.
#' The seed is stored in the element seed of the options list.
#' @param seed seed for the random number generator. Defaults to 123456
#' @param nmc.is nb of samples used when computing the likelihood through
#' importance sampling
#' @param nu.is number of degrees of freedom of the Student distribution used
#' for the estimation of the log-likelihood by Importance Sampling. Defaults to
#' 4
#' @param print.is when TRUE, a plot of the likelihood as a function of the
#' number of MCMC samples when computing the likelihood through importance
#' sampling is produced and updated every 500 samples. Defaults to FALSE
#' @param nbdisplay nb of iterations after which to display progress
#' @param displayProgress when TRUE, the convergence plots are plotted after
#' every nbdisplay iteration, and a dot is written in the terminal window to
#' indicate progress. When FALSE, plots are not shown and the algorithm runs
#' silently. Defaults to FALSE
#' @param nnodes.gq number of nodes to use for the Gaussian quadrature when
#' computing the likelihood with this method (defaults to 12)
#' @param nsd.gq span (in SD) over which to integrate when computing the
#' likelihood by Gaussian quadrature. Defaults to 4 (eg 4 times the SD)
#' @param maxim.maxiter Maximum number of iterations to use when maximising the
#' fixed effects in the algorithm. Defaults to 100
#' @param nb.sim number of simulations to perform to produce the VPC plots or
#' compute npde. Defaults to 1000
#' @param nb.simpred number of simulations used to compute mean predictions
#' (ypred element), taken as a random sample within the nb.sim simulations used
#' for npde
#' @param ipar.lmcmc number of iterations required to assume convergence for
#' the conditional estimates. Defaults to 50
#' @param ipar.rmcmc confidence interval for the conditional mean and variance.
#' Defaults to 0.95
#' @param print whether the results of the fit should be printed out. Defaults
#' to TRUE
#' @param save whether the results of the fit should be saved to a file.
#' Defaults to TRUE
#' @param save.graphs whether diagnostic graphs and individual graphs should be
#' saved to files. Defaults to TRUE
#' @param directory the directory in which to save the results. Defaults to
#' "newdir" in the current directory
#' @param warnings whether warnings should be output during the fit. Defaults to FALSE
#' 
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' 
#' @seealso \code{\link{SaemixData}},\code{\link{SaemixModel}},
#' \code{\link{SaemixObject}}, \code{\link{saemix}}
#' 
#' @references E Comets, A Lavenu, M Lavielle M (2017). Parameter estimation in nonlinear mixed effect models using saemix,
#' an R implementation of the SAEM algorithm. Journal of Statistical Software, 80(3):1-41.
#' 
#' E Kuhn, M Lavielle (2005). Maximum likelihood estimation in nonlinear mixed effects models. 
#' Computational Statistics and Data Analysis, 49(4):1020-1038.
#' 
#' E Comets, A Lavenu, M Lavielle (2011). SAEMIX, an R version of the SAEM algorithm. 20th meeting of the 
#' Population Approach Group in Europe, Athens, Greece, Abstr 2173.
#' 
#' B Karimi, M Lavielle , E Moulines E  (2019). f-SAEM: A fast Stochastic Approximation of the EM algorithm for nonlinear mixed effects models.
#' Computational Statistics & Data Analysis, 141:123-38
#' 
#' @keywords models
#' @examples
#' 
#' # All default options
#' saemix.options<-saemixControl()
#' 
#' # All default options, changing seed
#' saemix.options<-saemixControl(seed=632545)
#' 
#' @export saemixControl

saemixControl<-function(map=TRUE,fim=TRUE,ll.is=TRUE,ll.gq=FALSE,nbiter.saemix=c(300,100), nbiter.sa=NA, nb.chains=1,fix.seed=TRUE,seed=23456,nmc.is=5000,nu.is=4, print.is=FALSE, nbdisplay=100, displayProgress=FALSE, nbiter.burn=5,nbiter.map=5, nbiter.mcmc=c(2,2,2,0),proba.mcmc=0.4,stepsize.rw=0.4,rw.init=0.5,alpha.sa=0.97,  nnodes.gq=12,nsd.gq=4,maxim.maxiter=100,nb.sim=1000,nb.simpred=100, ipar.lmcmc=50,ipar.rmcmc=0.05, print=TRUE, save=TRUE, save.graphs=TRUE,directory="newdir",warnings=FALSE) {
  # ECO: Need to duplicate some code in the creation of the object (initialize from SaemixObject) so that arguments like nbiter.sa or the random seed get correctly assigned :-/ 
  if(fix.seed) seed<-seed else {
    rm(.Random.seed)
    runif(1)
    seed<-.Random.seed[5]
  }
  if(ipar.lmcmc<2) {
    ipar.lmcmc<-2
    cat("Value of L_MCMC too small, setting it to 2 (computation of the conditional means and variances of the individual parameters)\n")
  }
  if(is.na(nbiter.sa)) nbiter.sa<-nbiter.saemix[1]/2
  if(nbiter.sa>nbiter.saemix[1]) {
    if(warnings) message("The number of iterations for the simulated annealing should be lower or equal to the number of iterations in the first stage of the algorithm, setting it to K1=nbiter.saemix[1]. We advise setting it to nbiter.saemix[1]/2.\n")
    nbiter.sa<-nbiter.saemix[1]
  }
  list(map=map,fim=fim,ll.is=ll.is,ll.gq=ll.gq,nbiter.saemix=nbiter.saemix, nbiter.sa=nbiter.sa, nbiter.burn=nbiter.burn, nbiter.map=nbiter.map,nb.chains=nb.chains,
       fix.seed=fix.seed,seed=seed, nmc.is=nmc.is,nu.is=nu.is,print.is=print.is, nbdisplay=nbdisplay,displayProgress=displayProgress,
       print=print,save=save, save.graphs=save.graphs,directory=directory,warnings=warnings, nbiter.mcmc=nbiter.mcmc,proba.mcmc=proba.mcmc,
       stepsize.rw=stepsize.rw, rw.init=rw.init,alpha.sa=alpha.sa,nnodes.gq=nnodes.gq,nsd.gq=nsd.gq, maxim.maxiter=maxim.maxiter,
       nb.sim=nb.sim,nb.simpred=nb.simpred, ipar.lmcmc=ipar.lmcmc,ipar.rmcmc=ipar.rmcmc)
}

####################################################################################
####			saemixObject class - accesseur				####
####################################################################################

##' Get/set methods for SaemixObject object
##' 
##' Access slots of a SaemixObject object using the object\["slot"\] format
##' 
#' @param x object
#' @param i element to be replaced
#' @param j element to replace with
#' @param drop whether to drop unused dimensions
#' @keywords methods
#' @exportMethod [
#' @exportMethod [<-
#' @exportPattern "^[[:alpha:]]+"
##### @name [
##### @aliases [,SaemixObject-method
##### @docType methods
##### @rdname extract-methods

# Getteur
setMethod(
  f ="[",
  signature = "SaemixObject" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "data"={return(x@data)},
    "model"={return(x@model)},
    "results"={return(x@results)},
    "rep.data"={return(x@rep.data)},
    "sim.data"={return(x@sim.data)},
    "options"={return(x@options)},
    "prefs"={return(x@prefs)},
    stop("No such attribute\n")
   )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixObject" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "data"={x@data<-value},
    "model"={x@model<-value},
    "results"={x@results<-value},
    "rep.data"={x@rep.data<-value},
    "sim.data"={x@sim.data<-value},
    "options"={x@options<-value},
    "prefs"={x@prefs<-value},
    stop("No such attribute\n")
   )
   validObject(x)
   return(x)
  }
)


####################################################################################
####				Summary method for SaemixObject			####
####################################################################################
#' @rdname summary-methods
#' 
#' @exportMethod summary

setMethod("summary","SaemixObject",
  function(object, print=TRUE, ...) {
    if(object@results@status %in% c("empty")) {
      cat("Object of class SaemixObject, no fit performed yet.\n")
      return()
    }
    digits<-2;nsmall<-2
    if(print) {
    cat("----------------------------------------------------\n")
    cat("-----------------  Fixed effects  ------------------\n")
    cat("----------------------------------------------------\n")
    }
#    browser()
    if(length(object@results@se.fixed)==0) {
       if(object@model@modeltype=="structural") {
      tab<-data.frame(c(object@results@name.fixed, object@results@name.sigma[object@results@indx.res]), c(object@results@fixed.effects,object@results@respar[object@results@indx.res]))
        }else{
          tab<-data.frame(c(object@results@name.fixed), c(object@results@fixed.effects))
        }
      colnames(tab)<-c("Parameter","Estimate")
    } else {
       if(object@model@modeltype=="structural") {
      tab<-data.frame(c(object@results@name.fixed, object@results@name.sigma[object@results@indx.res]), c(object@results@fixed.effects,object@results@respar[object@results@indx.res]),c(object@results@se.fixed,object@results@se.respar[object@results@indx.res]), stringsAsFactors=FALSE)
        }else{
            tab<-data.frame(c(object@results@name.fixed), c(object@results@fixed.effects),c(object@results@se.fixed), stringsAsFactors=FALSE)
        }
      tab<-cbind(tab,100*abs(as.double(tab[,3])/as.double(tab[,2])))
      colnames(tab)<-c("Parameter","Estimate","SE","CV(%)")
      if(length(object@results@indx.cov)>0) {
      wstat<-as.double(tab[,2])/as.double(tab[,3])
      pval<-rep("-",length(wstat))
      pval[object@results@indx.cov]<-1-normcdf(abs(wstat[object@results@indx.cov]))
      tab<-cbind(tab,"p-value"=pval,stringsAsFactors=FALSE)
      }
    }
    tab.fix<-tab
    for(i in 2:dim(tab)[2]) {
     xcol<-as.double(as.character(tab[,i]))
     idx<-which(!is.na(xcol) & xcol!="-")
     tab[idx,i]<-format(xcol[idx],digits=digits,nsmall=nsmall)
    }
    if(print) {
    print(tab,quote=FALSE)
    cat("----------------------------------------------------\n")
    cat("-----------  Variance of random effects  -----------\n")
    cat("----------------------------------------------------\n")
#  cat("   ECO TODO: check if Omega or Omega2 (SD or variances) and can we choose ?\n")
    }
    if(length(object@results@se.omega)==0) {
      tab<-data.frame(object@results@name.random, diag(object@results@omega)[object@results@indx.omega])
      colnames(tab)<-c("Parameter","Estimate")
    } else {
      tab<-data.frame(object@results@name.random, diag(object@results@omega)[object@results@indx.omega], object@results@se.omega[object@results@indx.omega])
      tab<-cbind(tab,100*as.double(tab[,3])/as.double(tab[,2]))
      colnames(tab)<-c("Parameter","Estimate","SE","CV(%)")
    }
    tab.random<-tab
    for(i in 2:dim(tab)[2]) 
      tab[,i]<-format(as.double(as.character(tab[,i])),digits=digits,nsmall=nsmall)
    if(print) {
    print(tab,quote=FALSE)
    cat("----------------------------------------------------\n")
    cat("------  Correlation matrix of random effects  ------\n")
    cat("----------------------------------------------------\n")
    }
    tab<-cov2cor(object@results@omega[object@results@indx.omega, object@results@indx.omega,drop=FALSE])
    tab.corr<-tab
    for(i in 1:dim(tab)[2])
      tab[,i]<-format(as.double(as.character(tab[,i])),digits=digits,nsmall=nsmall)
    try(colnames(tab)<-rownames(tab)<-object@results@name.random)
    if(print) print(tab,quote=FALSE)
    l1<-rep(NA,3)
    tab.ll<-data.frame(Method=c("Linearisation","Importance Sampling","Gaussian Quadrature"),"-2xLL"=l1,AIC=l1,BIC=l1)
    if(length(object@results@ll.lin)>0 | length(object@results@ll.is)>0 | length(object@results@ll.gq)>0) {
    	if(print) {
    cat("----------------------------------------------------\n")
    cat("---------------  Statistical criteria  -------------\n")
    cat("----------------------------------------------------\n")
    	}
    if(length(object@results@ll.lin)>0) {
    	if(print) {
    cat("Likelihood computed by linearisation\n")
    cat("      -2LL=",(-2*object@results@ll.lin),"\n")
    cat("      AIC =",object@results@aic.lin,"\n")
    cat("      BIC =",object@results@bic.lin,"\n")
    	}
    tab.ll[1,2:4]<-c((-2*object@results@ll.lin),object@results@aic.lin, object@results@bic.lin)
    }
    if(length(object@results@ll.is)>0) {
    	if(print) {
    cat("\nLikelihood computed by importance sampling\n")
    cat("      -2LL=",(-2*object@results@ll.is),"\n")
    cat("      AIC =",object@results@aic.is,"\n")
    cat("      BIC =",object@results@bic.is,"\n")
    	}
    tab.ll[2,2:4]<-c((-2*object@results@ll.is),object@results@aic.is, object@results@bic.is)
    }  
    if(length(object@results@ll.gq)>0) {
    	if(print) {
    cat("\nLikelihood computed by Gaussian quadrature\n")
    cat("      -2LL=",(-2*object@results@ll.gq),"\n")
    cat("      AIC =",object@results@aic.gq,"\n")
    cat("      BIC =",object@results@bic.gq,"\n")
    	}
    tab.ll[3,2:4]<-c((-2*object@results@ll.gq),object@results@aic.gq, object@results@bic.gq)
    }
    if(print) cat("----------------------------------------------------\n")
    }
    tab<-data.frame(Id=unique(object@data@data[,object@data@name.group]), object@results@cond.mean.psi)
    try(colnames(tab)[-c(1)]<-object@model@name.modpar,silent=TRUE)
    npar<-length(object@results@name.fixed)
    coef<-list(fixed=tab.fix[1:npar,2],random=list(map.psi=object@results@map.psi, cond.mean.psi=tab))
    sigma<-tab.fix[-c(1:npar),2]

    res<-list(modeltype=object@model@modeltype,fixed.effects=tab.fix,sigma=sigma,random.effects=tab.random, correlation.matrix=tab.corr,logLik=tab.ll,coefficients=coef)
    if(length(object@results@fim)>0) res$FIM<-object@results@fim
    res$data<-list(N=object@data@N,nobs=list(ntot=object@data@ntot.obs, nind=object@data@nind.obs),data=object@data@data)
    if(length(object@results@ypred)>0 | length(object@results@ipred)>0  | length(object@results@ppred)>0 | length(object@results@icpred)>0) {
      res$fitted<-list(population=list(pop.param=object@results@ppred, pop.mean=object@results@ypred),individual=list(map.ipred=object@results@ipred, cond.ipred=object@results@icpred))
    }
     if(length(object@results@wres)>0 | length(object@results@iwres)>0  | length(object@results@icwres)>0 | length(object@results@pd)>0) {
      res$residuals<-list(population=list(wres=object@results@wres), individual=list(map.iwres=object@results@iwres,cond.iwres=object@results@icwres, pd=object@results@pd, npde=object@results@npde))
    }
   
    invisible(res)
 }
)

####################################################################################
####			Print and show methods for SaemixObject			####
####################################################################################

#' @rdname print-methods
#' @exportMethod print

setMethod("print","SaemixObject",
  function(x,nlines=10,...) {
    cat("Nonlinear mixed-effects model fit by the SAEM algorithm\n")
    cat("-----------------------------------\n")
    cat("----          Data             ----\n")
    cat("-----------------------------------\n")
    print(x@data,nlines=nlines)
    cat("-----------------------------------\n")
    cat("----          Model            ----\n")
    cat("-----------------------------------\n")
    print(x@model)
    cat("-----------------------------------\n")
    cat("----    Key algorithm options  ----\n")
    cat("-----------------------------------\n")
    if(x@options$map) cat("    Estimation of individual parameters (MAP)\n")
    if(x@options$fim) cat("    Estimation of standard errors and linearised log-likelihood\n")
    if(x@options$ll.is) cat("    Estimation of log-likelihood by importance sampling\n")
    if(x@options$ll.gq) cat("    Estimation of log-likelihood by gaussian quadrature\n")
    if(as.integer(x@options$map+x@options$fim+x@options$ll.is+x@options$ll.gq)==0) cat("    Algorithms: estimation only\n")
    st1<-paste(c("K1=","K2="),x@options$nbiter.saemix,sep="",collapse=", ")
    cat("    Number of iterations: ",st1,"\n")
    cat("    Number of chains: ",x@options$nb.chains,"\n")
    cat("    Seed: ",x@options$seed,"\n")
    if(x@options$ll.is) cat("    Number of MCMC iterations for IS: ",x@options$nmc.is,"\n")
    cat("    Simulations:\n")
    cat("        nb of simulated datasets used for npde: ",x@options$nb.sim,"\n")
    cat("        nb of simulated datasets used for VPC: ",x@options$nb.simpred,"\n")
    cat("    Input/output\n")
    cat("        save the results to a file: ",x@options$save,"\n")
    cat("        save the graphs to files: ",x@options$save.graphs,"\n")
    if(x@options$save | x@options$save.graphs) cat("        directory where results should be saved: ",x@options$directory,"\n")
    cat("----------------------------------------------------\n")
    cat("----                  Results                   ----\n")
    print(x@results)
  }
)

#' @rdname show-methods
#' 
#' @exportMethod show

setMethod("show","SaemixObject",
  function(object) {
#    cat("Object of class SaemixObject\n")
    cat("Nonlinear mixed-effects model fit by the SAEM algorithm\n")
    cat("-----------------------------------------\n")
    cat("----         Data and Model          ----\n")
    cat("-----------------------------------------\n")
#    show(object@data)
    cat("Data\n")
    cat("    Dataset",object@data@name.data,"\n")
    st1<-paste(object@data@name.response," ~ ",paste(object@data@name.predictors, collapse=" + ")," | ", object@data@name.group,sep="")
    cat("    Longitudinal data:",st1,"\n\n")
#    show(object@model)
    cat("Model:\n")
    if(length(object@model@description)>0 && nchar(object@model@description)>0) {
      cat("   ",object@model@description,"\n")}
    fix1<-ifelse(object@model@fixed.estim==1,""," [fixed]")
    cat("    ",object@model@nb.parameters,"parameters:", paste(object@model@name.modpar,fix1,sep=""),"\n")
    cat("     error model:",object@model@error.model,"\n")
    if(dim(object@model@covariate.model)[1]>0) {
      cat("     covariate model:\n")
      print(object@model@covariate.model) 
    } else cat("     No covariate\n")
    cat("\n")
    cat("Key options\n")
    if(object@options$map) cat("    Estimation of individual parameters (MAP)\n")
    if(object@options$fim) cat("    Estimation of standard errors and linearised log-likelihood\n")
    if(object@options$ll.is) cat("    Estimation of log-likelihood by importance sampling\n")
    if(object@options$ll.gq) cat("    Estimation of log-likelihood by gaussian quadrature\n")
    if(as.integer(object@options$map+object@options$fim+object@options$ll.is+object@options$ll.gq)==0) cat("    Algorithms: estimation only\n")
    st1<-paste(c("K1=","K2="),object@options$nbiter.saemix,sep="",collapse=", ")
    cat("    Number of iterations: ",st1,"\n")
    cat("    Number of chains: ",object@options$nb.chains,"\n")
    cat("    Seed: ",object@options$seed,"\n")
    if(object@options$ll.is) cat("    Number of MCMC iterations for IS: ",object@options$nmc.is,"\n")
    cat("    Input/output\n")
    if(object@options$save)
      cat("        save the results to a file: ",object@options$save,"\n") else 
      cat("        results not saved\n")
    if(object@options$save.graphs)
      cat("        save the graphs to files: ",object@options$save.graphs,"\n") else 
      cat("        no graphs\n")
    if(object@options$save | object@options$save.graphs) cat("        directory where results are saved: ",object@options$directory,"\n")
    if(FALSE) {
    if(length(object@rep.data)>0) irep<-1 else irep<-0
    if(length(object@sim.data)>0) isim<-1 else isim<-0
    if(irep>0 | isim>0) {
      cat("-----------------------------------\n")
      cat("----      Other components     ----\n")
      cat("-----------------------------------\n")
      if(irep>0) cat("    Replicated data on ",object@options$nb.chains," chains\n")
      if(isim>0) cat("    Simulated data, ",object@options$nb.sim," simulations\n")
    }
    cat("-----------------------------------\n")
  }
  if(length(object@results@fixed.effects)>0) {
    cat("----------------------------------------------------\n")
    cat("----                  Results                   ----\n")
    show(object@results)
  }}
)


#' @rdname showall-methods
#' 
#' @exportMethod showall

# Could be print, with only head of data
setMethod("showall","SaemixObject",
  function(object) {
#    cat("Object of class SaemixObject\n")
    cat("Nonlinear mixed-effects model fit by the SAEM algorithm\n")
    cat("-----------------------------------\n")
    cat("----          Data             ----\n")
    cat("-----------------------------------\n")
    showall(object@data)
    cat("-----------------------------------\n")
    cat("----          Model            ----\n")
    cat("-----------------------------------\n")
    show(object@model)
    cat("-----------------------------------\n")
    cat("----      Algorithm options    ----\n")
    cat("-----------------------------------\n")
    if(object@options$map) cat("    Estimation of individual parameters (MAP)\n")
    if(object@options$fim) cat("    Estimation of standard errors and linearised log-likelihood\n")
    if(object@options$ll.is) cat("    Estimation of log-likelihood by importance sampling\n")
    if(object@options$ll.gq) cat("    Estimation of log-likelihood by gaussian quadrature\n")
    if(as.integer(object@options$map+object@options$fim+object@options$ll.is+object@options$ll.gq)==0) cat("    Algorithms: estimation only\n")
    st1<-paste(c("K1=","K2="),object@options$nbiter.saemix,sep="",collapse=", ")
    cat("    Number of chains: ",object@options$nb.chains,"\n")
    cat("    Number of iterations: ",st1,"\n")
    cat("        nb of iterations for SA: ",object@options$nbiter.sa,"\n")
    cat("        nb of burning iterations: ",object@options$nbiter.burn,"\n")
    cat("        nb of map kernel iterations: ",object@options$nbiter.map,"\n")
    cat("    Seed:\n")
    cat("        setting a random seed: ",!object@options$fix.seed,"\n")
    cat("        seed for the random number generator: ",object@options$seed,"\n")
    if(object@options$ll.is) {
    	cat("    Estimation of LL by Importance Sampling:\n")
    	cat("        number of MCMC samples: ",object@options$nmc.is,"\n")
    	cat("        nu for IS: ",object@options$nu.is,"\n")
    	cat("        produce plots during the estimation of LL by IS: ",object@options$print.is,"\n")
    }
    cat("    Input/output\n")
    cat("        display progress during the estimation process: ",object@options$displayProgress,"\n")
    cat("        nb of iterations after which to display progress: ",object@options$nbdisplay,"\n")
    cat("        print out the results after the fit: ",object@options$print,"\n")
    cat("        save the results to a file: ",object@options$save,"\n")
    cat("        save the graphs to files: ",object@options$save.graphs,"\n")
    if(object@options$save | object@options$save.graphs) cat("        directory where results should be saved: ",object@options$directory,"\n")
    cat("        whether warnings should be output during the fit: ",object@options$warnings,"\n")
    cat("    SAEM algorithm\n")
    cat("        number of MCMC iterations for each kernel: ",object@options$nbiter.mcmc,"\n")
    cat("        probability of acceptance: ",object@options$proba.mcmc,"\n")
    cat("        : ",object@options$stepsize.rw,"\n")
    cat("        : ",object@options$rw.init,"\n")
    cat("        : ",object@options$alpha.sa,"\n")
    cat("        maximum nb of iterations for estimation of fixed effects: ",object@options$maxim.maxiter,"\n")
    if(object@options$ll.gq) {
    	cat("    Estimation of LL by Gaussian Quadrature:\n")
    	cat("        number of nodes: ",object@options$nnodes.gq,"\n")
    	cat("        width of integral: ",object@options$nsd.gq,"\n")
    }
    cat("    Simulations:\n")
    cat("        nb of simulated datasets used for npde: ",object@options$nb.sim,"\n")
    cat("        nb of simulated datasets used for VPC: ",object@options$nb.simpred,"\n")
    cat("    Estimation of individual parameters\n")
    cat("        nb of iterations: ",object@options$ipar.lmcmc,"\n")
    cat("        : ",object@options$ipar.rhomcmc,"\n")
    cat("        size of confidence interval: ",object@options$ipar.rmcmc,"\n")
    cat("-----------------------------------\n")
    if(length(object@results@fixed.effects)>0) {
      cat("----------------------------------------------------\n")
      cat("----                  Results                   ----\n")
      show(object@results)
    }
  }
)

####################################################################################
####			SaemixObject class - method to predict			####
####################################################################################

#' @rdname predict-methods
#' 
#' @param object an SaemixObject
#' @param newdata an optional dataframe for which predictions are desired. If newdata is given, it must contain the predictors needed for the model in object
#' @param type the type of predictions (ipred= individual, ppred=population predictions obtained with the population estimates, ypred=mean of the population predictions, icpred=conditional predictions). With newdata, individual parameters can be estimated if the new data contains observations; otherwise, predictions correspond to the population predictions ppred, and type is ignored.
#' @param se.fit whether the SE are to be taken into account in the model predictions
#' @param ... additional arguments passed on to fitted()
#' 
#' @return a vector or a dataframe (if more than one type) with the corresponding predictions for each observation in the dataframe
#' 
#' @exportMethod predict

# Default is to return individual predictions using MAP. If newdata is given and does not have response data, population predictions are returned instead.

setMethod(f="predict",
          signature="SaemixObject",
          def=function(object,newdata=NULL,type=c("ipred", "ypred", "ppred", "icpred"), se.fit=FALSE, ...) {
            namObj<-deparse(substitute(object))
            type<-match.arg(type)
            type<-intersect(type,c("ipred", "ypred", "ppred", "icpred"))
            if(length(type)>1) type<-type[1]
            if(length(type)==0) type<-"ipred"
            #            print(namObj)
            saemix.data<-object["data"]
            saemix.model<-object["model"]
            #    se.fit<-match.arg(se.fit) # doesn't work with logical type, change
            #    if(se.fit & object@options$warnings) message("Currently predict() does not handle argument se.fit=TRUE.\n")
            if(missing(newdata)) { # Return predictions from fitted object
              xpred<-fitted(object,type,...)
              if(length(xpred)==0) {
                opred<-saemix.predict(object)
                assign(namObj,opred,envir=parent.frame()) # update object invisibly with predictions - does not work, pb with environment...
                #                print(parent.frame())
                #                print(namObj)
                xpred<-fitted(opred,type,...)
              }
              #              print(xpred)
            } else {
              # If individual predictions are requested, check the data contains a column with observations
              if(type %in% c("ipred","icpred")) {
                if(length(grep(saemix.data["name.response"],colnames(newdata)))==0) {
                  if(object["options"]$warnings) cat("Observed",saemix.data["name.response"],"in the new subjects are required to estimate individual parameters. The population predictions will be computed instead.\n")
                  type<-"ppred"
                }
              }
              has.y<-grep(saemix.data["name.response"],colnames(newdata))
              if(length(has.y)==0 & saemix.model["modeltype"]=="likelihood") {
                if(object["options"]$warnings) cat("Please provide values of the response to obtain predictions for a model defined by loglikelihood\n")
                return()
              }
              listpred<-saemixPredictNewdata(object,newdata,type=type)
              xpred<-listpred$predictions[,type]
            }
            return(xpred)
          }
)

#' Compute model predictions after an saemix fit
#' 
#' In nonlinear mixed effect models, different types of predictions may be obtained, including individual predictions and population predictions.
#' This function takes an SaemixObject and adds any missing predictions for maximum a posteriori and conditional mean estimations of the individual
#' parameters, and for the different types of individual and population predictions for the response variable.
#' 
#' @param object an SaemixObject object
#' @param type the type of predictions (ipred= individual, ppred=population predictions obtained with the population estimates, 
#' ypred=mean of the population predictions, icpred=conditional predictions). 
#' By default, computes all the predictions and residuals, along with the corresponding parameter estimates
#' 
#' @return an updated SaemixObject object
#' 
#' @details This function is used internally by saemix to automatically compute a number of elements needed for diagnostic plots.
#' It is normally executed directly during a call to saemix() but can be called to add residuals
#' 
#' @keywords methods
#' @export saemix.predict

saemix.predict<-function(object, type=c("ipred", "ypred", "ppred", "icpred")) {
  type<-intersect(type, c("ipred", "ypred", "ppred", "icpred"))
  # Estimate individual parameters
  if(length(object["results"]["map.psi"])==0 & ("ipred" %in% type))
    object<-map.saemix(object)
  # en principe n'arrive jamais car on les calcule pdt le fit...
  if(length(object["results"]["cond.mean.phi"])==0 & ("icpred" %in% type))
    object<-conddist.saemix(object)
  saemix.res<-object["results"]
  
  # Population predictions using the population parameters [ f(mu) ] - always computed
  xind<-object["data"]["data"][,c(object["data"]["name.predictors"],object["data"]["name.cens"],object["data"]["name.mdv"],object["data"]["name.ytype"]),drop=FALSE]
  # If exponential model, this is the transformed data
  yobs<-object["data"]["data"][,object["data"]["name.response"]]
  index<-object["data"]["data"][,"index"]
  psiM<-transphi(saemix.res["mean.phi"],object["model"]["transform.par"])
  ppred<-object["model"]["model"](psiM,index,xind)
  saemix.res["ppred"]<-unname(ppred)
  if(length(saemix.res["predictions"])==0)
    saemix.res["predictions"]<-data.frame(ppred=ppred) # create dataframe for predictions if not yet available
  
  # Compute predictions and residuals
  if(length(intersect(c("ipred","icpred"),type))>0) {
    # Individual predictions
    pres<-saemix.res["respar"]
    if("ipred" %in% type) {
      ipred<-object["model"]["model"](saemix.res["map.psi"],index,xind) # Predictions with MAP
      ires<-yobs-ipred
      # Individual weighted residuals
      gpred<-error(ipred,pres,xind$ytype)
      iwres<-(yobs-ipred)/gpred
      saemix.res["ipred"]<-ipred
      saemix.res["ires"]<-ires
      saemix.res["iwres"]<-iwres
      saemix.res["predictions"]$ipred<-ipred
      saemix.res["predictions"]$ires<-ires
      saemix.res["predictions"]$iwres<-iwres
    } 
    if("icpred" %in% type) {
      psiM<-transphi(saemix.res["cond.mean.phi"],object["model"]["transform.par"])
      icond.pred<-object["model"]["model"](psiM,index,xind) # Predictions with Conditional mean
      gpred<-error(icond.pred,pres,xind$ytype)
      icwres<-(yobs-icond.pred)/gpred
      saemix.res["icpred"]<-icond.pred
      saemix.res["icwres"]<-icwres
      saemix.res["predictions"]$icpred<-icond.pred
      saemix.res["predictions"]$icwres<-icwres
    }
  } 
  
  # Population weighted residuals, npde: needs the individual variance-covariance matrix => use compute.sres to estimate these by simulations
  object["results"]<-saemix.res
  return(object)
}

####################################################################################
####			SaemixObject class - method to plot			####
####################################################################################
#' General plot function from SAEM
#' 
#' Several plots (selectable by the type argument) are currently available:
#' convergence plot, individual plots, predictions versus observations,
#' distribution plots, VPC, residual plots, mirror.
#' 
#' This is the generic plot function for an SaemixObject object, which
#' implements different graphs related to the algorithm (convergence plots,
#' likelihood estimation) as well as diagnostic graphs. A description is
#' provided in the PDF documentation. Arguments such as main, xlab, etc... that
#' can be given to the generic plot function may be used, and will be
#' interpreted according to the type of plot that is to be drawn.
#' 
#' A special argument plot.type can be set to determine the type of plot; it
#' can be one of: 
#' \describe{ 
#' \item{data:}{A spaghetti plot of the data,displaying the observed data y 
#' as a function of the regression variable (time for a PK application)} 
#' \item{convergence:}{For each parameter in the model, this plot shows the 
#' evolution of the parameter estimate versus the iteration number} 
#' \item{likelihood:}{Graph showing the evolution of the
#' log-likelihood during the estimation by importance sampling}
#' \item{observations.vs.predictions:}{Plot of the predictions computed with
#' the population parameters versus the observations (left), and plot of the
#' predictions computed with the individual parameters versus the observations (right)} 
#' \item{residuals.scatter:}{Scatterplot of the residuals versus the
#' predictor (top) and versus predictions (bottom), for weighted residuals
#' (population residuals, left), individual weighted residuals (middle) and npde (right).} 
#' \item{residuals.distribution:}{Distribution of the residuals,
#' plotted as histogram (top) and as a QQ-plot (bottom), for weighted residuals
#' (population residuals, left), individual weighted residuals (middle) and npde (right).} 
#' \item{individual.fit:}{Individual fits are obtained using the
#' individual parameters with the individual covariates}
#' \item{population.fit:}{Population fits are obtained using the population
#' parameters with the individual covariates} 
#' \item{both.fit:}{Individual fits, superposing fits obtained using the population 
#' parameters with the individual covariates (red) and using the individual parameters 
#' with the individual covariates (green)} 
#' \item{mirror:}{Mirror plots assessing the compatibility of simulated data compared to the original}
#' \item{marginal.distribution:}{Distribution of
#' the parameters (conditional on covariates when some are included in the
#' model). A histogram of individual parameter estimates can be overlayed on
#' the plot, but it should be noted that the histogram does not make sense when
#' there are covariates influencing the parameters and a warning will be
#' displayed} 
#' \item{random.effects:}{Boxplot of the random effects}
#' \item{correlations:}{Correlation between the random effects}
#' \item{parameters.vs.covariates:}{Plots of the estimates of the individual
#' parameters versus the covariates, using scatterplot for continuous
#' covariates, boxplot for categorical covariates}
#' \item{randeff.vs.covariates:}{Plots of the estimates of the random effects
#' versus the covariates, using scatterplot for continuous covariates, boxplot
#' for categorical covariates} 
#' \item{npde:}{Plots 4 graphs to evaluate the shape of the distribution of the 
#' normalised prediction distribution errors (npde)} 
#' \item{vpc:}{Visual Predictive Check, with options to include the
#' prediction intervals around the boundaries of the selected interval as well
#' as around the median (50th percentile of the simulated data).} 
#' } 
#' In addition, the following values for plot.type produce a series of plots:
#' \describe{ 
#' \item{reduced:}{ produces the following plots: plot of the data,
#' convergence plots, plot of the likelihood by importance sampling (if
#' computed), plots of observations versus predictions. This is the default
#' behaviour of the plot function applied to an SaemixObject object}
#' \item{full:}{ produces the following plots: plot of the data, convergence
#' plots, plot of the likelihood by importance sampling (if computed), plots of
#' observations versus predictions, scatterplots and distribution of residuals,
#' VPC, npde, boxplot of the random effects, distribution of the parameters,
#' correlations between random effects, plots of the relationships between
#' individually estimated parameters and covariates, plots of the relationships
#' between individually estimated random effects and covariates}
#' 
#' Each plot can be customised by modifying options, either through a list of
#' options set by the \code{\link{saemix.plot.setoptions}} function, or on the
#' fly by passing an option in the call to the plot (see examples). 
#' }
#' 
#' @aliases plot.saemix plot,SaemixObject plot
#' @aliases plotnpde
#' @param x an object returned by the \code{\link{saemix}} function
#' @param y empty
#' @param \dots optional arguments passed to the plots
#' @return None
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{SaemixObject}},\code{\link{saemix}},
#' \code{\link{saemix.plot.setoptions}}, \code{\link{saemix.plot.select}},
#' \code{\link{saemix.plot.data}}
#' @references E Comets, A Lavenu, M Lavielle M (2017). Parameter estimation in nonlinear mixed effect models using saemix,
#' an R implementation of the SAEM algorithm. Journal of Statistical Software, 80(3):1-41.
#' 
#' E Kuhn, M Lavielle (2005). Maximum likelihood estimation in nonlinear mixed effects models. 
#' Computational Statistics and Data Analysis, 49(4):1020-1038.
#' 
#' E Comets, A Lavenu, M Lavielle (2011). SAEMIX, an R version of the SAEM algorithm. 20th meeting of the 
#' Population Approach Group in Europe, Athens, Greece, Abstr 2173.
#' 
#' @keywords plot
#' @examples
#' 
#' data(theo.saemix)
#' 
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
#'   name.group=c("Id"),name.predictors=c("Dose","Time"),
#'   name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#'   units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
#' 
#' model1cpt<-function(psi,id,xidep) { 
#' 	  dose<-xidep[,1]
#' 	  tim<-xidep[,2]  
#' 	  ka<-psi[id,1]
#' 	  V<-psi[id,2]
#' 	  CL<-psi[id,3]
#' 	  k<-CL/V
#' 	  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
#' 	  return(ypred)
#' }
#' 
#' saemix.model<-saemixModel(model=model1cpt,
#'   description="One-compartment model with first-order absorption", 
#'   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,
#'   dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),
#'   covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
#'   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#'   omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")
#' 
#' saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
#' 
#' # Not run (strict time constraints for CRAN)
#' # saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' 
#' # Set of default plots
#' # plot(saemix.fit)
#' 
#' # Data
#' # plot(saemix.fit,plot.type="data")
#' 
#' # Convergence
#' # plot(saemix.fit,plot.type="convergence")
#' 
#' # Individual plot for subject 1, smoothed
#' # plot(saemix.fit,plot.type="individual.fit",ilist=1,smooth=TRUE)
#' 
#' # Individual plot for subject 1 to 12, with ask set to TRUE 
#' # (the system will pause before a new graph is produced)
#' # plot(saemix.fit,plot.type="individual.fit",ilist=c(1:12),ask=TRUE)
#' 
#' # Diagnostic plot: observations versus population predictions
#' # par(mfrow=c(1,1))
#' # plot(saemix.fit,plot.type="observations.vs.predictions",level=0,new=FALSE)
#' 
#' # LL by Importance Sampling
#' # plot(saemix.fit,plot.type="likelihood")
#' 
#' # Scatter plot of residuals
#' # Data will be simulated to compute weighted residuals and npde
#' # the results shall be silently added to the object saemix.fit
#' # plot(saemix.fit,plot.type="residuals.scatter")
#' 
#' # Boxplot of random effects
#' # plot(saemix.fit,plot.type="random.effects")
#' 
#' # Relationships between parameters and covariates
#' # plot(saemix.fit,plot.type="parameters.vs.covariates")
#' 
#' # Relationships between parameters and covariates, on the same page
#' # par(mfrow=c(3,2))
#' # plot(saemix.fit,plot.type="parameters.vs.covariates",new=FALSE)
#' 
#' # VPC
#' # Not run (time constraints for CRAN)
#' # plot(saemix.fit,plot.type="vpc")
#' 
#' @exportMethod plot

setMethod(f="plot",
  signature="SaemixObject",
  def=function(x,y,...) {
    args1<-match.call(expand.dots=TRUE)
    i1<-match("plot.type",names(args1))
    if(!is.na(i1)) {
      plot.type<-as.character(args1[[i1]])
      plot.type<-plot.type[plot.type!="c"]
    } else plot.type<-"reduced"
#    cat("plot.type=",plot.type,"\n")
    if(plot.type[1]=="reduced") plot.type<-c("data","convergence","likelihood", "observations.vs.predictions")
    if(plot.type[1]=="full") plot.type<-c("data","convergence","likelihood", "observations.vs.predictions") #,"residuals.scatter","residuals.distribution","vpc")
    
    pltyp<-c("data","convergence","likelihood","individual.fit", "population.fit", "both.fit","observations.vs.predictions","residuals.scatter", "residuals.distribution","vpc","npde","random.effects","marginal.distribution", "correlations","parameters.vs.covariates","randeff.vs.covariates", "mirror")
    ifnd<-pmatch(plot.type,pltyp)
    if(sum(is.na(ifnd))>0) {
      if(x@options$warnings) message("The following plot types were not found or are ambiguous:", plot.type[is.na(ifnd)],"\n")
    }
    ifnd<-ifnd[!is.na(ifnd)]
    if(length(ifnd)==0) return()
    plot.type<-pltyp[ifnd]
    interactive<-x["prefs"]$interactive
    id.pred<-match(plot.type,c("observations.vs.predictions","individual.fit", "residuals.scatter","residuals.distribution"))
    if(x@prefs$which.poppred=="ppred") id.pred<-c(id.pred,match(plot.type, c("population.fit", "both.fit")))
    id.map<-match(plot.type,c("randeff.vs.covariates","parameters.vs.covariates"))
    id.sim<-match(plot.type, c("vpc","mirror"))
    id.res<-match(plot.type, c("npde","residuals.scatter", "residuals.distribution"))
    if(x@prefs$which.poppred=="ppred") id.sim<-c(id.sim,match(plot.type, c("population.fit", "both.fit")))
    id.pred<-id.pred[!is.na(id.pred)]
    id.sim<-id.sim[!is.na(id.sim)]
    id.map<-id.map[!is.na(id.map)]
    id.res<-id.res[!is.na(id.res)]
    namObj<-deparse(substitute(x))
#    cat(namObj,"\n")
    if(length(id.pred)>0) {
      if(length(x["results"]["ipred"])==0 | length(x["results"]["iwres"])) {
        boolpred<-TRUE
        if(interactive) {
      cok<-readline(prompt="Computations will be performed to obtain model predictions, proceed ? (y/Y) [default=yes] ")
        if(!cok %in% c("y","Y","yes","")) return()
        }
        if(boolpred) {
          x<-saemix.predict(x)
          assign(namObj,x,envir=parent.frame())
        }
      }
    }
    if(length(id.sim)>0) {
      if(length(x["sim.data"]["N"])==0 || x["sim.data"]["nsim"]==0) {
        boolpred<-TRUE
        if(interactive) {
        cok<-readline(prompt="Simulations will be performed. This might take a while, proceed ? (y/Y) [default=yes] ")
        if(!cok %in% c("y","Y","yes","")) return()
        } else {
        	if(x@options$warnings) message("Performing simulations under the model.\n")
        }
        if(boolpred) {
          x<-simulate(x)
          assign(namObj,x,envir=parent.frame())
        }
      }
    }
    if(length(id.res)>0) {
      if(length(x["results"]["npde"])==0) {
        boolpred<-TRUE
        if(interactive) {
        cok<-readline(prompt="Simulations will be performed to obtain residuals, VPC and npde. This might take a while, proceed ? (y/Y) [default=yes] ")
        if(!cok %in% c("y","Y","yes","")) return()
        }
        if(boolpred) {
          x<-compute.sres(x)
          assign(namObj,x,envir=parent.frame())
        }
      }
    }
    if(length(id.map)>0) {
      if(length(x["results"]["map.eta"])==0) {
        if(x@options$warnings) message("Computing ETA estimates and adding them to fitted object.\n")
      	x<-compute.eta.map(x)
        assign(namObj,x,envir=parent.frame())
      }
    }
    for(ipl in plot.type) {
      switch (EXPR=ipl,
    "data"={
       if(x@options$warnings) message("Plotting the data\n")
       saemix.plot.data(x,...)
    },
    "convergence"={
       if(x@options$warnings) message("Plotting convergence plots\n")
       saemix.plot.convergence(x,...)
    },
    "likelihood"={  
       if(x@options$warnings) message("Plotting the likelihood\n")
       saemix.plot.llis(x,...)
    },
    "observations.vs.predictions"={
       if(x@options$warnings) message("Plotting observations versus predictions\n")      
       saemix.plot.obsvspred(x,...)
    },
    "individual.fit"={
      if(x@options$warnings) message("Plotting individual fits\n")
      saemix.plot.fits(x,...)
    },
    "population.fit"={
      if(x@options$warnings) message("Plotting fits obtained with population predictions\n")
      saemix.plot.fits(x,level=0,...)
    },
    "both.fit"={
      if(x@options$warnings) message("Plotting the fits overlaying individual and population predictions\n")
      saemix.plot.fits(x,level=c(0,1),...)
    },
    "individual.fit"={
      if(x@options$warnings) message("Mirror plots\n")
      saemix.plot.mirror(x,...)
    },
    "residuals.scatter"={
      message("Please use npdeSaemix to obtain VPC and npde\n")
#      if(x@options$warnings) message("Plotting scatterplots of residuals\n")
#      saemix.plot.scatterresiduals(x,...)
    },
    "residuals.distribution"={
      message("Please use npdeSaemix to obtain VPC and npde\n")
#      if(x@options$warnings) message("Plotting the distribution of residuals\n")
#      saemix.plot.distribresiduals(x,...)
    },
    "random.effects"={
      saemix.plot.randeff(x,...)
    },
    "correlations"={
      if(length(x@model@indx.omega>1)) saemix.plot.correlations(x,...)
    },
    "parameters.vs.covariates"={
      if(length(x@data@name.covariates)==0) {
        if(x@options$warnings) message("No covariates in the dataset\n")
        return()
      } else saemix.plot.parcov(x,...)
    },
    "randeff.vs.covariates"={
      if(length(x@data@name.covariates)==0) {
        if(x@options$warnings) message("No covariates in the dataset\n")
        return()
      } else saemix.plot.randeffcov(x,...)
    },
    "marginal.distribution"={
      saemix.plot.distpsi(x,...)
    },
    "vpc"={
      if(x@options$warnings) message("Direct call to VPC will soon be deprecated, please use npdeSaemix for VPC\n")
      saemix.plot.vpc(x,...)
    },
    "npde"={
      message("Please use npdeSaemix to obtain VPC and npde\n")
#      if(x@options$warnings) message("Plotting npde\n")
#      saemix.plot.npde(x,...)
    },
    cat("Plot ",ipl," not implemented yet\n")
     )
   }
  }
)

####################################################################################
####			Likelihood and tests		####
####################################################################################

#' Extract likelihood from a saemixObject resulting from a call to saemix
#'
#' The likelihood in saemix can be computed by one of three methods: linearisation (linearisation of the model), importance sampling (stochastic integration) and gaussian quadrature (numerical integration). The linearised likelihood is obtained as a byproduct of the computation of the Fisher Information Matrix (argument FIM=TRUE in the options given to the saemix function).
#' If no method argument is given, this function will attempt to extract the likelihood computed by importance sampling (method="is"), unless the object contains the likelihood computed by linearisation, in which case the function will extract this component instead.
#' If the requested likelihood is not present in the object, it will be computed and aded to the object before returning.
#' 
#' BIC.covariate implements the computation of the BIC from Delattre et al. 2014.
#'
#' @name logLik
#' @aliases AIC.SaemixObject BIC.SaemixObject logLik.SaemixObject BIC.covariate
#'
#' @param object name of an SaemixObject object
#' @param method character string, one of c("is","lin","gq"), to select one of the available approximations to the log-likelihood (is: Importance Sampling; lin: linearisation and gq: Gaussian Quadrature). See documentation for details
#' @param ... additional arguments
#' @param k numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC
#' @return Returns the selected statistical criterion (log-likelihood, AIC, BIC) extracted from the SaemixObject, computed with the 'method' argument if given (defaults to IS).
#' @references E Comets, A Lavenu, M Lavielle M (2017). Parameter estimation in nonlinear mixed effect models using saemix,
#' an R implementation of the SAEM algorithm. Journal of Statistical Software, 80(3):1-41.
#' 
#' E Kuhn, M Lavielle (2005). Maximum likelihood estimation in nonlinear mixed effects models. 
#' Computational Statistics and Data Analysis, 49(4):1020-1038.
#' 
#' E Comets, A Lavenu, M Lavielle (2011). SAEMIX, an R version of the SAEM algorithm. 20th meeting of the 
#' Population Approach Group in Europe, Athens, Greece, Abstr 2173.
#' 
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}
#' @author Audrey Lavenu
#' @author Marc Lavielle
#' @author Maud Delattre
#' @seealso \code{\link{AIC}},\code{\link{BIC}}, \code{\link{saemixControl}}, \code{\link{saemix}}
#' 
#' @docType methods
#' @keywords methods
#' @export

## log-likelihood for SaemixObject objects (changes adding stop() by Johannes)
logLik.SaemixObject<- function(object, method=c("is","lin","gq"), ...) {
  method <- match.arg(method)
  
  if(method=="gq" & length(object@results@ll.gq)==0) {
    #    object<-llis.saemix(object)
    #    assign(namObj,object,envir=parent.frame())
    stop("The log-likelihood by Gaussian Quadrature has not yet been computed.")
  }
  if(method=="is" & length(object@results@ll.is)==0) {
    #    object<-llgq.saemix(object)
    #    assign(namObj,object,envir=parent.frame())
    stop("The log-likelihood by Importance Sampling has not yet been computed.")
  }
  if(method=="lin" & length(object@results@ll.lin)==0) {
    #    object<-fim.saemix(object)
    #    assign(namObj,object,envir=parent.frame())
    stop("The log-likelihood by linearisation has not yet been computed.")
  }
  val<-switch(method,is=object@results@ll.is,lin=object@results@ll.lin, gq=object@results@ll.gq)
  attr(val, "nall") <- object@data@N
  attr(val, "nobs") <- object@data@ntot.obs
  attr(val, "df") <- object@results@npar.est
  class(val) <- "logLik"
  val
}

#' @export 
#' @rdname logLik

AIC.SaemixObject<-function(object, method=c("is","lin","gq"), ..., k=2) {
  method <- match.arg(method)

  if(method=="is" & length(object@results@ll.is)==0) {
    stop("The log-likelihood by Importance Sampling has not yet been computed.")
  }
  if(method=="gq" & length(object@results@ll.gq)==0) {
    stop("The log-likelihood by Gaussian Quadrature has not yet been computed.")
  }
  if(method=="lin" & length(object@results@ll.lin)==0) {
    stop("The log-likelihood by linearisation has not yet been computed.")
  }
  val<-switch(method,is=object@results@aic.is,lin=object@results@aic.lin, gq=object@results@aic.gq)
  val
}

#' @export
#' @rdname logLik

BIC.SaemixObject<-function(object, method=c("is","lin","gq"), ...) {
  method <- match.arg(method)
#  -2 * as.numeric(object) + attr(object, "df") * log(nobs(object))

  if(method=="is" & length(object@results@ll.is)==0) {
    stop("The log-likelihood by Importance Sampling has not yet been computed.")
  }
  if(method=="gq" & length(object@results@ll.gq)==0) {
    stop("The log-likelihood by Gaussian Quadrature has not yet been computed.")
  }
  if(method=="lin" & length(object@results@ll.lin)==0) {
    stop("The log-likelihood by linearisation has not yet been computed.")
  }
  val<-switch(method,is=object@results@bic.is,lin=object@results@bic.lin, gq=object@results@bic.gq)
  val
}

#' @export
#' @rdname logLik

BIC.covariate<-function(object, method=c("is","lin","gq"), ...) {
  method <- match.arg(method)
  #  -2 * as.numeric(object) + attr(object, "df") * log(nobs(object))

  if(method=="is" & length(object@results@ll.is)==0) {
    stop("The log-likelihood by Importance Sampling has not yet been computed.")
  }
  if(method=="gq" & length(object@results@ll.gq)==0) {
    stop("The log-likelihood by Gaussian Quadrature has not yet been computed.")
  }
  if(method=="lin" & length(object@results@ll.lin)==0) {
    stop("The log-likelihood by linearisation has not yet been computed.")
  }
  val<-switch(method,is=object@results@bic.covariate.is,lin=object@results@bic.covariate.lin, gq=object@results@bic.covariate.gq)
  val
}

# Log-likelihood ratio test, given 2 models
# ECO TODO

####################################################################################
####			saemixObject class - accesseurs parametres		####
####################################################################################
#' @rdname psi-methods
#' @aliases psi
#' @exportMethod psi

# Extract individual parameter estimates (psi_i)
setMethod("psi","SaemixObject",
  function(object,type=c("mode","mean")) {
    type <- match.arg(type)
    namObj<-deparse(substitute(object))
    if(type=="mode") { # MAP
      if(length(object@results@map.psi)==0) {
        object<-map.saemix(object)
        assign(namObj,object,envir=parent.frame())
      }
      psi<-object@results@map.psi
    } else { # conditional means
      if(length(object@results@cond.mean.psi)==0) {
        object<-conddist.saemix(object)
        assign(namObj,object,envir=parent.frame())
      }
      psi<-object@results@cond.mean.psi
    }
    return(psi)
  }
)

#' @rdname psi-methods
#' @exportMethod phi

# Extract individual parameter estimates on non-transformed scale (phi_i)
setMethod("phi","SaemixObject",
  function(object,type=c("mode","mean")) {
    type <- match.arg(type)
    namObj<-deparse(substitute(object))
    if(type=="mode") { # MAP
      if(length(object@results@map.phi)==0) {
        object<-map.saemix(object)
        assign(namObj,object,envir=parent.frame())
      }
      phi<-object@results@map.phi
    } else { # conditional means
      if(length(object@results@cond.mean.phi)==0) {
        object<-conddist.saemix(object)
        assign(namObj,object,envir=parent.frame())
      }
      phi<-object@results@cond.mean.phi
    }
    return(phi)
  }
)

# Extract individual estimates of random effects (eta_i)
#' @rdname psi-methods
#' @exportMethod eta

setMethod("eta","SaemixObject",
  function(object,type=c("mode","mean")) {
    type <- match.arg(type)
    namObj<-deparse(substitute(object))
#    cat("Nom objet",namObj,"\n")
    if(type=="mode") { # MAP
      if(length(object@results@map.eta)==0) {
        object<-compute.eta.map(object)
        assign(namObj,object,envir=parent.frame())
      }
      eta<-object@results@map.eta
    } else { # conditional means
      if(length(object@results@cond.mean.eta)==0) {
        object<-conddist.saemix(object)
        assign(namObj,object,envir=parent.frame())
      }
      eta<-object@results@cond.mean.eta
    }
    return(eta)
  }
)

#' Extract coefficients from a saemix fit
#' 
#' @name coef.saemix
#' 
#' @param object a SaemixObject
#' @param ... further arguments to be passed to or from other methods
#' @return a list with 3 components:
#' \describe{
#' \item{fixed}{fixed effects}
#' \item{population}{a list of population parameters with two elements, a matrix containing the untransformed parameters psi and a matrix containing the transformed parameters phi}
#' \item{individual}{a list of individual parameters with two elements, a matrix containing the untransformed parameters psi and a matrix containing the transformed parameters phi}
#' }
#' 
#' @aliases coef coef.SaemixObject
#' @aliases coef,SaemixObject coef,SaemixObject-method
#' @export

coef.SaemixObject<-function(object, ...) {
#    if(missing(level)) level<-1
    pfix<-object@results@fixed.effects[object@results@indx.fix]
    names(pfix)<-object@results@name.fixed[object@results@indx.fix]
#c(object@results@fixed.effects,object@name.sigma[object@indx.res])
#    names(pfix)<-c(object@results@name.fixed,object@name.sigma[object@indx.res])    
    pop.phi<-object@results@mean.phi
    pop.psi<-transphi(pop.phi,object@model@transform.par)
    ind.psi<-list(map=object@results@map.psi, cond=object@results@cond.mean.psi)
    ind.phi<-list(map=object@results@map.phi, cond=object@results@cond.mean.phi)
    eta<-list(map=object@results@map.eta,cond=object@results@cond.mean.eta)
    colnames(pop.phi)<-colnames(pop.psi)<-colnames(ind.phi$map)<-names(pfix)
    colnames(ind.psi$map)<-colnames(ind.phi$cond)<-colnames(ind.psi$cond)<-names(pfix)
    coef<-list(fixed=pfix,population=list(psi=pop.psi,phi=pop.phi), individual=list(psi=ind.psi,phi=ind.phi,eta=eta))
    return(coef)
  }

#' @rdname resid.saemix
#' @export 

# Extract residuals and fitted
resid.SaemixObject<-function (object, type = c("ires", "wres", "npde", "pd", "iwres", "icwres"), ...) {
  type <- match.arg(type)
  resid(object@results, type=type, ...)
}

#' @rdname fitted.saemix
#' @export 

fitted.SaemixObject<-function (object, type = c("ipred", "ypred", "ppred", "icpred"), ...) 
          {
            type <- match.arg(type)
            fitted(object@results, type=type, ...)
}


####################################################################################
####            SaemixObject class - extract variance-covariance matrix         ####
####################################################################################

#' @rdname vcov
#' @export

vcov.SaemixObject<-function(object, ...) {
  vcov(object@results)
}

####################################################################################
####			saemixObject class - Replace the data component in an saemixObject 		####
####################################################################################

#' Replace the data element in a SaemixObject object
#' 
#' Returns an SaemixObject object where the data object has been replaced by the data provided in a dataframe
#'
#' @name replaceData
#' @aliases replaceData-methods replaceData.saemixObject
#' 
#' @param saemixObject an SaemixObject object
#' @param newdata a dataframe containing data 
#' @return an object of class \code{"\linkS4class{SaemixObject}"}. The population parameters are retained but all the predictions, individual parameters and statistical criteria are removed. The function attempts to extract the elements entering the statistical model (subject id, predictors, covariates and response).
#' @examples 
#' # TODO
#' @keywords methods
#' @export 

replaceData.saemixObject<-function(saemixObject, newdata) {
  # Takes a saemixObject fit, changes the data object and removes all results pertaining to individual parameters and LL
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  orig.data<-saemixObject["data"]
  ipb<-numeric(0)
  if(is.na(match(orig.data["name.group"],colnames(newdata)))) ipb<-c(ipb,1) else ipb<-c(ipb,0)
  if(sum(is.na(match(orig.data["name.predictors"],colnames(newdata))))>0) ipb<-c(ipb,1) else ipb<-c(ipb,0)
  if(is.na(match(orig.data["name.response"],colnames(newdata)))) ipb<-c(ipb,1) else ipb<-c(ipb,0)
  if(sum(is.na(match(rownames(saemixObject["model"]["covariate.model"]),colnames(newdata))))>0) ipb<-c(ipb,1) else ipb<-c(ipb,0)
  if(ipb[1]>0) cat("The new dataset needs to include a group column named",orig.data["name.group"],"\n")
  if(ipb[2]>0) cat("The new dataset needs to include the following predictors:",orig.data["name.predictors"],"\n")
  if(sum(ipb[1:2])>0) return()  
  iflag<-0
  if(ipb[3]>0) {
    cat("No response column named",orig.data["name.response"],"given. Only population predictions can be obtained for the new data.\n")
    iflag<-1
  }
  if(ipb[4]>0) cat("The following covariates are present in the model but not in the new dataset",orig.data["name.covariates"],"\nTheir value will be set to the median value (for continuous) or the mode (for categorial) of the covariate in the original dataset.")
  if(ipb[3]>0) {
    newdata[orig.data["name.response"]]<-0
  }
  if(ipb[4]>0) {
    idx<-setdiff(rownames(saemixObject["model"]["covariate.model"]),colnames(newdata))
    ocov<-orig.data["data"][!duplicated(orig.data["data"][orig.data["name.group"]]), idx,drop=F]
    for(icov in idx) {
      if(typeof(ocov[icov])=="double") newdata[icov]<-median(ocov[,icov]) else {
        newdata[icov]<-getmode(ocov[,icov])
      }
    }
  }
  tempname<-tempfile()
  write.table(newdata,tempname,quote=F,col.names=T, sep=";")
  saemix.newdata<-saemixData(name.data=tempname, name.group=orig.data["name.group"], name.response =orig.data["name.response"], name.predictors=orig.data["name.predictors"], name.covariates=rownames(saemixObject["model"]["covariate.model"]), units=orig.data["units"], name.X=orig.data["name.X"],verbose=FALSE, sep=";")
  if(iflag==1) saemix.newdata["data"][,orig.data["name.response"]]<-NA
  
  saemix.newObj<-saemixObject
  saemix.newObj["data"]<-saemix.newdata
  saemix.newObj["results"]["cond.mean.phi"] <-matrix(nrow=0,ncol=0)
  saemix.newObj["results"]["cond.mean.psi"] <-matrix(nrow=0,ncol=0)
  saemix.newObj["results"]["cond.var.phi"] <-matrix(nrow=0,ncol=0)
  saemix.newObj["results"]["cond.mean.eta"] <-matrix(nrow=0,ncol=0)
  saemix.newObj["results"]["cond.shrinkage"] <-numeric(0)
  saemix.newObj["results"]["mean.phi"] <-matrix(nrow=0,ncol=0)
  saemix.newObj["results"]["map.psi"] <-data.frame()
  saemix.newObj["results"]["map.phi"] <-data.frame()
  saemix.newObj["results"]["map.eta"] <-matrix(nrow=0,ncol=0)
  saemix.newObj["results"]["map.shrinkage"] <-numeric(0)
  saemix.newObj["results"]["phi"] <-matrix(nrow=0,ncol=0)
  saemix.newObj["results"]["psi.samp"] <-array(dim=c(0,0,0),data=0)
  saemix.newObj["results"]["phi.samp"] <-array(dim=c(0,0,0),data=0)
  saemix.newObj["results"]["phi.samp.var"] <-array(dim=c(0,0,0),data=0)
  saemix.newObj["results"]["ll.lin"] <-numeric(0)
  saemix.newObj["results"]["aic.lin"] <-numeric(0)
  saemix.newObj["results"]["bic.lin"] <-numeric(0)
  saemix.newObj["results"]["bic.covariate.lin"] <-numeric(0)
  saemix.newObj["results"]["ll.gq"] <-numeric(0)
  saemix.newObj["results"]["aic.gq"] <-numeric(0)
  saemix.newObj["results"]["bic.gq"] <-numeric(0)
  saemix.newObj["results"]["bic.covariate.gq"] <-numeric(0)
  saemix.newObj["results"]["ll.is"] <-numeric(0)
  saemix.newObj["results"]["aic.is"] <-numeric(0)
  saemix.newObj["results"]["bic.is"] <-numeric(0)
  saemix.newObj["results"]["bic.covariate.is"] <-numeric(0)
  saemix.newObj["results"]["predictions"] <-data.frame()
  saemix.newObj["results"]["ypred"] <-numeric(0)
  saemix.newObj["results"]["ppred"] <-numeric(0)
  saemix.newObj["results"]["ipred"] <-numeric(0)
  saemix.newObj["results"]["icpred"] <-numeric(0)
  saemix.newObj["results"]["ires"] <-numeric(0)
  saemix.newObj["results"]["iwres"] <-numeric(0)
  saemix.newObj["results"]["icwres"] <-numeric(0)
  saemix.newObj["results"]["wres"] <-numeric(0)
  saemix.newObj["results"]["npde"] <-numeric(0)
  saemix.newObj["results"]["pd"] <-numeric(0)
  return(saemix.newObj)
}


####################################################################################

#' Create saemix objects with only data filled in
#' 
#' Create saemix objects either with empty results or with parameters set by the user. 
#' This is an internal function used as a preliminary step to obtain predictions for new data
#' and is not intended to be used directly.
#'
#' @name createSaemixObject
#' @aliases createSaemixObject.empty createSaemixObject.initial
#' 
#' @param data an saemixData object
#' @param model an saemixModel object
#' @param control a list of options (if empty, will be set to the default values by \code{saemixControl})
#' 
#' @return an object of class \code{"\linkS4class{SaemixObject}"}. 
#' 
#' @details with createSaemixObject.empty, the data component is set to the data object, the model component is set to the model object, and the result component is empty
#' 
#' @details with createSaemixObject.initial, the data and model are set as with createSaemixObject.empty, 
#' but the population parameter estimates are used to initialise the result component
#' as in the initialisation step of the algorithm (initialiseMainAlgo)
#' 
#' @examples 
#' # TODO
#' @keywords methods
#' @export 


createSaemixObject.empty<-function(model,data,control=list()) {
  if(!is(model,"SaemixModel")) {
    message("Please provide a valid model object (see the help page for SaemixModel)\n")
    return("Need a valid model object")
  }
  if(!is(data,"SaemixData")) {
    message("Please provide a valid data object (see the help page for SaemixData)\n")
    return("Need a valid data object")
  }
  
  saemixObject<-new(Class="SaemixObject",data=data,model=model,options=control)
  #  saemixObject<-new(Class="SaemixObject",data=saemix.data, model=saemix.model,options=saemix.options)
  opt.warn<-getOption("warn")
  if(!saemixObject["options"]$warnings) options(warn=-1)
  
  saemix.data<-saemixObject["data"]
  saemix.data@ocov<-saemix.data@ocov[saemix.data@data[,"mdv"]==0,,drop=FALSE]
  saemix.data@data<-saemix.data@data[saemix.data@data[,"mdv"]==0,]
  saemix.data@ntot.obs<-dim(saemix.data@data)[1]
  saemixObject["data"]<-saemix.data
  
  return(saemixObject)
}

createSaemixObject.initial<-function(model,data,control=list()) {
  # Checking validity of input
  if(!is(model,"SaemixModel")) {
    message("Please provide a valid model object (see the help page for SaemixModel)\n")
    return("Need a valid model object")
  }
  if(!is(data,"SaemixData")) {
    message("Please provide a valid data object (see the help page for SaemixData)\n")
    return("Need a valid data object")
  }
  
  # Creating saemixObject (empty results)
  saemixObject<-new(Class="SaemixObject",data=data,model=model,options=control)
  opt.warn<-getOption("warn")
  if(!saemixObject["options"]$warnings) options(warn=-1)
  
  saemix.options<-saemixObject["options"]
  saemix.model<-saemixObject["model"]
  saemix.data<-saemixObject["data"]
  saemix.data@ocov<-saemix.data@ocov[saemix.data@data[,"mdv"]==0,,drop=FALSE]
  saemix.data@data<-saemix.data@data[saemix.data@data[,"mdv"]==0,]
  saemix.data@ntot.obs<-dim(saemix.data@data)[1]
  
  # Initialising results component to initial estimates
  xinit<-initialiseMainAlgo(saemix.data,saemix.model,saemix.options)
  saemix.model<-xinit$saemix.model
  Uargs<-xinit$Uargs
  varList<-xinit$varList
  
  if(saemix.model["modeltype"]=="structural"){
    xres1<-new(Class="SaemixRes",modeltype="structural",status="initial",
           name.fixed=saemix.model["name.fixed"], name.random=saemix.model["name.random"],name.sigma=saemix.model["name.sigma"],
             fixed.effects=saemix.model@psi0[saemix.model@betaest.model==1],
             fixed.psi=xinit$fixedpsi.ini,
             betaC=xinit$betas[xinit$Uargs$indx.betaC],betas=xinit$betas,
             omega=varList$omega,respar=varList$pres,MCOV=varList$MCOV)
  } else{
    xres1<-new(Class="SaemixRes",modeltype=saemix.model["modeltype"],status="initial",
           name.fixed=saemix.model["name.fixed"], name.random=saemix.model["name.random"],name.sigma=saemix.model["name.sigma"],
           fixed.effects=saemix.model@psi0[saemix.model@betaest.model==1],
           fixed.psi=xinit$fixedpsi.ini,
           betaC=xinit$betas[xinit$Uargs$indx.betaC],betas=xinit$betas,
           omega=varList$omega,MCOV=varList$MCOV)
  }

  xres1@indx.cov<-saemix.model@indx.cov
  xres1@indx.res<-saemix.model@indx.res
  xres1@indx.fix<-saemix.model@indx.fix
  xres1@indx.omega<-saemix.model@indx.omega
  saemixObject["results"]<-xres1
  
  return(saemixObject)
}

####################################################################################

