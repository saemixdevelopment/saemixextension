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
  if(length(nbiter.mcmc)<4) {
    length(nbiter.mcmc)<-4
    nbiter.mcmc[is.na(nbiter.mcmc)]<-2
    nbiter.mcmc[4]<-0
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
