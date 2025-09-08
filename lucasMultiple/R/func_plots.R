#######################	   Default plot options (list)	 ########################

#' Function setting the default options for the plots in SAEM
#' 
#' This function can be used to create a list containing the default options
#' and arguments used by the plot functions.
#' 
#' A more detailed description of the options set via these lists is provided
#' in the PDF documentation. The "replace" functions are helper functions used
#' within the plot functions. saemix.plot.setoptions has more available options
#' than saemix.data.setoptions since it applies to all possible plots while the
#' latter only applies to data. 
#' @description 
#' \itemize{
#' \item ablinecol Color of the lines added to the plots (default: "DarkRed")
#' \item ablinelty Type of the lines added to the plots. Defaults to 2 (dashed line)
#' \item ablinelwd Width of the lines added to the plots (default: 2)
#' \item ask A logical value. If TRUE, users will be prompted before each new plot. Defaults to FALSE
#' \item cex A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default. Defaults to 1 (no magnification)
#' \item cex.axis Magnification to be used for axis annotation relative to the current setting of 'cex'. Defaults to 1 (no magnification)
#' \item cex.main Magnification to be used for main titles relative to the current setting of 'cex'. Defaults to 1 (no magnification)
#' \item cex.lab Magnification to be used for x and y labels relative to the current setting of 'cex'. Defaults to 1 (no magnification)
#' \item col.fillmed For the VPC plots: color filling the prediction interval for the median. Defaults to "pink"
#'\item col.fillpi For the VPC plots: color filling the prediction interval for the limits of the prediction interval. Defaults to "slategray1"
#' \item col.lmed For the VPC plots: color of the line showing the median of the simulated data. Defaults to "indianred4"
#' \item col.lobs For the VPC plots: color of the lines showing the median, 2.5 and 97.5th percentiles (for a 95% prediction interval). Defaults to "steelblue4"
#' \item col.lpi For the VPC plots: color of the line showing the boundaries of the prediction intervals. Defaults to "slategray4"
#' \item col.obs For the VPC plots: color used to plot the observations. Defaults to "steelblue4"
#' \item cov.name Name of the covariate to be used in the plots. Defaults to the first covariate in the model
#' \item cov.value Value of the covariate to be used in the plots. Defaults to NA, indicating that the median value of the covariate (for continuous covariates) or the reference category (for categorical covariates) will be used 
#' \item ilist List of indices of subjects to be included in the individual plots (defaults to all subjects)
#' \item indiv.par  a string, giving the type of the individual estimates ("map"= conditional mode, "eap"=conditional mean). Defaults to conditional mode
#' \item lcol Main line color (default: black)
#' \item line.smooth Type of smoothing when a smoothed line is used in the plot ("m": mean value, "l": linear regression; "s": natural splines). Several options may be combined, for instance "ls" will add both a linear regression line and a line representing the fit of a natural spline. Defaults to "s"
#' \item lty Line type. Defaults to 1, corresponding to a straight line
#' \item lty.lmed For the VPC plots: type of the line showing the median of the simulated data. Defaults to 2 (dashed)
#' \item lty.obs For the VPC plots: type of the line showing the observed data. Defaults to 1
#' \item lty.lpi For the VPC plots: type of the line showing the boundaries of the simulated data. Defaults to 2 (dashed)
#' \item lwd Line width (default: 1)
#' \item lwd.lmed For the VPC plots: thickness of the line showing the median of the simulated data. Defaults to 2
#' \item lwd.obs For the VPC plots: thickness of the line showing the median and boundaries of the observed data. Defaults to 2
#' \item lwd.lpi For the VPC plots: thickness of the line showing the boundaries of the simulated data. Defaults to 1
#' \item par.name Name of the parameter to be used in the plots. Defaults to the first parameter in the model
#' \item pch Symbol type. Defaults to 20, corresponding to small dots
#' \item pcol Main symbol color (default: black)
#' \item range Range (expressed in number of SD) over which to plot the marginal distribution. Defaults to 4, so that the random effects for the marginal distribution is taken over the range \[-4 SD; 4 SD\]
#' \item res.plot Type of residual plot ("res.vs.x": scatterplot versus X, "res.vs.pred": scatterplot versus predictions, "hist": histogram, "qqplot": QQ-plot) (default: "res.vs.x")
#'\item smooth When TRUE, smoothed lines are added in the plots of predictions versus observations (default: FALSE)
#' \item tit Title of the graph (default: none) 
#' \item type Type of the plot (as in the \emph{R} plot function. Defaults to "b", so that both lines and symbols are shown
#' \item units Name of the predictor used in the plots (X). Defaults to the name of the first predictor in the model (saemix.data$names$predictors\[1\])
#' \item vpc.bin Number of binning intervals when plotting the VPC (the (vpc.bin-1) breakpoints are taken as the empirical quantiles of the X data). Defaults to 10
#' \item vpc.interval Size of the prediction intervals.Defaults to 0.95 for the 95\% prediction interval
#' \item vpc.obs Should the observations be overlayed on the VPC plot. Defaults to TRUE
#' \item vpc.pi Should prediction bands be computed around the median and the bounds of the prediction intervals for the VPC. Defaults to TRUE
#' \item xlab Label for the X-axis. Defaults to the name of the X predictor followed by the unit in bracket (eg "Time (hr)")
#' \item xlim Range for the X-axis. Defaults to NA, indicating that the range is to be set by the plot function
#' \item xlog A logical value. If TRUE, a logarithmic scale is in use. Defaults to FALSE
#' \item xname Name of the predictor used in the plots (X)
#' \item ylab Label for the Y-axis. Defaults to the name of the response followed by the unit in bracket (eg "Concentration (mg/L)" (Default: none)
#' \item ylim Range for the Y-axis. Defaults to NA, indicating that the range is to be set by the plot function
#' \item ylog A logical value. If TRUE, a logarithmic scale is in use. Defaults to FALSE
#' }
#' 
#' Plotting a SaemixData object also allows the following options: 
#' \describe{ 
#' \item{individual}{ if TRUE, plots separate plots for each individual, otherwise plots a spaghetti plot of all the data. Defaults to FALSE} 
#' \item{limit}{ for individual plots, plots only a limited number of subjets (nmax). Defaults to TRUE} 
#' \item{nmax}{ for individual plots, when limit is TRUE, the maximum number of plots to produce. Defaults to 12}
#' \item{sample}{ for individual plots, if TRUE, randomly samples nmax different subjects to plot. Defaults to FALSE (the first nmax subjects are used in the plots)} 
#' } 
#' 
#' @aliases saemix.plot.setoptions saemix.data.setoptions replace.plot.options replace.data.options
#' 
#' @param saemixObject an object returned by the \code{\link{saemix}} function
#' 
#' @return A list containing the options set at their default value. This list
#' can be stored in an object and its elements modified to provide suitable graphs.
#' 
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}, Audrey Lavenu, Marc Lavielle.
#' @seealso \code{\link{SaemixObject}},\code{\link{saemix}},
#' \code{\link{saemix.plot.data}}, \code{\link{saemix.plot.convergence}},
#' \code{\link{saemix.plot.llis}}, \code{\link{saemix.plot.randeff}},
#' \code{\link{saemix.plot.obsvspred}}, \code{\link{saemix.plot.fits}},
#' \code{\link{saemix.plot.parcov}}, \code{\link{saemix.plot.distpsi}},
#' \code{\link{saemix.plot.scatterresiduals}}, \code{\link{saemix.plot.vpc}},
#' \code{\link{saemix.plot.mirror}}
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
#' # Theophylline example, after a call to fit.saemix (see examples)
#' # Not run
#' # sopt<-saemix.plot.setoptions(saemix.fit)
#' # sopt$ask<-TRUE
#' 
#' @export saemix.plot.setoptions

saemix.plot.setoptions<-function(saemixObject) {
# setting default plot options
  plot.opt<-list(
# Options for plot types
    ilist=c(1:saemixObject["data"]["N"]),
    level=0:1,
    smooth=FALSE,
    line.smooth="s",
    indiv.par="map",			# type of individual parameters
    which.par="all",			# which parameters to plot 
    which.cov="all",			# which covariates to plot 
    which.resplot=c("res.vs.x","res.vs.pred","dist.qqplot","dist.hist"), # which type of residual plots
    which.pres=c("wres","npde"),	# which population weighted residuals
    which.poppred=c("ppred"),		# which population predictions to use (ypred=E(f(theta_i)), ppred=f(population parameters))
    indiv.histo=FALSE,			# whether to include an histogram of estimated individual parameters
    cov.value=rep(NA,length(saemixObject["model"]["name.cov"])),
# General graphical options
    new=TRUE,				# whether a new page should be called
    ask=FALSE,				# whether the program should ask before creating a new page
    interactive=FALSE,			# whether the user should be prompted before computing predictions or performing simulations for VPC, npde and wres
    mfrow=c(),				# page layout (if empty, defaults to the default layout for each graph type)
    main="",				# title
    xlab="",
    ylab="",
    col="black",
    pch=20,
    lty=1,
    lwd=1,
    xlim=c(),
    ylim=c(),
    xlog=FALSE,
    ylog=FALSE,
    type="b",
    cex=1,
    cex.axis=1,
    cex.lab=1,
    cex.main=1,
    obs.pch=20,
    pred.pch=20,
    obs.col="black",
    obs.lty=1,
    obs.lwd=1,
    obs.pch=20,
    ipred.col="black",
    ipred.lty=2,
    ipred.lwd=1,
    ipred.pch=20,
    ppred.col="black",
    ppred.lty=3,
    ppred.lwd=1,
    ppred.pch=20,
    pcol="black",
    lcol="black",
    fillcol="lightblue1",
    ablinecol="DarkRed",
    ablinelty=2,
    ablinelwd=2,
# 
    range=3,
    col.fillmed="pink",
    col.fillpi="slategray1",
    col.lmed="indianred4",
    col.lpi="slategray4",
    col.pobs="steelblue4",
    col.lobs="steelblue4",
    lty.lmed=2,
    lty.lpi=2,
    lwd.lmed=2,
    lwd.lpi=1,
    lwd.lobs=2,
    lty.lobs=1,
# Options for VPC plot
    vpc.method="equal",			# method (one of "equal"=same nb of points in each interval, "width"=equally spaced intervals (on the log-scale if xlog=TRUE), "user"=user-defined breaks, "optimal"=Marc's optimal binning algorithm); for "user", the breaks must be specified in vpc.breaks (otherwise defaults back to "equal"), while for the other methods the number of bins must be specified in vpc.bin
    vpc.breaks=NULL,			# user-defined breaks
    vpc.bin=10,				# nb of bins
    vpc.beta=0.2,			# value of beta used to compute the variance-based criterion (Jopt,beta(I)) in the clustering algorithm
    vpc.lambda=0.3,			# value of lambda used in the penalised criterion to select the number of bins (if vpc.bin=NULL)
    vpc.interval=0.95,
    vpc.pi=TRUE,
    vpc.obs=TRUE)
    
     if(is.null(plot.opt$name.X)) {
        if(length(saemixObject["data"]["name.X"])>0) plot.opt$name.X<-saemixObject["data"]["name.X"] else plot.opt$name.X<-saemixObject["data"]["name.predictors"][1]
    }
    plot.opt$xlab<-paste(plot.opt$name.X," (",saemixObject["data"]["units"]$x,")", sep="")
     if(length(saemixObject["data"]["name.response"])>0)
    plot.opt$ylab<-paste(saemixObject["data"]["name.response"]," (", saemixObject["data"]["units"]$y,")",sep="")
   return(plot.opt)
}

#################    Function to supersede default plot options	 ##################

replace.plot.options<-function(plot.opt,...) {
  args1<-match.call(expand.dots=TRUE)
  if(length(args1)>2) {
# General arguments: col, pch
    i1<-match("col",names(args1))
    if(!is.na(i1)) {
      plot.opt$col<-eval(args1[[i1]])
      plot.opt$obs.col<-eval(args1[[i1]])
      plot.opt$ipred.col<-eval(args1[[i1]])
      plot.opt$ppred.col<-eval(args1[[i1]])
      plot.opt$pcol<-eval(args1[[i1]])
      plot.opt$lcol<-eval(args1[[i1]])
    }
    i1<-match("pch",names(args1))
    if(!is.na(i1)) {
      plot.opt$pch<-eval(args1[[i1]])
      plot.opt$obs.pch<-eval(args1[[i1]])
      plot.opt$ipred.pch<-eval(args1[[i1]])
      plot.opt$ppred.pch<-eval(args1[[i1]])
    }
# Other arguments
    for(i in 3:length(args1)) {
      if(match(names(args1)[i],names(plot.opt),nomatch=0)>0)    
#    plot.opt[[names(args1)[i]]]<-args1[[i]] else {
    plot.opt[[names(args1)[i]]]<-eval(args1[[i]]) else {
      if(names(args1)[i]!="plot.type") {
        if(plot.opt$interactive) cat("Argument",names(args1)[i],"not available, check spelling.\n")
      }
    }
   }
  }
  return(plot.opt)
}

#####################################################################################
###########################		Plots		#############################
#####################################################################################
###############################	   Wrapper functions  #############################

#' Plots of the results obtained by SAEM
#' 
#' Several plots (selectable by the type argument) are currently available:
#' convergence plot, individual plots, predictions versus observations,
#' distribution plots, residual plots, VPC.
#' 
#' This function plots different graphs related to the algorithm (convergence
#' plots, likelihood estimation) as well as diagnostic graphs. A description is
#' provided in the PDF documentation. 
#' \describe{ 
#' \item{data}{A spaghetti plot of the data, displaying the observed data y as a function of the regression
#' variable (eg time for a PK application)} 
#' \item{convergence}{For each parameter in the model, this plot shows the evolution of the parameter
#' estimate versus the iteration number} 
#' \item{likelihood}{Estimation of the likelihood estimated by importance sampling, as a function of the number of
#' MCMC samples} 
#' \item{individual.fit}{Individual fits, using the individual parameters with the individual covariates} 
#' \item{population.fit}{Individual fits, using the population parameters with the individual covariates}
#' \item{both.fit}{Individual fits, using the population parameters with the individual covariates and the individual parameters with the individual covariates} 
#' \item{observations.vs.predictions}{Plot of the predictions computed with the population parameters versus the observations (left), and plot of the predictions computed with the individual parameters versus the observations (right)}
#' \item{residuals.scatter}{Scatterplot of standardised residuals versus the X predictor and versus predictions. These plots are shown for individual and population residuals, as well as for npde when they are available} 
#' \item{residuals.distribution}{Distribution of standardised residuals, using histograms and QQ-plot. These plots are shown for
#' individual and population residuals, as well as for npde when they are available} 
#' \item{random.effects}{Boxplot of the random effects}
#' \item{correlations}{Correlation between the random effects, with a smoothing spline} 
#' \item{parameters.versus.covariates}{Plots of the estimate of the individual parameters versus the covariates, using scatterplot for continuous covariates, boxplot for categorical covariates}
#' \item{randeff.versus.covariates}{Plots of the estimate of the individual random effects versus the covariates, using scatterplot for continuous covariates, boxplot for categorical covariates}
#' \item{marginal.distribution}{Distribution of each parameter in the model (conditional on covariates when some are included in the model)}
#' \item{npde}{Plot of npde as in package npde (deprecated in 3.0, please use \code{\link{npdeSaemix}} instead)} 
#' \item{vpc}{Visual Predictive Check (we suggest to use \code{\link{npdeSaemix}} instead)} 
#' }
#' 
#' @param saemixObject an object returned by the \code{\link{saemix}} function
#' @param data if TRUE, produce a plot of the data. Defaults to FALSE
#' @param convergence if TRUE, produce a convergence plot. Defaults to FALSE
#' @param likelihood if TRUE, produce a plot of the estimation of the LL by
#' importance sampling. Defaults to FALSE
#' @param individual.fit if TRUE, produce individual fits with individual
#' estimates. Defaults to FALSE
#' @param population.fit if TRUE, produce individual fits with population
#' estimates. Defaults to FALSE
#' @param both.fit if TRUE, produce individual fits with both individual and
#' population estimates. Defaults to FALSE
#' @param observations.vs.predictions if TRUE, produce a plot of observations
#' versus predictions. Defaults to FALSE
#' @param residuals.scatter if TRUE, produce scatterplots of residuals versus
#' predictor and predictions. Defaults to FALSE
#' @param residuals.distribution if TRUE, produce plots of the distribution of
#' residuals. Defaults to FALSE
#' @param random.effects if TRUE, produce boxplots of the random effects.
#' Defaults to FALSE
#' @param correlations if TRUE, produce a matrix plot showing the correlation
#' between random effects. Defaults to FALSE
#' @param parameters.vs.covariates if TRUE, produce plots of the relationships
#' between parameters and covariates, using the Empirical Bayes Estimates of
#' individual parameters. Defaults to FALSE
#' @param randeff.vs.covariates if TRUE, produce plots of the relationships
#' between random effects and covariates, using the Empirical Bayes Estimates
#' of individual random effects. Defaults to FALSE
#' @param marginal.distribution if TRUE, produce plots of the marginal
#' distribution of the random effects. Defaults to FALSE
#' @param vpc if TRUE, produce Visual Predictive Check plots. Defaults to FALSE (we suggest to use \code{\link{npdeSaemix}} instead)
#' @param npde if TRUE, produce plots of the npde. Defaults to FALSE (deprecated in 3.0, please use \code{\link{npdeSaemix}} instead)
#' @param \dots optional arguments passed to the plots
#' 
#' @return None
#' 
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{SaemixObject}},\code{\link{saemix}},
#' \code{\link{default.saemix.plots}}, \code{\link{saemix.plot.setoptions}},
#' \code{\link{saemix.plot.data}}, \code{\link{saemix.plot.convergence}},
#' \code{\link{saemix.plot.llis}}, \code{\link{saemix.plot.randeff}},
#' \code{\link{saemix.plot.obsvspred}}, \code{\link{saemix.plot.fits}},
#' \code{\link{saemix.plot.parcov}}, \code{\link{saemix.plot.randeffcov}},
#' \code{\link{saemix.plot.distpsi}}, \code{\link{npdeSaemix}}
#' \code{\link{saemix.plot.scatterresiduals}},
#' \code{\link{saemix.plot.distribresiduals}}, \code{\link{saemix.plot.vpc}}
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
#' # saemix.plot.select(saemix.fit,data=TRUE,main="Spaghetti plot of data")
#' 
#' # Putting several graphs on the same plot
#' # par(mfrow=c(2,2))
#' # saemix.plot.select(saemix.fit,data=TRUE,vpc=TRUE,observations.vs.predictions=TRUE, new=FALSE)
#' 
#' @export saemix.plot.select 
##### @exportMethod plot


saemix.plot.select<-function(saemixObject,data=FALSE,convergence=FALSE, likelihood=FALSE,individual.fit=FALSE,population.fit=FALSE,both.fit=FALSE, observations.vs.predictions=FALSE,residuals.scatter=FALSE, residuals.distribution=FALSE,random.effects=FALSE,correlations=FALSE, parameters.vs.covariates=FALSE,randeff.vs.covariates=FALSE, marginal.distribution=FALSE,vpc=FALSE,npde=FALSE,...) {
# Function selecting which plots are to be drawn
  namObj<-deparse(substitute(saemixObject))
  interactive<-saemixObject["prefs"]$interactive
  boolsim<-boolpred<-boolres<-FALSE
  if(vpc) {
    if(length(saemixObject["sim.data"]["nsim"])==0) boolsim<-TRUE
  }
  if(individual.fit | both.fit | observations.vs.predictions) {
    if(length(saemixObject["results"]["ipred"])==0) boolpred<-TRUE
  }
  if(population.fit | both.fit | observations.vs.predictions) {
    if(saemixObject["prefs"]$which.poppred=="ppred" & length(saemixObject["results"]["ppred"])==0) boolpred<-TRUE
    if(saemixObject["prefs"]$which.poppred=="ypred" & length(saemixObject["results"]["ypred"])==0) boolres<-TRUE
  }
  if(residuals.scatter | residuals.distribution) {
      if(length(saemixObject["results"]["npde"])==0) boolres<-TRUE
  }
  if(boolsim & !boolres & interactive) {
    cok<-readline(prompt="Simulations will be performed to obtain residuals, VPC and npde. This might take a while, proceed ? (y/Y) [default=yes] ")
    if(!cok %in% c("y","Y","yes","")) boolsim<-FALSE 
  }
  if(boolres & interactive) {
    cok<-readline(prompt="Simulations will be performed to obtain residuals, VPC and npde. This might take a while, proceed ? (y/Y) [default=yes] ")
    if(!cok %in% c("y","Y","yes","")) boolres<-FALSE 
  }
  if(boolpred & interactive) {
    cok<-readline(prompt="Computations will be performed to obtain model predictions, proceed ? (y/Y) [default=yes] ")
    if(!cok %in% c("y","Y","yes","")) boolpred<-FALSE 
  }
  if(boolsim & !boolres) {
    saemixObject<-simulate(saemixObject)
    assign(namObj,saemixObject,envir=parent.frame())
  }
  if(boolpred) {
    saemixObject<-saemix.predict(saemixObject)
    assign(namObj,saemixObject,envir=parent.frame())
  }
    if(boolres) {
    saemixObject<-compute.sres(saemixObject)
    assign(namObj,saemixObject,envir=parent.frame())
  }
  if(parameters.vs.covariates) {
    if(length(saemixObject["results"]["map.psi"])==0) {
      saemixObject<-map.saemix(saemixObject)
      assign(namObj,saemixObject,envir=parent.frame())
    }
  }
# ECO TODO: replace with partial matching
  if(data) plot(saemixObject,plot.type="data",...)
  if(convergence) plot(saemixObject,plot.type="convergence",...)
  if(likelihood) plot(saemixObject,plot.type="likelihood",...)
  if(observations.vs.predictions) plot(saemixObject,plot.type="observations.vs.predictions", ...)
  if(individual.fit) plot(saemixObject,plot.type="individual.fit",...)
  if(population.fit) plot(saemixObject,plot.type="population.fit",...)
  if(both.fit) plot(saemixObject,plot.type="both.fit",...)
  if(residuals.scatter) plot(saemixObject,plot.type="residuals.scatter",...)
  if(residuals.distribution) plot(saemixObject,plot.type="residuals.distribution",...)
  if(random.effects) plot(saemixObject,plot.type="random.effects",...)
  if(correlations) plot(saemixObject,plot.type="correlations",...)
  if(parameters.vs.covariates) plot(saemixObject,plot.type="parameters.vs.covariates", ...)
  if(randeff.vs.covariates) plot(saemixObject,plot.type="randeff.vs.covariates",...)
  if(marginal.distribution) plot(saemixObject,plot.type="marginal.distribution",...)
  if(vpc) {
    message("Direct call to VPC will soon be deprecated, please use npdeSaemix for VPC\n")
    plot(saemixObject,plot.type="vpc",...)
    #message("Please use npdeSaemix to obtain VPC and npde\n")
  }
  if(npde) message("Please use npdeSaemix to obtain npde\n")
#  message("Please use npdeSaemix to obtain VPC and npde\n")
    # plot(saemixObject,plot.type="npde",...)
}

#### Metalevel


#' Wrapper functions to produce certain sets of default plots
#' 
#' These functions produce default sets of plots, corresponding to diagnostic
#' or individual fits.
#' 
#' These functions are wrapper functions designed to produce default sets of
#' plots to help the user assess their model fits.
#' 
#' @name default.saemix.plots
#' @aliases default.saemix.plots basic.gof advanced.gof covariate.fits individual.fits
#' @param saemixObject an object returned by the \code{\link{saemix}} function
#' @param \dots optional arguments passed to the plots
#' 
#' @return Depending on the type argument, the following plots are produced:
#' \itemize{
#' \item default.saemix.plots by default, the following plots are
#' produced: a plot of the data, convergence plots, plot of the likelihood by
#' importance sampling (if computed), plots of observations versus predictions,
#' scatterplots and distribution of residuals, boxplot of the random effects,
#' correlations between random effects, distribution of the parameters, VPC
#' \item basic.gof basic goodness-of-fit plots: convergence plots, plot of
#' the likelihood by importance sampling (if computed), plots of observations
#' versus predictions
#' \item advanced.gof advanced goodness-of-fit plots:
#' scatterplots and distribution of residuals, VPC,...
#' \item covariate.fits
#' plots of all estimated parameters versus all covariates in the dataset
#' \item individual.fits plots of individual predictions (line) overlayed on
#' individual observations (dots) for all subjects in the dataset
#' }
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{saemix}}, \code{\link{saemix.plot.data}},
#' \code{\link{saemix.plot.setoptions}}, \code{\link{plot.saemix}}
#' @references E Comets, A Lavenu, M Lavielle M (2017). Parameter estimation in nonlinear mixed effect models using saemix,
#' an R implementation of the SAEM algorithm. Journal of Statistical Software, 80(3):1-41.
#' 
#' E Kuhn, M Lavielle (2005). Maximum likelihood estimation in nonlinear mixed effects models. 
#' Computational Statistics and Data Analysis, 49(4):1020-1038.
#' 
#' E Comets, A Lavenu, M Lavielle (2011). SAEMIX, an R version of the SAEM algorithm. 20th meeting of the 
#' Population Approach Group in Europe, Athens, Greece, Abstr 2173.
#' 
#' @keywords plots
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
#' # Reducing the number of iterations due to time constraints for CRAN
#' saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE,nbiter.saemix=c(100,100))
#' 
#' saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' 
#' default.saemix.plots(saemix.fit)
#' 
#' # Not run (time constraints for CRAN)
#' # basic.gof(saemix.fit)
#' 
#' # Not run (time constraints for CRAN)
#' # advanced.gof(saemix.fit)
#' 
#' individual.fits(saemix.fit)
#' 
#' 
#' @export default.saemix.plots
#' @export basic.gof
#' @export advanced.gof
#' @export individual.fits
#' @export covariate.fits

####. Recognised arguments are: "which", allowing to select covariate fits, whether they should be produced with the EBE estimates of the parameters (the default) or with random effects (which="randeff")


default.saemix.plots<-function(saemixObject,...) {
# When plot(saemixObject) is called without plot.type  
  namObj<-deparse(substitute(saemixObject))
  if(length(saemixObject["results"]["ipred"])==0) {
    saemixObject<-saemix.predict(saemixObject)
    assign(namObj,saemixObject,envir=parent.frame())
  }
  if(length(saemixObject["results"]["npde"])==0) {
    saemixObject<-compute.sres(saemixObject)
    assign(namObj,saemixObject,envir=parent.frame())
  }
  saemix.plot.select(saemixObject,data=TRUE,convergence=TRUE,likelihood=TRUE, observations.vs.predictions=TRUE,residuals.scatter=FALSE, residuals.distribution=FALSE,random.effects=TRUE,correlations=TRUE, marginal.distribution=TRUE,vpc=TRUE,...)
#  saemix.plot.select(saemixObject,data=TRUE,convergence=TRUE,likelihood=TRUE, observations.vs.predictions=TRUE,residuals.scatter=TRUE, residuals.distribution=TRUE,random.effects=TRUE,correlations=TRUE, marginal.distribution=TRUE,vpc=TRUE,...)
}

basic.gof<-function(saemixObject,...) {
# Basic goodness of fit plots
  if(saemixObject@options$warnings) cat("Now producing basic goodness of fit plots\n")
  namObj<-deparse(substitute(saemixObject))
  if(length(saemixObject["results"]["ipred"])==0) {
    saemixObject<-saemix.predict(saemixObject)
    assign(namObj,saemixObject,envir=parent.frame())
  }
  saemix.plot.select(saemixObject,convergence=TRUE,likelihood=TRUE, observations.vs.predictions=TRUE, ...)
}

advanced.gof<-function(saemixObject,...) {
# Advanced goodness of fit plots
  if(saemixObject@options$warnings) cat("Now producing advanced goodness of fit plots\n")
  namObj<-deparse(substitute(saemixObject))
  if(length(saemixObject["results"]["ipred"])==0) {
    saemixObject<-saemix.predict(saemixObject)
    assign(namObj,saemixObject,envir=parent.frame())
  }
  if(length(saemixObject["results"]["npde"])==0) {
    saemixObject<-compute.sres(saemixObject)
    assign(namObj,saemixObject,envir=parent.frame())
  }
  saemix.plot.select(saemixObject, residuals.scatter=TRUE, residuals.distribution=TRUE, vpc=TRUE,...)
}

individual.fits<-function(saemixObject,...) {
# Individual plots
  if(saemixObject@options$warnings) cat("Now producing plots of individual fits\n")
  namObj<-deparse(substitute(saemixObject))
  if(length(saemixObject["results"]["ipred"])==0) {
    saemixObject<-saemix.predict(saemixObject)
    assign(namObj,saemixObject,envir=parent.frame())
  }
  plot(saemixObject,plot.type="individual.fit",...)
}

covariate.fits<-function(saemixObject,which="parameters",...) {
# Parameters or random effects versus covariates
  if(which=="parameters") {
    if(saemixObject@options$warnings) cat("Now producing plots of parameters versus covariates\n")
    plot(saemixObject,plot.type="parameters.vs.covariates",...)
  } else {
    if(saemixObject@options$warnings) cat("Now producing plots of random effects versus covariates\n")
    plot(saemixObject,plot.type="randeff.vs.covariates",...)
  }
}

###############################	   	Data	 #################################

# ECO FINISH THIS ONE (redo without using data part of object)

#' Functions implementing each type of plot in SAEM
#' 
#' Several plots (selectable by the type argument) are currently available:
#' convergence plot, individual plots, predictions versus observations,
#' distribution plots, VPC, residual plots, and mirror plots.
#' 
#' These functions implement plots different graphs related to the algorithm
#' (convergence plots, likelihood estimation) as well as diagnostic graphs. A
#' description is provided in the PDF documentation.
#' 
#' saemix.plot.parcov.aux, compute.sres and compute.eta.map are helper
#' functions, not intended to be called by the user directly.
#' 
#' By default, the following plots are produced:
#'  \describe{
#' \item{saemix.plot.data:}{A spaghetti plot of the data, displaying the
#' observed data y as a function of the regression variable (time for a PK
#' application)} 
#' \item{saemix.plot.convergence:}{For each parameter in the
#' model, this plot shows the evolution of the parameter estimate versus the
#' iteration number} 
#' \item{saemix.plot.llis:}{Graph showing the evolution of
#' the log-likelihood during the estimation by importance sampling}
#' \item{saemix.plot.obsvspred:}{Plot of the predictions computed with the
#' population parameters versus the observations (left), and plot of the
#' predictions computed with the individual parameters versus the observations
#' (right)} 
#' \item{saemix.plot.scatterresiduals:}{Scatterplot of the residuals
#' versus the predictor (top) and versus predictions (bottom), for weighted
#' residuals (population residuals, left), individual weighted residuals
#' (middle) and npde (right).}
#' \item{saemix.plot.distribresiduals:}{Distribution of the residuals, plotted
#' as histogram (top) and as a QQ-plot (bottom), for weighted residuals
#' (population residuals, left), individual weighted residuals (middle) and
#' npde (right).} 
#' \item{saemix.plot.fits:}{Model fits. Individual fits are
#' obtained using the individual parameters with the individual covariates.
#' Population fits are obtained using the population parameters with the
#' individual covariates (red) and the individual parameters with the
#' individual covariates (green). By default the individual plots are
#' displayed.} 
#' \item{saemix.plot.distpsi:}{Distribution of the parameters
#' (conditional on covariates when some are included in the model). A histogram
#' of individual parameter estimates can be overlayed on the plot, but it
#' should be noted that the histogram does not make sense when there are
#' covariates influencing the parameters and a warning will be displayed}
#' \item{saemix.plot.randeff:}{Boxplot of the random effects}
#' \item{saemix.plot.correlations:}{Correlation between the random effects}
#' \item{saemix.plot.parcov:}{Plots of the estimates of the individual
#' parameters versus the covariates, using scatterplot for continuous
#' covariates, boxplot for categorical covariates}
#' \item{saemix.plot.randeffcov:}{Plots of the estimates of the random effects
#' versus the covariates, using scatterplot for continuous covariates, boxplot
#' for categorical covariates} 
#' \item{saemix.plot.npde:}{Plots 4 graphs to
#' evaluate the shape of the distribution of the normalised prediction
#' distribution errors (npde)} 
#' \item{saemix.plot.vpc:}{Visual Predictive Check,
#' with options to include the prediction intervals around the boundaries of
#' the selected interval as well as around the median (50th percentile of the
#' simulated data). Several methods are available to define binning on the
#' X-axis (see methods in the PDF guide).} 
#' \item{saemix.plot.mirror:}{When simulated data is available in the object (component
#' sim.dat, which can be filled by a call to \code{\link{simulate}}), this function
#' plots the original data as spaghetti plot and compares it to several simulated
#' datasets under the fitted model.} 
#' }
#' 
#' Each plot can be customised by modifying options, either through a list of
#' options set by the \code{\link{saemix.plot.setoptions}} function, or on the
#' fly by passing an option in the call to the plot (see examples).
#' 
#' @aliases saemix.plot.data saemix.plot.convergence saemix.plot.llis
#' saemix.plot.obsvspred saemix.plot.distribresiduals
#' saemix.plot.scatterresiduals saemix.plot.fits saemix.plot.distpsi
#' saemix.plot.randeff saemix.plot.correlations saemix.plot.parcov
#' saemix.plot.randeffcov saemix.plot.npde saemix.plot.vpc
#' saemix.plot.mirror saemix.plot.parcov.aux compute.sres compute.eta.map 
#' 
#' @param saemixObject an object returned by the \code{\link{saemix}} function
#' @param \dots optional arguments passed to the plots
#' @return None
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{SaemixObject}},\code{\link{saemix}},
#' \code{\link{saemix.plot.setoptions}}, \code{\link{saemix.plot.select}},
#' \code{\link{plot.saemix}}
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
#' \donttest{
#' # Not run (strict time constraints for CRAN)
#' saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' 
#' # Simulate data and compute weighted residuals and npde
#' saemix.fit<-compute.sres(saemix.fit)
#' 
#' # Data
#' saemix.plot.data(saemix.fit)
#' 
#' # Convergence
#' saemix.plot.convergence(saemix.fit)
#' 
#' # Individual plot for subject 1, smoothed
#' saemix.plot.fits(saemix.fit,ilist=1,smooth=TRUE)
#' 
#' # Individual plot for subject 1 to 12, with ask set to TRUE 
#' # (the system will pause before a new graph is produced)
#' saemix.plot.fits(saemix.fit,ilist=c(1:12),ask=TRUE)
#' 
#' # Mirror plots (plots of simulated data compared to the original)
#' saemix.plot.mirror(saemix.fit)
#' 
#' # Diagnostic plot: observations versus population predictions
#' par(mfrow=c(1,1))
#' saemix.plot.obsvspred(saemix.fit,level=0,new=FALSE)
#' 
#' # LL by Importance Sampling
#' saemix.plot.llis(saemix.fit)
#' 
#' # Boxplot of random effects
#' saemix.plot.randeff(saemix.fit)
#' 
#' # Relationships between parameters and covariates
#' saemix.plot.parcov(saemix.fit)
#' 
#' # Relationships between parameters and covariates, on the same page
#' par(mfrow=c(3,2))
#' saemix.plot.parcov(saemix.fit,new=FALSE)
#' 
#' # Scatter plot of residuals
#' # Not run
#' # Works interactively but not in the contained environment of CRAN (it looks for a datafile 
#' # instesad of finding the dataset in the environment)
#' # saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' # npde.obj<-npdeSaemix(saemix.fit)
#' # plot(npde.obj)
#' # VPC, default options (10 bins, equal number of observations in each bin)
#' # plot(npde.obj, plot.type="vpc")
#' # plot(npde.obj, plot.type="covariates")
#' # plot(npde.obj, plot.type="cov.x.scatter")
#' # plot(npde.obj, plot.type="cov.ecdf")
#' # VPC, user-defined breaks for binning
#' # plot(npde.obj, plot.type="vpc", bin.method="user", bin.breaks=c(0.4,0.8,1.5,2.5,4,5.5,8,10,13))
#' 
#' }
#' @export saemix.plot.data

saemix.plot.data<-function(saemixObject,...) {
# Plot of the data as spaghetti plot
# options: change data point, line type, line color, lines plotted or not, points plotted or not...
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1
  plot.opt<-saemixObject["prefs"]
  plot.opt$new<-TRUE
  plot.opt$plot.type<-"l"
  plot.opt<-replace.plot.options(plot.opt,...)
  if(plot.opt$new) {
    mfrow=c(1,1)
    if(length(plot.opt$mfrow)>0) mfrow<-plot.opt$mfrow
  par(mfrow=mfrow,ask=plot.opt$ask)
  }
  plot(saemixObject["data"],plot.type=plot.opt$plot.type,...)
}


saemix.plot.mirror<-function(saemixObject, nplots=3,...) {
  # Plot of the data and nplots associated mirror plots (3 mirror plots by default)
  # Creates a new page to plot the graphs in
  # when main is used, it passes a title to all the plots; use main.origonly to pass a title only to
  # the first plot
  # options: change data point, line type, line color, lines plotted or not, points plotted or not...
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1
  userPlotOptions<-list(...)
  plot.opt<-saemixObject["prefs"]
  plot.opt$new<-FALSE
  plot.opt$main<-""
  plot.opt$plot.type<-"l"
  if(length(userPlotOptions)>0)
    plot.opt <- modifyList(plot.opt, userPlotOptions[intersect(names(userPlotOptions), names(plot.opt))])
  i1<-match("main.origonly",names(userPlotOptions))
  if(!is.na(i1)) plot.opt$main.origonly<-as.character(userPlotOptions[[i1]]) else plot.opt$main.origonly<-""
  
# Tried to use the do.call approach but it doesn't recognise the plot function to apply
#  list.args<-plot.opt
  
  if(saemixObject@sim.data@nsim>0) nplots<-min(saemixObject@sim.data@nsim, nplots)
  np<-nplots+1
  if(np>12) np<-12
  n1<-round(sqrt(np))
  n2<-ceiling(np/n1)
  par(mfrow=c(n1,n2),ask=plot.opt$ask)
  # data
  plot.opt2<-plot.opt
  if(plot.opt$main=="") plot.opt2$main<-"Original data"
  if(plot.opt2$main.origonly!="") plot.opt2$main<-plot.opt2$main.origonly
  plot(saemixObject["data"],plot.opt2)
#  print(plot.opt$new)
#  list.args$saemixObject<-saemixObject@data
#  do.call(plot, list.args)
  # simulated datasets
  if(saemixObject@sim.data@nsim>0) {
    irep1<-sort(sample(1:saemixObject@sim.data@nsim, nplots, replace=FALSE))
  } else {
    saemixObject<-simulate(saemixObject, nsim=nplots)
    irep1<-1:nplots
  }
#  list.args$saemixObject<-saemixObject@sim.data
  if(plot.opt$main=="") plot.opt2$main<-" "
  for(irep in irep1)
    plot(saemixObject["sim.data"],irep=irep, plot.opt2)
  #    do.call(plot, list.args)
}


#######################	   Convergence plots & LL	 ########################

saemix.plot.convergence<-function(saemixObject,niter=0,...) {
# Convergence plots for all the fixed effects, random effects and residual variability
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1
  plot.opt<-saemixObject["prefs"]
  plot.opt$xlab<-"Iteration"
  plot.opt<-replace.plot.options(plot.opt,...)
  change.ylab<-FALSE
  if(plot.opt$ylab!=saemixObject["prefs"]$ylab & length(plot.opt$which.par)==1) change.ylab<-TRUE
  change.main<-FALSE
  if(plot.opt$main!=saemixObject["prefs"]$main & length(plot.opt$which.par)==1) change.main<-TRUE
  K<-dim(saemixObject["results"]["allpar"])[1]
  if(niter==0) niter<-K
  if(plot.opt$which.par[1]=="all")
     np<-dim(saemixObject["results"]["allpar"])[2] else  
     np<-length(plot.opt$which.par)
  if(plot.opt$new) {
    if(length(plot.opt$mfrow)==0) {
    n1<-round(sqrt(np))
    n2<-ceiling(np/n1)
    if(n1>5 | n2>5) {
      n1<-3
      n2<-4
#      cat("Changing the plot layout\n")
    }
    par(mfrow=c(n1,n2),ask=plot.opt$ask)
  } else par(mfrow=plot.opt$mfrow,ask=plot.opt$ask)
  }
  if(plot.opt$which.par[1]=="all") { # default convergence plot, for all parameters
    for(j in 1:np) {
      laby<-"" #colnames(saemixObject["results"]["allpar"])[j]
      maintit<-colnames(saemixObject["results"]["allpar"])[j]
      plot(1:niter,saemixObject["results"]["allpar"][1:niter,j],type="l", xlab=plot.opt$xlab,ylab=laby, main=maintit,col=plot.opt$col,lty=plot.opt$lty, lwd=plot.opt$lwd,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
      abline(v=saemixObject["options"]$nbiter.saemix[1],col=plot.opt$ablinecol, lwd=plot.opt$ablinelwd)
    }
  } else {
      for(ipar in 1:length(plot.opt$which.par)) {
      j<-as.integer(plot.opt$which.par[ipar])
      if(is.na(j)) j<-which(colnames(saemixObject["results"]["allpar"])== plot.opt$which.par[ipar])
      if(length(j)>0) {
        laby<-""
        maintit<-colnames(saemixObject["results"]["allpar"])[j]
        if(change.ylab) laby<-plot.opt$ylab
        if(change.main) maintit<-plot.opt$main
        plot(1:niter,saemixObject["results"]["allpar"][1:niter,j],type="l", xlab=plot.opt$xlab,ylab=laby,main=maintit,col=plot.opt$col,lty=plot.opt$lty, lwd=plot.opt$lwd,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
      abline(v=saemixObject["options"]$nbiter.saemix[1],col=plot.opt$ablinecol, lwd=plot.opt$ablinelwd)
    }}
  }
}

saemix.plot.llis<-function(saemixObject,...) {
# Plot of the evolution of the log-likelihood by importance sampling
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1
    plot.opt<-saemixObject["prefs"]
    plot.opt$main<-"-2xLL by Importance Sampling"
    plot.opt$xlab<-"Iteration"
    plot.opt$ylab<-"-2 x LL"
    plot.opt<-replace.plot.options(plot.opt,...)
    MM<-100
    KM<-round(saemixObject["options"]$nmc.is/MM)
    kmin<-min(10,ceiling(KM/4))
    x1<-MM*c(kmin:KM)
    y1<-(-2)*saemixObject["results"]["LL"][kmin:KM]
    if(plot.opt$new) {
      if(length(plot.opt$mfrow)==0) mfrow=c(1,1) else mfrow<-plot.opt$mfrow
      par(mfrow=mfrow,ask=plot.opt$ask)
    }
    if(sum(!is.na(y1))) plot(x1,y1,type="l",xlab=plot.opt$xlab, ylab=plot.opt$ylab,main=plot.opt$main,col=plot.opt$col,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
}

#######################	   Basic GOF plots & residuals	 ########################

saemix.plot.obsvspred<-function(saemixObject,...) {
# Predictions versus observations
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1
  plot.opt<-saemixObject["prefs"]
  plot.opt$ylab<-"Observations"
  plot.opt$xlab<-"Predictions"
  plot.opt$main<-""
  plot.opt<-replace.plot.options(plot.opt,...)
  change.main<-FALSE
  if(plot.opt$main!=saemixObject["prefs"]$main) change.main<-TRUE
  if(plot.opt$new) {
    mfrow<-c(1,length(plot.opt$level))
    if(length(plot.opt$mfrow)>0) mfrow<-plot.opt$mfrow
    par(mfrow=mfrow,ask=plot.opt$ask)
  }
  if(saemixObject["model"]["error.model"]=="exponential")
    ydat<-saemixObject["data"]["yorig"] else ydat<-saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]
  if(length(grep(0,plot.opt$level))>0) {
    if(!change.main) main<-"Population predictions" else main<-plot.opt$main
    if(plot.opt$which.poppred=="ppred") xpl<-saemixObject["results"]["ppred"] else xpl<-saemixObject["results"]["ypred"]
    if(length(xpl)==length(ydat)) {
    plot(xpl,ydat,xlab=plot.opt$xlab, ylab=plot.opt$ylab,pch=plot.opt$pch, col=plot.opt$col,main=main,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    abline(0,1,col=plot.opt$ablinecol,lty=plot.opt$ablinelty, lwd=plot.opt$ablinelwd)
    }
     }
  if(length(grep(1,plot.opt$level))>0) {
    if(!change.main) main<-paste("Individual predictions", ifelse(plot.opt$indiv.par=="map","MAP","Cond mean"),sep=", ") else main<-plot.opt$main
    if(plot.opt$indiv.par=="map") xpl<-saemixObject["results"]["ipred"] else xpl<-saemixObject["results"]["icpred"]
    if(length(xpl)==length(ydat)) {
    plot(xpl,ydat,xlab=plot.opt$xlab, ylab=plot.opt$ylab,pch=plot.opt$pch, col=plot.opt$col,main=main,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    abline(0,1,col=plot.opt$ablinecol,lty=plot.opt$ablinelty,lwd=plot.opt$ablinelwd)
    }
   }
}

# Eco TODO: deprecate this function to use npde library instead
# current function only for wres (no PI)

saemix.plot.distribresiduals<-function(saemixObject,...) {
# Histogram and QQ-plot
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1
  plot.opt<-saemixObject["prefs"]
  plot.opt$main<-""
  plot.opt$level<-0:1
  plot.opt$smooth<-TRUE
  plot.opt$which.resplot<-c("dist.qqplot","dist.hist")
  plot.opt$which.pres<-c("wres") #,"npde")
  plot.opt<-replace.plot.options(plot.opt,...)
  change.main<-FALSE
  if(plot.opt$main!=saemixObject["prefs"]$main) change.main<-TRUE
  change.xlab<-FALSE
  if(plot.opt$xlab!=saemixObject["prefs"]$xlab) change.xlab<-TRUE
  change.ylab<-FALSE
  if(plot.opt$ylab!=saemixObject["prefs"]$ylab) change.ylab<-TRUE
  if(plot.opt$new) {
    if(length(plot.opt$mfrow)>0) mfrow<-plot.opt$mfrow else {
      ncol<-as.integer(1%in%plot.opt$level)+as.integer(0%in%plot.opt$level)* length(plot.opt$which.pres)
      mfrow=c(length(plot.opt$which.resplot),ncol)
    }
    par(mfrow=mfrow,ask=plot.opt$ask)
  }
  plot.ind<-FALSE
  if(1%in%plot.opt$level) {
    if(length(saemixObject["results"]["iwres"])==0) {
      if(saemixObject@options$warnings) cat("Please compute individual residuals first using predict().\n")
      return("No individual residuals")
    }
    plot.ind<-TRUE
    if(plot.opt$indiv.par=="map") {
      iwres<-saemixObject["results"]["iwres"]
    } else {
      iwres<-saemixObject["results"]["icwres"]
    }
  }
  plot.pop<-FALSE
  if(0%in%plot.opt$level) {
    if(length(saemixObject["results"]["wres"])==0 | length(saemixObject["results"]["npde"])==0) {
      if(saemixObject@options$warnings) cat("Please compute WRES and npde first by using compute.sres().\n")
      return("No WRES or npde")
    }
    plot.pop<-TRUE
    wres<-saemixObject["results"]["wres"]
    npde<-saemixObject["results"]["npde"]
  }
  if("dist.qqplot"%in%plot.opt$which.resplot) {
  if(plot.pop & "wres"%in%plot.opt$which.pres) {
    laby<-"Sample quantiles"
    if(change.ylab) laby<-plot.opt$ylab
    labx<-"Theoretical quantiles"
    if(change.xlab) labx<-plot.opt$xlab
    main<-"Population weighted residuals"
    if(change.main) main<-plot.opt$main
    qqnorm(wres,xlab=labx,ylab=laby,main=plot.opt$main, col=plot.opt$col)
    qqline(wres,lty=plot.opt$lty,col=plot.opt$col)
  }
  if(plot.ind) {
    laby<-"Sample quantiles"
    if(change.ylab) laby<-plot.opt$ylab
    labx<-"Theoretical quantiles"
    if(change.xlab) labx<-plot.opt$xlab
    main<-"Individual weighted residuals"
    if(change.main) main<-plot.opt$main
    qqnorm(iwres,xlab=labx,ylab=laby,main=plot.opt$main, col=plot.opt$col)
    qqline(iwres,lty=plot.opt$lty,col=plot.opt$col)
  }
  if(plot.pop & "npde"%in%plot.opt$which.pres) {
    laby<-"Sample quantiles"
    if(change.ylab) laby<-plot.opt$ylab
    labx<-"Theoretical quantiles"
    if(change.xlab) labx<-plot.opt$xlab
    main<-"NPDE"
    if(change.main) main<-plot.opt$main
    qqnorm(npde,xlab=labx,ylab=laby,main=plot.opt$main, col=plot.opt$col)
    qqline(npde,lty=plot.opt$lty,col=plot.opt$col)
  }
  }
  if("dist.hist"%in%plot.opt$which.resplot) {
  if(plot.pop & "wres"%in%plot.opt$which.pres) {
    labx<-"Population weighted residuals"
    if(change.xlab) labx<-plot.opt$xlab
    vec<-wres
    xh<-hist(vec,nclass=10,main=plot.opt$main, xlab=labx)
    if(plot.opt$smooth) {
      xpl<-min(vec)+c(0:100)/100*(max(vec)-min(vec))
      ypl<-dnorm(xpl)
      ypl<-ypl/max(ypl)*max(xh$counts)
      lines(xpl,ypl,lwd=2)
    }
  }
  if(plot.ind) {
    labx<-"Individual weighted residuals"
    if(change.xlab) labx<-plot.opt$xlab
    vec<-iwres
    xh<-hist(vec,nclass=10,main=plot.opt$main, xlab=labx)
    if(plot.opt$smooth) {
      xpl<-min(vec)+c(0:100)/100*(max(vec)-min(vec))
      ypl<-dnorm(xpl)
      ypl<-ypl/max(ypl)*max(xh$counts)
      lines(xpl,ypl,lwd=2)
    }
  }
  if(plot.pop & "npde"%in%plot.opt$which.pres) {
    labx<-"NPDE"
    if(change.xlab) labx<-plot.opt$xlab
    vec<-npde
    xh<-hist(vec,nclass=10,main=plot.opt$main, xlab=labx)
    if(plot.opt$smooth) {
      xpl<-min(vec)+c(0:100)/100*(max(vec)-min(vec))
      ypl<-dnorm(xpl)
      ypl<-ypl/max(ypl)*max(xh$counts)
      lines(xpl,ypl,lwd=2)
    }
  }
  }
}

# Eco TODO: deprecate this function to use npde library instead
# removed npde from list of residuals to plot
# current function should be only for wres (no PI)

saemix.plot.scatterresiduals<-function(saemixObject,...) {
# Graphs of residuals versus time and predictions
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1
  plot.opt<-saemixObject["prefs"]
  plot.opt$main<-""
  plot.opt$level<-0:1
  plot.opt$which.resplot<-c("res.vs.x","res.vs.pred")
  plot.opt$which.pres<-c("wres") # only wres, better plots with the npde package
  #  plot.opt$which.pres<-c("wres","npde")
  plot.opt<-replace.plot.options(plot.opt,...)
  change.main<-FALSE
  if(plot.opt$main!=saemixObject["prefs"]$main) change.main<-TRUE
  change.xlab<-FALSE
  if(plot.opt$xlab!=saemixObject["prefs"]$xlab) change.xlab<-TRUE
  change.ylab<-FALSE
  if(plot.opt$ylab!=saemixObject["prefs"]$ylab) change.ylab<-TRUE
  if(plot.opt$new) {
    if(length(plot.opt$mfrow)>0) mfrow<-plot.opt$mfrow else {
      ncol<-as.integer(1%in%plot.opt$level)+as.integer(0%in%plot.opt$level)* length(plot.opt$which.pres)
      mfrow=c(length(plot.opt$which.resplot),ncol)
    }
    par(mfrow=mfrow,ask=plot.opt$ask)
  }
  plot.ind<-FALSE
  if(1%in%plot.opt$level) {
    if(length(saemixObject["results"]["iwres"])==0) {
      if(saemixObject@options$warnings) cat("Please compute individual residuals first using predict().\n")
      return("No individual residuals")
    }
    plot.ind<-TRUE
    if(plot.opt$indiv.par=="map") {
      iwres<-saemixObject["results"]["iwres"]
      ipred<-saemixObject["results"]["ipred"]
    } else {
      iwres<-saemixObject["results"]["icwres"]
      ipred<-saemixObject["results"]["icpred"]
    }
  }
  plot.pop<-FALSE
  if(0%in%plot.opt$level) {
    if(length(saemixObject["results"]["wres"])==0 | length(saemixObject["results"]["npde"])==0) {
      if(saemixObject@options$warnings) cat("Please compute WRES and npde first by using compute.sres().\n")
      return("No WRES or npde")
    }
    plot.pop<-TRUE
    wres<-saemixObject["results"]["wres"]
    npde<-saemixObject["results"]["npde"]
    if(plot.opt$which.poppred=="ppred") ppred<-saemixObject["results"]["ppred"] else ppred<-saemixObject["results"]["ypred"]
  }
  if("res.vs.x"%in%plot.opt$which.resplot) {
  if(plot.pop & "wres"%in%plot.opt$which.pres) {
    laby<-"Population weighted residuals"
    if(change.ylab) laby<-plot.opt$ylab
    plot(saemixObject["data"]["data"][,saemixObject["data"]["name.X"]],wres, pch=plot.opt$pch, col=plot.opt$col,main=plot.opt$main,xlab=plot.opt$xlab,ylab=laby,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    abline(h=0,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=-1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
  }
  if(plot.ind) {
    laby<-"Individual weighted residuals"
    if(change.ylab) laby<-plot.opt$ylab
    plot(saemixObject["data"]["data"][,saemixObject["data"]["name.X"]],iwres, pch=plot.opt$pch,col=plot.opt$col,main=plot.opt$main,xlab=plot.opt$xlab,ylab=laby,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    abline(h=0,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=-1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
  }
  if(plot.pop & "npde"%in%plot.opt$which.pres) {
    laby<-"NPDE"
    if(change.ylab) laby<-plot.opt$ylab
    plot(saemixObject["data"]["data"][,saemixObject["data"]["name.X"]],npde, pch=plot.opt$pch,col=plot.opt$col,main=plot.opt$main,xlab=plot.opt$xlab,ylab=laby,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    abline(h=0,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=-1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
  }
  }
  if("res.vs.pred"%in%plot.opt$which.resplot) {
  if(plot.pop & "wres"%in%plot.opt$which.pres) {
    laby<-"Population weighted residuals"
    if(change.ylab) laby<-plot.opt$ylab
    labx<-"Population predictions"
    if(change.xlab) labx<-plot.opt$xlab
    plot(ppred,wres,pch=plot.opt$pch,col=plot.opt$col,main=plot.opt$main, xlab=labx,ylab=laby,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    abline(h=0,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=-1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
  }
  if(plot.ind) {
    laby<-"Individual weighted residuals"
    if(change.ylab) laby<-plot.opt$ylab
    labx<-"Individual predictions"
    if(change.xlab) labx<-plot.opt$xlab
    plot(ipred,iwres,pch=plot.opt$pch,col=plot.opt$col,main=plot.opt$main, xlab=labx,ylab=laby,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    abline(h=0,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=-1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
  }
  if(plot.pop & "npde"%in%plot.opt$which.pres) {
    laby<-"NPDE"
    if(change.ylab) laby<-plot.opt$ylab
    labx<-"Population predictions"
    if(change.xlab) labx<-plot.opt$xlab
    plot(ppred,npde,pch=plot.opt$pch,col=plot.opt$col,main=plot.opt$main, xlab=labx,ylab=laby,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    abline(h=0,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=-1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
  }
  }
}

#######################	  	Individual fits		 ########################

saemix.plot.fits<-function(saemixObject,...) {
# Plot of the model fits overlayed on the data
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1
  plot.opt<-saemixObject["prefs"]
  plot.opt$main<-""
  plot.opt$xlab<-paste(saemixObject["data"]["name.X"]," (",saemixObject["data"]["units"]$x,")",sep="")
  plot.opt$ylab<-paste(saemixObject["data"]["name.response"]," (",saemixObject["data"]["units"]$y,")",sep="")
  plot.opt$new<-TRUE
  plot.opt$ilist<-1:saemixObject["data"]["N"]  
  plot.opt$type<-"p"
  plot.opt$level<-c(1)
  plot.opt$ipred.lty<-1
  plot.opt$ppred.lty<-2
  plot.opt<-replace.plot.options(plot.opt,...)
  plot.opt$ilist<-plot.opt$ilist[plot.opt$ilist %in% 1:saemixObject["data"]["N"]]
  if(plot.opt$new) {
    if(length(plot.opt$mfrow)==0) {
    np<-length(plot.opt$ilist)
    if(np>12) np<-12
    n1<-round(sqrt(np))
    n2<-ceiling(np/n1)
    par(mfrow=c(n1,n2),ask=plot.opt$ask)
  } else par(mfrow=plot.opt$mfrow,ask=plot.opt$ask)
  }
  change.main<-FALSE
  if(plot.opt$main!=saemixObject["prefs"]$main) change.main<-TRUE
  indplot<-(length(grep(1,plot.opt$level))>0)
  popplot<-(length(grep(0,plot.opt$level))>0)
  if(indplot & plot.opt$smooth & length(saemixObject["results"]["map.psi"])==0) {
    if(saemixObject@options$warnings) cat("Individual parameter estimates should be computed to produce individual plots, conditional means will be used.\n")
  }
  if(indplot & !(plot.opt$smooth) & length(saemixObject["results"]["ipred"])==0) {
    if(saemixObject@options$warnings) cat("For graphs of predictions, please use predict first.\n") 
    return("Run predict first")
  }
  if(popplot & !(plot.opt$smooth) & length(saemixObject["results"]["ypred"])==0) {
    cat("For graphs of predictions, please use predict first.\n") 
    if(saemixObject@options$warnings) return("Run predict first")
  }
  logtyp<-""
  if(plot.opt$xlog) logtyp<-paste(logtyp,"x",sep="")
  if(plot.opt$ylog) logtyp<-paste(logtyp,"y",sep="")
  pl.line<-(length(plot.opt$level)>0)
  xind<-saemixObject["data"]["data"][,c(saemixObject["data"]["name.predictors"],saemixObject["data"]["name.cens"],saemixObject["data"]["name.mdv"],saemixObject["data"]["name.ytype"]),drop=FALSE]
  id<-saemixObject["data"]["data"][,saemixObject["data"]["name.group"]]
  if(saemixObject["model"]["error.model"]=="exponential")
    yobs<-saemixObject["data"]["yorig"] else yobs<-saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]
  for(i1 in plot.opt$ilist) {
    isuj<-unique(id)[i1]
    if((indplot | popplot) & plot.opt$smooth) {
# If smooth requested and either population and individual predictions
      xvec<-xind[id==isuj, saemixObject["data"]["name.X"]]
      xpred<-seq(min(xvec),max(xvec),length.out=100)
      if(dim(xind)[2]==1) xdep<-matrix(xpred,ncol=1) else {
        x1<-xind[id==isuj,]
# creating an expanded dataframe (warning: will not work with different occasions)
# ECO TODO change this when several occasions
        id1<-unlist(lapply(xpred,function(x,vec) max(which(vec<=x)),vec=xvec))
        xdep<-x1[id1,]
        xdep[,saemixObject["data"]["name.X"]]<-xpred
      }
      idx<-rep(i1,dim(xdep)[1])
      if(indplot) {
# ECO TODO change this when several occasions
        if(length(saemixObject["results"]["map.psi"])>0)
          ypred<-saemixObject["model"]["model"](saemixObject["results"]["map.psi"],idx,xdep)  else {
         psiM<-transphi(saemixObject["results"]["cond.mean.phi"], saemixObject["model"]["transform.par"])
         ypred<-saemixObject["model"]["model"](psiM,idx,xdep)
        }
      }
      if(popplot) {
       psiM<-transphi(saemixObject["results"]["mean.phi"], saemixObject["model"]["transform.par"])
       yppred<-saemixObject["model"]["model"](psiM,idx,xdep)
      }
    } else {
# else, use predictions at each observation time (no smooth)
      xpred<-xind[id==isuj,saemixObject["data"]["name.X"]]
      ypred<-saemixObject["results"]["ipred"][id==isuj]
      yppred<-saemixObject["results"]["ppred"][id==isuj]
    }
    vec<-yobs[id==isuj]
    if(indplot) vec<-c(vec,ypred)
    if(popplot) vec<-c(vec,yppred)
    if(length(plot.opt$ylim)>0) limy<-plot.opt$ylim else {
      if(!plot.opt$ylog) limy<-c(min(vec,na.rm=TRUE),max(vec,na.rm=TRUE)) else limy<-c(min(vec[!is.na(vec) & vec>0]),max(vec[!is.na(vec) & vec>0]))
    }
    main<-paste("Subject",isuj)
    if(change.main) main<-plot.opt$main
    plot(xind[id==isuj,saemixObject["data"]["name.X"]], yobs[id==isuj], xlab=plot.opt$xlab,ylab=plot.opt$ylab,log=logtyp,ylim=limy,type=plot.opt$type, main=main,pch=plot.opt$obs.pch,col=plot.opt$obs.col,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    if(pl.line) {
      if(indplot) lines(xpred,ypred,col=plot.opt$ipred.col, lty=plot.opt$ipred.lty,lwd=plot.opt$ipred.lwd)
      if(popplot) lines(xpred,yppred,col=plot.opt$ppred.col, lty=plot.opt$ppred.lty,lwd=plot.opt$ppred.lwd)
    }
  }
}

#######################	   Advanced GOF plots (VPC, npde) ########################
    
plotnpde<-function(xobs,npde,ypred) {
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1
  nclass<-10
    par(mfrow=c(2,2))
    qqnorm(sort(npde),xlab="Sample quantiles (npde)",ylab="Theoretical Quantiles", cex.lab=1.5,main="Q-Q plot versus N(0,1) for npde")
    qqline(sort(npde))
    #Histogram of npde, with N(0,1) superimposed on the plot
    xh<-hist(npde,nclass=nclass,xlab="npde",main="",cex.lab=1.5)
    xpl<-min(npde)+c(0:100)/100*(max(npde)-min(npde))
    ypl<-dnorm(xpl)
    ypl<-ypl/max(ypl)*max(xh$counts)
    lines(xpl,ypl,lwd=2)
    
    #residuals
    plot(xobs,npde,xlab="X",ylab="npde",cex.lab=1.5)
    abline(h=0,lty=2)
    x1<-qnorm(0.05)
    abline(h=x1,lty=3);abline(h=(-x1),lty=3)
    
    plot(ypred,npde,xlab="Predicted Y",ylab="npde",cex.lab=1.5)
    abline(h=0,lty=2)
    abline(h=x1,lty=3);abline(h=(-x1),lty=3)
}

# Eco TODO: Replaced by a call to the npde library in 3.0, remove in next version

saemix.plot.npde<-function(saemixObject,...) {
  if(length(saemixObject["results"]["npde"])==0) {
    if(saemixObject@options$warnings) cat("Please estimate the npde first\n")
    return("No npde")
  }
  plotnpde(saemixObject["data"]["data"][,saemixObject["data"]["name.X"]], saemixObject["results"]["npde"],saemixObject["results"]["ypred"])
  y<-testnpde(saemixObject["results"]["npde"])
  return(y)
}

saemix.plot.vpc<-function(saemixObject,npc=FALSE,...) {
  if(length(saemixObject["sim.data"]["nsim"])==0) {
    message("Please simulate data first, using the simulate function.\n") 
    return("No simulated data")
  }
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1
  
# Internal function
compute.vpc.pi<-function(ysim,xgrp,idrep,nbin,vpc.pi=0.95) {
  nsim<-length(unique(idrep))
  sim.pi.low<-sim.pi.med<-sim.pi.up<-matrix(0,nrow=nbin,ncol=nsim)
  alpha<-(1-vpc.pi)/2
  i0<-1
  for(irep in unique(idrep)) {
    ysim1<-ysim[idrep==irep]
    l1<-unlist(tapply(ysim1,xgrp,function(vec) quantile(vec,c(alpha,0.5,1-alpha))))
    l1<-matrix(l1,ncol=3,byrow=TRUE)
    sim.pi.low[,i0]<-l1[,1]
    sim.pi.med[,i0]<-l1[,2]
    sim.pi.up[,i0]<-l1[,3]
    i0<-i0+1
  }
  return(list(sim.pi.low=sim.pi.low,sim.pi.med=sim.pi.med,sim.pi.up=sim.pi.up))
}

  plot.opt<-saemixObject["prefs"]
  plot.opt$main<-"Visual Predictive Check"
  plot.opt<-replace.plot.options(plot.opt,...)
  if(plot.opt$new) {
    mfrow<-plot.opt$mfrow
    if(length(mfrow)==0) mfrow<-c(1,1)
    par(mfrow=mfrow,ask=plot.opt$ask)
  }
  logtyp<-""
  if(plot.opt$xlog) logtyp<-paste(logtyp,"x",sep="")
  if(plot.opt$ylog) logtyp<-paste(logtyp,"y",sep="")
  
  if(!is.na(pmatch(plot.opt$vpc.method,"optimal"))) {
    if(saemixObject@options$warnings) cat("Optimal binning not yet implemented, reverting to equal binning\n")
    plot.opt$vpc.method<-"equal"
  }
  if(!is.na(pmatch(plot.opt$vpc.method,"user")) & is.null(plot.opt$vpc.breaks)) {
    if(saemixObject@options$warnings) cat("User-defined method specified, but vpc.breaks is empty; reverting to equal binning\n")
    plot.opt$vpc.method<-"equal"
  }
  if(!is.na(pmatch(plot.opt$vpc.method,c("equal","width"))) & is.null(plot.opt$vpc.bin)) {
    plot.opt$vpc.bin<-10
  }

# Binning
  xvec<-saemixObject["data"]["data"][,saemixObject["data"]["name.X"]]
  ydat<-saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]
  ysim<-saemixObject["sim.data"]["datasim"]$ysim
  nbin<-plot.opt$vpc.bin
  alpha<-(1-plot.opt$vpc.interval)/2
# ECO TODO: implement the optimal binning algorithm of Marc
#  if(is.na(plot.opt$vpc.bin)) {
#  } else { # binning by quantiles
#    bnds<-unique(quantile(xvec,seq(0,1,length.out=nbin),type=8))
#    xgrp<-findInterval(xvec,bnds)
    if(!is.na(pmatch(plot.opt$vpc.method,"user"))) {
      bnds<-plot.opt$vpc.breaks
      if(min(bnds)>=min(xvec)) bnds<-c(min(xvec)-1,bnds)
      if(max(bnds)<max(xvec)) bnds<-c(bnds,max(xvec))
    }
    if(!is.na(pmatch(plot.opt$vpc.method,"equal"))) {
      xvec2<-xvec;xvec2[xvec2==min(xvec)]<-min(xvec)-1
      bnds<-unique(quantile(xvec2,(0:nbin)/nbin,type=8))
    }
    if(!is.na(pmatch(plot.opt$vpc.method,"width"))) {
      if(plot.opt$xlog) xvec2<-log(xvec) else xvec2<-xvec
      bnds<-seq(min(xvec2),max(xvec2),length.out=(nbin+1))
      if(plot.opt$xlog) bnds<-exp(bnds)
      bnds[1]<-bnds[1]-1
    }
    if(!is.na(pmatch(plot.opt$vpc.method,c("equal","width","user")))) {
      xgrp<-factor(cut(xvec,bnds,include.lowest=F))
      if(!is.na(pmatch(plot.opt$vpc.method,"equal")) & length(unique(xvec))<=nbin)
        xgrp<-match(xvec,sort(unique(xvec)))
    } else {
      
    }
    nbin<-length(unique(xgrp))
    xpl<-tapply(xvec,xgrp,mean)
    if(!is.na(pmatch(plot.opt$vpc.method,c("equal","width","user")))) {
      tab<-cbind(Interval=names(xpl),Centered.On=format(xpl,digits=2))
      row.names(tab)<-1:dim(tab)[1]
      xnam<-switch(EXPR=plot.opt$vpc.method,equal="binning by quantiles on X", width="equal sized intervals",user="user-defined bins")
      if(saemixObject@options$warnings) {
        cat("Method used for VPC:",xnam,", dividing into the following intervals\n")
        print(tab,quote=F)
      }
    }
# Observed data
    ypl<-tapply(ydat,xgrp,mean)
    obs.bnd<-cbind(tapply(ydat,xgrp,quantile,alpha),tapply(ydat,xgrp,mean), tapply(ydat,xgrp,quantile,1-alpha))
#  }
  if(plot.opt$vpc.pi) {
#    idsim<-saemixObject["sim.data"]["datasim"]$idsim
    idrep<-saemixObject["sim.data"]["datasim"]$irep
    isamp<-sample(1:saemixObject["options"]$nb.sim, saemixObject["options"]$nb.simpred)
    idx<-match(idrep,isamp,nomatch=0)>0
    sbnd<-compute.vpc.pi(ysim[idx],xgrp,idrep[idx],nbin,0.95)
    pi.low<-apply(sbnd$sim.pi.low,1,quantile,c(0.025,0.5,0.975))
    pi.med<-apply(sbnd$sim.pi.med,1,quantile,c(0.025,0.5,0.975))
    pi.up<-apply(sbnd$sim.pi.up,1,quantile,c(0.025,0.5,0.975))
    vec1<-c(pi.low,obs.bnd[,1])
    vec2<-c(obs.bnd[,3],pi.up)
    if(length(grep("y",logtyp))>0) {
      vec1<-vec1[vec1>0]
      vec2<-vec2[vec2>0]
    }
    limy<-c(min(vec1),max(vec2))
    xvec1<-xvec
    if(length(grep("x",logtyp))>0) xvec1<-xvec1[xvec1>0]
    limx<-c(min(xvec1),max(xvec1))

    plot(xpl,ypl,type="n",xlim=limx,ylim=limy,xlab=plot.opt$xlab, ylab=plot.opt$ylab,main=plot.opt$main,log=logtyp,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    polygon(c(xpl,rev(xpl)),c(pi.low[1,],rev(pi.low[3,])), col=plot.opt$col.fillpi,lty=plot.opt$lty.lpi, border=plot.opt$col.lpi)
    polygon(c(xpl,rev(xpl)),c(pi.up[1,],rev(pi.up[3,])), col=plot.opt$col.fillpi,lty=plot.opt$lty.lpi, border=plot.opt$col.lpi)
    polygon(c(xpl,rev(xpl)),c(pi.med[1,],rev(pi.med[3,])), col=plot.opt$col.fillmed,lty=plot.opt$lty.lmed, border=plot.opt$col.lmed)
    lines(xpl,pi.low[2,],lty=plot.opt$lty.lpi, col=plot.opt$col.lpi,lwd=plot.opt$lwd.lpi)
    lines(xpl,pi.med[2,],lty=plot.opt$lty.lmed, col=plot.opt$col.lmed,lwd=plot.opt$lwd.lmed)
    lines(xpl,pi.up[2,],lty=plot.opt$lty.lpi, col=plot.opt$col.lpi,lwd=plot.opt$lwd.lpi)
    lines(xpl,obs.bnd[,2],lty=plot.opt$lty.lobs, col=plot.opt$col.lmed,lwd=plot.opt$lwd.lobs)
    for (icol in c(1,3)) lines(xpl,obs.bnd[,icol],lty=plot.opt$lty.lobs, col=plot.opt$col.lobs,lwd=plot.opt$lwd.lobs)
    if(plot.opt$vpc.obs)
      points(xvec,ydat,pch=plot.opt$pch,col=plot.opt$col.pobs)
  } else {
# Simulated data
    nsim<-length(ysim)/length(ydat)
    id.grp<-rep(xgrp,nsim)
    sim.bnd<-cbind(tapply(ysim,id.grp,quantile,alpha),tapply(ysim,id.grp,mean), tapply(ysim,id.grp,quantile,1-alpha))
    vec1<-c(obs.bnd[,1],sim.bnd[,1])
    vec2<-c(obs.bnd[,3],sim.bnd[,3])
    if(length(grep("y",logtyp))>0) {
      vec1<-vec1[vec1>0]
      vec2<-vec2[vec2>0]
    }
    limy<-c(min(vec1),max(vec2))
    xvec1<-xvec
    if(length(grep("x",logtyp))>0) xvec1<-xvec1[xvec1>0]
    limx<-c(min(xvec1),max(xvec1))
    plot(xpl,ypl,type="n",xlim=limx,ylim=limy,xlab=plot.opt$xlab, ylab=plot.opt$ylab,main=plot.opt$main,log=logtyp,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    polygon(c(xpl,rev(xpl)),c(sim.bnd[,3],rev(sim.bnd[,1])), col=plot.opt$fillcol,lty=plot.opt$ablinelty, border=plot.opt$ablinecol)
    lines(xpl,sim.bnd[,2],lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    lines(xpl,obs.bnd[,2],lty=plot.opt$lty, col=plot.opt$lcol,lwd=plot.opt$lwd)
    for (icol in c(1,3)) lines(xpl,obs.bnd[,icol],lty=plot.opt$lty, col=plot.opt$lcol,lwd=plot.opt$lwd)
    if(plot.opt$vpc.obs)
      points(xvec,ydat,pch=plot.opt$pch,col=plot.opt$pcol)  
  }
  npc.stat<-c()
  if(npc==TRUE) {
    # ECO TODO: compute NPC - interpolation ? 
  }
  return(npc=npc.stat)
}

#######################	   Distribution of random effects ########################
saemix.plot.correlations<-function(saemixObject,...) {
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1
  plot.opt<-saemixObject["prefs"]
  plot.opt$which.par<-"all"
  plot.opt$main<-"Correlations between random effects"
  plot.opt<-replace.plot.options(plot.opt,...)
  if(plot.opt$which.par[1]=="all") plot.opt$which.par<-saemixObject["model"]["name.modpar"]
  plist<-match(plot.opt$which.par,saemixObject["model"]["name.modpar"])
  plist<-plist[!is.na(match(plist,saemixObject["model"]["indx.omega"]))]
  if(plot.opt$indiv.par=="map" & length(saemixObject["results"]["map.psi"])) {
    indiv.par<-saemixObject["results"]["map.psi"]
  } else {
    if(plot.opt$indiv.par=="map" & (saemixObject@options$warnings) ) cat("No MAP estimates, using the conditional means for individual parameters.\n")
    indiv.par<-transphi(saemixObject["results"]["cond.mean.phi"], saemixObject["model"]["transform.par"])
  }
  labs<-saemixObject["model"]["name.modpar"][plist]
  pairs(indiv.par[,plist],labels=labs,panel=panel.smooth,main=plot.opt$main, pch=plot.opt$pch,col=plot.opt$col)
}

saemix.plot.randeff<-function(saemixObject,...) {
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1
  plot.opt<-saemixObject["prefs"]
  plot.opt$which.par<-"all"
  plot.opt$main<-""
  plot.opt<-replace.plot.options(plot.opt,...)
  if(plot.opt$which.par[1]=="all") plot.opt$which.par<-saemixObject["model"]["name.modpar"]
  change.xlab<-FALSE
  if(plot.opt$xlab!=saemixObject["prefs"]$xlab & length(plot.opt$which.par)==1) change.xlab<-TRUE
  change.main<-FALSE
  if(plot.opt$main!=saemixObject["prefs"]$main) change.main<-TRUE
  if(plot.opt$new) {
    if(length(plot.opt$mfrow)==0) {
    np<-length(plot.opt$which.par)
    n1<-round(sqrt(np))
    n2<-ceiling(np/n1)
    if(n1>5 | n2>5) {
      n1<-3
      n2<-4
#      cat("Changing the plot layout\n")
    }
    par(mfrow=c(n1,n2),ask=plot.opt$ask)
  } else par(mfrow=plot.opt$mfrow,ask=plot.opt$ask)
  }
#  if(length(grep("l",plot.opt$line.smooth))>0)
  if(plot.opt$smooth)
    pl.psi<-transpsi(matrix(saemixObject["results"]["fixed.effects"],nrow=1), saemixObject["model"]["transform.par"])
  plist<-match(plot.opt$which.par,saemixObject["model"]["name.modpar"])
  for(ipar in plist) {
    tit<-saemixObject["model"]["name.modpar"][ipar]
    if(change.main) tit<-plot.opt$main    
    labx<-""
    if(change.xlab) labx<-plot.opt$xlab
    boxplot(saemixObject["results"]["cond.mean.phi"][,ipar],xlab=labx,main=tit,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
#    if(length(grep("l",plot.opt$line.smooth))>0)
    if(plot.opt$smooth)
      abline(h=pl.psi[1,ipar],lty=plot.opt$lty,lwd=plot.opt$ablinelwd, col=plot.opt$ablinecol)
  }
}

saemix.plot.distpsi<-function(saemixObject,...) {
# Plots the distribution of the model parameters conditional on covariates 
# plot.opt$cov.value: value of the covariates used
# Adds an histogram of individual parameter estimates if plot.opt$indiv.histo==TRUE
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1
  plot.opt<-saemixObject["prefs"]
  plot.opt$which.par<-"all"
  plot.opt$main<-""
  plot.opt<-replace.plot.options(plot.opt,...)
  if(plot.opt$which.par[1]=="all") plot.opt$which.par<-saemixObject["model"]["name.modpar"]
  change.xlab<-FALSE
  if(plot.opt$xlab!=saemixObject["prefs"]$xlab & length(plot.opt$which.par)==1) change.xlab<-TRUE
  change.ylab<-FALSE
  if(plot.opt$ylab!=saemixObject["prefs"]$ylab & length(plot.opt$which.par)==1) change.ylab<-TRUE
  change.main<-FALSE
  if(plot.opt$main!=saemixObject["prefs"]$main) change.main<-TRUE
  nampar<-saemixObject["model"]["name.modpar"]
  plist<-match(plot.opt$which.par,nampar)
  if(plot.opt$new) {
    if(length(plot.opt$mfrow)==0) {
    np<-length(plist)
    n1<-round(sqrt(np))
    n2<-ceiling(np/n1)
    if(n1>5 | n2>5) {
      n1<-3
      n2<-4
#      cat("Changing the plot layout\n")
    }
    par(mfrow=c(n1,n2),ask=plot.opt$ask)
  } else par(mfrow=plot.opt$mfrow,ask=plot.opt$ask)
  }
  if(plot.opt$indiv.histo) {
    if(plot.opt$indiv.par=="map" & length(saemixObject["results"]["map.psi"])) {
      indiv.par<-saemixObject["results"]["map.psi"]
    } else {
      if(plot.opt$indiv.par=="map" & (saemixObject@options$warnings) ) cat("No MAP estimates, using the conditional means for individual parameters.\n")
      indiv.par<-transphi(saemixObject["results"]["cond.mean.phi"], saemixObject["model"]["transform.par"])
    }
  }
  xpl<-c(1:100)/100*4
  xpl<-c(rev(-xpl),0,xpl)
  mpar<-saemixObject["results"]["betas"][saemixObject["model"]["indx.fix"]]
  for(ipar in plist) {
    labx<-nampar[ipar]
    if(change.xlab) labx<-plot.opt$xlab
    tit<-plot.opt$main
    if(colSums(saemixObject["model"]["covariate.model"])[ipar]>0) {
      idcov<-which(saemixObject["model"]["covariate.model"][,ipar]==1)
      for(icov in idcov) {
        xcov<-plot.opt$cov.value[icov]
        nlev<-length(unique(saemixObject["model"]["Mcovariates"][,(icov+1)]))
# covariable binaire
        if(is.na(xcov)) {
        if(nlev==2) 
          xcov<-min(saemixObject["model"]["Mcovariates"][,(icov+1)]) else 
# covariable continue
          xcov<-median(saemixObject["model"]["Mcovariates"][,(icov+1)])
        }
# ECO TODO securiser + dans le cas binaire, faire les 2 distributions si xcov=="all"
# Verifier le code ci-dessous
        ig1<-grep(saemixObject["model"]["name.cov"][idcov], saemixObject["model"]["name.fixed"])
        ig2<-grep(nampar[ipar],saemixObject["model"]["name.fixed"])
        iig<-c(ig1,ig2)
   mpar[ipar]<-mpar[ipar]+xcov*saemixObject["results"]["betas"][iig[duplicated(iig)]]
        if(!change.main) {
        if(tit!="") sep1<-"-" else sep1<-""
	xunit<-saemixObject["data"]["units"]$covariates[icov]
        tit<-paste(tit,sep1,saemixObject["model"]["name.cov"][icov],"=",xcov, ifelse(xunit=="-","",xunit),sep="")
        }
      }
      if(length(idcov)>0 & plot.opt$indiv.histo & (saemixObject@options$warnings) ) cat("Warning: histograms of individual parameter estimates do not make sense since covariates enter the model for parameter",nampar[ipar],"\n")
    }
    xpl1<-mpar[ipar]+xpl*sqrt(diag(saemixObject["results"]["omega"]))[ipar]
    xpl2<-transphi(matrix(xpl1,ncol=1),saemixObject["model"]["transform.par"][ipar])
    if(saemixObject["model"]["transform.par"][ipar]==2) {
      ypl<-pnorm(xpl)*derivphi(matrix(xpl1,ncol=1), saemixObject["model"]["transform.par"][ipar])
    } else
      ypl<-dnorm(xpl)*derivphi(matrix(xpl1,ncol=1), saemixObject["model"]["transform.par"][ipar])
    if(plot.opt$indiv.histo) {
      vec<-c(indiv.par[,(ipar+1)],xpl2)
      limx<-c(min(vec),max(vec))
    } else limx<-c(min(xpl2),max(xpl2))
    if(limx[1]<0) limx[1]<-limx[1]*1.05 else limx[1]<-limx[1]*0.95
    if(limx[2]>0) limx[2]<-limx[2]*1.05 else limx[2]<-limx[2]*0.95
    if(plot.opt$indiv.histo) {
      laby<-"Counts"
      if(change.ylab) laby<-plot.opt$ylab
      h1<-hist(indiv.par[,(ipar+1)],xlim=limx,main=tit,xlab=labx,ylab=laby, col=plot.opt$fillcol)
      ypl<-ypl/max(ypl)*max(h1$counts)
      lines(xpl2,ypl,lty=plot.opt$lty,col=plot.opt$lcol,lwd=plot.opt$lwd)
    } else {
      laby<-"Density"
      if(change.ylab) laby<-plot.opt$ylab
      plot(xpl2,ypl,type="l",xlab=labx,ylab=laby,xlim=limx, main=tit,lty=plot.opt$lty,col=plot.opt$lcol,lwd=plot.opt$lwd,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
      }
  }
}

#######################	   Parameters versus covariates  ##########################

saemix.plot.parcov<-function(saemixObject,...) {
# non-user level function
# parameters versus covariates
  saemix.plot.parcov.aux(saemixObject,partype="p",...)  
}

saemix.plot.randeffcov<-function(saemixObject,...) {
# non-user level function
# random effects versus covariates
  saemix.plot.parcov.aux(saemixObject,partype="r",...)  
}

saemix.plot.parcov.aux<-function(saemixObject,partype="p",...) {
# Plot of parameters (parytype="p") or random effects ("r") versus covariates
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1
  plot.opt<-saemixObject["prefs"]
  plot.opt$which.par<-"all"
  plot.opt$which.cov<-"all"
  plot.opt$main<-""
  plot.opt<-replace.plot.options(plot.opt,...)
  change.xlab<-FALSE
  if(plot.opt$xlab!=saemixObject["prefs"]$xlab) change.xlab<-TRUE
  change.ylab<-FALSE
  if(plot.opt$ylab!=saemixObject["prefs"]$ylab) change.ylab<-TRUE
  nampar<-saemixObject["model"]["name.modpar"]
  namcov<-saemixObject["data"]["name.covariates"]
  if(plot.opt$which.par[1]=="all") plot.opt$which.par<-nampar
  if(plot.opt$which.cov[1]=="all") plot.opt$which.cov<-namcov
  if(!is.integer(plot.opt$which.par)) plist<-match(plot.opt$which.par,nampar)
  plist<-plist[!is.na(plist)]
  if(!is.integer(plot.opt$which.cov)) clist<-match(plot.opt$which.cov,namcov)
  clist<-clist[!is.na(clist)]
  if(length(plist)==0) {
    if(saemixObject@options$warnings) cat("Cannot find parameter",plot.opt$which.par,", please check parameter names\n")
    return("Parameter not found")
  }
  if(length(clist)==0) {
    if(saemixObject@options$warnings) cat("Cannot find covariates",plot.opt$which.cov,", please check covariate names\n")
    return("Covariate not found")
  }
  replot<-FALSE
  mfrow<-plot.opt$mfrow
  if(plot.opt$new) {
    if(length(plot.opt$mfrow)==0) {
      if(length(plist)>1 & length(clist)>1) replot<-TRUE
      if(length(clist)>1) np<-length(clist) else np<-length(plist)   
      n1<-round(sqrt(np))
      n2<-ceiling(np/n1)
      mfrow<-c(n1,n2)
    }
    if(!replot) par(mfrow=mfrow,ask=plot.opt$ask)
  }
  # ECO removed column Id from all parameter dataframes
  if(partype=="r") { # random effects versus covariates
  if(tolower(plot.opt$indiv.par)=="map") {
    if(length(saemixObject["results"]["map.eta"])==0) {
      if(saemixObject@options$warnings) cat("Computing ETA estimates and adding them to fitted object.\n")
      saemixObject<-compute.eta.map(saemixObject)
    }
    param<-saemixObject["results"]["map.eta"]
  } else 
    param<-saemixObject["results"]["cond.mean.phi"]
  } else { # parameters versus covariates
    if(tolower(plot.opt$indiv.par)=="map") 
      param<-saemixObject["results"]["map.psi"] else 
        param<-transphi(saemixObject["results"]["cond.mean.phi"], saemixObject["model"]["transform.par"])
  }

# ECO: will not work with IOV  
  id<-saemixObject["data"]["data"][,saemixObject["data"]["name.group"]]
  idlist<-unique(id)
  matcov<-saemixObject["data"]["data"][match(idlist,id), saemixObject["data"]["name.covariates"],drop=FALSE]
  for(ipar in plist) {
    if(replot) par(mfrow=mfrow,ask=plot.opt$ask) # new page for each parameter (only if several covariates & several parameters, & plot.new==TRUE)
    xpar<-param[,ipar]
    laby<-nampar[ipar] 
    if(partype=="r") laby<-paste("ETA(",laby,")",sep="")
    if(change.ylab) laby<-plot.opt$ylab
    for(icov in clist) {
      covar<-matcov[,icov]
      labx<-saemixObject["data"]["name.covariates"][icov]
      if(change.xlab) labx<-plot.opt$xlab
      if(length(unique(covar))<=2) {
        boxplot(xpar~covar,xlab=labx,ylab=laby,main=plot.opt$main,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
        if (length(grep("m",plot.opt$line.smooth))>0) {
        y1<-saemixObject["results"]["fixed.psi"][ipar]
        abline(h=y1,lty=plot.opt$ablinelty,col=plot.opt$ablinecol, lwd=plot.opt$ablinelwd)
        }
       } else {
         plot(covar,xpar,xlab=labx,ylab=laby,main=plot.opt$main,pch=plot.opt$pch, col=plot.opt$col,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
         if (length(grep("l",plot.opt$line.smooth))>0) {
          y1<-lm(xpar~covar)
          abline(y1,lty=plot.opt$ablinelty,col=plot.opt$ablinecol, lwd=plot.opt$ablinelwd)
         }
         if (length(grep("s",plot.opt$line.smooth))>0) {
          lines(lowess(covar[!is.na(covar)],xpar[!is.na(covar)]), lty=plot.opt$ablinelty,col=plot.opt$ablinecol, lwd=plot.opt$ablinelwd)
         }
         if (length(grep("m",plot.opt$line.smooth))>0) {
          y1<-saemixObject["results"]["fixed.psi"][ipar]
          abline(h=y1,lty=plot.opt$ablinelty,col=plot.opt$ablinecol,lwd=plot.opt$lwd)
       }
     }
    }
  }
#  return()
}
