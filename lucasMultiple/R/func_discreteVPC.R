###########################	VPC for non Gaussian data models		#############################

#' VPC for non Gaussian data models
#' 
#' This function provides VPC plots for non Gaussian data models (work in progress)
#' 
#' @param object an saemixObject object returned by the \code{\link{saemix}} function.
#' The object must include simulated data under the empirical design, using the model and 
#' estimated parameters from a fit, produced via the \code{\link{simulateDiscreteSaemix}} function.
#' @param outcome type of outcome (valid types are "TTE", "binary", "categorical", "count")
#' @param verbose whether to print messages (defaults to FALSE)
#' @param \dots additional arguments, used to pass graphical options (to be implemented, currently not available)
#' 
#' @details This function is a very rough first attempt at automatically creating VPC plots for 
#' models defined through their log-likelihood (categorical, count, or TTE data). It makes use of the
#' new element simulate.function in the model component of the object
#' - for TTE data, a KM-VPC plot will be produced
#' - for count, categorical and binary data, a plot showing the proportion of each score/category across time will be shown
#' along with the corresponding prediction intervals from the model
#' These plots can be stratified over a covariate in the data set (currently only categorical covariates) 
#' by passing an argument which.cov='name' to the call
#' 
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>
#' 
#' @seealso \code{\link{SaemixObject}}, \code{\link{saemix}},
#' \code{\link{saemix.plot.vpc}}, \code{\link{simulateDiscreteSaemix}}
#' 
#' @references Brendel, K, Comets, E, Laffont, C, Laveille, C, Mentre, F.
#' Metrics for external model evaluation with an application to the population
#' pharmacokinetics of gliclazide, Pharmaceutical Research 23 (2006),
#' 2036-2049.
#' 
#' Holford, N. The Visual Predictive Check: superiority to standard diagnostic
#' (Rorschach) plots (Abstract 738), in: 14th Meeting of the Population
#' Approach Group in Europe, Pamplona, Spain, 2005.
#' 
#' Ron Keizer, tutorials on VPC TODO
#' 
#' @keywords plot
#' 
#' @examples
#' data(lung.saemix)
#' 
#' saemix.data<-saemixData(name.data=lung.saemix,header=TRUE,name.group=c("id"),
#' name.predictors=c("time","status","cens"),name.response=c("status"),
#' name.covariates=c("age", "sex", "ph.ecog", "ph.karno", "pat.karno", "wt.loss","meal.cal"),
#' units=list(x="days",y="",covariates=c("yr","","-","%","%","cal","pounds")))
#' 
#' weibulltte.model<-function(psi,id,xidep) {
#'   T<-xidep[,1]
#'   y<-xidep[,2] # events (1=event, 0=no event)
#'   cens<-which(xidep[,3]==1) # censoring times (subject specific)
#'   init <- which(T==0)
#'   lambda <- psi[id,1] # Parameters of the Weibull model
#'   beta <- psi[id,2]
#'   Nj <- length(T)
#'   ind <- setdiff(1:Nj, append(init,cens)) # indices of events
#'   hazard <- (beta/lambda)*(T/lambda)^(beta-1) # H'
#'   H <- (T/lambda)^beta # H
#'   logpdf <- rep(0,Nj) # ln(l(T=0))=0
#'   logpdf[cens] <- -H[cens] + H[cens-1] # ln(l(T=censoring time))
#'   logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind]) # ln(l(T=event time))
#'   return(logpdf)
#' }
#' 
#' simulateWeibullTTE <- function(psi,id,xidep) {
#'   T<-xidep[,1]
#'   y<-xidep[,2] # events (1=event, 0=no event)
#'   cens<-which(xidep[,3]==1) # censoring times (subject specific)
#'   init <- which(T==0)
#'   lambda <- psi[,1] # Parameters of the Weibull model
#'   beta <- psi[,2]
#'   tevent<-T
#'   Vj<-runif(dim(psi)[1])
#'   tsim<-lambda*(-log(Vj))^(1/beta) # nsuj events
#'   tevent[T>0]<-tsim
#'   tevent[tevent[cens]>T[cens]] <- T[tevent[cens]>T[cens]]
#'   return(tevent)
#'   }
#' saemix.model<-saemixModel(model=weibulltte.model,description="time model",modeltype="likelihood",
#'        simulate.function = simulateWeibullTTE,
#'                 psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("lambda","beta"))),
#'                 transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE))
#' saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)
#' 
#' \donttest{
#' tte.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' simtte.fit <- simulateDiscreteSaemix(tte.fit, nsim=100)
#' gpl <- discreteVPC(simtte.fit, outcome="TTE")
#' plot(gpl)
#' }
#' 
#' @aliases discreteVPCcount discreteVPCcat discreteVPC.aux
#' 
#' @export 

discreteVPC <- function(object, outcome="categorical", verbose=FALSE, ...) {
  # object
  ## an object resulting from a call to saemix
  ## must have a sim.data object containing data simulated under the model
  # outcome: type of outcome (valid types are TTE, binary, categorical, count)
  # verbose: whether to print messages
  # ...: to pass additional plot options (such as title, colours, etc...) which will supersede the options in the prefs slot of object
  if(!is(object,"SaemixObject")) {
    message("Please provide a valid object \n")
    return()
  }
  outcome<-tolower(outcome)
  if(outcome=="event") outcome<-"tte"
  if(is.na(match(outcome,c("continuous","tte","count","categorical","binary")))) {
    if(verbose) message("Please specify a valid outcome type (TTE, count, categorical, binary). For continuous data, VPC can be obtained via the npde package loaded along with saemix.\n")
    return()
  }
  if(object@sim.data@nsim==0) {
    if(verbose) message("Please run simulations under the model using simulateDiscreteSaemix with a matching simulation function in the model component of the object.\n")
    return()
  } else {
    if(object@sim.data@nsim<100 & verbose) message(paste0("Number of simulations in the sim.data component (",object@sim.data@nsim,") is too low, consider increasing\n"))
  }
  if(outcome=="tte") {
    # TODO: set ngrid from plot options
    xplot <- discreteVPCTTE(object, ngrid=200, verbose=verbose, ...)
  }
  if(outcome=="count")
    xplot <- discreteVPCcount(object, verbose=verbose, ...)
  if(outcome %in% c("binary","categorical")) {
    if(outcome=="binary") {
      max.cat<-2
    } else {
      max.cat<-length(unique(object@data@data[,object@data@name.response]))
    }
    xplot <- discreteVPCcat(object, max.cat=max.cat, verbose=verbose, ...)
  }
  
  return(xplot)
}

###########################	VPC for TTE		#############################
# TODO: clarify the use of cens !!! Currently used as a censoring indicator but may be automatically filled  (ok if =0 when automatically filled, probably)
# TODO: RTTE (more than 1 event, do a VPC for first, second,... until maxevent)
# TODO: passing plot options + extending current graphical options to match those of npde (through npdeControl)

#' VPC for time-to-event models
#' 
#' This function provides VPC plots for time-to-event data models (work in progress)
#' 
#' @param object an saemixObject object returned by the \code{\link{saemix}} function.
#' The object must include simulated data under the empirical design,
#' using the model and estimated parameters from a fit, produced via the 
#' \code{\link{simulateDiscreteSaemix}} function
#' @param ngrid number of grid points on the X-axis to extrapolate the KM-VPC
#' @param interpolation.method method to use for the interpolation of the KM for the simulated datasets. Available methods are 
#' "step": the value of the survival function for a given grid point is set to the value of the last time
#' "lin": a linear approximation is used between two consecutive times (defaults to "step")
#' @param verbose whether to print messages (defaults to FALSE)
#' @param \dots additional arguments, used to pass graphical options (to be implemented, currently not available)
#' 
#' @details add details on TTE VPC, RTTE VPC, etc...
#' 
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>
#' 
#' @seealso \code{\link{SaemixObject}}, \code{\link{saemix}},
#' \code{\link{saemix.plot.vpc}}, \code{\link{simulateDiscreteSaemix}}
#' 
#' Tutorials on TTE-VPC TODO
#' 
#' @keywords plot
#' 
#' @aliases interpol.locf interpol.lin
#' 
#' @importFrom scales trans_format  math_format pretty_breaks trans_breaks
#' @export 

discreteVPCTTE <- function(object, ngrid=200, interpolation.method="step", verbose=FALSE, ...) {
  # object: a saemixObject including simulated data
  # ngrid: number of grid points for interpolation of the PI
  # verbose: whether to print messages
  # ...: to pass additional plot options (such as title, colours, etc...)
  
  # Creates:
  # obsdat: dataframe with the observation data
  ## id: subject id
  ## time: event time
  ## event: 0 for no event, 1 for event
  ## cens: 0 for observed event/no event, 1 for censored (optional, if not present all 1 for events are assumed to be observed events and a 0 at the last observation time indicates a censored event)
  # simdat:
  ## irep: simulation number
  ## id: subject id (replicated)
  ## time: simulated times
  obsdat <- data.frame(id=object@data@data[,object@data@name.group], time=object@data@data[,object@data@name.X], event=object@data@data[,object@data@name.response])
  if(length(object@data@name.cens)>0)
    obsdat$cens <- object@data@data[,object@data@name.cens]
  simdat <- data.frame(irep=object@sim.data@datasim$irep, id=object@sim.data@datasim$idsim, time=object@sim.data@datasim$ysim)
  # Observed KM
  event.obs<-obsdat[obsdat$time>0,]
  t1<-table(event.obs$time) # nb of subjects dropping off at different times (event or censoring)
  tab.obs <- data.frame(tobs=as.numeric(names(t1)),nk=c(t1))
  if(match("cens", colnames(obsdat)) && sum(obsdat$cens)>0) { # cens indicates censoring (cens=1) versus event (cens=0)
    t2<-table(event.obs$time[event.obs$cens==0]) # nb of subjects with an event 
    tab.obs$dk<-0
    idx<-match(names(t2), names(t1))
    tab.obs$dk[idx]<-t2
    idx.cens <- which(tab.obs$nk != tab.obs$dk)
  } else { # all events are observed (no cens or all cens=0)
    tab.obs$dk<-tab.obs$nk
    idx.cens <- c()
  }
  nk.obs <- c(sum(tab.obs$nk), sum(tab.obs$nk)-cumsum(tab.obs$nk))
  tab.obs$nk <- nk.obs[1:(length(nk.obs)-1)]
  tab.obs$sk <- (1-tab.obs$dk/tab.obs$nk)
  tab.obs$km <- cumprod(tab.obs$sk)
  
  # Simulated KM
  ymax<-max(obsdat$time)
  ytime<- seq(0,ymax,length.out=ngrid)
  
  nsim<-object@sim.data@nsim
  km.sim<-NULL
  for (isim in 1:nsim) {
    xsim<-simdat[simdat$irep==isim,]
    tsim <- xsim$time[xsim$time>0]
    t1 <- table(tsim)
    tab1<-data.frame(time=as.numeric(names(t1)),dk=c(t1))
    nk.sim <- c(sum(tab1[,2]), sum(tab1[,2])-cumsum(tab1[,2]))
    tab1$nk <- nk.sim[1:(length(nk.sim)-1)]
    tab1$sk <- (1-tab1$dk/tab1$nk)
    tab1$km <- cumprod(tab1$sk)
    if(interpolation.method=="lin") 
      km <- sapply(ytime, interpol.lin, tab1) else 
        km <- sapply(ytime, interpol.locf, tab1)
    km.sim<-cbind(km.sim, km)
  }
  
  # Quantiles 
  alpha <- object@prefs$vpc.interval
  tab.quant<-t(apply(km.sim, 1, quantile, c((1-alpha)/2, 0.5, 1-(1-alpha)/2)))
  tab.quant <- cbind(time=ytime,tab.quant)
  colnames(tab.quant)<-c("time","lower","median","upper")
  tab.quant <- as.data.frame(tab.quant)
  
  # Plot options - TODO
  plot.opt<-object@prefs
  plot.opt$bands <- TRUE # plot PI
  plot.opt$alpha.bands <-0.3
  plot.opt$col.pcens <- "steelblue3"
  plot.opt$breaks.x <- plot.opt$breaks.y <- 10
  plot.opt$plot.censTTE <- TRUE
  
  plot.opt$ylab <- "Survival (%)"
  # Future: stratify over eg covariates or groups ?
  numberCategories <- 1
  
  # Creating plot 
  ## TODO: add options
  ## TODO: add censored data (eg red dots to indicate censored events in the KM plot)
  xplot <- ggplot(data=tab.quant, aes(x=.data$time, y=.data$median)) +
    # theme of the ggplot template - from npde, to integrate here
    # theme(plot.title = element_text(hjust = 0.5, size = plot.opt$size.sub),
    #       axis.title.x = element_text(size = plot.opt$size.xlab),
    #       axis.title.y = element_text(size = plot.opt$size.ylab),
    #       axis.text.x = element_text(size=plot.opt$size.text.x),
    #       axis.text.y = element_text(size=plot.opt$size.text.y),
    #       axis.line.x = element_line(color=ifelse(plot.opt$xaxt==TRUE,"black","white")),
    #       axis.line.y = element_line(color=ifelse(plot.opt$yaxt==TRUE,"black","white")),
    #       panel.background=element_rect("white"),
    #       panel.grid.major.x = element_line(ifelse(plot.opt$grid==TRUE,"grey80","white"),linetype = plot.opt$lty.grid),
    #       panel.grid.minor.x = element_line(ifelse(plot.opt$grid==TRUE,"grey80","white"),linetype = plot.opt$lty.grid),
    #       panel.grid.major.y = element_line(ifelse(plot.opt$grid==TRUE,"grey80","white"),linetype = plot.opt$lty.grid),
    #       panel.grid.minor.y = element_line(ifelse(plot.opt$grid==TRUE,"grey80","white"),linetype = plot.opt$lty.grid))  +
    
    # coordinates x-y
    # coord_cartesian(xlim=x.limits, ylim=y.limits) +
    { if ( plot.opt$bands == TRUE )
      geom_ribbon(aes(ymin = .data$lower, ymax = .data$upper), fill = plot.opt$fillcol, alpha = plot.opt$alpha.bands) } +
    geom_line(linetype = plot.opt$lty.lpi,colour = plot.opt$col.lpi, size = plot.opt$lwd.lpi)+
    geom_line(aes(y=.data$lower), colour = plot.opt$col.lpi, size = plot.opt$lwd.lpi, alpha = plot.opt$alpha.bands)+
    geom_line(aes(y=.data$upper), colour = plot.opt$col.lpi, size = plot.opt$lwd.lpi, alpha = plot.opt$alpha.bands)+
    geom_line(data=tab.obs,aes(x=.data$tobs, y = .data$km),  linetype = plot.opt$lty.lobs,colour = plot.opt$col.lobs,size = plot.opt$lwd.lobs)+
    # Censored events
    {if(length(idx.cens)>0 & plot.opt$plot.censTTE)
      geom_point(data=tab.obs[idx.cens,],aes(x=.data$tobs, y=.data$km), colour=plot.opt$col.pcens)} +
    # x-y log-scales
    { if (plot.opt$xlog == FALSE)  scale_x_continuous(plot.opt$xlab, scales::pretty_breaks(n = plot.opt$breaks.x))
    } +
    { if (plot.opt$ylog == FALSE)  scale_y_continuous(plot.opt$ylab, scales::pretty_breaks(n = plot.opt$breaks.y))
    } +
    { if (plot.opt$xlog == TRUE)
      scale_x_log10(plot.opt$xlab,breaks = scales::trans_breaks("log10", function(x) 10 ^ x),
                    labels = scales::trans_format("log10", scales::math_format(10 ^ .x)))
    } +
    { if (plot.opt$ylog == TRUE)
      scale_y_log10(plot.opt$ylab, breaks = scales::trans_breaks("log10", function(x) 10 ^ x),
                    labels = scales::trans_format("log10", scales::math_format(10 ^ .x)))
    } +
    # if log scales plot logticks
    { if (plot.opt$xlog == TRUE) annotation_logticks(sides = "b")} +
    { if (plot.opt$ylog == TRUE) annotation_logticks(sides = "l")} +
    
    {if (plot.opt$main!="") ggtitle(plot.opt$main)} +
    
    # facet wrap over covariate categories
#    facet_wrap(.~factor(category), nrow=1, scales=plot.opt$scales) +
    
    {if(numberCategories==1)
      theme(strip.background = element_blank(), strip.text.x = element_blank())
    }
  return(xplot)
}


###########################	VPC for Count data		#############################
# TODO: plot options (colours, line types, ...)

discreteVPCcount <- function(object, max.cat=10, breaks=NULL, catlabel=NULL, verbose=FALSE, ...) {
  # object: a saemixObject including simulated data
  # ngrid: number of grid points for interpolation of the PI
  # verbose: whether to print messages
  # ...: to pass additional plot options (such as title, colours, etc...)
  
  # Creates:
  # obsdat: dataframe with the observation data
  ## id: subject id
  ## time: event time
  ## event: 0 for no event, 1 for event
  ## cens: 0 for observed event/no event, 1 for censored (optional, if not present all 1 for events are assumed to be observed events and a 0 at the last observation time indicates a censored event)
  # simdat:
  ## irep: simulation number
  ## id: subject id (replicated)
  ## y: simulated counts

  # if(plot.opt$new) {
  #   mfrow<-plot.opt$mfrow
  #   if(length(mfrow)==0) mfrow<-c(1,1)
  #   par(mfrow=mfrow,ask=plot.opt$ask)
  # }
  object@prefs$ylab <- "Proportion of counts (-)"
  x<-discreteVPC.aux(object, max.cat=max.cat, breaks=breaks, verbose=verbose, ...)
  xtab<-x$xtab
  stab<-x$stab
  plot.opt<-x$plot.opt
  cgroups<-unique(xtab$group)
  # xgroups <- sort(unique(xtab$x.group))
  # sgroups <- sort(unique(xtab$y.group))
  
  # ggplot(data=stab, aes(x=x, ymin=lower, ymax=upper, group=as.factor(group), fill=as.factor(group))) + geom_ribbon(alpha=0.2) + 
  #   geom_line(data=stab, aes(x=x, y=median, group=as.factor(group), colour=as.factor(group)), linetype="dashed") + 
  #   geom_line(data=xtab, aes(x=x, y=freq, group=as.factor(group), colour=as.factor(group))) + 
  #   xlab(plot.opt$xlab) + ylab("Proportion of counts") + 
  #   {if(length(cgroups)>1) guides(fill=guide_legend(title='Gender'), colour=guide_legend(title='Gender'))} +
  #   facet_wrap(.~score, ncol=4)
  if(length(plot.opt$mfrow)==2) {
    nrow1 <- plot.opt$mfrow[1]
    ncol1 <- plot.opt$mfrow[2]
  } else {
    nrow1<-ncol1<-NULL
  }
  # ggplot
  plot.counts <- ggplot(data=xtab, aes(x=.data$x.group, y=.data$freq, group=as.factor(.data$group), colour=as.factor(.data$group))) + geom_line() + 
    geom_line(data=stab, aes(x=.data$x.group, y=.data$median, group=as.factor(.data$group), colour=as.factor(.data$group)), linetype="dashed") + 
    geom_ribbon(data=stab, aes(x=.data$x.group, ymin=.data$lower, ymax=.data$upper, group=as.factor(.data$group), fill=as.factor(.data$group)), alpha=0.2, linetype=0) + 
    xlab(plot.opt$xlab) + ylab(plot.opt$ylab) + 
    {if(length(cgroups)==1) theme(legend.position="none")} +
    {if(length(cgroups)>1) guides(fill=guide_legend(title=plot.opt$which.cov), colour=guide_legend(title=plot.opt$which.cov))} +
    facet_wrap(.~.data$y.group, nrow=nrow1, ncol=ncol1)
  return(plot.counts)
}

###########################	VPC for Categorical data		#############################

discreteVPCcat <- function(object, max.cat=10, breaks=NULL, catlabel=NULL, verbose=FALSE, ...) {
  # object: a saemixObject including simulated data
  # ngrid: number of grid points for interpolation of the PI
  # verbose: whether to print messages
  # ...: to pass additional plot options (such as title, colours, etc...)
  
  # Creates:
  # obsdat: dataframe with the observation data
  ## id: subject id
  ## time: event time
  ## event: 0 for no event, 1 for event
  ## cens: 0 for observed event/no event, 1 for censored (optional, if not present all 1 for events are assumed to be observed events and a 0 at the last observation time indicates a censored event)
  # simdat:
  ## irep: simulation number
  ## id: subject id (replicated)
  ## y: simulated counts
  
  # if(plot.opt$new) {
  #   mfrow<-plot.opt$mfrow
  #   if(length(mfrow)==0) mfrow<-c(1,1)
  #   par(mfrow=mfrow,ask=plot.opt$ask)
  # }
  if(!is.null(breaks))
    max.cat<-length(breaks)
  if(max.cat==2) {
    object@prefs$ylab <- "Proportion of events (-)"
  } else {
    object@prefs$ylab <- "Proportion of category (-)"
  }
  x<-discreteVPC.aux(object, max.cat=max.cat, breaks=breaks, catlabel=catlabel, verbose=verbose, ...)
  xtab<-x$xtab
  stab<-x$stab
  plot.opt<-x$plot.opt
  cgroups<-unique(xtab$group)
  # xgroups <- sort(unique(xtab$x.group))
  # sgroups <- sort(unique(xtab$y.group))
  
  # ggplot(data=stab, aes(x=x, ymin=lower, ymax=upper, group=as.factor(group), fill=as.factor(group))) + geom_ribbon(alpha=0.2) + 
  #   geom_line(data=stab, aes(x=x, y=median, group=as.factor(group), colour=as.factor(group)), linetype="dashed") + 
  #   geom_line(data=xtab, aes(x=x, y=freq, group=as.factor(group), colour=as.factor(group))) + 
  #   xlab(plot.opt$xlab) + ylab("Proportion of counts") + 
  #   {if(length(cgroups)>1) guides(fill=guide_legend(title='Gender'), colour=guide_legend(title='Gender'))} +
  #   facet_wrap(.~score, ncol=4)
  if(length(plot.opt$mfrow)==2) {
    nrow1 <- plot.opt$mfrow[1]
    ncol1 <- plot.opt$mfrow[2]
  } else {
    nrow1<-ncol1<-NULL
  }
  plot.cat <- ggplot(data=xtab, aes(x=.data$x.group, y=.data$freq, group=as.factor(.data$group), colour=as.factor(.data$group))) + geom_line() + 
    geom_line(data=stab, aes(x=.data$x.group, y=.data$median, group=as.factor(.data$group), colour=as.factor(.data$group)), linetype="dashed") + 
    geom_ribbon(data=stab, aes(x=.data$x.group, ymin=.data$lower, ymax=.data$upper, group=as.factor(.data$group), fill=as.factor(.data$group)), alpha=0.2, linetype=0) + 
    xlab(plot.opt$xlab) + ylab(plot.opt$ylab) + 
    {if(length(cgroups)==1) theme(legend.position="none")} +
    {if(length(cgroups)>1) guides(fill=guide_legend(title=plot.opt$which.cov), colour=guide_legend(title=plot.opt$which.cov))} +
    facet_wrap(.~.data$y.group, nrow=nrow1, ncol=ncol1)
  return(plot.cat)
}


###########################	Auxiliary functions		#############################

# Interpolation on a grid - step function or linear extrapolation
interpol.locf <- function(x, y) {
  ykm<-y$km[y[,1]<=x]
  nminus <- length(ykm)
  if(nminus>0) return(ykm[nminus]) else return(1)
}
interpol.lin <- function(x, y) {
  # y must contain time and km
  nminus <- length(y$km[y$time<=x])
  if(nminus==0) return(1)
  if(nminus==dim(y)[1]) return(y$km[nminus])
  y$km[nminus] + (x-y$time[nminus])*(y$km[nminus+1]-y$km[nminus])/(y$time[nminus+1]-y$time[nminus])
}

# Binning and computing PI
discreteVPC.aux <- function(object, max.cat=10, breaks=NULL, catlabel=NULL, verbose=FALSE, ...) {
  # object: a saemixObject including simulated data
  # ngrid: number of grid points for interpolation of the PI
  # verbose: whether to print messages
  # ...: to pass additional plot options (such as title, colours, etc...)
  
  # Creates:
  # obsdat: dataframe with the observation data
  ## id: subject id
  ## time: event time
  ## event: 0 for no event, 1 for event
  ## cens: 0 for observed event/no event, 1 for censored (optional, if not present all 1 for events are assumed to be observed events and a 0 at the last observation time indicates a censored event)
  # simdat:
  ## irep: simulation number
  ## id: subject id (replicated)
  ## y: simulated counts
  
  # Returns
  ## xtab: a dataframe with columns
  ### group: covariate/stratification group
  ### x.group: independent variable (usually time), grouped
  ### y.group: dependent variable (category, score), grouped
  ### freq: frequency of y.group at the corresponding combination of group and x.group
  ### nobs: number of observations at the corresponding combination of group and x.group
  ## stab:
  ### group, x.group, y.group: same as for xtab
  ### lower, median, upper: lower, median and upper values of the PI (upper and lower correspond to the PI defined by the plot.opt$vpc.interval)
  ## updated plot.opt
  ydat <- object@data
  ysim <- object@sim.data@datasim$ysim
  obsdat<-data.frame(id=ydat@data[,ydat@name.group], x=ydat@data[,ydat@name.X], y=ydat@data[,ydat@name.response])
  #  if(length(ydat@name.covariates)>0) obsdat<-cbind(obsdat, covariate.group=ydat@data[,ydat@name.covariates,drop=FALSE])
  nsim<-length(ysim)/dim(ydat@data)[1]
  #  simdat <- data.frame(irep=object@sim.data@datasim$irep, id=object@sim.data@datasim$idsim, y=object@sim.data@datasim$ysim)
  
  plot.opt<-object["prefs"]
  plot.opt$main<-"Visual Predictive Check"
  plot.opt<-replace.plot.options(plot.opt,...)
  logtyp<-""
  if(plot.opt$xlog) logtyp<-paste(logtyp,"x",sep="")
  if(plot.opt$ylog) logtyp<-paste(logtyp,"y",sep="")
  
  plot.opt$bin.number <- plot.opt$vpc.bin # TODO: change and harmonise with npde
  if(is.null(plot.opt$bin.method)) plot.opt$bin.method<- "equal"
  if(plot.opt$bin.method=="user" & is.null(plot.opt$bin.breaks)) {
    if(verbose) message("When using the bin.method='user' option, please provide the breaks in bin.breaks\n")
    plot.opt$bin.method<- "equal"
  }
  # Create covariate groups - TODO
  if(length(plot.opt$which.cov)!=1 || length(grep(plot.opt$which.cov, ydat@name.covariates))==0)
    obsdat$covariate.group <- 'all' else obsdat$covariate.group <- ydat@data[,plot.opt$which.cov]
  
  # Binning on X
  if(length(unique(obsdat$x))>plot.opt$vpc.bin) {
    xbnd <- xbinning(obsdat$x, plot.opt, verbose=verbose)
    #    obsdat$x.group<-xbnd$xgrp
    obsdat$x.group<-xbnd$xcent[xbnd$xgrp]
  } else {
    xgrp <- sort(unique(obsdat$x))
    #    obsdat$x.group<-match(obsdat$x, xgrp)
    obsdat$x.group <- obsdat$x
    #    xbreaks<-c(sort(unique(obsdat$x))-0.01, max(obsdat$x)+0.01)
  }
  
  # Binning on score
  if(length(unique(obsdat$y))<=max.cat & is.null(breaks)) {
    max.cat<-length(unique(obsdat$y))
    breaks <- sort(unique(obsdat$y))
    obsdat$score.group<-obsdat$y
    sim.score<-ysim
  } else {
    #    mybreaks <- c(0:9, 16, 25, 80)
    if(is.null(breaks)) {
      if(is.null(max.cat) || max.cat<=1 | max.cat>20) {
        if(verbose) message("Please give a valid number of categories to bin scores (between 2 and 20)\n")
        max.cat<-10
      }
      breaks<-unique(quantile(obsdat$y, seq(0,1, length.out=(max.cat+1))))
      breaks[1]<-breaks[1]-0.01
      #      print(breaks)
      if(max(ysim)>breaks[length(breaks)])  breaks[length(breaks)]<-max(ysim)
      xgrp <- cut(obsdat$y, breaks=breaks, include.lowest=TRUE, right=FALSE)
      sim.score<-cut(ysim, breaks=breaks, include.lowest = TRUE, right=FALSE)
    } else {
      breaks<-sort(unique(breaks))
      if(min(breaks)>min(obsdat$y)) breaks<-c(min(obsdat$y), breaks)
      if(max(breaks)<max(obsdat$y, ysim)) breaks<-c(breaks,max(obsdat$y, ysim))
      breaks[1]<-breaks[1]-0.01 # make sure all observations are included
      xgrp <- cut(obsdat$y, breaks=breaks, include.lowest=TRUE, right=FALSE)
      sim.score<-cut(ysim, breaks=breaks, include.lowest = TRUE, right=FALSE)
    }
    obsdat$score.group<-xgrp
  }
  # Computing frequency of each score at each time in each covariate group - observed data
  # Computing corresponding PI - simulated data
  xtab<-stab<-NULL
  cgroups <- unique(obsdat$covariate.group)
  xgroups <- sort(unique(obsdat$x.group))
  sgroups <- sort(unique(obsdat$score.group))
  #  print(sgroups)
  #  print(xgroups)
  simdat <- data.frame(irep=rep(1:nsim, each=dim(ydat@data)[1]), x.group=rep(obsdat$x.group, nsim), covariate.group=rep(obsdat$covariate.group, nsim), score.group=sim.score)
  alpha <- (1-plot.opt$vpc.interval)/2
  for(igroup in cgroups) {
    tab1<-table(obsdat$x.group[obsdat$covariate.group==igroup],obsdat$score.group[obsdat$covariate.group==igroup])
    ncol<-colSums(tab1)
    nobs<-rowSums(tab1)
    tab1<-tab1[,(ncol>0),drop=FALSE]
    freqtab <- tab1/nobs
    xtab <- rbind(xtab,
                  data.frame(group=igroup, x.group=rep(xgroups,length(sgroups)), y.group=rep(sgroups, each=length(xgroups)), freq=c(freqtab), nobs=rep(nobs, length(sgroups))))
    # array with 3 dimensions for the simulated data
    tab1<-table(simdat$x.group[simdat$covariate.group==igroup],simdat$score.group[simdat$covariate.group==igroup], simdat$irep[simdat$covariate.group==igroup])
    for(irep in unique(simdat$irep)) {
      tab1[,,irep] <- tab1[,,irep]/nobs # frequencies
    }
    quant<-apply(tab1, c(1,2), quantile, c(alpha, 0.5, 1-alpha))
    #    print(quant)
    for(iscr in 1:length(sgroups)) {
      #      cat("score group",sgroups[iscr],"\n")
      #      print(quant[,,iscr])
      stab<-rbind(stab,data.frame(group=igroup, x.group=xgroups, y.group=sgroups[iscr], t(quant[,,match(sgroups[iscr],dimnames(quant)[[3]])])))
    }
  }
  colnames(stab)[4:6]<-c("lower","median","upper")
  # ggplot(data=xtab, aes(x=x, y=freq, group=group, colour=as.factor(group))) + geom_line() + facet_wrap(.~score, ncol=4)
  stab$freq<-0
  # More pleasing labels
  if(!is.null(catlabel) && length(catlabel)==length(levels(xtab$y.group)))
    levels(xtab$y.group)<-levels(stab$y.group)<-catlabel else {# try to define nicer names...
      if(max.cat==length(unique(xtab$y.group)) & max.cat==length(unique(ysim))) {
        catlabel1<-sort(unique(obsdat$y))
      } else {
        catlabel<-gsub("[","",levels(xtab$y.group), fixed=TRUE)
        catlabel<-gsub(")","",catlabel, fixed=TRUE)
        catlabel<-gsub("]","",catlabel, fixed=TRUE)
        catlabel<-matrix(as.numeric(unlist(strsplit(catlabel,",",fixed=TRUE))), ncol=2, byrow=T)
        for(i in 1:2) catlabel[,i]<-round(catlabel[,i])
        if(catlabel[1,1]<0) catlabel[1,1]<-0
        catlabel1<-catlabel[,1]
        for(i in 1:length(catlabel1)) {
          if(catlabel[i,2]>(catlabel[i,1]+1)) catlabel1[i]<-paste0("[",catlabel[i,1],"-",catlabel[i,2]-1,"]")
        }
      }
      levels(xtab$y.group)<-levels(stab$y.group)<-catlabel1
    }
  
  return(list(stab=stab, xtab=xtab, plot.opt=plot.opt))
}

