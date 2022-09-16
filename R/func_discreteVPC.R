###########################	VPC for non Gaussian data models		#############################

#' VPC for non Gaussian data models
#' 
#' This function provides VPC plots for non Gaussian data models (work in progress)
#' 
#' @param object an saemixObject object returned by the \code{\link{saemix}} function.
#' The object must include simulated data under the empirical design,
#' using the model and estimated parameters from a fit, produced via the 
#' \code{\link{simulateDiscreteSaemix}} function
#' @param outcome type of outcome (valid types are "TTE", "binary", "categorical", "count")
#' @param verbose whether to print messages (defaults to FALSE)
#' @param \dots additional arguments, used to pass graphical options (to be implemented, currently not available)
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
#' @export 

discreteVPC <- function(saemixObject, outcome="TTE", verbose=FALSE, ...) {
  # saemixObject
  ## an SaemixObject resulting from a call to saemix
  ## must have a sim.data object containing data simulated under the model
  # outcome: type of outcome (valid types are TTE, binary, categorical, count)
  # verbose: whether to print messages
  # ...: to pass additional plot options (such as title, colours, etc...) which will supersede the options in the prefs slot of saemixObject
  
  if(!is(saemixObject,"SaemixObject")) {
    message("Please provide a valid SaemixObject \n")
    return()
  }
  if(is.na(match(outcome,c("TTE","count","categorical","binary")))) {
    if(verbose) message("Please specify a valid outcome type (TTE, count, categorical, binary). For continuous data, VPC can be obtained via the npde package loaded along with saemix.\n")
    return()
  }
  if(saemixObject@sim.data@nsim==0) {
    if(verbose) message("Please run simulations under the model using simulateDiscreteSaemix with a matching simulation function in the model component of the object.\n")
    return()
  } else {
    if(saemixObject@sim.data@nsim<100 & verbose) message(paste0("Number of simulations in the sim.data component (",saemixObject@sim.data@nsim,") is too low, consider increasing\n"))
  }
  if(outcome=="TTE") {
    # TODO: set ngrid from plot options
    xplot <- discreteVPCTTE(saemixObject, ngrid=200, verbose=verbose, ...)
  }
  if(outcome=="count")
    xplot <- discreteVPCcount(saemixObject, verbose=verbose, ...)
  if(outcome %in% c("binary","categorical"))
    xplot <- discreteVPCcat(saemixObject, verbose=verbose, ...)
  
  return(xplot)
}

###########################	VPC for TTE		#############################
# TODO: clarify the use of cens !!! Currently used as a censoring indicator but may be automatically filled  (ok if =0 when automatically filled, probably)
# passing plot options + extending current graphical options to match those of npde (through npdeControl)

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

discreteVPCTTE <- function(saemixObject, ngrid=200, interpolation.method="step", verbose=FALSE, ...) {
  # saemixObject
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
  obsdat <- data.frame(id=saemixObject@data@data[,saemixObject@data@name.group], time=saemixObject@data@data[,saemixObject@data@name.X], event=saemixObject@data@data[,saemixObject@data@name.response])
  if(length(saemixObject@data@name.cens)>0)
    obsdat$cens <- saemixObject@data@data[,saemixObject@data@name.cens]
  simdat <- data.frame(irep=saemixObject@sim.data@datasim$irep, id=saemixObject@sim.data@datasim$idsim, time=saemixObject@sim.data@datasim$ysim)
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
  
  nsim<-saemixObject@sim.data@nsim
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
  alpha <- saemixObject@prefs$vpc.interval
  tab.quant<-t(apply(km.sim, 1, quantile, c((1-alpha)/2, 0.5, 1-(1-alpha)/2)))
  tab.quant <- cbind(time=ytime,tab.quant)
  colnames(tab.quant)<-c("time","lower","median","upper")
  tab.quant <- as.data.frame(tab.quant)
  
  # Plot options - TODO
  plot.opt<-saemixObject@prefs
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
  xplot <- ggplot(data=tab.quant, aes(x=time, y=median)) +
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
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = plot.opt$fillcol, alpha = plot.opt$alpha.bands) } +
    geom_line(linetype = plot.opt$lty.lpi,colour = plot.opt$col.lpi, size = plot.opt$lwd.lpi)+
    geom_line(aes(y=lower), colour = plot.opt$col.lpi, size = plot.opt$lwd.lpi, alpha = plot.opt$alpha.bands)+
    geom_line(aes(y=upper), colour = plot.opt$col.lpi, size = plot.opt$lwd.lpi, alpha = plot.opt$alpha.bands)+
    geom_line(data=tab.obs,aes(x=tobs, y = km),  linetype = plot.opt$lty.lobs,colour = plot.opt$col.lobs,size = plot.opt$lwd.lobs)+
    # Censored events
    {if(length(idx.cens)>0 & plot.opt$plot.censTTE)
      geom_point(data=tab.obs[idx.cens,],aes(x=tobs, y=km), colour=plot.opt$col.pcens)} +
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


###########################	VPC for Categorical data		#############################


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

