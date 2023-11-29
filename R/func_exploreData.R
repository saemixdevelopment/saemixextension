# fix "no visible binding for global variable '.x'" NOTE (from scales package)
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".x"))

###########################	Visualising non Gaussian data		#############################

#' Plot non Gaussian data
#' 
#' This function provides exploration plots for non Gaussian longitudinal data (work in progress, doesn't work yet for RTTE)
#' 
#' @param object an SaemixData object returned by the \code{\link{saemixData}} function. 
#' For plotDiscreteDataElement, an SaemixObject object returned by the \code{\link{saemix}} function
#' @param outcome type of outcome (valid types are "TTE", "binary", "categorical", "count")
#' @param verbose whether to print messages (defaults to FALSE)
#' @param \dots additional arguments, used to pass graphical options (to be implemented, currently not available)
#' 
#' @details This function is a very rough first attempt at automatically creating plots to explore
#' discrete longitudinal data.
#' - for TTE data, a KM plot will be produced
#' - for count, categorical and binary data, a plot showing the proportion of each score/category across time will be shown
#' These plots can be stratified over a covariate in the data set (currently only categorical covariates)
#' by passing an argument which.cov='name' to the call
#' #' 
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
#' # Time-to-event data
#' data(lung.saemix)
#' 
#' saemix.data<-saemixData(name.data=lung.saemix,header=TRUE,name.group=c("id"),
#' name.predictors=c("time","status","cens"),name.response=c("status"),
#' name.covariates=c("age", "sex", "ph.ecog", "ph.karno", "pat.karno", "wt.loss","meal.cal"),
#' units=list(x="days",y="",covariates=c("yr","","-","%","%","cal","pounds")))
#' 
#' # Plots a KM survival plot
#' plotDiscreteData(saemix.data, outcome="TTE")
#' # Plots a KM survival plot, stratified by sex
#' plotDiscreteData(saemix.data, outcome="TTE", which.cov="sex")
#' 
#' # Count data
#' data(rapi.saemix)
#' saemix.data<-saemixData(name.data=rapi.saemix, name.group=c("id"),
#'                  name.predictors=c("time","rapi"),name.response=c("rapi"),
#'                  name.covariates=c("gender"),units=list(x="months",y="",covariates=c("")))
#' 
#' # Plots a histogram of the counts
#' plotDiscreteData(saemix.data, outcome="count")
#' 
#' @importFrom utils modifyList 
#' 
#' @aliases exploreDataTTE exploreDataCat exploreDataCountHist plotDiscreteData.aux
#' 
#' @export 

plotDiscreteData <- function(object, outcome="continuous", verbose=FALSE, ...) {
  # object
  ## an SaemixData object
  # outcome: type of outcome (valid types are TTE, binary, categorical, count)
  # verbose: whether to print messages
  # ...: to pass additional plot options (such as title, colours, etc...) which will supersede the options in the prefs slot of object
  if(!is(object,"SaemixData")) {
    message("Please provide a valid SaemixData object (created by a call to saemixData()) \n")
    return()
  }
  outcome<-tolower(outcome)
  if(outcome=="event") outcome<-"tte"
  if(is.na(match(outcome,c("continuous","tte","count","categorical","binary")))) {
    if(verbose) message("Please specify a valid outcome type (continuous (default), tte, count, categorical, binary).")
    return()
  }
  userPlotOptions  = list(...)
  if(!is.na(pmatch(outcome,"continuous"))) # TODO: switch to ggplot
    xplot<-plot(object)
  if(outcome=="tte") {
    xplot <- exploreDataTTE(object, verbose=verbose, ...)
  }
  if(outcome=="count") {
    # different options; if hist=TRUE, plot a histogram
#    xplot <- exploreDataCountHist(object, verbose=verbose, ...)
    # if trend=TRUE, plot trend of counts with time TODO
    xplot <- exploreDataCat(object, verbose=verbose, ...)
  }
  if(outcome %in% c("binary","categorical")) {
    if(outcome=="binary") {
      max.cat<-2
    } else {
      max.cat<-length(unique(object@data[,object@name.response]))
    }
    xplot <- exploreDataCat(object, max.cat=max.cat, verbose=verbose, ...)
  }
  
  return(xplot)
}

#' @rdname plotDiscreteData
#' 
#' @param mirror if TRUE, plots a mirror plot of the same type as the data (the object must include simulated data)
#' @param irep number of the replication to use in the mirror plot

plotDiscreteDataElement <- function(object, outcome="categorical", mirror=FALSE, irep=1, verbose=FALSE, ...) {
  # object
  ## an object resulting from a call to saemix
  ## if mirror=FALSE, plot the data 
  ## if mirror=TRUE, plot mirror plots (must have a sim.data object containing data simulated under the model)
  # outcome: type of outcome (valid types are TTE, binary, categorical, count)
  # irep: for mirror plots, replication to plot
  # verbose: whether to print messages
  # ...: to pass additional plot options (such as title, colours, etc...) which will supersede the options in the prefs slot of object
  if(!is(object,"SaemixObject")) {
    message("Please provide a valid SaemixObject object (created by a call to saemix()) \n")
    return()
  }
  objectData <- object@data
  outcome<-tolower(outcome)
  if(outcome=="event") outcome<-"tte"
  if(mirror) {
    if(object@sim.data@nsim==0) {
      message("Please simulate data before trying to plot it \n")
      return()
    }
    irep<-try(trunc(irep))
    if(irep>object@sim.data@nsim | irep<1) irep<-1
    if(outcome %in% c("tte","rtte")) {
     if(outcome=="tte") 
       objectData@data[objectData@name.X]<-object@sim.data@datasim$ysim[object@sim.data@datasim$irep==irep] else {# RTTE... more complicated
         origdata <- objectData@data
         covdata <- origdata[!duplicated(origdata$index), c(objectData@name.group, objectData@name.cens, objectData@name.occ,objectData@name.covariates)]
         yevent <- object@sim.data@datasim$ysim[object@sim.data@datasim$irep==irep]
         yresp <- as.integer(yevent>0) # 1=event, 0=initial
         idx<-which(yevent==0) # time=0
         idx1<-c(idx[-c(1)],length(yevent)+1)
         nind.obs<-c(idx1[1]-1,diff(idx))
         idx<-c(idx-1,length(yresp)) # last time before time=0, remove -1
         idx<-idx[idx>0]
         yresp[idx]<-0 # last event of each subject is censored
         simdata <- data.frame(index=rep(1:objectData@N, times=nind.obs), time=yevent, status=yresp)
         colnames(simdata)[2:3]<-c(objectData@name.X, objectData@name.response)
         simdata<-cbind(simdata, covdata[rep(1:objectData@N, times=nind.obs),])
         simdata[,objectData@name.cens]<-yresp
         objectData@data <- simdata # Replicating covariate lines for the original data and replacing time and status/cens columns by the censoring indicator (0=no event, 1=event)
         objectData@nind.obs <- nind.obs
         objectData@ntot.obs<- dim(simdata)[1]
     } 
    } else 
        objectData@data[objectData@name.response] <- object@sim.data@datasim$ysim[object@sim.data@datasim$irep==irep]
  }
  if(is.na(match(outcome,c("continuous","tte","count","categorical","binary","rtte")))) {
#    if(verbose) message("Please specify a valid outcome type (continuous (default), tte, count, categorical, binary). RTTE are currently not supported.")
    if(verbose) message("Please specify a valid outcome type (continuous (default), tte, rtte, count, categorical, binary).")
    return()
  }
  userPlotOptions  = list(...)
  if(!is.na(pmatch(outcome,"continuous"))) # TODO: switch to ggplot
    xplot<-plot(objectData)
  if(outcome %in% c("tte","rtte")) {
    xplot <- exploreDataTTE(objectData, verbose=verbose, ...)
  }
  if(outcome=="count") {
    # different options; if hist=TRUE, plot a histogram
    #    xplot <- exploreDataCountHist(object, verbose=verbose, ...)
    # if trend=TRUE, plot trend of counts with time TODO
    xplot <- exploreDataCat(objectData, verbose=verbose, ...)
  }
  if(outcome %in% c("binary","categorical")) {
    if(outcome=="binary") {
      max.cat<-2
    } else {
      max.cat<-length(unique(objectData@data[,objectData@name.response]))
    }
    xplot <- exploreDataCat(objectData, max.cat=max.cat, verbose=verbose, ...)
  }
  return(xplot)
}


###########################	KM plot for TTE data		#############################
# TODO: RTTE (more than 1 event, do a VPC for first, second,... until maxevent)
# TODO: passing plot options + extending current graphical options to match those of npde (through npdeControl)

exploreDataTTE <- function(object, verbose=FALSE, ...) {
  # object: SaemixData object
  # verbose: whether to print messages
  # ...: to pass additional plot options (such as title, colours, etc...)

  # Creates:
  # obsdat: dataframe with the observation data
  ## id: subject id
  ## time: event time
  ## event: 0 for no event, 1 for event
  ## cens: 0 for observed event/no event, 1 for censored (optional, if not present all 1 for events are assumed to be observed events and a 0 at the last observation time indicates a censored event)
  
  # Plot options - TODO
  userPlotOptions  = list(...)
  plot.opt<-saemix.data.setoptions(object)
  plot.opt <- modifyList( plot.opt, userPlotOptions[ intersect( names( userPlotOptions ), names( plot.opt ) ) ] )
  plot.opt$bands <- FALSE # plot PI (from KM)
  plot.opt$alpha.bands <-0.3
  plot.opt$col.pcens <- "steelblue3"
  plot.opt$breaks.x <- plot.opt$breaks.y <- 10
  plot.opt$plot.censTTE <- TRUE
  plot.opt$lty.lobs<-plot.opt$lty
  plot.opt$col.lobs<-plot.opt$col
  plot.opt$lwd.lobs<-plot.opt$lwd
  plot.opt$ylab <- "Survival (%)"
  
  obsdat <- data.frame(id=object@data[,object@name.group], time=object@data[,object@name.X], event=object@data[,object@name.response])
  if(length(object@name.cens)>0)
    obsdat$cens <- object@data[,object@name.cens]
  if(length(plot.opt$which.cov)!=1 || length(grep(plot.opt$which.cov, object@name.covariates))==0)
    obsdat$covariate.group <- 'all' else obsdat$covariate.group <- object@data[,plot.opt$which.cov]
  cgroups<-unique(obsdat$covariate.group)
  
  # Observed KM for each covariate group
  tab.obs1<-tab.cens1<-NULL
  has.cens<-FALSE
  for(igroup in cgroups) {
    obsdat1<-obsdat[obsdat$covariate.group==igroup,]
    event.obs<-obsdat1[obsdat1$time>0,]
    t1<-table(event.obs$time) # nb of subjects dropping off at different times (event or censoring)
    tab.obs <- data.frame(tobs=as.numeric(names(t1)),nk=c(t1))
    if(match("cens", colnames(obsdat1)) && sum(obsdat1$cens)>0) { # cens indicates censoring (cens=1) versus event (cens=0)
      t2<-table(event.obs$time[event.obs$cens==0]) # nb of subjects with an event 
      tab.obs$dk<-0
      idx<-match(names(t2), names(t1))
      tab.obs$dk[idx]<-t2
      idx.cens <- which(tab.obs$nk != tab.obs$dk)
      has.cens<-TRUE
    } else { # all events are observed (no cens or all cens=0)
      tab.obs$dk<-tab.obs$nk
      idx.cens<-c()
    }
    nk.obs <- c(sum(tab.obs$nk), sum(tab.obs$nk)-cumsum(tab.obs$nk))
    tab.obs$nk <- nk.obs[1:(length(nk.obs)-1)]
    tab.obs$sk <- (1-tab.obs$dk/tab.obs$nk)
    tab.obs$km <- cumprod(tab.obs$sk)
    if(length(idx.cens)>0) {
      tab.cens<-tab.obs[idx.cens,]
      tab.cens1<-rbind(tab.cens1, cbind(covariate.group=igroup,tab.cens))
    }
    tab.obs1<-rbind(tab.obs1, cbind(covariate.group=igroup,tab.obs))
  }
  tab.obs<-tab.obs1

  # Creating plot 
  ## TODO: add options
  ## TODO: add censored data (eg red dots to indicate censored events in the KM plot)
  xplot <- ggplot(data=tab.obs, aes(x=.data$tobs, y=.data$km, group=as.factor(.data$covariate.group), colour=as.factor(.data$covariate.group))) + 
    geom_line(linetype = plot.opt$lty.lobs, size = plot.opt$lwd.lobs) +
    # Censored events
    {if(has.cens & plot.opt$plot.censTTE)
      geom_point(data=tab.cens1, aes(x=.data$tobs, y=.data$km), colour=plot.opt$col.pcens)} +
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
    {if(length(cgroups)==1) theme(legend.position="none")} +
    {if(length(cgroups)>1) guides(fill=guide_legend(title=plot.opt$which.cov), colour=guide_legend(title=plot.opt$which.cov))}
  return(xplot)
}

###########################	Histogram for count data		#############################
# Current version
# TODO: time trend, using ggstream to produce streamgraphs (with options to stack them and sum them to 1)

exploreDataCountHist <- function(object, verbose=FALSE, ...) {
  # object: SaemixData object
  # verbose: whether to print messages
  # ...: to pass additional plot options (such as title, colours, etc...)
  
  # Creates:
  # obsdat: dataframe with the observation data
  ## id: subject id
  ## x: independent variable used for the plot (eg time)
  ## y: count at time x
  ## covariates in the model if present
  obsdat <- data.frame(id=object@data[,object@name.group], x=object@data[,object@name.X], y=object@data[,object@name.response])
  if(length(object@name.covariates)>0)
    obsdat<-cbind(obsdat, object@data[,object@name.covariates, drop=FALSE]) # retain column names
  
  # Plot options - TODO
  userPlotOptions  = list(...)
  plot.opt<-saemix.data.setoptions(object)
  plot.opt$bands <- FALSE # plot PI (from KM)
  plot.opt$alpha.bands <-0.3
  plot.opt$col.pcens <- "steelblue3"
  plot.opt$breaks.x <- plot.opt$breaks.y <- 10
  plot.opt$plot.censTTE <- TRUE
  plot.opt$lty.lobs<-plot.opt$lty
  plot.opt$col.lobs<-plot.opt$col
  plot.opt$lwd.lobs<-plot.opt$lwd
  plot.opt$hist.freq <- TRUE
  plot.opt <- modifyList( plot.opt, userPlotOptions[ intersect( names( userPlotOptions ), names( plot.opt ) ) ] )
  plot.opt$ylab <- ifelse(plot.opt$hist.freq, "Frequency","Density")
  if ("xlim" %in% names(plot.opt) & length(plot.opt$xlim)==2)
    x.limits = c(plot.opt$xlim[1],plot.opt$xlim[2])  else
      x.limits = c(min(obsdat$y,na.rm=TRUE),max(obsdat$y,na.rm=TRUE))
  
  # Future: stratify over eg covariates or groups ?
  numberCategories <- 1
  
  # Creating plot 
  ## TODO: add options
  xplot <- ggplot(data=obsdat, aes(x=.data$y)) + geom_histogram() +
    # x log-scale (no sense in y log-scale)
    { if (plot.opt$xlog == FALSE)  scale_x_continuous(plot.opt$xlab, limits = x.limits,  scales::pretty_breaks(n = plot.opt$breaks.x))
    } +
    { if (plot.opt$xlog == TRUE)
      scale_x_log10(plot.opt$xlab,breaks = scales::trans_breaks("log10", function(x) 10 ^ x),
                    limits = x.limits, labels = scales::trans_format("log10", scales::math_format(10 ^ .x))) 
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

###########################	Proportion of each category across X, categorical data		#############################
# TODO: passing plot options + extending current graphical options to match those of npde (through npdeControl)

exploreDataCat <- function(object, max.cat=10, breaks=NULL, catlabel=NULL, verbose=FALSE, ...) {
  # object: SaemixData object
  # verbose: whether to print messages
  # ...: to pass additional plot options (such as title, colours, etc...)
  x<-plotDiscreteData.aux(object, max.cat=max.cat, breaks=breaks, catlabel=catlabel, verbose=verbose, ...)
  xtab<-x$xtab
  plot.opt<-x$plot.opt
  if(length(plot.opt$mfrow)==2) {
    nrow1 <- plot.opt$mfrow[1]
    ncol1 <- plot.opt$mfrow[2]
  } else {
    nrow1<-ncol1<-NULL
  }
  cgroups<-unique(xtab$group)
  plot.cat <- ggplot(data=xtab, aes(x=.data$x.group, y=.data$freq, group=as.factor(.data$group), colour=as.factor(.data$group))) + geom_line() + 
    xlab(plot.opt$xlab) + ylab(plot.opt$ylab) + 
    {if(length(cgroups)==1) theme(legend.position="none")} +
    {if(length(cgroups)>1) guides(fill=guide_legend(title=plot.opt$which.cov), colour=guide_legend(title=plot.opt$which.cov))} +
    facet_wrap(.~.data$y.group, nrow=nrow1, ncol=ncol1)
  return(plot.cat)

}

###########################	Binning on X-axis and on scores		#############################

plotDiscreteData.aux <- function(object, max.cat=10, breaks=NULL, catlabel=NULL, verbose=FALSE, ...) {
  # object: a SaemixData object including simulated data
  # max.cat, breaks: used to bin the categories of the response
  # verbose: whether to print messages
  # ...: to pass additional plot options (such as title, colours, etc...)
  
  # Creates:
  # obsdat: dataframe with the observation data
  ## id: subject id
  ## time: event time
  ## event: 0 for no event, 1 for event
  ## cens: 0 for observed event/no event, 1 for censored (optional, if not present all 1 for events are assumed to be observed events and a 0 at the last observation time indicates a censored event)
  
  # Returns
  ## xtab: a dataframe with columns
  ### group: covariate/stratification group
  ### x.group: independent variable (usually time), grouped
  ### y.group: dependent variable (category, score), grouped
  ### freq: frequency of y.group at the corresponding combination of group and x.group
  ### nobs: number of observations at the corresponding combination of group and x.group
  ## updated plot.opt
  ydat <- object
  obsdat<-data.frame(id=ydat@data[,ydat@name.group], x=ydat@data[,ydat@name.X], y=ydat@data[,ydat@name.response])

  plot.opt<-saemix.data.setoptions(object)
  plot.opt$bin.method <- "equal"
  plot.opt$bin.number <- 10
  plot.opt<-replace.plot.options(plot.opt,...)
  logtyp<-""
  if(plot.opt$xlog) logtyp<-paste(logtyp,"x",sep="")
  if(plot.opt$ylog) logtyp<-paste(logtyp,"y",sep="")
  
  if(is.null(plot.opt$bin.method)) plot.opt$bin.method<- "equal"
  if(plot.opt$bin.method=="user" & is.null(plot.opt$bin.breaks)) {
    if(verbose) message("When using the bin.method='user' option, please provide the breaks in bin.breaks\n")
    plot.opt$bin.method<- "equal"
  }
  # Create covariate groups - TODO
  if(length(plot.opt$which.cov)!=1 || length(grep(plot.opt$which.cov, ydat@name.covariates))==0)
    obsdat$covariate.group <- 'all' else obsdat$covariate.group <- ydat@data[,plot.opt$which.cov]
  
  # Binning on X
  if(length(unique(obsdat$x))>plot.opt$bin.number) {
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
      xgrp <- cut(obsdat$y, breaks=breaks, include.lowest=TRUE, right=FALSE)
    } else {
      breaks<-sort(unique(breaks))
      if(min(breaks)>min(obsdat$y)) breaks<-c(min(obsdat$y), breaks)
      breaks[1]<-breaks[1]-0.01 # make sure all observations are included
      xgrp <- cut(obsdat$y, breaks=breaks, include.lowest=TRUE, right=FALSE)
    }
    obsdat$score.group<-xgrp
  }
  # Computing frequency of each score at each time in each covariate group - observed data
  # Computing corresponding PI - simulated data
  xtab<-NULL
  cgroups <- unique(obsdat$covariate.group)
  xgroups <- sort(unique(obsdat$x.group))
  sgroups <- sort(unique(obsdat$score.group))
  #  print(sgroups)
  #  print(xgroups)
  for(igroup in cgroups) {
    tab1<-table(obsdat$x.group[obsdat$covariate.group==igroup],obsdat$score.group[obsdat$covariate.group==igroup])
    ncol<-colSums(tab1)
    tab1<-tab1[,(ncol>0),drop=FALSE]
    nobs<-rowSums(tab1)
    freqtab <- tab1/nobs
    xtab <- rbind(xtab,
                  data.frame(group=igroup, x.group=rep(xgroups,length(sgroups)), y.group=rep(sgroups, each=length(xgroups)), freq=c(freqtab), nobs=rep(nobs, length(sgroups))))
  }
  # More pleasing labels
  if(!is.null(catlabel) && length(catlabel)==length(levels(xtab$y.group)))
    levels(xtab$y.group)<-catlabel else {# try to define nicer names...
      if(max.cat==length(unique(xtab$y.group))) {
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
      levels(xtab$y.group)<-catlabel1
    }
  return(list(xtab=xtab, plot.opt=plot.opt))
}

######################	Binning a vector of values according to different methods ################################

#' Internal functions used to produce prediction intervals (from the npde package)
#'
#' Functions used by plot functions to define the boundaries of the bins on the X-axis
#'
#' These functions are normally not called by the end-user.
#'
#' @param xvec a vector of values for the X-axis of the plot
#' @param plot.opt graphical options
#' @param verbose boolean (defaults to FALSE). If TRUE, a table showing how the binning was performed
#'
#' @return a list with 3 elements, xgrp (the number of the bin associated with each element of xvec), xcent (a named vector containing the mean of the elements of xvec contained in each bin; the name of the bin is the interval), and xgroup (a vector with the group associated to each element of xvec after binning)
#' If verbose is TRUE, a table showing the bins is shown, giving the interval of xvec associated with each bin, the mean value
#' of xvec in each bin, and the number of observations
#'
#' @name xbinning
#' 
## #' @aliases
#' @author Emmanuelle Comets <emmanuelle.comets@@bichat.inserm.fr>
#' @seealso \code{\link{npde}}, \code{\link{autonpde}}
#' @references K. Brendel, E. Comets, C. Laffont, C. Laveille, and F.
#' Mentre. Metrics for external model evaluation with an application to the
#' population pharmacokinetics of gliclazide. \emph{Pharmaceutical Research},
#' 23:2036--49, 2006.
#'
#' @importFrom mclust Mclust
#' @export

# Binning the X data
xbinning<-function(xvec, plot.opt, verbose=FALSE) {
  # plot.opt should contain the following items
  # bin.method: binning method (one of "optimal","width","user","equal")
  # bin.number: number of bins (default: 10, maximum: 20)
  # bin.breaks
  # bin.extreme
  # xlog
  
  nbin<-plot.opt$bin.number
  bnds <- plot.opt$bin.breaks
  
  xvec1<-xvec
  xvec<-xvec[!is.na(xvec)]
  if(is.na(pmatch(plot.opt$bin.method,c("optimal","width","user","equal")))) {
    if(verbose) cat("Binning method",plot.opt$bin.method,"not found, reverting to equal binning\n")
    plot.opt$bin.method<-"equal"
  }
  if(!is.na(pmatch(plot.opt$bin.method,"optimal"))) {
    if(!("mclust"%in%.packages(all.available = TRUE))) {
      if(verbose) cat("mclust library not installed, reverting to equal binning\n")
      plot.opt$bin.method<-"equal"
    } else {
      #			require(mclust)
      if(is.null(plot.opt$bin.number) || is.na(plot.opt$bin.number)) plot.opt$bin.number<-10
    }
  }
  if(!is.na(pmatch(plot.opt$bin.method,"user")) & is.null(plot.opt$bin.breaks)) {
    if(verbose) cat("User-defined method specified, but bin.breaks is empty; reverting to equal binning\n")
    plot.opt$bin.method<-"equal"
  }
  if(!is.na(pmatch(plot.opt$bin.method,c("equal","width"))) & is.null(plot.opt$bin.number)) {
    nbin<-length(unique(xvec))
    if(nbin>20) nbin<-20
    plot.opt$bin.number<-nbin
  }
  if(is.na(pmatch(plot.opt$bin.method,c("optimal","user"))) && length(unique(xvec))<=nbin) {
    xgrp<-match(xvec,sort(unique(xvec)))
    xpl<-tapply(xvec,xgrp,mean)
  } else {
    if(!is.na(pmatch(plot.opt$bin.method,"user"))) {
      bnds<-plot.opt$bin.breaks
      if(min(bnds)>=min(xvec)) bnds<-c(min(xvec)*(1-sign(min(xvec))*0.001),bnds)
      if(max(bnds)<max(xvec)) bnds<-c(bnds,max(xvec))
    }
    if(!is.na(pmatch(plot.opt$bin.method,"equal"))) {
      xvec2<-xvec;xvec2[xvec2==min(xvec)]<-min(xvec)-1
      if(!is.null(plot.opt$bin.extreme) & length(plot.opt$bin.extreme)==2) {
        xq<-plot.opt$bin.extreme
        xquant<-c(0,seq(xq[1],xq[2],length.out=(nbin-1)),1)
      } else xquant<-(0:nbin)/nbin
      bnds<-unique(quantile(xvec2,xquant,type=8))
    }
    if(!is.na(pmatch(plot.opt$bin.method,"width"))) {
      if(plot.opt$xlog) xvec2<-log(xvec[xvec>0]) else xvec2<-xvec
      if(!is.null(plot.opt$bin.extreme) & length(plot.opt$bin.extreme)==2) {
        xq<-plot.opt$bin.extreme
        xq1<-quantile(xvec2,xq,type=8)
        bnds<-c(min(xvec2),seq(xq1[1],xq1[2],length.out=(nbin-1)),max(xvec2))
        bnds<-sort(unique(bnds))
      } else bnds<-seq(min(xvec2),max(xvec2),length.out=(nbin+1))
      if(plot.opt$xlog) {
        bnds<-exp(bnds)
        bnds[length(bnds)]<-bnds[length(bnds)]*(1+sign(bnds[length(bnds)])*0.001)
        if(sum(xvec<=0)>0) bnds<-c(min(xvec),bnds)
      }
      bnds[1]<-bnds[1]*(1-sign(bnds[1])*0.001)
    }
    if(!is.na(pmatch(plot.opt$bin.method,"optimal"))) {
      yfit<-mclust::Mclust(xvec,G=((nbin-5):(nbin+5)))
      xgrp<-yfit$classification
      xpl<-yfit$parameters$mean
      xpl<-xpl[match(names(table(xgrp)),names(xpl))]
      minx<-tapply(xvec,xgrp,min)
      maxx<-tapply(xvec,xgrp,max)
      bnds <- c(minx[1],(minx[-c(1)]+maxx[-length(maxx)])/2,maxx[length(maxx)])
      names(xpl)<-paste("[",tapply(xvec,xgrp,min),"-",tapply(xvec,xgrp,max),"]", sep="")
    } else {
      xgrp<-factor(cut(xvec,bnds,include.lowest=F))
      xpl<-tapply(xvec,xgrp,mean)
    }
  }
  
  nbin<-length(unique(xgrp))
  npl<-tapply(xvec,xgrp,length)
  tab<-cbind(Interval=names(xpl),Centered.On=format(xpl,digits=2),Nb.obs=npl)
  row.names(tab)<-1:dim(tab)[1]
  
  if(verbose) {
    xnam<-switch(EXPR=plot.opt$bin.method,equal="by quantiles on X", width="equal sized intervals",user="user-defined bins",optimal="clustering algorithm")
    cat("Method used for binning:",xnam,", dividing into the following",nbin,"intervals\n")
    print(tab,quote=F)
  }
  
  xgrp2<-rep(NA,length(xvec1))
  xgrp2[!is.na(xvec1)]<-xgrp
  
  # ------------------------------------------------------------------------------------------------
  # Take the min and max values of xpl to extend prediction band to min/max X
  # start to 0 instead of the minimal value of bnds (can give negative values for Time in x-axis)
  # is.null test for non numerical covariate
  if (!is.null(plot.opt$vpc.extend) && plot.opt$vpc.extend) { # extend vpc bands to min/max
    xpl.extended = xpl
    xpl.extended[1] = min(xvec)
    xpl.extended[nbin] = max(xvec)
    xpl = xpl.extended
  }
  
  return(list(xgrp=xgrp2,xcent=xpl,xbound=bnds))
}
