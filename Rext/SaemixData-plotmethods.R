#' @include aaa_generics.R
#' @include SaemixData.R
#' @include SaemixData-methods.R
#' @include SaemixData-printmethods.R
NULL

####################################################################################
####			SaemixData class - method to plot			####
####################################################################################
### #' @export
### #' @docType methods
### #' @rdname plot-methods
### #' @aliases plot,SaemixData 


saemix.data.setoptions<-function(saemix.data) {
  # setting default plot options
  plot.opt<-list(
    # General graphical options
    new=TRUE,				# whether a new page should be called
    ask=FALSE,				# whether the program should ask before creating a new page
    ilist=c(1:saemix.data["N"]),
    separate=FALSE,	# if TRUE, plots individual subjects (a la nlme), if FALSE plots a single plot with all the subjects
    # Options for individual plots
    nmax=12,					  # maximum number of subjects
    limit=TRUE,					# limit to nmax plots
    sample=FALSE,				# if FALSE=use the (nmax) first subjects; TRUE=randomly sample (nmax) subjects from the dataset
    interactive=FALSE,  # whether the program should prompt the user for the number of subjects to plot in the individual plots if this number exceeds nmax
    which.cov="none",  # whether to split over covariates
    # Layout and plots options
    mfrow=c(),				# page layout (if empty, defaults to the default layout for each graph type)
    main=" ",				# title
    xlab=" ",
    ylab=" ",
    units=saemix.data@units,
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
    cex.main=1)
  
  if(is.null(saemix.data@name.X))
    plot.opt$name.X<-saemix.data["name.predictors"][1] else plot.opt$name.X<-saemix.data@name.X
    plot.opt$xlab<-paste(plot.opt$name.X," (",ifelse(plot.opt$units$x=="","-",plot.opt$units$x),")", sep="")
    if(length(saemix.data["name.response"])>0)
      plot.opt$ylab<-paste(saemix.data["name.response"]," (",ifelse(plot.opt$units$y=="","-",plot.opt$units$y),")", sep="")
    return(plot.opt)
}

replace.data.options<-function(plot.opt,...) {
  args1<-match.call(expand.dots=TRUE)
  # These arguments are used by other functions and may be passed on via "...", so we want to ignore them. Other arguments not in list will raise warnings
  legacy<-c("plot.type","individual")
  if(length(args1)>2) {
    # Other arguments
    for(i in 3:length(args1)) {
      if(match(names(args1)[i],names(plot.opt),nomatch=0)>0) {
        #    plot.opt[[names(args1)[i]]]<-args1[[i]] else {
        if(!is.null(eval(args1[[i]]))) plot.opt[[names(args1)[i]]]<-eval(args1[[i]])
      } else {
        if(is.na(match(names(args1)[i],legacy))) message(paste("Argument",names(args1)[i],"not available, check spelling"))
      }
    }
  }
  return(plot.opt)
}

#' Plot of longitudinal data 
#' 
#' This function will plot a longitudinal dataframe contained in an SaemixData object. By default it produces a spaghetti plot, but arguments can be passed on to modify this behaviour. 
#' 
## #' @name plot-SaemixData
#' 
#' @param x an SaemixData object or an SaemixSimData object
#' @param y unused, present for compatibility with base plot function
#' @param ... additional arguments to be passed on to plot (titles, legends, ...)
#' 
#' @aliases plot,SaemixData-methods 
#' @aliases plot-SaemixData
#' @aliases plot,SaemixData plot,SaemixData,ANY-method
#' @keywords plot
### #' @docType methods
#' @rdname plot-SaemixData
#' 
#' @import ggplot2 grid gridExtra
#' @importFrom graphics plot
#' @importFrom rlang is_missing
#' @method plot SaemixData
#' @export 


# Plot the data, either as points or as lines grouped by x@name.group
plot.SaemixData<-function(x,y,...) {
  if(length(x@data)==0) {
    message("No data to plot.\n")
    return("Missing data")
  }
  # Eco: commented, otherwise was resetting the graphical layout on exit and preventing the graphs to be set on the same page
  #    oldpar <- par(no.readonly = TRUE)    # code line i
  #    on.exit(par(oldpar))            # code line i + 1
  
  # User-defined options
  userPlotOptions<-list(...)
  if(!is_missing(y) && is.list(y)) {
    userPlotOptions<-c(y,userPlotOptions)
  }
  i1<-match("individual",names(userPlotOptions))
  if(!is.na(i1)) {
    individual<-as.logical(eval(userPlotOptions[[i1]]))
  } else individual<-FALSE
  i1<-match("type",names(userPlotOptions))
  if(!is.na(i1)) {
    plot.type<-as.character(userPlotOptions[[i1]])
    plot.type<-plot.type[plot.type!="c"]
  } else plot.type<-c()
  if(length(plot.type)==0) plot.type<-ifelse(individual,"b","l")
  
  # Default options for data plot
  plot.opt <- saemix.data.setoptions(x)
  plot.opt$xlab<-paste(x@name.X," (",x@units$x,")",sep="")
  plot.opt$ylab<-paste(x@name.response," (",x@units$y,")",sep="")
  plot.opt$type<-ifelse(individual,"b","l")
  
  # Replace default options by options passed explicitly
  if(length(userPlotOptions)>0)
    plot.opt <- modifyList(plot.opt, userPlotOptions[intersect(names(userPlotOptions), names(plot.opt))])
  
  #plot.opt<-replace.data.options(plot.opt,...)
  logtyp<-paste(ifelse(plot.opt$xlog,"x",""),ifelse(plot.opt$ylog,"y",""),sep="")
  if(individual) { # separate plots subject per subject
    if(length(plot.opt$ilist)>plot.opt$nmax & plot.opt$limit) {
      if(plot.opt$interactive) {
        x1<-readline(prompt=paste("The number of subjects may be too large to be plotted. Should I plot only",plot.opt$nmax,"subjects ? (Y/n) \n"))
        if(tolower(x1)=="y") {
          plot.opt$limit<-TRUE
          plot.opt$ilist<-plot.opt$ilist[1:plot.opt$nmax]
          if(plot.opt$sample) plot.opt$ilist<-sort(sample(plot.opt$ilist, plot.opt$nmax)) else plot.opt$ilist<-plot.opt$ilist[1:plot.opt$nmax]
          if(!plot.opt$ask) {
            x1<-readline(prompt="Stop after each page of plot ? (Y/n) \n")
            if(tolower(x1)=="y") plot.opt$ask<-TRUE
          }
        }
      } else {
        if(plot.opt$interactive) {
          cat("The number of subjects is too large, I will plot only")
          if(plot.opt$sample) cat(" the data for",plot.opt$nmax,"subjects sampled randomly;") else cat(" only the data for the first",plot.opt$nmax,"subjects;")
          cat(" use limit=FALSE in the call to plot to force plotting all the subjects.\n")
        }
        if(plot.opt$sample) plot.opt$ilist<-sort(sample(plot.opt$ilist, plot.opt$nmax)) else plot.opt$ilist<-plot.opt$ilist[1:plot.opt$nmax]
      }
    } # end of test on length(ilist)
    if(plot.opt$new) {
      if(length(plot.opt$mfrow)==0) {
        np<-length(plot.opt$ilist)
        if(np>12) np<-12
        n1<-round(sqrt(np))
        n2<-ceiling(np/n1)
        par(mfrow=c(n1,n2),ask=plot.opt$ask)
      } else par(mfrow=plot.opt$mfrow,ask=plot.opt$ask)
    }
    xind<-x["data"][,x["name.predictors"], drop=FALSE]
    id<-x["data"][,"index"]
    yobs<-x["data"][,x["name.response"]]
    for(isuj in plot.opt$ilist) {
      if(plot.opt$main=="") main<-paste("Subject",isuj) else main<-plot.opt$main
      plot(xind[id==isuj,x@name.X],yobs[id==isuj],type=plot.type, xlab=plot.opt$xlab,ylab=plot.opt$ylab,col=plot.opt$col,pch=plot.opt$pch,log=logtyp, xlim=plot.opt$xlim,ylim=plot.opt$ylim,main=main,cex=plot.opt$cex, cex.axis=plot.opt$cex.axis,cex.lab=plot.opt$cex.lab,lty=plot.opt$lty, lwd=plot.opt$lwd)
    }
  } else {	# One plot for all the data
    if(plot.opt$new) par(mfrow=c(1,1))
    if(plot.type=="p" | plot.type=="b") {
      plot(x@data[,x@name.X],x@data[,x@name.response],xlab=plot.opt$xlab, ylab=plot.opt$ylab,col=plot.opt$col,pch=plot.opt$pch,log=logtyp,xlim=plot.opt$xlim, ylim=plot.opt$ylim,main=plot.opt$main,cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab) }
    if(plot.type=="l") {
      plot(x@data[,x@name.X],x@data[,x@name.response],xlab=plot.opt$xlab, ylab=plot.opt$ylab,col=plot.opt$col,lty=plot.opt$lty,lwd=plot.opt$lwd,type="n", log=logtyp,xlim=plot.opt$xlim,ylim=plot.opt$ylim,main=plot.opt$main, cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab)
    }
    if(plot.type=="l" | plot.type=="b") {
      for(isuj in unique(x@data[,x@name.group])) {
        lines(x@data[x@data[,x@name.group]==isuj,x@name.X], x@data[x@data[,x@name.group]==isuj,x@name.response],col=plot.opt$col, lty=plot.opt$lty,lwd=plot.opt$lwd)
      }
    }
  }
}


#' When applied to an SaemixSimData object, mirror plots are produced which help assess whether the simulated data has similar features when compared to the original data.
#'
#' @name plot-SaemixData
#' 
#' @param irep which replicate datasets to use in the mirror plot (defaults to -1, causing a random simulated dataset to be sampled from the nsim
#' simulated datasets)
#' @param prediction if TRUE, plot the predictions without residual variability (ypred instead of ysim). Defaults to FALSE.
#' 
#' @aliases plot,SaemixSimData-method plot,SaemixSimData plot,SaemixSimData,ANY-method
#' 
#' @details this function can also be used to visualise the predictions for simulated values of the individual parameters,
#' using the ypred element instead of the ysim element normally used here
#' 
### #' @docType methods
#' @rdname plot-SaemixData
#' @importFrom graphics plot
#' @method plot SaemixSimData
#' @export 

# Check simulations using mirror plots
plot.SaemixSimData<-function(x, y, irep=-1, prediction=FALSE, ...) {
  #    oldpar <- par(no.readonly = TRUE)    # code line i
  #    on.exit(par(oldpar))            # code line i + 1
  # User-defined options
  userPlotOptions<-list(...)
  if(!is_missing(y) && is.list(y)) {
    userPlotOptions<-c(y,userPlotOptions)
  }
  i1<-match("interactive",names(userPlotOptions))
  if(!is.na(i1)) interactive<-as.logical(userPlotOptions[[i1]]) else interactive<-FALSE
  
  if(prediction) yplot <- x@datasim$ypred else yplot<-x@datasim$ysim
  if(length(yplot)==0) {
    message("No simulated data to plot\n")
    return()
  }
  
  i1<-match("warnings",names(userPlotOptions))
  if(!is.na(i1)) printwarnings<-as.logical(userPlotOptions[[i1]]) else printwarnings<-FALSE
  if(dim(x@datasim)[1]==0) {
    if(interactive | printwarnings) message("No simulated data.\n")} else {  
      i1<-match("type",names(userPlotOptions))
      if(!is.na(i1)) {
        plot.type<-as.character(userPlotOptions[[i1]])
        plot.type<-plot.type[plot.type!="c"]
      } else plot.type<-"l"
      # Default options for data plot
      plot.opt <- saemix.data.setoptions(x)
      plot.opt$new<-FALSE
      plot.opt$plot.type<-"b"
      plot.opt$xlab<-paste(x@name.X," (",x@units$x,")",sep="")
      plot.opt$ylab<-paste(x@name.response," (",x@units$y,")",sep="")
      # Replace default options by options passed explicitly
      if(length(userPlotOptions)>0)
        plot.opt <- modifyList(plot.opt, userPlotOptions[intersect(names(userPlotOptions), names(plot.opt))])
      logtyp<-paste(ifelse(plot.opt$xlog,"x",""),ifelse(plot.opt$ylog,"y",""),sep="")
      
      if(irep[1]<0) irep<-sample(unique(x@nsim),1)
      for(irep1 in irep) {
        if(plot.opt$main==" ") tit<-paste("Mirror plot (replication ",irep1,")",sep="") else tit<-plot.opt$main
        tab<-data.frame(id=x@data[,x@name.group],x=x@data[,x@name.X], y=yplot[x@datasim$irep==irep1])
        if(plot.type=="p" | plot.type=="b") {
          plot(tab[,"x"],tab[,"y"],xlab=plot.opt$xlab, ylab=plot.opt$ylab, col=plot.opt$col,pch=plot.opt$pch,log=logtyp,xlim=plot.opt$xlim, ylim=plot.opt$ylim,main=tit,cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab) }
        if(plot.type=="l") {
          plot(tab[,"x"],tab[,"y"],type="n",xlab=plot.opt$xlab, ylab=plot.opt$ylab,col=plot.opt$col,lty=plot.opt$lty,lwd=plot.opt$lwd, log=logtyp,xlim=plot.opt$xlim,ylim=plot.opt$ylim,main=tit, cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab)
        }
        if(plot.type=="l" | plot.type=="b") {
          for(isuj in unique(tab[,"id"])) {
            lines(tab[tab[,"id"]==isuj,"x"],tab[tab[,"id"]==isuj,"y"])
          }
        }
        
      }
    }
}
