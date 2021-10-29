# Setting some options by default and intersecting with options set by user - howto
library(rlang)

afunction<-function(x, y, ...) {
  #  userPlotOptions<-match.call(expand.dots=TRUE)
  userPlotOptions<-list(...)
  if(!is_missing(y) && is.list(y)) {
    userPlotOptions<-c(y,userPlotOptions)
  }
  print(x)
  if(!is_missing(y)) print(y) else message("y is missing")
  cat("User options:\n")
  print(userPlotOptions)
}

afunction(x=c(3,4))

afunction(x=c(3,4),list(x=6, col="blue", y=1:6))

afunction(x=c(3,4),list(x=6, col="blue", y=1:6))

afunction(x=c(3,4),pch=6, col="blue", isuj=1:6)




# Fiddling with plots
# for some reason calling the plot() functions resets par(mfrow) even though plot.opt$new=FALSE, which is a pain

saemix.plot.mirrorplot<-function(saemixObject, nplots=2,...) {
  # Plot of the data and nplots associated mirror plots (2 mirror plots by default)
  # Creates a new page
  # options: change data point, line type, line color, lines plotted or not, points plotted or not...
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1
  
  plot.opt<-saemixObject["prefs"]
  plot.opt$new<-FALSE
  plot.opt$plot.type<-"l"
  plot.opt<-replace.plot.options(plot.opt,...)
  # Tried to use the do.call approach but it doesn't recognise the plot function to apply
  #  list.args<-plot.opt
  
  if(saemixObject@sim.data@nsim>0) nplots<-min(saemixObject@sim.data@nsim, nplots)
  np<-nplots+1
  if(np>12) np<-12
  n1<-round(sqrt(np))
  n2<-ceiling(np/n1)
  par(mfrow=c(n1,n2),ask=plot.opt$ask)
  # data
  plot(saemixObject["data"],plot.opt)
  #  list.args$saemixObject<-saemixObject@data
  #  do.call(plot, list.args)
  # simulated datasets
  if(saemixObject@sim.data@nsim>0) {
    irep1<-sort(sample(1:saemixObject@sim.data@nsim, nplots, replace=FALSE))
  } else {
    saemixObject<-simulate(saemixObject, nsim=nplots)
    irep1<-1:nplots
  }
  list.args$saemixObject<-saemixObject@sim.data
  for(irep in irep1)
    plot(saemixObject["sim.data"],irep=irep)
  #    do.call(plot, list.args)
}

############################ plot functions as S4 methods
# Plot the data, either as points or as lines grouped by x@name.group
setMethod("plot","SaemixData",
          function(x,y,...) {
            if(length(x@data)==0) {
              message("No data to plot.\n")
              return("Missing data")
            }
            # Eco: commented, otherwise was resetting the graphical layout on exit and preventing the graphs to be set on the same page
            #    oldpar <- par(no.readonly = TRUE)    # code line i
            #    on.exit(par(oldpar))            # code line i + 1
            
            args1<-match.call(expand.dots=TRUE)
            i1<-match("individual",names(args1))
            if(!is.na(i1)) {
              individual<-as.logical(eval(args1[[i1]]))
            } else individual<-FALSE
            i1<-match("type",names(args1))
            if(!is.na(i1)) {
              plot.type<-as.character(args1[[i1]])
              plot.type<-plot.type[plot.type!="c"]
            } else plot.type<-c()
            if(length(plot.type)==0) plot.type<-ifelse(individual,"b","l")
            plot.opt<-saemix.data.setoptions(x)
            mainkeep<-plot.opt$main
            plot.opt$new<-TRUE
            plot.opt$xlab<-paste(x@name.X," (",x@units$x,")",sep="")
            plot.opt$ylab<-paste(x@name.response," (",x@units$y,")",sep="")
            plot.opt$type<-ifelse(individual,"b","l")
            plot.opt<-replace.data.options(plot.opt,...)
            print(plot.opt$new)
            change.main<-FALSE
            if(plot.opt$main!=mainkeep) change.main<-TRUE
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
                if(!change.main) main<-paste("Subject",isuj) else main<-plot.opt$main
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
)


# Check for mirror plots
setMethod("plot","SaemixSimData",
          function(x,y,irep=-1,...) {
            #    oldpar <- par(no.readonly = TRUE)    # code line i
            #    on.exit(par(oldpar))            # code line i + 1
            args1<-match.call(expand.dots=TRUE)
            i1<-match("type",names(args1))
            if(!is.na(i1)) {
              plot.type<-as.character(args1[[i1]])
              plot.type<-plot.type[plot.type!="c"]
            } else plot.type<-"l"
            plot.opt<-saemix.data.setoptions(x)
            plot.opt$new<-TRUE
            plot.opt$plot.type<-"b"
            plot.opt$xlab<-paste(x@name.X," (",x@units$x,")",sep="")
            plot.opt$ylab<-paste(x@name.response," (",x@units$y,")",sep="")
            plot.opt<-replace.data.options(plot.opt,...)
            logtyp<-paste(ifelse(plot.opt$xlog,"x",""),ifelse(plot.opt$ylog,"y",""),sep="")
            if(dim(x@datasim)[1]==0) {
              if(plot.opt$interactive & plot.opt$warnings) message("No simulated data.\n")} else {
                if(irep<0) irep<-sample(unique(x@nsim),1)
                tit<-paste("Mirror plot (replication ",irep,")",sep="")
                tab<-data.frame(id=x@data[,x@name.group],x=x@data[,x@name.X], y=x@datasim$ysim[x@datasim$irep==irep])
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
)