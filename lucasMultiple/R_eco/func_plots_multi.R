
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
  ytype <- saemix.data["data"][,"ytype"]
  ydat<-saemixObject["data"]["yorig"]
  idx.exp<-which(saemixObject["model"]["error.model"]=="exponential")
  if(length(idx.exp)>0) {
    ydat[ytype %in% idx.exp]<-saemixObject["data"]["data"][ytype %in% idx.exp,saemixObject["data"]["name.response"]]
  }
  if(length(grep(0,plot.opt$level))>0) {
    if(!change.main) main<-"Population predictions" else main<-plot.opt$main
    if(plot.opt$which.poppred=="ppred") xpl<-saemixObject["results"]["ppred"] else xpl<-saemixObject["results"]["ypred"]
    if(length(xpl)==length(ydat)) {
      for(itype in unique(ytype)) {
        plot(xpl[ytype==itype],ydat[ytype==itype],xlab=paste0(plot.opt$xlab," (outcome ",itype,")"), 
             ylab=paste0(plot.opt$ylab," (outcome ",itype,")"),pch=plot.opt$pch, col=plot.opt$col,
             main=main,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
        abline(0,1,col=plot.opt$ablinecol,lty=plot.opt$ablinelty, lwd=plot.opt$ablinelwd)
      }
    }
  }
  if(length(grep(1,plot.opt$level))>0) {
    if(!change.main) main<-paste("Individual predictions", ifelse(plot.opt$indiv.par=="map","MAP","Cond mean"),sep=", ") else main<-plot.opt$main
    if(plot.opt$indiv.par=="map") xpl<-saemixObject["results"]["ipred"] else xpl<-saemixObject["results"]["icpred"]
    if(length(xpl)==length(ydat)) {
      for(itype in unique(ytype)) {
        plot(xpl[ytype==itype],ydat[ytype==itype],xlab=paste0(plot.opt$xlab," (outcome ",itype,")"), 
             ylab=paste0(plot.opt$ylab," (outcome ",itype,")"),pch=plot.opt$pch, col=plot.opt$col,
             main=main,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
        abline(0,1,col=plot.opt$ablinecol,lty=plot.opt$ablinelty, lwd=plot.opt$ablinelwd)
      }
    }
  }
}

