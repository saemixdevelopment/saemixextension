# Grouping data by time and score
yfit <- ysim.hurdle1
ydat <- yfit@data
ysim <- yfit@sim.data@datasim$ysim
nsim<-length(ysim)/dim(ydat@data)[1]
obsmat<-data.frame(id=ydat@data[,ydat@name.group], x=ydat@data[,ydat@name.X], y=ydat@data[,ydat@name.response], covariate.group=ydat@data[,ydat@name.covariates])

# Regrouping times - not needed here as everyone has the same times

# Regrouping scores - observed data
mybreaks <- c(0:9, 16, 25, 80)
x <- cut(obsmat$y, breaks=mybreaks, include.lowest = TRUE)
obsmat$score.group <- x

# R-base
data1<-NULL
for(ilev in sort(unique(obsmat$x))) {
  for(igroup in 0:1) {
    xtab<-obsmat[obsmat$x==ilev & obsmat$covariate.group==igroup,]
    x1<-data.frame(table(xtab$score.group))
    x1[,2]<-x1[,2]/sum(x1[,2])
    data2<-cbind(x=ilev, group=igroup, x1)
    data1<-rbind(data1, data2)
  }
}
data1 <- data1  %>%
  mutate(gender=ifelse(group==0,"Men","Women")) # here
colnames(data1)[3]<-"Score"
ggplot(data=data1, aes(x=x, y=Freq, group=gender, colour=gender)) + geom_line() + facet_wrap(.~Score, ncol=4)

# With tidyverse
counting.scores <- obsmat %>%
  group_by(x, covariate.group) %>%
  count(score.group)
number.samples <- obsmat %>%
  group_by(x, covariate.group) %>%
  summarise(n=n())
freq.scores <- number.samples %>%
  left_join(counting.scores, 
            by = c("x","covariate.group")) %>%
  mutate(freq=n.y/n.x)

ggplot(data=freq.scores, aes(x=x, y=freq, group=covariate.group, colour=as.factor(covariate.group))) + geom_line() + facet_wrap(.~score.group, ncol=4)


# Regrouping scores - simulated data
ysim.tab <- data.frame(irep=rep(1:nsim, each=dim(ydat@data)[1]), x=rep(obsmat$x, nsim), covariate.group=rep(obsmat$covariate.group, nsim), score.group=cut(ysim, breaks=mybreaks, include.lowest = TRUE))
sim.scores <- ysim.tab %>%
  group_by(irep, x, covariate.group) %>%
  count(score.group)
simfreq.scores <- number.samples %>%
  left_join(sim.scores, 
            by = c("x","covariate.group")) %>%
  mutate(freq=n.y/n.x)

simfreq.bands <- simfreq.scores %>%
  group_by(x, covariate.group, score.group) %>%
  summarise(lower=quantile(freq, c(0.05)), median=quantile(freq, c(0.5)), upper=quantile(freq, c(0.95)), freq=mean(freq)) 

ggplot(data=freq.scores, aes(x=x, y=freq, group=covariate.group, colour=as.factor(covariate.group))) + geom_line() + 
  geom_line(data=simfreq.bands, aes(x=x, y=median, group=covariate.group, colour=as.factor(covariate.group)), linetype="dashed") + 
  geom_ribbon(data=simfreq.bands, aes(ymin=lower, ymax=upper,  group=covariate.group, fill=as.factor(covariate.group)), alpha=0.2) + 
  facet_wrap(.~score.group, ncol=4)


countVPC <- function(obs, ysim, verbose=FALSE) {
  # obs: dataframe with the observed data, containing 3 columns in the following order
  ## subject ID
  ## predictor (x-axis)
  ## response (y-axis)
  # ysim: vector of simulated responses (the length of sim is assumed to be a multiple of the size of obs (nrow(obs)))
  # verbose: set to TRUE for error messages
  nsim <- length(ysim)/dim(obs)[1]
  if(!(as.integer(nsim)==nsim)) {
    if(verbose) message(paste("ysim must contain simulations of the response column in obs, please check dimensions between ysim of length",length(ysim),"and obs of size ",dim(obs)[1],"\n"))
    return()
  }
  
}

function(groups=NULL, verbose=FALSE) {
  # One plot produced per group
  if(is.null(groups)) {
    nunique<-sort(unique(obs[,3]))
    if(length(nunique)>12) {
      groups <- 12
      xcat<-cut(obs[,3], groups)
    } else xcat <- obs[,3]
  } else {
    if(is.numeric(groups)) 
      xcat<-cut(obs[,3], groups)
  }
  obs[,3]<-xcat
  countVPC(obs, sim, verbose=verbose)
}



aux.npdeplot.hist<-function(obsmat,  plot.opt, distrib="norm", nclass=10, sim.ypl=NULL) {
  # input
  ## obsmat: matrix with the data to plot, with columns
  ### x; variable to plot
  ### category: covariate category (if "all", overall plot)
  ## plot.opt: list of graphical options
  ## distrib: reference distribution plot (one of norm, unif)
  ## nclass: number of classes for the histogram
  ## sim.ypl: if given, a vector of simulated data for the variable to plot
  nrep<-100
  nameCovariate = plot.opt$which.cov    # nom de la covariable
  namesCategories = sort(unique(obsmat$category ))   # catégories pour la covariable
  if(is.null(namesCategories)) namesCategories<-c("all")
  numberCategories =  length(namesCategories)  # nombre de catégories pour la covariable
  
  if(plot.opt$approx.pi && is.null(sim.ypl)) plot.opt$approx.pi<-TRUE else nrep<-length(sim.ypl)/dim(obsmat)[1]
  # -----------------------------------------------------------------------------------
  # Observed histogram
  xhist<-hist(obsmat$x,breaks=nclass,plot=FALSE) # breaks overall
  nB<-length(xhist$breaks)
  xwidth<-unique(diff(xhist$breaks))
  obshist<-NULL
  for(icat in namesCategories) {
    obsmat.cov<-obsmat[obsmat$category==icat,]
    x1<-hist(obsmat.cov$x, breaks=xhist$breaks, plot=FALSE)
    obshist.cov<-data.frame(name=xhist$mids, value=x1$counts)
    zecat<-as.character(icat)
    obshist.cov<-data.frame(obshist.cov, category=zecat, stringsAsFactors = FALSE)
    obshist<-rbind(obshist,obshist.cov)
  }
  obsmat$category<-factor(obsmat$category, levels=namesCategories, ordered=TRUE)
  
  # -----------------------------------------------------------------------------------
  # PI for histogram
  pimat<-NULL
  if(plot.opt$bands) {
    alpha<-plot.opt$pi.size
    if(alpha>0.5) alpha<-1-alpha
    for(icat in namesCategories) {
      obsmat.cov<-obsmat[obsmat$category==icat,]
      ndat<-dim(obsmat.cov)[1]
      if(plot.opt$approx.pi) {
        sim.ypl<-switch(distrib,norm=rnorm(ndat*nrep),unif=runif(ndat*nrep)) # PI computed with 100 replicate
      }
      sim.ypl<-matrix(sim.ypl,nrow=ndat)
      tmat<-matrix(nrow=length(xhist$breaks)-1,ncol=nrep)
      for(j in 1:nrep) {
        xvec<-cut(sim.ypl[,j],breaks=xhist$breaks,include.lowest=TRUE, ordered_result=TRUE)
        tmat[,j]<-table(xvec)
      }
      row.names(tmat)<-names(table(xvec))
      bnds<-apply(tmat,1,quantile,c(alpha/2,0.5,1-alpha/2))
      bnds<-t(bnds)
      pimat.cov<-data.frame(x=xhist$mids, lower=bnds[,1], median=bnds[,2], upper=bnds[,3]) # using geom_crossbar
      zecat<-as.character(icat)
      pimat.cov<-data.frame(pimat.cov, category=zecat, stringsAsFactors = FALSE)
      pimat<-rbind(pimat, pimat.cov)
      # x1<-rep(xhist$breaks,each=2) # using geom_ribbon
      # x1<-x1[-c(1,length(x1))]
      # pimat.cov<-data.frame(x=x1, lower=rep(bnds[,1],each=2), median=rep(bnds[,2],each=2), upper=rep(bnds[,3],each=2))
      # zecat<-as.character(icat)
      # pimat.cov<-data.frame(pimat.cov, category=zecat, stringsAsFactors = FALSE)
      #      pimat<-rbind(pimat, pimat.cov)
    }
    pimat$category<-factor(pimat$category, levels=namesCategories, ordered=TRUE)
  }
  
  # -----------------------------------------------------------------------------------
  # Plot, facetting by covariate category
  
  if(is.null(plot.opt$xlim)) plot.opt$xlim<-c(min(xhist$breaks,na.rm=TRUE), max(xhist$breaks,na.rm=TRUE))
  if(is.null(plot.opt$ylim)) plot.opt$ylim<-c(0, max(c(obshist$value,pimat$upper),na.rm=TRUE))
  
  p <- ggplot(obsmat, aes(group=category)) +
    theme(plot.title = element_text(hjust = 0.5, size = plot.opt$size.sub),
          axis.title.x = element_text(size = plot.opt$size.xlab),
          axis.title.y = element_text(size = plot.opt$size.ylab),
          axis.text.x = element_text(size=plot.opt$size.text.x, color = ifelse(plot.opt$xaxt==TRUE,"black","white")),
          axis.text.y = element_text(size=plot.opt$size.text.y, color = ifelse(plot.opt$yaxt==TRUE,"black","white")),
          axis.line.x = element_line(color=ifelse(plot.opt$xaxt==TRUE,"black","white")),
          axis.line.y = element_line(color=ifelse(plot.opt$yaxt==TRUE,"black","white")),
          panel.background=element_rect("white"),
          panel.grid.major.x = element_line(ifelse(plot.opt$grid==TRUE,"grey80","white"),linetype = plot.opt$lty.grid),
          panel.grid.major.y = element_line(ifelse(plot.opt$grid==TRUE,"grey80","white"),linetype = plot.opt$lty.grid))+
    
    # coordinates x-y
    coord_cartesian(xlim=plot.opt$xlim, ylim=plot.opt$ylim) +
    
    # Plot bands
    {  if(plot.opt$bands==TRUE)
      geom_crossbar(data=pimat, aes(x=.data$x, y=median, ymin=lower, ymax=upper),
                    width=diff(xhist$mids)[1],
                    colour = plot.opt$col.ther,
                    fill = plot.opt$fill.bands,
                    alpha = plot.opt$alpha.bands,
                    linetype =  plot.opt$lty.ther,
                    size = plot.opt$lwd.ther )} +
    
    # Plot observed histograms
    geom_bar(data=obshist, aes(x=name, y=value),
             stat="identity",
             width = diff(xhist$mids)[1],
             colour = plot.opt$col.lobs,
             fill = plot.opt$fill,
             alpha = plot.opt$alpha,
             linetype =  plot.opt$lty,
             size = plot.opt$lwd ) +
    # x-y logscales
    { if (plot.opt$xlog == FALSE)
      scale_x_continuous(plot.opt$xlab,
                         scales::pretty_breaks(n = plot.opt$breaks.x))
    } +
    {if (plot.opt$ylog == FALSE)
      scale_y_continuous(plot.opt$ylab,
                         scales::pretty_breaks(n = plot.opt$breaks.y))
    } +
    { if (plot.opt$xlog == TRUE)
      scale_x_log10(plot.opt$xlab,
                    breaks = scales::trans_breaks("log10", function(x) 10 ^ x),labels = scales::trans_format("log10", scales::math_format(10 ^ .x)))
    } +
    {if (plot.opt$ylog == TRUE)
      scale_y_log10(plot.opt$ylab,
                    breaks = scales::trans_breaks("log10", function(x) 10 ^ x),labels = scales::trans_format("log10", scales::math_format(10 ^ .x)))
    } +
    #if log scales plot logticks
    { if (plot.opt$xlog == TRUE) annotation_logticks(sides = "b")} +
    { if (plot.opt$ylog == TRUE) annotation_logticks(sides = "l")} +
    
    # facet wrap over covariate categories
    facet_wrap(.~factor(category), nrow=1) +
    {if(numberCategories==1)
      theme(strip.background = element_blank(), strip.text.x = element_blank())
    } +
    {if (plot.opt$main!="") ggtitle(plot.opt$main)}
  
  return(p)
} # End function histogram
