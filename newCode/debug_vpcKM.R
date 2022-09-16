############################################################################
# To run, run saemix3_tteModel.R first
source("/home/eco/work/saemix/saemixextension/paperSaemix3/saemix3_tteModel.R")

# Plot functions for discrete VPC
source("/home/eco/work/saemix/saemixextension/R/func_discreteVPC.R")
xpl1 <- discreteVPCTTE(simtte.fit)
plot(xpl1)

xpl2 <- discreteVPC(simtte.fit, outcome="TTE")
plot(xpl2)

############################################################################
# Compute KM

# obs: dataframe with the observation data
## id: subject id
## time: event time
## event: 0 for no event, 1 for event
## cens: 0 for observed event/no event, 1 for censored (optional, if not present all 1 for events are assumed to be observed events and a 0 at the last observation time indicates a censored event)
# sim:
## irep: simulation number
## id: subject id (replicated)
## time: simulated times
obsdat <- data.frame(id=lung.saemix$id, time=lung.saemix$time, event=lung.saemix$status, cens=lung.saemix$cens)
simdat <- data.frame(irep=ysim.tte@sim.data@datasim$irep,id=ysim.tte@sim.data@datasim$idsim, time=ysim.tte@sim.data@datasim$ysim)

# Observed KM
event.obs<-obsdat[obsdat$time>0,]
t1<-table(event.obs$time) # nb of subjects dropping off at different times (event or censoring)
tab.obs <- data.frame(tobs=as.numeric(names(t1)),nk=c(t1))
if(match("cens", colnames(obsdat))) { # cens indicates censoring (cens=1) versus event (cens=0)
  t2<-table(event.obs$time[event.obs$cens==0]) # nb of subjects with an event 
  tab.obs$dk<-0
  idx<-match(names(t2), names(t1))
  tab.obs$dk[idx]<-t2
} else { # all events are observed
  tab.obs$dk<-tab.obs$nk
}
nk.obs <- c(sum(tab.obs$nk), sum(tab.obs$nk)-cumsum(tab.obs$nk))
tab.obs$nk <- nk.obs[1:(length(nk.obs)-1)]
tab.obs$sk <- (1-tab.obs$dk/tab.obs$nk)
tab.obs$km <- cumprod(tab.obs$sk)

# Check compared to survfit, pretty close (maybe need to plot not lines() but something stepwise to make it exactly similar)
lung.surv<-lung.saemix[lung.saemix$time>0,]
lung.surv$status<-lung.surv$status+1
Surv(lung.surv$time, lung.surv$status) # 1=censored, 2=dead
yfit1 <- survfit(Surv(time, status) ~ 1, data = lung.surv)
plot(yfit1)
lines(tab.obs$tobs, tab.obs$km, col="Red")

# KM for simulated data
# If no censored column, add it [no, simulated data only has observed events...]
# if(is.na(match("cens",colnames(simdat)))) {
#   simdat$cens<-as.integer(simdat$time>0)
# }
if(is.na(match("event",colnames(simdat)))) {
  simdat$event<-as.integer(simdat$time>0)
}

ngrid<-200 # nb of grid points
ymax<-max(obsdat$time)
ytime<- seq(0,ymax,length.out=ngrid)

isim<-1
xsim<-simdat[simdat$irep==isim,]
tsim <- xsim$time[xsim$time>0]
t1 <- table(tsim)
tab1<-data.frame(time=as.numeric(names(t1)),dk=c(t1))
nk.sim <- c(sum(tab1[,2]), sum(tab1[,2])-cumsum(tab1[,2]))
tab1$nk <- nk.sim[1:(length(nk.sim)-1)]
tab1$sk <- (1-tab1$dk/tab1$nk)
tab1$km <- cumprod(tab1$sk)

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
  if(nminus==dim(y)[1]) return(ykm[nminus])
  y$km[nminus] + (x-y$time[nminus])*(y$km[nminus+1]-y$km[nminus])/(y$time[nminus+1]-y$time[nminus])
}

interpol.locf(km.sim[2,1], tab1)
interpol.lin(km.sim[2,1], tab1)

interpol.locf(8, tab1)
interpol.lin(8, tab1)

km.locf <- sapply(ytime, interpol.locf, tab1)
km.lin <- sapply(ytime, interpol.lin, tab1)
# Check with survfit
yfit1 <- survfit(Surv(time, event) ~ 1, data = xsim)
plot(yfit1)
lines(tab1$time, tab1$km, col="Red")
lines(ytime, km.locf, col="Blue") # interpolated
lines(ytime, km.lin, col="green", lty=2) # interpolated

# Loop
nsim<-length(unique(simdat$irep))
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
  km <- sapply(ytime, interpol.locf, tab1)
  km.sim<-cbind(km.sim, km)
}
tab.quant<-t(apply(km.sim, 1, quantile, c(0.025, 0.5, 0.975)))
tab.quant <- cbind(time=ytime,tab.quant)


# Quantiles by simulations compared to standard KM interval
plot(yfit1)
for(icol in 1:3) 
  lines(tab.quant[,1], tab.quant[,(icol+1)], lty=2, col="Red")

# 
