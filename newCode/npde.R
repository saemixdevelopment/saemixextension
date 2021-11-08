# Alternate call to library(npde) to compute the npde

library(npde)

npdeSaemix<-function(saemixObject, nsim=1000) {
  if(saemixObject@results@status!="fitted") {
    if(saemixObject@options$warnings) message("Please fit the model first\n")
    return()
  }
  namobs<-saemixObject@data@data
  if(saemixObject@sim.data@nsim==0) {
    if(saemixObject@options$warnings) message("Simulating from the model to compute npde, nsim=",nsim,"\n")
    saemixObject<-simulate(saemixObject, nsim=nsim)
  }
  namsim<-data.frame(idsim=saemixObject@sim.data@datasim[,"idsim"],
                     xsim=rep(namobs[,saemixObject@data@name.X],saemixObject@sim.data@nsim),
                     ysim=saemixObject@sim.data@datasim[,"ysim"])
  npdeObject<-autonpde(namobs, namsim, iid=saemixObject@data@name.group, ix=saemixObject@data@name.X, 
           iy=saemixObject@data@name.response, icens="cens", imdv="mdv", icov=saemixObject@data@name.covariates,
           units=list(x=saemixObject@data@units$x, y=saemixObject@data@units$y))
  
  return(npdeObject)
}

xnpde<-npdeSaemix(theo.fit1)

plot(xnpde, plot.type="vpc")

plot(xnpde, plot.type="default")

plot(xnpde, plot.type="cov.x.scatter")

# Check size if cens=1 or mdv=1 for some subjects