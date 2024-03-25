# simulated.annealing: if FALSE, deactivate simulated annealing. By default, use nbiter.sa from the options slot of saemixObject, or the value passed through saemix.options.
# CI: one of "final", "initial", "adaptive" (defaults to adaptive)

saemix.bootstrap2<-function(saemixObject, method="conditional", nboot=200, nsamp=100, saemix.options=NULL, simulated.annealing=TRUE, CI="adaptive") {
  if(method!="case") {
    if(saemixObject@model@modeltype!="structural" & (is.null(body(saemixObject@model@simulate.function)) | length(formals(saemixObject@model@simulate.function))!=3)) {
      if(saemixObject@options$warnings) message("A simulation function needs to be provided for non Gaussian models to obtain bootstrap distributions by other methods than the Case bootstrap. This function needs to have the same structure as the model function and return simulated values based on the same model. \nPlease provide a simulation function the simulate.function slot of the model or use method='case' for Case bootstrap. \nExiting bootstrap.")
      return(NULL)
    }
  }
  if(method=="residual" | method=="conditional") {
    saemixObject<-saemix.predict(saemixObject) # estimate individual parameters and compute residuals (currently iwres are needed also for conditional but need to modify this in a further extension to cNP ECO TODO)
  }
  if(!(CI %in% c("final", "initial", "adaptive"))) {
    if(saemixObject@options$warnings) cat("Option",CI,"not recognised, please use one of final, initial, or adaptive (starting from final, testing various changes and trying initial last if all else fails).")
    CI<-"adaptive"
  }
  if(method=="conditional") {
    ndone <- dim(saemixObject@results@phi.samp)
    if(!is.null(ndone)) ndone<-ndone[3] else ndone<-0
    if(ndone<nsamp) {
      if(saemixObject@options$warnings) message("Not enough samples in the object, sampling from the conditional distribution\n")
      saemixObject<-conddist.saemix(saemixObject, nsamp=nsamp) # estimate conditional distributions and sample residuals
    }
    eta.sampc<-centerDist.NPcond(saemixObject, nsamp=nsamp) # Center eta samples from the conditional distribution, to avoid doing this repeatedly
  }
  if(is.null(saemix.options)) {
    #      saemix.options<-list(directory="current",fix.seed=FALSE,map=FALSE,ll.is=FALSE,displayProgress=FALSE,save.graphs=FALSE,print=FALSE)
    saemix.options<-saemixObject["options"]
    saemix.options$directory<-"current"
    saemix.options$fix.seed<-FALSE
    saemix.options$map<-FALSE   # Only parameter estimates are required for bootstrap
    saemix.options$fim<-FALSE
    saemix.options$displayProgress<-FALSE 
    saemix.options$save.graphs<-FALSE
    saemix.options$save<-FALSE
    saemix.options$ll.is<-FALSE
    saemix.options$print<-FALSE
  }
  verbose <- saemix.options$warnings
  if(!simulated.annealing) saemix.options$nbiter.sa <- 5 # need at least a few iterations of SA for the burn-in
  if(saemix.options$nbiter.sa<5) saemix.options$nbiter.sa<-5
  if(saemixObject@model@modeltype=="structural") idx.eps<-saemixObject@model@indx.res else idx.eps<-integer(0)
  idx.iiv<-saemixObject@model@indx.omega
  idx.rho<-which(saemixObject@model@covariance.model[lower.tri(saemixObject@model@covariance.model)]==1)
  bootstrap.distribution<-failed.runs<-data.frame()
  nelements <- length(saemixObject@results@fixed.effects)+length(idx.iiv)+length(idx.rho)+length(idx.eps)
  # Starting point: estimates from the fit 
  model.boot<-saemixObject["model"]
  if(CI %in% c("final","adaptive") ) {
    model.boot@psi0 <- model.boot["betaest.model"]
    model.boot@psi0[model.boot["betaest.model"]==1]<-saemixObject@results@fixed.effects
  }
  for(iboot in 1:nboot) {
    if(method=="case")  
      data.boot <- dataGen.case(saemixObject)
    if(method=="residual")
      data.boot <- dataGen.NP(saemixObject,conditional=FALSE)
    if(method=="conditional")
      data.boot <- dataGen.NP(saemixObject, nsamp=nsamp,eta.sampc=eta.sampc, conditional=TRUE)
    if(method=="parametric")
      data.boot <- dataGen.Par(saemixObject)
    fit.boot<-try(saemix(model.boot, data.boot, saemix.options))
    if(is(fit.boot,"try-error") & CI=="adaptive") { # try changing the initial parameters (needed for categorical models)
      model.boot2<-model.boot
      model.boot2@psi0[model.boot2["betaest.model"]==1]<-saemixObject@results@fixed.effects*.5
      fit.boot<-try(saemix(model.boot2, data.boot, saemix.options))
      if(is(fit.boot,"try-error")) {
        model.boot2@psi0[model.boot2["betaest.model"]==1]<-saemixObject@results@fixed.effects*.2
        fit.boot<-try(saemix(model.boot2, data.boot, saemix.options))
        if(is(fit.boot,"try-error")) {
          model.boot2@psi0[1,]<-saemixObject@model@psi0[,1]
          fit.boot<-try(saemix(model.boot2, data.boot, saemix.options))
        }
      }
    }
    if(is(fit.boot,"try-error")) {
      l1<-c(iboot,rep(NA,nelements))
      failed.runs <- rbind(failed.runs, c(iboot, fit.boot))
    } else {
      res<-fit.boot@results
      l1<-c(iboot,res@fixed.effects, diag(res@omega)[idx.iiv])
      if(length(idx.rho)>0) l1<-c(l1,res@omega[lower.tri(res@omega)][idx.rho])
      if(length(idx.eps)>0) l1<-c(l1, res@respar[idx.eps])
      if(length(res@ll.lin)>0) l1<-c(l1, res@ll.lin)
      
    }
    #    l1<-c(iboot,res@fixed.effects, diag(res@omega)[idx.iiv],res@omega[lower.tri(res@omega)][idx.rho],res@respar[idx.eps], res@se.fixed, res@se.omega[idx.iiv],res@se.cov[lower.tri(res@se.cov)][idx.rho], res@se.respar[idx.eps],res@ll.lin)
    bootstrap.distribution<-rbind(bootstrap.distribution,l1) 
  }
  # Names
  nampar<-colnames(saemixObject@model@covariance.model)
  namcol<-c(saemixObject@results@name.fixed, saemixObject@results@name.random)
  if(length(idx.rho)>0) {
    for(i in 1:(length(nampar)-1)) {
      for(j in (i+1):length(nampar)) {
        if(saemixObject@model@covariance.model[i,j]==1) {
          namcol<-c(namcol,paste("cov.",nampar[i],nampar[j],sep=""))
        }
      }
    }
  }
  if(length(idx.eps)>0) namcol<-c(namcol,saemixObject@model@name.sigma[idx.eps])
  if(length(res@ll.lin)>0) namcol<-c(namcol,"LL.lin")
  #  namcol<-c(saemixObject@results@name.fixed,saemixObject@results@name.random,namcol, saemixObject@results@name.sigma[saemixObject@results@indx.res])
  #  colnames(bootstrap.distribution)<-c("Replicate",namcol,paste("SE",namcol,sep="."),"LL.lin")
  colnames(bootstrap.distribution)<-c("Replicate",namcol)
  if(verbose && dim(failed.runs)[1]>0) {
    cat(dim(failed.runs)[1],"failed:\n")
    print(head(failed.runs))
  }
  return(bootstrap.distribution)
}
