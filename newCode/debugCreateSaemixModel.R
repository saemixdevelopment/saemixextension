# Setup
if(FALSE) {
  saemixDir<-"/home/eco/work/saemix/saemixextension"
  progDir<-file.path(saemixDir,"R")
  progDirExt<-file.path(saemixDir,"Rext")
  source(file.path(progDir,"aaa_generics.R"))
  source(file.path(progDirExt,"SaemixOutcome.R"))
  source(file.path(progDirExt,"SaemixData.R"))
  source(file.path(progDirExt,"SaemixCovariateModel.R"))
  source(file.path(progDirExt,"SaemixCovariate.R"))
  source(file.path(progDirExt,"SaemixVarLevel.R"))
  source(file.path(progDirExt,"SaemixParameter.R"))
  source(file.path(progDirExt,"SaemixParameter-methods.R"))
  source(file.path(progDirExt,"SaemixModel.R"))
  source(file.path(progDirExt,"SaemixModel-methods.R"))
}

## Parameter model (with variability and covariates)
lpar <- list(ka=saemixPar(mu.start=2),
             cl=saemixPar(mu.start=20, covariate=c(age=contCov(), sex=binCov(), wt=contCov(beta=0.75, beta.fix=1))),
             vd=saemixPar(mu.start=10, covariate=c(wt=contCov(beta=1, beta.fix=1)), rho.param=c("cl","ka")),
             ic50=saemixPar(mu.start=2, distribution="normal", omega.level=c("id","occ"), omega.start = c(1,0.5), rho.param=list(c("cl"))),
             imax=saemixPar(distribution="logit", omega.level=c()))
lpar5 <- list(ka=lognormalPar(mu.start=2),
             cl=lognormalPar(mu.start=20, covariate=c(age=contCov(), sex=binCov(), wt=contCov(beta=0.75, beta.fix=1))),
             vd=lognormalPar(mu.start=10, covariate=c(wt=contCov(beta=1, beta.fix=1)), rho.param=c("cl","ka")),
             ic50=normalPar(mu.start=2, omega.level=c("id","occ"), omega.start = c(1,0.5), rho.param=list(c("cl"))),
             imax=logitPar())
# imax=logitPar(omega.fix=1, omega.start=c(0.2))
             
## Just the parameter names
lpar2 <- c("ka","cl","vd","ic50","imax") # imax will be considered as lognormal
## Just the parameter types
lpar3 <- c(ka="lognormal",cl="lognormal",vd="lognormal",ic50="lognormal",imax="logit") 
## A mix of names, types and named types
lpar4 <- c("ka","cl","vd","lognormal",imax="logit") 

## Outcome
### list of outcomes objects
lout<-list(conc=saemixOutcome(unit="mg/L", model="combined1", start=c(0.5, 0.3), fix=c(1,0)), 
           effect=saemixOutcome(unit="%", model="proportional"),
           pain=saemixOutcome(type="categorical", levels=1:3)) # 3rd response to be ignored
### lists defined by continousOutcome() or discreteOutcome()
lout2<-list(conc=continuousOutcome(unit="mg/L",model="combined1", start=c(0.5, 0.3), fix=c(1,0)), 
            effect=continuousOutcome(unit="%", model="proportional"),
            pain=discreteOutcome(levels=1:3),discreteOutcome(type="event"))
### Just the parameter names+types or names alone (defaults to continuous outcome)
lout3<-c(conc="continuous", effect="continuous", pain="categorical")
lout4<-c("conc","effect","pain") # pain will be considered as continuous
lout5<-c(conc="continuous","continuous",pain="categorical")

### If given as a vector, automatically converted to a list
if(FALSE) {
  lout5<-c(conc=saemixOutcome(unit="mg/L", model="combined1", start=c(0.5, 0.3), fix=c(1,0)), 
           effect=saemixOutcome(unit="%", model="proportional"),
           pain=saemixOutcome(type="categorical", levels=1:3)) # 3rd response to be ignored
  is(lout5, "list")  
}

## Mixed - gets automatically converted to a list
lout6<-c(conc=saemixOutcome(unit="mg/L", model="combined1", start=c(0.5, 0.3), fix=c(1,0)), 
         conc="continuous","categorical","effect")

## Model
model1cptdirect<-function(psi,id,xidep) { 
  tim<-xidep[,1]
  dose<-xidep[,2]
  ytype<-xidep$ytype
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  ic50<-psi[id, 4]
  imax<-psi[id, 5]
  k<-CL/V
  ypk<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  ypd<-100*(1-imax*ypk/(ypk+ic50))
  ypk[ytype==2]<-ypd[ytype==2]
  # need here to change to lpdf for discrete data models
#  yresp<-as.integer(ypd>25)+as.integer(ypd>75)
#  ypk[ytype==3]<-yresp[ytype==3]
  return(ypk)
}

######################################################################
# Parameter

## Parameter given as just the names
parameter<-lpar2
convertArg2Parameter(lpar)

parameter<-lpar3
convertArg2Parameter(parameter)

## Parameter given as a list of SaemixPar objects
parameter<-lpar
convertArg2Parameter(parameter)

## Mixed list
parameter<-lpar4
convertArg2Parameter(parameter)


xmat1<-diag(6)
xmat1[2,4]<-xmat1[4,2]<-xmat1[3,4]<-xmat1[4,3]<-xmat1[5,6]<-xmat1[6,5]<-1
completeBlocks(xmat1)

# Create a SaemixVarLevel object using lpar
getVarianceModel(lpar)

# same, output the individual components as a list
getVarianceModel(lpar, output ="list")

######################################################################
# Outcome
outcome<-lout
is(outcome,"list")

## Outcome given as a vector with just the names
outcome<-lout4
is(outcome,"list")
saemix.outcome<-vector(mode="list", length=length(outcome))
if(is.null(names(outcome))) names(saemix.outcome)<-outcome
for(i in 1:length(outcome)) saemix.outcome[[i]]<-saemixOutcome(name=outcome[i],type="continuous")

## Outcome given as a vector of names+type
outcome<-lout3
is(outcome,"list")
saemix.outcome<-vector(mode="list", length=length(outcome))
if(is.null(names(outcome))) names(saemix.outcome)<-outcome else names(saemix.outcome)<-names(outcome)
for(i in 1:length(outcome)) saemix.outcome[[i]]<-saemixOutcome(name=names(outcome)[i],type=outcome[i])

## Outcome given as a list of SaemixOutcome objects
outcome<-lout
saemix.outcome<-outcome

## Outcome given as a list of continuousOutcome() and discreteOutcome()
outcome<-lout2
is(outcome,"list")
saemix.outcome<-vector(mode="list", length=length(outcome))
if(is.null(names(outcome))) names(saemix.outcome)<-outcome else names(saemix.outcome)<-names(outcome)
for(i in 1:length(outcome)) {
  if(names(saemix.outcome)[i]=="") {
    if(is.null(outcome[[i]]$name) || outcome[[i]]$name=="") names(saemix.outcome)[i]<-paste0("out",i)
  }
}
for(i in 1:length(outcome))
  saemix.outcome[[i]] <- createSaemixOutcome(outcome[[i]], name=names(saemix.outcome)[i])

## Mixed list
outcome<-lout5

smx.out1 <- convertArg2Outcome(lout)
smx.out2 <- convertArg2Outcome(lout2)
smx.out3 <- convertArg2Outcome(lout3)
smx.out4 <- convertArg2Outcome(lout4)
smx.out5 <- convertArg2Outcome(lout5)
smx.out6 <- convertArg2Outcome(lout6)


######################################################################

xobj<-saemixModel(model=model1cptdirect, outcome=lout, parameter=lpar, verbose=TRUE)

checkNested<-function(parameter, var.level) {
  # Check variability levels are nested in the same order for all parameters
  
  # Need a second check when data is coupled with a model
}
