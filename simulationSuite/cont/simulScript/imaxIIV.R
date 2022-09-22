# Folders
simulDir <-"/home/eco/work/saemix/saemixextension"

# Scenario
scenario <- "imaxRichComb"
datDir <- file.path(simulDir,"simulationSuite","cont", scenario,"data")
resDir <- file.path(simulDir,"simulationSuite","cont", scenario,"results")

# Create directories
dir.create(file.path(simulDir,"simulationSuite","cont", scenario))
dir.create(file.path(simulDir,"simulationSuite","cont", scenario, "data"))
dir.create(file.path(simulDir,"simulationSuite","cont", scenario, "results"))

# Parameters
nsim<-100
namsimdat<-"pdimax"
parpop<-c(100,0.7,500,2)
parcov<-(-0.5) # ED50=500 in group 0, 300 in group 1
nampar<-c("E0","Imax","ED50","gamma")
omega<-diag(c(0.01,0,0.09,0)) # 10% IIV on E0 and 30% IIV on ID50, no correlation
respar<-c(5,0.2) # additive and proportional error terms

# settings: n=6, N=100, 2 groups
tdose <- c(0, 200, 400, 500, 700, 1000)
nsuj <- 100

sigma<-respar
if(length(grep("Add", scenario))>0) sigma<-respar[1]
if(length(grep("Prop", scenario))>0) sigma<-respar[2]
pvrai <- c(parpop[1:3], parcov, parpop[4], omega[1,1],omega[3,3], respar)

sigmoidImax<-function(psi,id,xidep) {
  # input:
  #   psi : matrix of parameters (3 columns, E0, Emax, E50)
  #   id : vector of indices 
  #   xidep : dependent variables (same nb of rows as length of id)
  # returns:
  #   a vector of predictions of length equal to length of id
  dose<-xidep[,1]
  e0<-psi[id,1]
  imax<-psi[id,2]
  ed50<-psi[id,3]
  gamma<-psi[id,4]
  f<-e0*(1-imax*(dose**gamma)/((ed50**gamma)+(dose**gamma)))
  return(f)
}

# Simulation
if(length(grep("Prop", scenario))>0) set.seed(2009227419)
if(length(grep("Add", scenario))>0) set.seed(43509722)
if(length(grep("Comb", scenario))>0) set.seed(6185678)

for(isim in 1:nsim) {
  xidep<-data.frame(dose=rep(tdose, nsuj))
  psi1 <- data.frame(e0=parpop[1]*exp(rnorm(nsuj, mean=0, sd=sqrt(omega[1,1]))), imax=parpop[2], 
                     ed50=parpop[3]*exp(rnorm(nsuj, mean=0, sd=sqrt(omega[3,3]))), gamma=parpop[4])
  psi1[51:100,3]<-psi1[51:100, 3]*exp(parcov[1])
  id1<-rep(c(1:nsuj),each=length(tdose))
  
  fpred <- sigmoidImax(psi1, id1, xidep)
  if(length(grep("Add", scenario))>0) fpred<-fpred+rnorm(length(fpred), mean=0, sd=respar[1])
  if(length(grep("Prop", scenario))>0) fpred<-fpred*(1+rnorm(length(fpred), mean=0, sd=respar[2]))
  if(length(grep("Comb", scenario))>0) {
    gpred<-sqrt(respar[1]**2 + (respar[2]*fpred)**2)
    fpred<-fpred+gpred*rnorm(length(fpred), mean=0, sd=1)
  }
  xdat<-data.frame(id=id1, dose=xidep[,1], eff=fpred, group=0)
  xdat$group[xdat$id>50]<-1
  namfich<-paste('data_',namsimdat,isim,".tab",sep="")
  write.table(xdat, file.path(datDir,namfich), quote=F, row.names = FALSE)
}

# plot(xidep[,1], fpred)


