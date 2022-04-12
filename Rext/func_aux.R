###########################	Computational fcts	#############################

#' @rdname saemix.internal
#' 
# #' @aliases cutoff cutoff.eps cutoff.max cutoff.res
#' @aliases normcdf norminv
#' @aliases transpsi transphi derivphi dtransphi
# #' @aliases compute.Uy compute.LLy conditional.distribution trnd.mlx gammarnd.mlx tpdf.mlx
# #' @aliases conditional.distribution_c conditional.distribution_d
#' 
#' @keywords internal

# Moved these functions to SaemixOutcome.R (define error models)
# cutoff<-function(x,seuil=.Machine$double.xmin) {x[x<seuil]<-seuil; return(x)}
# cutoff.max<-function(x) max(x,.Machine$double.xmin)
# cutoff.eps<-function(x) max(x,.Machine$double.eps)
# cutoff.res<-function(x,ares,bres) max(ares+bres*abs(x),.Machine$double.xmin)

# Inverse of the normal cumulative distribution fct: using erfcinv from ?pnorm
norminv<-function(x,mu=0,sigma=1)  mu-sigma*qnorm(x,lower.tail=FALSE)

# Truncated gaussian distribution (verifie par rapport a definition de erf/matlab)
normcdf<-function(x,mu=0,sigma=1)
  cutoff(pnorm(-(x-mu)/sigma,lower.tail=FALSE),1e-30)

transpsi<-function(psi,tr) {
  phi<-psi
  #  if(is.null(dim(psi))) phi<-as.matrix(t(phi),nrow=1)
  # ECO TODO: pourquoi ce test ??? Dans le cas ou psi est un vecteur ?
  i1<-which(tr==1) # log-normal
  phi[,i1]<-log(phi[,i1])
  i2<-which(tr==2) # probit
  phi[,i2]<-norminv(phi[,i2])
  i3<-which(tr==3) # logit
  phi[,i3]<-log(phi[,i3]/(1-phi[,i3]))
  if(is.null(dim(psi))) phi<-c(phi)
  return(phi)
}

transphi<-function(phi,tr) {
  psi<-phi
  #  if(is.null(dim(psi))) psi<-as.matrix(t(psi),nrow=1)
  i1<-which(tr==1) # log-normal
  psi[,i1]<-exp(psi[,i1])
  i2<-which(tr==2) # probit
  psi[,i2]<-normcdf(psi[,i2])
  i3<-which(tr==3) # logit
  psi[,i3]<-1/(1+exp(-psi[,i3]))
  if(is.null(dim(phi))) psi<-c(psi)
  return(psi)
}
derivphi<-function(phi,tr) {
  # Fonction calculant ??? (only used to plot parameter distributions/histograms of individual parameters)
  psi<-phi # identite
  i1<-which(tr==1) # log-normal
  psi[,i1]<-1/exp(phi[,i1])
  i2<-which(tr==2) # probit
  psi[,i2]<-1/(sqrt(2*pi))*exp(-(phi[,i2]**2)/2)
  i3<-which(tr==3) # logit
  psi[,i3]<-2+exp(phi[,i3])+exp(-phi[,i3])
  if(is.null(dim(phi))) psi<-c(psi)
  return(psi)
}

dtransphi<-function(phi,tr) {
  # Fonction computing the derivative of h, used to compute the Fisher information matrix
  psi<-phi
  if(is.null(dim(phi))) {
    dpsi<-as.matrix(t(rep(1,length(phi))))
    psi<-as.matrix(t(phi),nrow=1)
  } else 
    dpsi<-matrix(1,dim(phi)[1],dim(phi)[2])
  i1<-which(tr==1) # log-normal
  dpsi[,i1]<-exp(psi[,i1])
  i2<-which(tr==2) # probit
  dpsi[,i2]<-1/dnorm(qnorm(dpsi[,i2]))   # derivee de la fonction probit, dqnorm <- function(p) 1/dnorm(qnorm(p))
  i3<-which(tr==3) # logit
  dpsi[,i3]<-1/(2+exp(-psi[,i3])+exp(psi[,i3]))
  if(is.null(dim(phi))) dpsi<-c(dpsi)
  return(dpsi)
}

###########################	Computing the acceptance ratio	#############################
# par: parameters to be optimised
# y: observations
# f: model predictions
# error.function: error function 

ssq<-function(par,y,f,error.function) { # Sum of squares; need to put par first as these parameters are optimised by optim
  g<-error.function(f,par)
  e<-sum(((y-f)**2/g**2)+2*log(g))
  return(e)
}

sumSquare<-function(par,y,f,error.function) { # Sum of squares; need to put par first as these parameters are optimised by optim
  g<-error.function(f,par)
  e<-sum(((y-f)**2/g**2))
  return(e)
}

# Input
## phiM: matrix of parameters
## Dargs: list, here we use the elements transform.par, model (structural model), IdM, XM, yM, ind.ioM (a vector of 0/1 indicating which elements of the square matrix DYF have data), list of outcomes (type, error.model, error.function)
## error.parameters: list of the error parameters for the different outcomes in the model (NULL for non-continuous outcome)
## DYF: a matrix of dimension nmax=max(n_i) times (N*nb.chains); the column for subject i has the last nmax-n_i elements set to 0, so that colSums sums on the observations for that subject 
# Output 
## U: vector with either sum(-LL) (continuous models, minus the constant sum(log(1/sqrt(2*pi)))) or sum(-logpdf) (likelihood models)

compute.LLy<-function(phiM, Dargs, DYF, error.parameters) {
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-Dargs$model(psiM,Dargs$IdM,Dargs$XM)
  lpred<-fpred
  for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
  for(iout in 1:length(Dargs$outcome)) {
    idx1<-which(Dargs$XM$ytype==iout)
    # print(summary(fpred[idx1]))
    if(Dargs$outcome[[iout]]@type=="continuous") {
      if(Dargs$outcome[[iout]]@error.model=="exponential") fpred[idx1]<-log(cutoff(fpred[idx1]))
      gpred<-Dargs$outcome[[iout]]@error.function(fpred[idx1], error.parameters[[iout]])
      lpred[idx1]<-0.5*((Dargs$yM[idx1]-fpred[idx1])/gpred)**2+log(gpred)
      # print(summary(gpred))
    } else lpred[idx1]<- (-fpred[idx1])
  }
  DYF[Dargs$ind.ioM]<-lpred
  U<-colSums(DYF)
  return(U)
}

########################################################### Remove
# Deprecated in version 4.0 - multiple outcomes
error<-function(f,ab,etype) { # etype: error model
  g<-f
  for(ityp in sort(unique(etype))) {
    g[etype==ityp]<-error.typ(f[etype==ityp],ab[((ityp-1)*2+1):(ityp*2)])
  }
  return(g)
}
error.typ<-function(f,ab) {
  #  g<-cutoff(ab[1]+ab[2]*abs(f))
  g<-cutoff(sqrt(ab[1]^2+ab[2]^2*f^2))  # Johannes 02/21
  return(g)
}
