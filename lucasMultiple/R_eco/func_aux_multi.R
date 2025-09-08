###########################	Computational fcts	#############################

#' @rdname saemix.internal
#' 
#' @aliases cutoff cutoff.eps cutoff.max cutoff.res
#' @aliases normcdf norminv
#' @aliases error error.typ ssq
#' @aliases transpsi transphi derivphi dtransphi
#' @aliases compute.Uy compute.LLy conditional.distribution trnd.mlx gammarnd.mlx tpdf.mlx
#' @aliases conditional.distribution_c conditional.distribution_d
#' 
#' @keywords internal

# Moved these functions to SaemixOutcome.R (define error models) in Rext
cutoff<-function(x,seuil=.Machine$double.xmin) {x[x<seuil]<-seuil; return(x)}
cutoff.max<-function(x) max(x,.Machine$double.xmin)
cutoff.eps<-function(x) max(x,.Machine$double.eps)
cutoff.res<-function(x,ares,bres) max(ares+bres*abs(x),.Machine$double.xmin)

# Inverse of the normal cumulative distribution fct: using erfcinv from ?pnorm
norminv<-function(x,mu=0,sigma=1)  mu-sigma*qnorm(x,lower.tail=FALSE)

# Truncated gaussian distribution (verifie par rapport a definition de erf/matlab)
normcdf<-function(x,mu=0,sigma=1)
  cutoff(pnorm(-(x-mu)/sigma,lower.tail=FALSE),1e-30)

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

# ssq<-function(ab,y,f) { # Sum of squares
# 	g<-abs(ab[1]+ab[2]*f)
# 	e<-sum(((y-f)/g)**2+2*log(g))
# 	return(e)
# }
# ssq<-function(ab,y,f,ytype) { # Sum of squares
#   g<-f
#   for(ityp in sort(unique(ytype))) {
#     g[ytype==ityp]<-(ab[((ityp-1)*2+1)]+f[ytype==ityp]*ab[(ityp*2)])**2
#   }
#   e<-sum(((y-f)**2/g)+log(g))
#   return(e)
# }

ssq<-function(ab,y,f,etype) { # Sum of squares; need to put ab first as these parameters are optimised by optim
  g<-(error(f,ab,etype))
  e<-sum(((y-f)**2/g**2)+2*log(g))
  return(e)
}

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


compute.Uy<-function(b0,phiM,pres,args,Dargs,DYF) {
# Attention, DYF variable locale non modifiee en dehors
  args$MCOV0[args$j0.covariate]<-b0
  phi0<-args$COV0 %*% args$MCOV0
  phiM[,args$i0.omega2]<-do.call(rbind,rep(list(phi0),args$nchains))
  psiM<-transphi(phiM,Dargs$transform.par)
  if (Dargs$modeltype=="structural"){
  fpred<-Dargs$structural.model(psiM,Dargs$IdM,Dargs$XM)
  for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
  gpred<-error(fpred,pres,Dargs$XM$ytype)
  DYF[args$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
  } else {
    fpred<-Dargs$structural.model(psiM,Dargs$IdM,Dargs$XM)
    for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
    DYF[args$ind.ioM]<- -fpred
  }
  U<-sum(DYF)
  return(U)
}

compute.LLy<-function(phiM,args,Dargs,DYF,pres) {
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-Dargs$structural.model(psiM,Dargs$IdM,Dargs$XM)
  for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
  if (Dargs$modeltype=="structural"){
    gpred<-error(fpred,pres,Dargs$XM$ytype)
    DYF[args$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
  } else {
    DYF[args$ind.ioM]<- -fpred
  }
  U<-colSums(DYF)
  return(U)
}



conditional.distribution_c<-function(phi1,phii,idi,xi,yi,mphi,idx,iomega,trpar,model,pres,err) {
  phii[idx]<-phi1
  psii<-transphi(matrix(phii,nrow=1),trpar)
  if(is.null(dim(psii))) psii<-matrix(psii,nrow=1)
  fi<-model(psii,idi,xi)
#  if(err=="exponential") # Reverted this bit to the previous version to avoid a compiler error, not sure why it was changed...
#    fi<-log(cutoff(fi))
#  if (!(is.null(pres)) && pres[1] == pres) { # package compile throws an error when comparing a vector of length 2 (pres) to a vector of length 1
  # gi <- cutoff(pres[1])
  # } 
  # else{
  #   gi<-error(fi,pres) #    cutoff((pres[1]+pres[2]*abs(fi)))
  # }
  ind.exp<-which(err=="exponential")
  for(ityp in ind.exp) 
    fi[xi$ytype==ityp]<-log(cutoff(fi[xi$ytype==ityp]))
  gi<-error(fi,pres,xi$ytype)      #    cutoff((pres[1]+pres[2]*abs(fi)))
  Uy<-sum(0.5*((yi-fi)/gi)**2+log(gi))
  dphi<-phi1-mphi
  Uphi<-0.5*sum(dphi*(dphi%*%iomega))
  return(Uy+Uphi)
}

conditional.distribution_d<-function(phi1,phii,idi,xi,yi,mphi,idx,iomega,trpar,model) {
  phii[idx]<-phi1
  psii<-transphi(matrix(phii,nrow=1),trpar)
  if(is.null(dim(psii))) psii<-matrix(psii,nrow=1)
  fi<-model(psii,idi,xi)
  Uy <- sum(-fi)
  dphi<- phi1-mphi
  Uphi<- 0.5*sum(dphi*(dphi%*%iomega))
  return(Uy+Uphi)
}


trnd.mlx<-function(v,n,m) {
  r<-rnorm(n*m)*sqrt(v/2/gammarnd.mlx(v/2,n,m))
  return(r=matrix(r,nrow=n,ncol=m))
}

gammarnd.mlx<-function(a,n,m) {
  nm<-n*m
  y0 <- log(a)-1/sqrt(a)
  c <- a - exp(y0)
  b <- ceiling(nm*(1.7 + 0.6*(a<2)))
  y <- log(runif(b))*sign(runif(b)-0.5)/c + log(a)
  f <- a*y-exp(y) - (a*y0 - exp(y0))
  g <- c*(abs((y0-log(a))) - abs(y-log(a)))
  reject <- ((log(runif(b)) + g) > f)
  y<-y[!reject]
  if(length(y)>=nm) x<-exp(y[1:nm]) else 
    x<-c(exp(y),gammarnd.mlx(a,(nm-length(y)),1))
#  x<-matrix(x,nrow=n,ncol=m) # not useful ?
  return(x)
}

tpdf.mlx<-function(x,v) {
# TPDF_MLX  Probability density function for Student's T distribution

    term<-exp(lgamma((v + 1) / 2) - lgamma(v/2))
    return(term/(sqrt(v*pi)*(1+(x**2)/v)**((v+1)/2)))
}
