###########################	Computational fcts	#############################

#' @rdname saemix.internal
#' 
#' @aliases compute.Uy compute.LLy 
#' 
#' @keywords internal

compute.LLy.multi<-function(phiM,args,Dargs,DYF,pres) {
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-Dargs$structural.model(psiM,Dargs$IdM,Dargs$XM)
  ytype = Dargs[["XM"]][["ytype"]]
  for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
  gpred<-error(fpred,pres,Dargs$XM$ytype)
  for(itype in 1:length(Dargs$modeltype)) {
    if (Dargs$modeltype[itype]=="structural"){
      DYF[args$ind.ioM][ytype==itype]<-0.5*((Dargs$yM[ytype==itype]-fpred[ytype==itype])/gpred[ytype==itype])**2+log(gpred[ytype==itype])
    } else {
      DYF[args$ind.ioM][ytype==itype]<- -fpred[ytype==itype]
    }
  }
  U<-colSums(DYF)
  return(U)
}

# juste pour le calcul des SE, même fonctions mais moyennée sur les phiM 
compute.LLy.multi2<-function(phiM,args,Dargs,DYF,pres) {
  DYF = DYF[,1:Dargs$N]
  args$ind.ioM =  args$ind.ioM[1:length(Dargs$yobs)]
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-Dargs$structural.model(psiM,Dargs$IdM[1:length(Dargs$yobs)],Dargs$XM[1:length(Dargs$yobs),])
  ytype = Dargs[["XM"]][["ytype"]][1:length(Dargs$yobs)]
  for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
  gpred<-error(fpred,pres,Dargs$XM$ytype[1:length(Dargs$yobs)])
  for(itype in 1:length(Dargs$modeltype)) {
    if (Dargs$modeltype[itype]=="structural"){
      DYF[args$ind.ioM][ytype==itype]<-0.5*((Dargs$yM[1:length(Dargs$yobs)][ytype==itype]-fpred[ytype==itype])/gpred[ytype==itype])**2+log(gpred[ytype==itype])
    } else {
      DYF[args$ind.ioM][ytype==itype]<- -fpred[ytype==itype]
    }
  }
  U<-colSums(DYF)
  return(U)
}

compute.se_res_add<-function(phiM,args,Dargs,DYF,pres) {
  sigma = pres[1]
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-Dargs$structural.model(psiM,Dargs$IdM[1:length(Dargs$yobs)],Dargs$XM[1:length(Dargs$yobs),])
  ni = sapply(unique(Dargs$IdM[1:length(Dargs$yobs)]),function(i) length(which(Dargs$IdM[1:length(Dargs$yobs)]==i))) - length(Dargs$modeltype[Dargs$modeltype=="likelihood"])
  err = (Dargs$yobs-fpred)^2
  err2 = sapply(unique(Dargs$IdM[1:length(Dargs$yobs)]), function(i) sum(err[Dargs$IdM[1:length(Dargs$yobs)]==i & Dargs$XM$ytype[1:length(Dargs$yobs)]==1]/(2*sigma**2)))
  f = -ni*log(sigma) - err2
  return(f)
}

compute.LLy.multi_selog<-function(phiM,args,Dargs,DYF,pres,j,coef,l,d) {
  DYF = DYF[,1:Dargs$N]
  args$ind.ioM =  args$ind.ioM[1:length(Dargs$yobs)]
  psiM<-transphi(phiM,Dargs$transform.par)
  psiM[,j] =  psiM[,j]+coef[l]*d
  fpred<-Dargs$structural.model(psiM,Dargs$IdM[1:length(Dargs$yobs)],Dargs$XM[1:length(Dargs$yobs),])
  ytype = Dargs[["XM"]][["ytype"]][1:length(Dargs$yobs)]
  for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
  gpred<-error(fpred,pres,Dargs$XM$ytype[1:length(Dargs$yobs)])
  for(itype in 1:length(Dargs$modeltype)) {
    if (Dargs$modeltype[itype]=="structural"){
      DYF[args$ind.ioM][ytype==itype]<-0.5*((Dargs$yM[1:length(Dargs$yobs)][ytype==itype]-fpred[ytype==itype])/gpred[ytype==itype])**2+log(gpred[ytype==itype])
    } else {
      DYF[args$ind.ioM][ytype==itype]<- -fpred[ytype==itype]
    }
  }
  U<-colSums(DYF)
  return(U)
}

compute.Uy.multi<-function(b0,phiM,pres,args,Dargs,DYF) {
  # Attention, DYF variable locale non modifiee en dehors
  ytype = Dargs[["XM"]][["ytype"]]
  args$MCOV0[args$j0.covariate]<-b0
  phi0<-args$COV0 %*% args$MCOV0
  phiM[,args$i0.omega2]<-do.call(rbind,rep(list(phi0),args$nchains))
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-Dargs$structural.model(psiM,Dargs$IdM,Dargs$XM)
  for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
  gpred<-error(fpred,pres,Dargs$XM$ytype)
  for(itype in 1:length(Dargs$modeltype)) {
    if (Dargs$modeltype[itype]=="structural"){
      DYF[args$ind.ioM][ytype==itype]<-0.5*((Dargs$yM[ytype==itype]-fpred[ytype==itype])/gpred[ytype==itype])**2+log(gpred[ytype==itype])
    } else {
      DYF[args$ind.ioM][ytype==itype]<- -fpred[ytype==itype]
    }
  }
  U<-sum(DYF)
  return(U)
}



compute.Uy.multi.lasso<-function(b0,phiM,pres,args,Dargs,DYF) {
  # Attention, DYF variable locale non modifiee en dehors
  lambda = 0.3
  ytype = Dargs[["XM"]][["ytype"]]
  args$MCOV0[args$j0.covariate]<-b0
  phi0<-args$COV0 %*% args$MCOV0
  phiM[,args$i0.omega2]<-do.call(rbind,rep(list(phi0),args$nchains))
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-Dargs$structural.model(psiM,Dargs$IdM,Dargs$XM)
  for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))
  gpred<-error(fpred,pres,Dargs$XM$ytype)
  for(itype in 1:length(Dargs$modeltype)) {
    if (Dargs$modeltype[itype]=="structural"){
      DYF[args$ind.ioM][ytype==itype]<-0.5*((Dargs$yM[ytype==itype]-fpred[ytype==itype])/gpred[ytype==itype])**2+log(gpred[ytype==itype])
    } else {
      DYF[args$ind.ioM][ytype==itype]<- -fpred[ytype==itype]
    }
  }
  U<-sum(DYF)+lambda*sum(abs(psiM[,8,9,10]))
  return(U)
}

###########################	Computational fcts	#############################

#' @rdname saemix.internal
#' 
#' @aliases cutoff cutoff.eps cutoff.max cutoff.res
#' @aliases normcdf norminv
#' @aliases error error.typ ssq
#' @aliases transpsi transphi derivphi dtransphi
#' @aliases trnd.mlx gammarnd.mlx tpdf.mlx
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
  g<-cutoff(abs(ab[1]+ab[2]*f))
  #g<-cutoff(sqrt(ab[1]^2+ab[2]^2*f^2))  # Johannes 02/21
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

###########################	Functions for npde	#############################

kurtosis<-function (x) 
{
  #from Snedecor and Cochran, p 80
  x<-x[!is.na(x)]
  m4<-sum((x - mean(x))^4)
  m2<-sum((x - mean(x))^2)
  kurt<-m4*length(x)/(m2**2)-3
  return(kurtosis=kurt)
}
skewness<-function (x) 
{
  #from Snedecor and Cochran, p 79
  x<-x[!is.na(x)]
  m3<-sum((x - mean(x))^3)
  m2<-sum((x - mean(x))^2)
  skew<-m3/(m2*sqrt(m2/length(x)))
  return(skewness=skew)
}

#' Tests for normalised prediction distribution errors
#' 
#' Performs tests for the normalised prediction distribution errors returned by
#' \code{npde}
#' 
#' Given a vector of normalised prediction distribution errors (npde), this
#' function compares the npde to the standardised normal distribution N(0,1)
#' using a Wilcoxon test of the mean, a Fisher test of the variance, and a
#' Shapiro-Wilks test for normality. A global test is also reported.
#' 
#' The helper functions \code{kurtosis} and \code{skewness} are called to
#' compute the kurtosis and skewness of the distribution of the npde.
#' 
#' @aliases testnpde kurtosis skewness
#' @param npde the vector of prediction distribution errors
#' @return a list containing 4 components: \item{Wilcoxon test of
#' mean=0}{compares the mean of the npde to 0 using a Wilcoxon test}
#' \item{variance test }{compares the variance of the npde to 1 using a Fisher
#' test} \item{SW test of normality}{compares the npde to the normal
#' distribution using a Shapiro-Wilks test} \item{global test }{an adjusted
#' p-value corresponding to the minimum of the 3 previous p-values multiplied
#' by the number of tests (3), or 1 if this p-value is larger than 1.}
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>
#' @seealso \code{\link{saemix}}, \code{\link{saemix.plot.npde}}
#' @references K. Brendel, E. Comets, C. Laffont, C. Laveille, and F. Mentr\'e.
#' Metrics for external model evaluation with an application to the population
#' pharmacokinetics of gliclazide. \emph{Pharmaceutical Research}, 23:2036--49,
#' 2006.
#' @keywords models
#' @export testnpde
testnpde<-function(npde) 
{
  cat("---------------------------------------------\n")
  cat("Distribution of npde:\n")
  sev<-var(npde)*sqrt(2/(length(npde)-1))
  sem<-sd(npde)/sqrt(length(npde))
  cat("           mean=",format(mean(npde),digits=4),"  (SE=",format(sem,digits=2),")\n")
  cat("       variance=",format(var(npde),digits=4),"  (SE=",format(sev,digits=2),")\n")
  cat("       skewness=",format(skewness(npde),digits=4),"\n")
  cat("       kurtosis=",format(kurtosis(npde),digits=4),"\n")
  cat("---------------------------------------------\n\n")
  myres<-rep(0,4)
  y<-wilcox.test(npde)
  myres[1]<-y$p.val
  y<-shapiro.test(npde)
  myres[3]<-y$p.val
  
  # test de variance pour 1 ?chantillon
  # chi=s2*(n-1)/sigma0 et test de H0={s=sigma0} vs chi2 ? n-1 df
  semp<-sd(npde)
  n1<-length(npde)
  chi<-(semp**2)*(n1-1)
  y<-2*min(pchisq(chi,n1-1),1-pchisq(chi,n1-1))
  myres[2]<-y
  xcal<-3*min(myres[1:3])
  myres[4]<-min(1,xcal)
  names(myres)<-c("  Wilcoxon signed rank test ","  Fisher variance test      ",
                  "  SW test of normality      ","Global adjusted p-value     ")
  cat("Statistical tests\n")
  for(i in 1:4) {
    cat(names(myres)[i],": ")
    #if (myres[i]<1) 
    cat(format(myres[i],digits=3)) 
    #else cat(myres[i])
    if(as.numeric(myres[i])<0.1 & as.numeric(myres[i])>=0.05) cat(" .")
    if(as.numeric(myres[i])<0.05) cat(" *")
    if(as.numeric(myres[i])<0.01) cat("*")
    if(as.numeric(myres[i])<0.001) cat("*")
    cat("\n")
  }
  cat("---\n")
  cat("Signif. codes: '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 \n")
  cat("---------------------------------------------\n")
  return(myres)
}

