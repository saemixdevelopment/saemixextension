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

###########################	Test for npde	#############################
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
