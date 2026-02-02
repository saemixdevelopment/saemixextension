# Table 8.6, Agresti 1996
party <- factor(rep(c("Rep","Dem"), c(407, 428)), 
                levels=c("Rep","Dem"))  
rpi <- c(30, 46, 148, 84, 99) # cell counts
dpi <- c(80, 81, 171, 41, 55) # cell counts
ideology <- c("Very Liberal","Slightly Liberal","Moderate","Slightly Conservative","Very Conservative")
pol.ideology <- factor(c(rep(ideology, rpi), 
                         rep(ideology, dpi)), levels = ideology)
dat <- data.frame(party,pol.ideology)

# fit proportional odds model
library(MASS)
pom <- polr(pol.ideology ~ party, data=dat)


# fit multinomial models
library(nnet)
mlm <- multinom(pol.ideology ~ party, data=dat)

M1 <- logLik(pom)
M2 <- logLik(mlm)
(G <- -2*(M1[1] - M2[1]))
pchisq(G,3,lower.tail = FALSE)


# Removed from .Rmd on categorical data, previous attempt at covariate modelling.

# We can then fit different models. Here we considered covariates on the two parameters with interindividual variability, first testing all covariates then reducing to a model with Age2 influencing $\theta_1$ and treatment affecting the slope $\beta$, which had the lowest BICc, as shown using the *compare.saemix()* function.

#```{r kneeCovariateModel}
# Fitting
covmodel3<-covmodel2<-covmodel1<-matrix(data=0,ncol=5,nrow=4)
covmodel1[,1]<-1
covmodel1[,5]<-1
covmodel2[3,5]<-covmodel2[4,1]<-1
covmodel3[,1]<-1

saemix.model.cov1<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",simulate.function=simulateOrdinal,
                               psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5,byrow=TRUE,dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))),
                               transform.par=c(0,1,1,1,1),omega.init=diag(rep(1,5)), covariance.model = diag(c(1,0,0,0,1)),
                               covariate.model = covmodel1, verbose=FALSE)
saemix.model.cov2<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",simulate.function=simulateOrdinal,
                               psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5,byrow=TRUE,dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))),
                               transform.par=c(0,1,1,1,1),omega.init=diag(rep(1,5)), covariance.model = diag(c(1,0,0,0,1)),
                               covariate.model = covmodel2, verbose=FALSE)

ord.fit.cov1<-saemix(saemix.model.cov1,ordknee.data,saemix.options)
ord.fit.cov2<-saemix(saemix.model.cov2,ordknee.data,saemix.options)
BIC(ord.fit)
BIC(ord.fit.cov1)
BIC(ord.fit.cov2)
summary(ord.fit.cov2)

# Comparing the 3 covariate models - model with Age2 on alp1 and treatment on beta best, but not as good as model from the stepwise algorithm
compare.saemix(ord.fit, ord.fit.cov1, ord.fit.cov2, ord.fit.cov)
#compare.saemix(ord.fit.cov1, ord.fit.cov2, ord.fit.cov)

# Model fit in Tutz et al. (book Regression data): unclear (book not available and corresponding sections not available online)
# the models in the vignette don't seem to be longitudinal, more like different analyses of dichotomised data (Time never explicitely taken into account)
# one vignette ()
# similar model as vignette ordinal-knee1.pdf but this analysis only fits R4 (last visit and not the longitudinal profiles)
if(FALSE) {
  saemix.model.cst<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
                                simulate.function=simulateOrdinal, psi0=matrix(c(0,0.2, 0.6, 3, 0),ncol=5, byrow=TRUE, 
                                                                               dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))), transform.par=c(0,1,1,1,1), fixed.estim = c(rep(1,4),0),
                                covariate.model=covmodel3,
                                omega.init=diag(c(100, 1, 1, 1, 1)), covariance.model = diag(c(1,0,0,0,0)), verbose=FALSE)
  
  # Fitting
  saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, fim=FALSE, nb.chains=10, nbiter.saemix=c(600,100), print=FALSE)
  #saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, nb.chains=10, fim=FALSE)
  
  ord.fit.cst<-saemix(saemix.model.cst,ordknee.data,saemix.options)
  summary(ord.fit.cst)
  # estimates not that close to polr3 (caution: signs of covariates are inverted in polr) or vglm in VGAM (both consistent and both different from our results)
  # but these only fit R4 (last visit and not the longitudinal profiles) so not sure the comparison is valid anyway
}
#```