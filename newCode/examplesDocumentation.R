library(catdata)

########################################################################
# Binary data - Toenail 
# library(HSAUR3)
# data(toenail)
# toe<-transform(toenail,y=as.integer(toenail$outcome=="moderate or severe"))
# saemix.data<-saemixData(name.data=toe,name.group=c("patientID"),name.predictors=c("time","y"), name.response="y",
#                        name.covariates=c("treatment"),name.X=c("time"))

library(prLogistic)
data(Toenail)
if(FALSE) {
  toenail.saemix<-Toenail
  colnames(toenail.saemix)<-c("id","y","treatment","time","visit")
  toenail.saemix<-toenail.saemix[,c(1,4,2,3,5)]
  write.table(toenail.saemix,file.path(datDir,"toenail.saemix.tab"), quote=F, row.names=F)
}
toe<-toenail.saemix

saemix.data<-saemixData(name.data=toenail.saemix,name.group=c("id"),name.predictors=c("time","y"), name.response="y",
                        name.covariates=c("treatment"),name.X=c("time"))

binary.model<-function(psi,id,xidep) {
  tim<-xidep[,1]
  y<-xidep[,2]
  inter<-psi[id,1]
  slope<-psi[id,2]
  logit<-inter+slope*tim
  pevent<-exp(logit)/(1+exp(logit))
  logpdf<-rep(0,length(tim))
  P.obs = (y==0)*(1-pevent)+(y==1)*pevent
  logpdf <- log(P.obs)
  return(logpdf)
}

saemix.model<-saemixModel(model=binary.model,description="Binary model",
                          modeltype="likelihood",
                          psi0=matrix(c(0,-.5,0,0.5),ncol=2,byrow=TRUE,dimnames=list(NULL,c("theta1","theta2"))),
                          transform.par=c(0,0), covariate.model=c(0,1),covariance.model=matrix(c(1,0,0,1),ncol=2))

saemix.options<-list(seed=1234567,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, nb.chains=10, fim=FALSE)

binary.fit<-saemix(saemix.model,saemix.data,saemix.options)

# Some graphs
barplot(table(toenail.saemix$y,toenail.saemix$visit),beside=T)

barplot(table(toenail.saemix$y,toenail.saemix$visit),beside=F)


## Binomial model in saemix for the toenail data - complete data 
## Fit only to the subjects with complete observations over the 7 visits :

# Graph of raw data
nobs<-tapply(toenail$visit,toenail$patientID,length)
isuj.nomiss<-names(nobs)[nobs==7]
toe2<-toe[toenail$patientID %in% isuj.nomiss,]

saemix.data2<-saemixData(name.data=toe2,name.group=c("patientID"),name.predictors=c("time","y"), name.response="y",
                        name.covariates=c("treatment"),name.X=c("time"))

binary.fit2<-saemix(saemix.model,saemix.data2,saemix.options)

########################################################################
# Count data - Encephalitis
if(FALSE) {
  data(encephalitis)
  encephalitis$BAV<-encephalitis$country
  encephalitis$BAV[encephalitis$BAV==2]<-0
  
  # log-linear Poisson model, with dependence on country and year
  enc1 <- glm(count ~ year+I(year^2)+BAV+year*BAV, family = poisson, data=encephalitis)
  summary(enc1)
  
  enc2 <- glm(count ~ year+BAV, family = poisson, data=encephalitis)
  summary(enc2)
}

# homi <- read.table("http://www.stat.ufl.edu/~aa/glm/data/Homicides.dat",  header = TRUE)

## Salamander data from glmmTMB
library(glmmTMB)
data(Salamanders)
summary(Salamanders)

barplot(Salamanders$DOY, Salamanders$count)

table(Salamanders$spp, Salamanders$count)

ecl<-Salamanders[Salamanders$spp=="EC-L",]
twospecies<-Salamanders[Salamanders$spp %in% c("EC-L","DF"),]
twospecies$spp <- as.character(twospecies$spp)

saemix.data<-saemixData(name.data=ecl, name.group=c("site"),
                        name.predictors=c("DOY","count"),name.response=c("count"),
                        name.covariates=c("mined","cover", "DOP", "Wtemp"),
                        units=list(x="day (scaled)",y="",covariates=c("","","","")))

countmodel.poisson<-function(psi,id,xidep) { 
  y<-xidep[,2]
  lambda<-psi[id,1]
  logp <- -lambda + y*log(lambda) - log(factorial(y))
  return(logp)
}

saemix.model<-saemixModel(model=countmodel.poisson,description="count model Poisson",modeltype="likelihood",   
                          psi0=matrix(c(0.5),ncol=1,byrow=TRUE,dimnames=list(NULL, c("lambda"))), 
                          transform.par=c(1)) #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE))

saemix.model.cov<-saemixModel(model=countmodel.poisson,description="count model Poisson",modeltype="likelihood",   
                          psi0=matrix(c(0.5),ncol=1,byrow=TRUE,dimnames=list(NULL, c("lambda"))), 
                          covariate.model = matrix(c(1,rep(0,3)), ncol=1, byrow=T),
                          transform.par=c(1)) #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)
poisson.fit<-saemix(saemix.model,saemix.data,saemix.options)
poisson.fit.cov<-saemix(saemix.model.cov,saemix.data,saemix.options)
plot(poisson.fit, plot.type="convergence")

## ZIP Poisson model
countmodel.zip<-function(psi,id,xidep) {
  y<-xidep[,2]
  lambda<-psi[id,1]
  p0<-psi[id,2]
  logp <- log(1-p0) -lambda + y*log(lambda) - log(factorial(y))
  logp0 <- log(p0+(1-p0)*exp(-lambda))
  logp[y==0]<-logp0[y==0]
  return(logp)
}

saemix.model.zip<-saemixModel(model=countmodel.zip,description="count model ZIP",modeltype="likelihood",   
                          psi0=matrix(c(0.5,0.2),ncol=2,byrow=TRUE,dimnames=list(NULL, c("lambda","p0"))), 
                          transform.par=c(1,3), #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE),
                          covariance.model=matrix(c(1,0,0,0),ncol=2,byrow=TRUE))
saemix.model.zipcov<-saemixModel(model=countmodel.zip,description="count model ZIP",modeltype="likelihood",   
                              psi0=matrix(c(0.5,0.2),ncol=2,byrow=TRUE,dimnames=list(NULL, c("lambda","p0"))), 
                              transform.par=c(1,3), #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE),
                              covariate.model = matrix(c(1,rep(0,7)), ncol=2, byrow=T),
                              covariance.model=matrix(c(1,0,0,0),ncol=2,byrow=TRUE))
zippoisson.fit<-saemix(saemix.model.zip,saemix.data,saemix.options)
zippoisson.fitcov<-saemix(saemix.model.zipcov,saemix.data,saemix.options)
plot(zippoisson.fit, plot.type="convergence")

BIC(poisson.fit)
BIC(zippoisson.fit)

if(FALSE) {
  salamander.saemix<-ecl
  salamander.saemix$mined<-as.integer(salamander.saemix$mined=="yes")
  write.table(salamander.saemix, file.path(saemixDir, "testbelhal", "salamander.saemix.tab"), quote=F, row.names=F)
  
  salamander.saemix<-Salamanders
  species<-as.character(unique(Salamanders$spp))
  for(isp in species) {
    xcov<-as.integer(salamander.saemix$spp==isp)
    salamander.saemix<-cbind(salamander.saemix, sp=xcov)
    colnames(salamander.saemix)[dim(salamander.saemix)[2]]<-isp
  }
  salamander.saemix$mined<-as.integer(salamander.saemix$mined=="yes")
  salamander.saemix$spp<-NULL
  salamander.saemix<-salamander.saemix[,c("site","DOY","count","mined","cover","DOP","Wtemp","sample",colnames(salamander.saemix)[9:15])]
  write.table(salamander.saemix, file.path(saemixDir, "testbelhal", "sal.saemix.tab"), quote=F, row.names=F)
  
  saemix.data<-saemixData(name.data=salamander.saemix, name.group=c("site"),
                          name.predictors=c("DOY","count"),name.response=c("count"),
                          name.covariates=c("mined","PR","DM","EC-A","EC-L","DES-L","DF"), # reference class will be GP
                          units=list(x="day (scaled)",y=""))
  saemix.model.zipcov<-saemixModel(model=countmodel.zip,description="count model ZIP",modeltype="likelihood",   
                                   psi0=matrix(c(0.5,0.2),ncol=2,byrow=TRUE,dimnames=list(NULL, c("lambda","p0"))), 
                                   transform.par=c(1,3), #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE),
                                   covariate.model = matrix(c(rep(1,2*length(saemix.data@name.covariates))), ncol=2, byrow=T),
                                   covariance.model=matrix(c(1,0,0,0),ncol=2,byrow=TRUE))
  zippoisson.fitcov<-try(saemix(saemix.model.zipcov,saemix.data,saemix.options)) # fails with Lapack problem => really need to debug this...
  
}

# Same, with the 2 most frequent species, using spp as a covariate - fails

saemix.data<-saemixData(name.data=twospecies, name.group=c("site"),
                        name.predictors=c("DOY","count"),name.response=c("count"),
                        name.covariates=c("spp","mined","cover", "DOP", "Wtemp"),
                        units=list(x="day (scaled)",y="",covariates=c("","","","","")))

countmodel.poisson<-function(psi,id,xidep) { 
  y<-xidep[,2]
  lambda<-psi[id,1]
  logp <- -lambda + y*log(lambda) - log(factorial(y))
  return(logp)
}

saemix.model<-saemixModel(model=countmodel.poisson,description="count model Poisson",modeltype="likelihood",   
                          psi0=matrix(c(1),ncol=1,byrow=TRUE,dimnames=list(NULL, c("lambda"))), 
                          covariate.model = matrix(c(1,rep(0,4)),ncol=1, byrow=T),
                          transform.par=c(1)) #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE, fim=FALSE)
poisson.fit<-try(saemix(saemix.model,saemix.data,saemix.options)) # fails, pb with solving Lapack

saemix.model.zip<-saemixModel(model=countmodel.zip,description="count model ZIP",modeltype="likelihood",   
                              psi0=matrix(c(0.5,0.2),ncol=2,byrow=TRUE,dimnames=list(NULL, c("lambda","p0"))), 
                              transform.par=c(1,3), #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE),
                              covariate.model = matrix(c(1,rep(0,9)),ncol=2, byrow=T),
                              covariance.model=matrix(c(1,0,0,0),ncol=2,byrow=TRUE))
zippoisson.fit<-try(saemix(saemix.model.zip,saemix.data,saemix.options))

## Fit with glmmTMB
zipm1 = glmmTMB(count~spp + mined + (1|site), Salamanders, family="poisson")
summary(zipm1)
zipm2 = glmmTMB(count~spp * mined + (1|site), Salamanders, family="poisson")
summary(zipm2)

zipm3 = glmmTMB(count~spp * mined + (1|site), zi=~spp * mined, Salamanders, family="poisson")
summary(zipm3)

### 

## nest data from glmmTMB

library(glmmTMB)
data(Owls)
Owls$time<-Owls$ArrivalTime-min(Owls$ArrivalTime)
Owls<-Owls[order(Owls$Nest, Owls$ArrivalTime),]
Owls[!duplicated(Owls$Nest),]
colnames(Owls)[5]<-"NbNego"

owls.saemix<-Owls

?# saemix
saemix.data<-saemixData(name.data=owls.saemix,header=TRUE,name.group=c("Nest"),
                        name.predictors=c("time","NbNego"),name.response=c("NbNego"),
                        name.covariates=c("FoodTreatment", "SexParent", "NegPerChick", "BroodSize"),
                        units=list(x="min",y="",covariates=c("","","","")))

#Basic Poisson model
countmodel.poisson<-function(psi,id,xidep) { 
  y<-xidep[,1]
  lambda<-psi[id,1]
  dummy<-psi[id,2]
  logp <- -lambda + y*log(lambda) - factorial(y)
  return(logp)
}

saemix.model<-saemixModel(model=countmodel.poisson,description="count model",modeltype="likelihood",   
                          psi0=matrix(c(0.5,1),ncol=2,byrow=TRUE,dimnames=list(NULL, c("lambda","dummy"))), 
                          transform.par=c(1,1), #omega.init=matrix(c(0.5,0,0,0.3),ncol=2,byrow=TRUE),
                          covariance.model=matrix(c(1,0,0,0),ncol=2,byrow=TRUE),
                          fixed.estim=c(1,0))

##Generalized Poisson model
countmodel.genpoisson<-function(psi,id,xidep) {
  y<-xidep[,1]
  delta<-psi[id,1]
  lambda<-psi[id,2]
  logp <- -lambda
  pos.ind <- which(y>0)
  logp[pos.ind] <- log(lambda) + (y-1)*log(lambda+y*delta) - (lambda+y*delta) - factorial(y)
  return(logp)
}

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)
poisson.fit<-saemix(saemix.model,saemix.data,saemix.options)
plot(poisson.fit, plot.type="convergence")

# Grouseticks Data on red grouse ticks from Elston et al. 2001 (Number of ticks on the heads of red grouse chicks sampled in the field (grouseticks) and an aggregated version (grouseticks_agg); see original source for more details)

library(lme4)
data("grouseticks")
form <- TICKS~YEAR+HEIGHT+(1|BROOD)+(1|INDEX)+(1|LOCATION)
full_mod1  <- glmer(form, family="poisson",data=grouseticks)


########################################################################
# Categorical data - Knee
data(knee)

# Transformation to long format and dichotomisation of response
# knee <- reshape(knee, direction="long", varying=list(5:8), v.names="R",timevar="Time")
knee2<-cbind(knee[,c(1:5)], time=0)
xtim<-c(0,3,7,10)
for(icol in 6:8) {
  xprov<-cbind(knee[,c(1:4,icol)], time=xtim[icol-4])
  colnames(xprov)<-colnames(knee2)
  knee2<-rbind(knee2, xprov)
}
knee2<-knee2[order(knee2$N, knee2$time),c(1,6,5,2:4)]
colnames(knee2)[3]<-"y"
knee2$RD<-as.integer(knee2$y>2)
knee2$Age <- knee2$Age - 30
knee2$treatment<-knee2$Th-1
knee2$Age2 <- knee2$Age^2

knee.saemix<-knee2[,-c(4,7)]
colnames(knee.saemix)[1]<-"id"
barplot(table(knee.saemix$y,knee.saemix$time),beside=T)

if(FALSE) 
  write.table(knee.saemix, file.path(datDir, "knee.saemix.tab"), row.names=F, quote=F)

saemix.data<-saemixData(name.data=knee.saemix,name.group=c("id"),
                        name.predictors=c("y", "time"), name.X=c("time"),
                        name.covariates = c("Age","Sex","treatment"))

ordinal.model<-function(psi,id,xidep) {
  y<-xidep[,1]
  time<-xidep[,2]
  alp1<-psi[id,1]
  alp2<-psi[id,2]
  alp3<-psi[id,3]
  alp4<-psi[id,4]
  beta<-psi[id,5]
  
  logit1<-alp1 + beta*time
  logit2<-logit1+alp2
  logit3<-logit2+alp3
  logit4<-logit3+alp4
  pge1<-exp(logit1)/(1+exp(logit1))
  pge2<-exp(logit2)/(1+exp(logit2))
  pge3<-exp(logit3)/(1+exp(logit3))
  pge4<-exp(logit4)/(1+exp(logit4))
  logpdf<-rep(0,length(y))
  P.obs = (y==1)*pge1+(y==2)*(pge2 - pge1)+(y==3)*(pge3 - pge2)+(y==4)*(pge4 - pge3)+(y==5)*(1 - pge4)
  logpdf <- log(P.obs)
  
  return(logpdf)
}
covmodel3<-covmodel2<-covmodel<-matrix(data=0,ncol=5,nrow=3)
covmodel[1:2,1]<-1
covmodel[,5]<-1
covmodel2[1,1]<-covmodel2[3,5]<-1
covmodel3[1,1]<-1

saemix.model<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
                          psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5,byrow=TRUE,dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))),
                          transform.par=c(0,1,1,1,1),omega.init=diag(rep(1,5)), covariance.model = diag(c(1,0,0,0,1)))

saemix.model.cov<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
                          psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5,byrow=TRUE,dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))),
                          transform.par=c(0,1,1,1,1),omega.init=diag(rep(1,5)), covariance.model = diag(c(1,0,0,0,1)),
                          covariate.model = covmodel)
saemix.model.cov2<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
                              psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5,byrow=TRUE,dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))),
                              transform.par=c(0,1,1,1,1),omega.init=diag(rep(1,5)), covariance.model = diag(c(1,0,0,0,1)),
                              covariate.model = covmodel2)
saemix.model.cov3<-saemixModel(model=ordinal.model,description="Ordinal categorical model",modeltype="likelihood",
                               psi0=matrix(c(0,0.2, 0.6, 3, 0.2),ncol=5,byrow=TRUE,dimnames=list(NULL,c("alp1","alp2","alp3","alp4","beta"))),
                               transform.par=c(0,1,1,1,1),omega.init=diag(rep(1,5)), covariance.model = diag(c(1,0,0,0,1)),
                               covariate.model = covmodel3)

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, nb.chains=10)

ord.fit<-saemix(saemix.model,saemix.data,saemix.options)
ord.fit.cov<-saemix(saemix.model.cov,saemix.data,saemix.options)
ord.fit.cov2<-saemix(saemix.model.cov2,saemix.data,saemix.options)
ord.fit.cov3<-saemix(saemix.model.cov3,saemix.data,saemix.options)
BIC(ord.fit)
BIC(ord.fit.cov)
BIC(ord.fit.cov2)
BIC(ord.fit.cov3)

#######################################
# Analysis of dichotomised response :-/
if(FALSE) { # install glmmML, same factors found significant
  library(glmmML)
  kneeGHQ <- glmmML(RD ~ as.factor(Th) + as.factor(Sex) + Age + Age2, data=knee2,
                    family=binomial(), method="ghq", n.points=20, cluster=N)
  summary(kneeGHQ)
}

kneePQL <- glmmPQL(RD ~ as.factor(Th) + as.factor(Sex) + Age + Age2, data=knee2,
                    random = ~ 1|N, family=binomial())
summary(kneePQL)

kneePQL2 <- glmmPQL(RD ~ as.factor(Th) + Age2, data=knee2,
                   random = ~ 1|N, family=binomial())
summary(kneePQL2)

########################################################################
# TTE - Lung cancer

weibulltte.model<-function(psi,id,xidep) {
  T<-xidep[,1]
  y<-xidep[,2] # events (1=event, 0=no event)
  cens<-which(xidep[,3]==1) # censoring times (subject specific)
  init <- which(T==0)
  lambda <- psi[id,1] # Parameters of the Weibull model
  beta <- psi[id,2]
  Nj <- length(T)
  
  ind <- setdiff(1:Nj, append(init,cens)) # indices of events
  hazard <- (beta/lambda)*(T/lambda)^(beta-1) # H'
  H <- (T/lambda)^beta # H
  logpdf <- rep(0,Nj) # ln(l(T=0))=0
  logpdf[cens] <- -H[cens] + H[cens-1] # ln(l(T=censoring time))
  logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind]) # ln(l(T=event time))
  return(logpdf)
}

# TTE - PBC
data("pbc2.id", package = "JM")

pbc<-pbc2.id
pbc$cens<-as.integer(pbc$status!="dead") # censored=1, non-censored=0
pbc$status<-as.integer(pbc$status=="dead")  # dead=1, alive=0
pbc2<-pbc
pbc2$years<-0
pbc2$status<-pbc2$cens<-0
pbc2<-rbind(pbc2, pbc)
colnames(pbc2)[colnames(pbc2)=="years"]<-"time"
pbc.saemix<-pbc2[order(pbc2$id, pbc2$time),]

saemix.data<-saemixData(name.data=pbc.saemix,header=TRUE,name.group=c("id"),
                        name.predictors=c("time","status","cens"),name.response=c("status"),
                        name.covariates=c("age", "sex", "drug", "ascites", "hepatomegaly", "spiders","alkaline"),
                        units=list(x="days",y="",covariates=c("yr","","-","","","","")))

exptte.model<-function(psi,id,xidep) {
  T<-xidep[,1]
  y<-xidep[,2] # events (1=event, 0=no event)
  cens<-which(xidep[,3]==1) # censoring times (subject specific)
  init <- which(T==0)
  beta <- psi[id,1]
  Nj <- length(T)
  
  ind <- setdiff(1:Nj, append(init,cens)) # indices of events
  hazard <- 1/beta # H'
  H <- T/beta # H
  logpdf <- rep(0,Nj) # ln(l(T=0))=0
  logpdf[cens] <- -H[cens] + H[cens-1] # ln(l(T=censoring time))
  logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind]) # ln(l(T=event time))
  return(logpdf)
}

saemix.Weibmodel<-saemixModel(model=weibulltte.model,description="time model",modeltype="likelihood",
                          psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL, c("lambda","beta"))),
                          transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE))
saemix.Expmodel<-saemixModel(model=exptte.model,description="time model",modeltype="likelihood",
                              psi0=matrix(c(2,1),ncol=2,byrow=TRUE,dimnames=list(NULL, c("beta","dummy"))),
                              transform.par=c(1,0), fixed.estim = c(1,0),
                             covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE))

saemix.WeibCov<-saemixModel(model=weibulltte.model,description="time model",modeltype="likelihood",
                          psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("lambda","beta"))),
                          covariate.model=matrix(c(1,0,1,0,1,0,rep(0,8)),ncol=2, byrow=TRUE),
                          transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE))
saemix.ExpCov<-saemixModel(model=exptte.model,description="time model",modeltype="likelihood",
                           psi0=matrix(c(2,1),ncol=2,byrow=TRUE,dimnames=list(NULL, c("beta","dummy"))),
                           transform.par=c(1,0), fixed.estim = c(1,0),
                           covariate.model=matrix(c(1,0,1,0,1,0,rep(0,8)),ncol=2, byrow=TRUE),
                           covariance.model=matrix(c(1,0,0,0),ncol=2, byrow=TRUE))

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, displayProgress=FALSE)
tte.Weibfit<-saemix(saemix.Weibmodel,saemix.data,saemix.options)
tte.Weibcov<-saemix(saemix.WeibCov,saemix.data,saemix.options)
tte.Expfit<-saemix(saemix.Expmodel,saemix.data,saemix.options)
tte.Expcov<-saemix(saemix.ExpCov,saemix.data,saemix.options)
plot(tte.fit, plot.type="convergence")

ypred<-predict(tte.fit)

if(FALSE) 
  write.table(pbc.saemix, file.path(datDir, "pbc.saemix.tab"), quote=F, row.names=F)

## Call:
## survreg(formula = Surv(years, status2) ~ drug + sex + age, data = pbc2.id)
##                  Value Std. Error     z      p
## (Intercept)    4.13862    0.48260  8.58 <2e-16
## drugD-penicil  0.13046    0.15576  0.84  0.402
## sexfemale      0.44472    0.20147  2.21  0.027
## age           -0.03874    0.00793 -4.89  1e-06
## Log(scale)    -0.10223    0.07423 -1.38  0.168

# Exponential, JM
##                 Value Std. Error     z       p
## (Intercept)    4.3282     0.5159  8.39 < 2e-16
## drugD-penicil  0.1438     0.1721  0.84    0.40
## sexfemale      0.4816     0.2217  2.17    0.03
## age           -0.0420     0.0084 -5.00 5.7e-07

########################################################################
# RTTE - Gaucher ? (few events...)

########################################################################
# RTTE - Simulated data

# For repeated time-to-events, a simulation algorithm is to simulate TTE repeatedly in an individual until the time to censoring has been reached (Penichou et al. 2014). Here the simulation for each successive event is perfomed by simulating a random variable in $\mathcal{U}[0,1]$ and using the inverse of the cumulative density function to generate the corresponding time to event.

simul.rtte.unif<-function(psi) { # xidep, id not important, we only use psi
  censoringtime <- 3
  maxevents <- 30
  lambda <- psi[,1]
  beta <- psi[,2]
  simdat<-NULL
  N<-nrow(psi)
  for(i in 1:N) {
    eventTimes<-c(0)
    T<-0
    Vj<-runif(1)
    #    T <- (-log(Vj)*lambda[i])^(beta[i])
    T<-lambda[i]*(-log(Vj))^(1/beta[i])
    nev<-0
    while (T < censoringtime & nev<maxevents){
      eventTimes <- c(eventTimes, T)  
      nev<-nev+1
      Vj<-runif(1)
      #      T <- T+(-log(Vj)*lambda[i])^(beta[i])
      #      T<-(-log(Vj)*lambda[i] + T^(1/beta[i]))^(beta[i])
      T<-lambda[i]*(-log(Vj) + (T/lambda[i])^(beta[i]))^(1/beta[i])
    }
    if(nev==maxevents) {
      message("Reached maximum number of events\n")
    }
    eventTimes<-c(eventTimes, censoringtime)
    cens<-rep(1,length(eventTimes))
    cens[1]<-cens[length(cens)]<-0
    simdat<-rbind(simdat,
                  data.frame(id=i, T=eventTimes, status=cens))
  }
  return(simdat)
}

# Subjects
set.seed(12345)
param<-c(2, 1.5, 0.5)
# param<-c(4, 1.2, 0.3)
omega<-c(0.25,0.25)
nsuj<-200
risk<-rep(0,nsuj)
risk[(nsuj/2+1):nsuj]<-1
psiM<-data.frame(lambda=param[1]*exp(rnorm(nsuj,sd=omega[1])), beta=param[2]*exp(param[3]*risk+rnorm(nsuj,sd=omega[2])))
simdat <- simul.rtte.unif(psiM)
simdat$risk<-as.integer(simdat$id>(nsuj/2))

if(FALSE) { # Check simulated parameters
  summary(psiM)
  apply(log(psiM),2,sd)
}

if(FALSE) {
  xtim<-seq(0:20)
  rtte.data<-data.frame(id=rep(1:nsuj,each=length(xtim)),time=rep(xtim,nsuj))
  preds <- simul.rtte.rweib(psiM, rtte.data$id, rtte.data[,c("time")])
  
  par(mfrow=c(1,2))
  hist(preds)
  hist(simdat$T[simdat$T>0])
  
  rtte.data$tlat<-preds
  rtte.data$status<-as.integer(rtte.data$tlat<3)
  
  # Remove duplicated censored times
  dat1<-NULL
  for(i in 1:nsuj) {
    idat<-rtte.data[rtte.data$id==i,]
    idat<-idat[!duplicated(idat$tlat),,drop=FALSE]
    dat1<-rbind(dat1, c(i,0,0), idat[,-c(2)])
  }
  rtte.data<-dat1
  table(tapply(rtte.data$id, rtte.data$id, length))
}

if(FALSE)
  write.table(simdat,file.path(ecoDir,"simulatedRTTE.csv"), quote=F, row.names=F)
saemix.data<-saemixData(name.data=simdat, name.group=c("id"), name.predictors=c("T"), name.response="status", name.covariates="risk")

### Fit RTTE
rtte.model<-function(psi,id,xidep) {
  T<-xidep[,1]
  N <- nrow(psi) # nb of subjects
  Nj <- length(T) # nb of events (including 0 and censoring times)
  # censoringtime = 6
  censoringtime = max(T) # same censoring for everyone
  lambda <- psi[id,1]
  beta <- psi[id,2]
  tinit <- which(T==0) # indices of beginning of observation period
  tcens <- which(T==censoringtime) # indices of censored events 
  tevent <- setdiff(1:Nj, append(tinit,tcens)) # indices of non-censored event times
  hazard <- (beta/lambda)*(T/lambda)^(beta-1)
  H <- (T/lambda)^beta
  logpdf <- rep(0,Nj)
  logpdf[tcens] <- -H[tcens] + H[tcens-1]
  logpdf[tevent] <- -H[tevent] + H[tevent-1] + log(hazard[tevent])
  return(logpdf)
}

saemix.model.base<-saemixModel(model=rtte.model,description="Repeated TTE model",modeltype="likelihood",
                               psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("lambda","beta"))),
                               transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, byrow=TRUE))
saemix.model<-saemixModel(model=rtte.model,description="Repeated TTE model",modeltype="likelihood",
                          psi0=matrix(c(1,2),ncol=2,byrow=TRUE,dimnames=list(NULL,  c("lambda","beta"))),
                          transform.par=c(1,1),covariate.model=matrix(c(0,1),ncol=2),
                          covariance.model=matrix(c(1,0,0,1),ncol=2, byrow=TRUE))
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE, fim=FALSE, displayProgress=FALSE)
rtte.fit<-saemix(saemix.model,saemix.data,saemix.options)
plot(rtte.fit, plot.type="convergence")
