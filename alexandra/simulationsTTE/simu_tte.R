### simulation single event selon un modèle de Weibull avec une covariable en 2 modalités 
# modele de la forme : (k/l)*(t/l)^(k-1)*exp(beta)

# paramètres utilisés pour la simulation 
N = 1000 # nb de sujets
t.max = 30 # censure administrative

## covariable 
x=c(1,2) # 2 cat d'age
prob=c(1/2,1/2)
sim_age=data.frame(id=1:N,age=sample(x,N,replace = T,prob))
sim_age$coeff = ifelse(sim_age$age==1,0,0.7)

k=1.5
l=50
DC=c()
delta=rep(1,N)
for(i in 1:N){
  h<-function(time){
    (k/l)*((time/l)^(k-1))*exp(sim_age$coeff[sim_age$id==i])
  }
  
  H<-function(time){
    vec<-rep(0, length(time))
    for (j in 1:length(time)) {
      int=integrate(Vectorize(h), lower=0, upper=time[j])$value
      vec[j]<-ifelse(int<100000000, int,100000000)
    }
    vec
  }
  
  #fdr de f
  fdr<-function(time) {
    vec<-rep(0, length(time))
    for (i in 1:length(time)) {
      vec[i]<-1-exp(-H(time[i]))}
    vec
  }
  
  #On inverse la fdr ()
  fdr.inv <- function(y){
    vec<-rep(0, length(y))
    for (i in 1:length(y)) {
      vec[i]<-uniroot(function(x){fdr(x)-y[i]},interval=c(0,10000), tol=0.1)$root
    }
    vec
  }
  
  u=runif(1,0,1)
  DCtemp=fdr.inv(u)
  if(DCtemp>30){
    DCtemp=30
    delta[i]=0
  }
  
  DC[i]=DCtemp
}
tte.data = data.frame(id=1:N,time = DC, ev = delta, age = sim_age$age)

table(delta)
plot(survfit(Surv(DC,delta)~1))

# mettre la table au bon format : on rajoute une ligne par individu t=0 ev=0
tte.data2<-data.frame(id=tte.data$id, time=0, ev=0, age = tte.data$age)
tte.data<-rbind(tte.data, tte.data2)
tte.data <- tte.data[order(tte.data$id, tte.data$time),]
tte.data$age = as.factor(tte.data$age)

################################################################################################################################
##########################################      SAEMIX      ####################################################################
################################################################################################################################

TTEmodel<-function(psi,id,xidep) {
  T<-xidep[,1] # time of the event
  idevent<-xidep[,2]  # events (0=no event (censored), 1=event: death)
  cens<-which(T==max(T))  # censoring time=30
  init <- which(T==0)
  k <- psi[id,1]
  l <- psi[id,2]
  beta <- psi[id,3]
  Nj <- length(T)
  ind <- setdiff(1:Nj, append(init,cens)) # indices of events
  haz <- (k/l)*((T/l)^(k-1))*exp(beta)
  H <- ((T/l)^k)*exp(beta)
  
  logpdf <- rep(0,Nj)
  logpdf[cens] <- -H[cens] + H[cens-1]
  logpdf[ind] <- -H[ind] + H[ind-1] + log(haz[ind])
  return(logpdf)
}

tte.data = read.table("C:/Users/AlexandraLAVALLEY/Documents/Code_saemix/essais/simus_tte/table_weib1.txt",header = T)

saemix.data<-saemixData(name.data=tte.data, name.group=c("id"), name.predictors=c("time","ev"), name.covariates = c("age"), name.response="ev")

saemix.model<-saemixModel(model=TTEmodel,description="TTE model with baseline Weibull risk and one covariate",modeltype="likelihood",
                          psi0=matrix(c(5, 40, 0),ncol=3,byrow=TRUE,dimnames=list(NULL, c("k","l","beta"))),
                          transform.par=c(1,1,0),covariate.model = matrix(c(0,0,1),ncol = 3, byrow = T), covariance.model=matrix(c(1,0,0,0,0,0,0,0,0),ncol=3, byrow=TRUE),fixed.estim = c(1,1,0))

saemix.options<-list(seed=12345,save=FALSE,save.graphs=FALSE, fim=FALSE, displayProgress=FALSE)

