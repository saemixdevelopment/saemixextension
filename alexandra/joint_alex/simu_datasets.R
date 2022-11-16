####### SIMU JOINT 1 mod longi + 1 TTE ############
library(survival)
## simu data

set.seed(1996)
N=1000

longi = data.frame(id=NA,obs=NA,time=NA)
b0=15
b1 = 0.3
omega_b0 = 0.5
omega_b1 = 0.1
sigma_a = 0.1
h0=0.01
alpha=0.10


rint = rnorm(N,mean=b0 ,sd=omega_b0) # vector des random intercepts
reff = rnorm(N,mean=b1 ,sd=omega_b1) # vector des random slopes
tt = 0:30	# design pour les observations longitudinales
longi = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt)))),time = rep(tt,N))
longi$obs = rint[longi$id] + reff[longi$id]*longi$time + rnorm(nrow(longi),sd=sigma_a) # on remplit avec les observations

gettimes = function(r1,r2){ # Fonction qui renvoit un temps d'évènement pour un patient avec intercept r1 et slope r2
  x = rexp(1)
  H = function(tt) h0*exp(alpha*r1)/alpha/r2*(exp(alpha*r2*tt)-1) - x
  z = try(uniroot(H,interval=c(0,30))$root,silent=T)
  ifelse(class(z)=="try-error",Inf,z)
}
tte.data = data.frame(id=1:N,time = mapply(gettimes,rint,reff), ev = 1) # dataset des données de survie
tte.data$ev[tte.data$time>30]=0	# Les temps > 30 seront censurées
tte.data$time[tte.data$time>30]=30	# Date de censure pour les temps >30
longi = longi[longi$time<=sapply(longi$id,function(x) tte.data$time[tte.data$id==x]),] # on supprime les observations postérieures au temps d'évènement

table(tte.data$ev)
summary(tte.data$time)
plot(survfit(Surv(time,ev)~1, data = tte.data))

tte.data2<-data.frame(id=tte.data$id, time=0, ev=0)
ttedata<-rbind(tte.data, tte.data2)
ttedata <- ttedata[order(ttedata$id, ttedata$time),]

names(ttedata)[3]="obs"

longi$ytype=1
ttedata$ytype=2

alldata = rbind(longi,ttedata)

#write.csv(alldata,file="C:/Users/AlexandraLAVALLEY/Documents/Code_saemix/essais/joint_tte2.csv",row.names = F)





###### SIMU 2 COMPETING EVENTS ########
library(survival)
library(cmprsk)

N=1000

res<-data.frame(id=rep(NA,N), temps=rep(NA,N), ev = rep(NA,N))


T = 1000

p=0.3

xx = seq(0,T,length.out=10000)
h = xx[2]-xx[1]

# on simule age en 4 categories avec les patterns de l'application SOFA
x=c(1,2,3) # 3 cat d'age
prob=c(1565/4050,1333/4050,1152/4050)
sim_age=data.frame(id=1:N,age=sample(x,N,replace = T,prob))
sim_age$coeff_1 = ifelse(sim_age$age==1,0,ifelse(sim_age$age==2,0.3,0.5))
sim_age$coeff_2 = ifelse(sim_age$age==1,0,ifelse(sim_age$age==2,-0.2,-0.3))

for(i in 1:N){
  #h1 = function(t) (p*exp(-t)*exp(sim_age$coeff_1[i]))/(1-p*(1-exp(-t))) on spécifie directement F1 dans la simu
  #F1 = 1-exp(-cumsum(Vectorize(h1)(xx))*h)
  F1 = 1-(1-p*(1-exp(-xx*h)))^(exp(sim_age$coeff_1[i]))
  F2 = (1-rev(F1)[1])*(1-exp(-xx*h))^(exp(sim_age$coeff_2[i]))
  
  p1 = rev(F1)[1]
  p2 = 1-p1
  
  u = runif(1)
  event = ifelse(u<p1,1,2)
  
  if(event==1) (t=(which(F1/p1>runif(1))[1])*h) #ceiling 
  if(event==2) (t=(which(F2/p2>runif(1))[1])*h) #ceiling 
  
  res[i,1]<-i
  res[i,2]<-round(t,2)
  res[i,3]<-event
}

plot(cuminc(res$temps,res$ev,cencode=3),xlim=c(0,30))
table(res$ev)
summary(res$temps)

d=30
data_comp<-data.frame(id=rep(res$id,each=4),time=0,obs=0,obs_id=rep(c(1,2),each=2,times=length(res$id)),age=NA)
for (i in res$id){
  if (res$temps[res$id==i]<=d & res$ev[res$id==i]==1){
    data_comp$time[data_comp$id==i & data_comp$obs_id==1][2]<-res$temps[res$id==i]
    data_comp$obs[data_comp$id==i & data_comp$obs_id==1][2]<-1
    data_comp$time[data_comp$id==i & data_comp$obs_id==2][2]<-d
  }
  if (res$temps[which(res$id==i)]<=d & res$ev[which(res$id==i)]==2){
    data_comp$time[data_comp$id==i & data_comp$obs_id==2][2]<-res$temps[res$id==i]
    data_comp$obs[data_comp$id==i & data_comp$obs_id==2][2]<-1
    data_comp$time[data_comp$id==i & data_comp$obs_id==1][2]<-d
  }
  if (res$temps[which(res$id==i)]>d){
    data_comp$time[data_comp$id==i][2]<-d
    data_comp$time[data_comp$id==i][4]<-d
  }
  data_comp$age[data_comp$id==i] = sim_age$age[sim_age$id==i]
}

write.csv(data_comp,file="C:/Users/AlexandraLAVALLEY/Documents/Code_saemix/essais/data_2comp.csv",row.names = F)


##### SIMU JOINT 1 mod longi + 2 COMPETING EVENTS #####

b0.pop<-15
b1.pop<- -0.1
alpha1<- 0.15
alpha2=0

omega_b0<-3
omega_b1<-0.1

sigma_a<-0.5

p1=0.04
g1 = 0.07

N=1000

rint = rnorm(N,mean=b0.pop ,sd=omega_b0) # vector des random intercepts
reff = rnorm(N,mean=b1.pop ,sd=omega_b1) # vector des random slopes
tt = 0:30	# design pour les observations longitudinales
longi = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt)))),time = rep(tt,N))
longi$obs = rint[longi$id] + reff[longi$id]*longi$time + rnorm(nrow(longi),sd=sigma_a) # on remplit avec les observations

gettimes = function(r1,r2){ # Fonction qui renvoit un temps d'évènement pour un patient avec intercept r1 et slope r2
  
  t = seq(0,1000,length.out=10000)
  pas = t[2]-t[1]
  h = (p1*g1*exp(-g1*t)/(1-p1*(1-exp(-g1*t))))*exp(alpha1*(r1+r2*t))
  F1 = 1-exp(-cumsum(h)*pas)
  F2 = (1-rev(F1)[1])*(1-exp(-t*pas))
  
  p1 = rev(F1)[1]
  p2 = 1-p1
  
  u = runif(1)
  event = ifelse(u<p1,1,2)
  
  if(event==1) z=(which(F1/p1>runif(1))[1])*pas 
  if(event==2) z=(which(F2/p2>runif(1))[1])*pas 
  c(z,event)
}

a=mapply(gettimes,rint,reff)

tte.data = data.frame(id=1:N,time = a[1,], obs = a[2,]) # dataset des données de survie
tte.data$obs[tte.data$time>30]=0	# Les temps > 30 seront censurées
tte.data$time[tte.data$time>30]=30	# Date de censure pour les temps >30
table(tte.data$obs)
summary(tte.data$time)
plot(cuminc(tte.data$time,tte.data$obs,cencode = 0))

longi = longi[longi$time<=sapply(longi$id,function(x) tte.data$time[tte.data$id==x]),] # on supprime les observations postérieures au temps d'évènement

DATA_SURV<-data.frame(id=rep(tte.data$id,each=4),time=0,obs=0,ytype=rep(c(2,3),each=2,times=length(tte.data$id)))
for (i in tte.data$id){
  if (tte.data$time[tte.data$id==i]<=30 & tte.data$obs[tte.data$id==i]==1){
    DATA_SURV$time[DATA_SURV$id==i & DATA_SURV$ytype==2][2]<-tte.data$time[tte.data$id==i]
    DATA_SURV$obs[DATA_SURV$id==i & DATA_SURV$ytype==2][2]<-1
    DATA_SURV$time[DATA_SURV$id==i & DATA_SURV$ytype==3][2]<-30
  }
  if (tte.data$time[which(tte.data$id==i)]<=30 & tte.data$obs[which(tte.data$id==i)]==2){
    DATA_SURV$time[DATA_SURV$id==i & DATA_SURV$ytype==3][2]<-tte.data$time[tte.data$id==i] 
    DATA_SURV$obs[DATA_SURV$id==i & DATA_SURV$ytype==3][2]<-1
    DATA_SURV$time[DATA_SURV$id==i & DATA_SURV$ytype==2][2]<-30
  }
  if (tte.data$obs[which(tte.data$id==i)]==0){
    DATA_SURV$time[DATA_SURV$id==i][2]<-30
    DATA_SURV$time[DATA_SURV$id==i][4]<-30
  }
}


longi$ytype=1

alldata = rbind(longi,DATA_SURV)

write.csv(alldata, "C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension/alexandra/joint_alex/datas/joint_comp.csv",row.names = F)





##### SIMU JOINT 3 mod longi + 1 TTE #####
N=1000

longi1 = data.frame(id=NA,obs=NA,time=NA)
b01=15
b11 = 0.3
omega_b01 = 0.5
omega_b11 = 0.1
sigma_a1 = 0.5
alpha1=0.10

longi2 = data.frame(id=NA,obs=NA,time=NA)
b02 = 7
b12 = -0.1
omega_b02 = 1
omega_b12 = 0.1
sigma_a2 = 0.1
alpha2 = -0.20

longi3 = data.frame(id=NA,obs=NA,time=NA)
b03 = 30
b13 = 0.8
omega_b03 = 4
omega_b13 = 0.6
sigma_a3 = 1
alpha3 = 0.15

h0 = 0.00005

rint1 = rnorm(N,mean=b01 ,sd=omega_b01) # vector des random intercepts
reff1 = rnorm(N,mean=b11 ,sd=omega_b11) # vector des random slopes
tt = 0:30	# design pour les observations longitudinales
longi1 = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt)))),time = rep(tt,N))
longi1$obs = rint1[longi1$id] + reff1[longi1$id]*longi1$time + rnorm(nrow(longi1),sd=sigma_a1) 

rint2 = rnorm(N,mean=b02 ,sd=omega_b02) # vector des random intercepts
reff2 = rnorm(N,mean=b12 ,sd=omega_b12) # vector des random slopes
longi2 = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt)))),time = rep(tt,N))
longi2$obs = rint2[longi2$id] + reff2[longi2$id]*longi2$time + rnorm(nrow(longi2),sd=sigma_a2)

rint3 = rnorm(N,mean=b03 ,sd=omega_b03) # vector des random intercepts
reff3 = rnorm(N,mean=b13 ,sd=omega_b13) # vector des random slopes
longi3 = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt)))),time = rep(tt,N))
longi3$obs = rint3[longi3$id] + reff3[longi3$id]*longi3$time + rnorm(nrow(longi3),sd=sigma_a3)

gettimes = function(r1,r2,r3,r4,r5,r6){ # Fonction qui renvoit un temps d'évènement pour un patient avec intercept r1 et slope r2
  x = rexp(1)
  t = seq(0,1000,length.out=10000)
  pas = t[2]-t[1]
  h = h0*exp(alpha1*(r1+r2*t)+alpha2*(r3+r4*t)+alpha3*(r5+r6*t))
  H = cumsum(h)*pas
  z=(which(H>x)[1])*pas
  z=ifelse(is.na(z)==T,Inf,z)
  z
}
tte.data = data.frame(id=1:N,time = mapply(gettimes,rint1,reff1,rint2,reff2,rint3,reff3), obs = 1) # dataset des données de survie
tte.data$obs[tte.data$time>30]=0	# Les temps > 30 seront censurées
tte.data$time[tte.data$time>30]=30	# Date de censure pour les temps >30

longi1 = longi1[longi1$time<=sapply(longi1$id,function(x) tte.data$time[tte.data$id==x]),] # on supprime les observations postérieures au temps d'évènement
longi2 = longi2[longi2$time<=sapply(longi2$id,function(x) tte.data$time[tte.data$id==x]),] # on supprime les observations postérieures au temps d'évènement
longi3 = longi3[longi3$time<=sapply(longi3$id,function(x) tte.data$time[tte.data$id==x]),] # on supprime les observations postérieures au temps d'évènement

table(tte.data$obs)
summary(tte.data$time)
plot(survfit(Surv(time,ev)~1, data = tte.data))

tte.data2<-data.frame(id=tte.data$id, time=0, obs=0)
ttedata<-rbind(tte.data, tte.data2)
ttedata <- ttedata[order(ttedata$id, ttedata$time),]

longi1$ytype=1
longi2$ytype=2
longi3$ytype=3
ttedata$ytype=4

alldata = rbind(longi1,longi2,longi3,ttedata)

write.csv(alldata,file="C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension/alexandra/joint_alex/datas/joint_multilongi_1tte.csv",row.names = F)
