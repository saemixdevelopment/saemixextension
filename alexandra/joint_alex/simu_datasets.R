####### SIMU 1 mod lin pour calcul SE STOCHASTIQUES ############
library(survival)
## simu data

set.seed(1996)
N=100

### erreur additive 
longi = data.frame(id=NA,obs=NA,time=NA)
b0=15
b1 = 0.3
omega_b0 = 3
omega_b1 = 0.2
sigma_a = 2

rint = rnorm(N,mean=b0 ,sd=omega_b0) # vector des random intercepts
reff = rnorm(N,mean=b1 ,sd=omega_b1) # vector des random slopes
tt = 0:10	# design pour les observations longitudinales
longi = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt)))),time = rep(tt,N))
longi$obs = rint[longi$id] + reff[longi$id]*longi$time + rnorm(nrow(longi),sd=sigma_a) # on remplit avec les observations

write.csv(longi,file="C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension/alexandra/joint_alex/datas/lin_mixed.csv",row.names = F)

## erreur proportionnelle 

longi = data.frame(id=NA,obs=NA,time=NA)
b0=15
b1 = 0.3
omega_b0 = 0.5
omega_b1 = 0.1
sigma_b = 0.3
h0=0.01
alpha=0.10


rint = rnorm(N,mean=b0 ,sd=omega_b0) # vector des random intercepts
reff = rnorm(N,mean=b1 ,sd=omega_b1) # vector des random slopes
tt = 0:30	# design pour les observations longitudinales
longi = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt)))),time = rep(tt,N))
longi$obs = rint[longi$id] + reff[longi$id]*longi$time + (rint[longi$id] + reff[longi$id]*longi$time)*rnorm(nrow(longi),sd=sigma_b) # on remplit avec les observations

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



####### SIMU 1 mod lin avec un param sans iiv pour calcul SE #######
set.seed(1996)
N=100

### erreur additive 
longi = data.frame(id=NA,obs=NA,time=NA)
b0=15
b1 = 0.3
omega_b0 = 0
omega_b1 = 0.9
sigma_a = 2

rint = rnorm(N,mean=b0 ,sd=omega_b0) # vector des random intercepts
reff = rnorm(N,mean=b1 ,sd=omega_b1) # vector des random slopes
tt = 0:10	# design pour les observations longitudinales
longi = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt)))),time = rep(tt,N))
longi$obs = rint[longi$id] + reff[longi$id]*longi$time + rnorm(nrow(longi),sd=sigma_a) # on remplit avec les observations

write.csv(longi,file="C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension/alexandra/joint_alex/datas/lin_mixed_without_var.csv",row.names = F)


####### SIMU JOINT 1 mod longi + 1 TTE ############
library(survival)
## simu data

set.seed(1996)
N=100

### erreur additive 
longi = data.frame(id=NA,obs=NA,time=NA)
b0=15
b1 = 0.3
omega_b0 = 0.5
omega_b1 = 0.1
sigma_a = 1
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
tte.data = data.frame(id=1:N,time = mapply(gettimes,rint,reff), obs = 1) # dataset des données de survie
tte.data$obs[tte.data$time>30]=0	# Les temps > 30 seront censurées
tte.data$time[tte.data$time>30]=30	# Date de censure pour les temps >30
longi = longi[longi$time<=sapply(longi$id,function(x) tte.data$time[tte.data$id==x]),] # on supprime les observations postérieures au temps d'évènement

table(tte.data$obs)
summary(tte.data$time)
plot(survfit(Surv(time,ev)~1, data = tte.data))

tte.data2<-data.frame(id=tte.data$id, time=0, obs=0)
ttedata2<-rbind(tte.data, tte.data2)
ttedata2 <- ttedata2[order(ttedata2$id, ttedata2$time),]

longi$ytype=1
tte.data$ytype=2
ttedata2$ytype=2

# base saemix 
alldata_saem = rbind(longi,tte.data)
alldata_mono = rbind(longi,ttedata2)

write.csv(alldata_saem,file="C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension/alexandra/joint_alex/datas/joint_tte.csv",row.names = F)
write.csv(alldata_mono,file="C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension/alexandra/joint_alex/datas/joint_tte_monol.csv",row.names = F)

## erreur proportionnelle 

longi = data.frame(id=NA,obs=NA,time=NA)
b0=15
b1 = 0.3
omega_b0 = 0.5
omega_b1 = 0.1
sigma_b = 0.3
h0=0.01
alpha=0.10


rint = rnorm(N,mean=b0 ,sd=omega_b0) # vector des random intercepts
reff = rnorm(N,mean=b1 ,sd=omega_b1) # vector des random slopes
tt = 0:30	# design pour les observations longitudinales
longi = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt)))),time = rep(tt,N))
longi$obs = rint[longi$id] + reff[longi$id]*longi$time + (rint[longi$id] + reff[longi$id]*longi$time)*rnorm(nrow(longi),sd=sigma_b) # on remplit avec les observations

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
b1.pop<- 0.1
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





##### SIMU JOINT 1 mod longi + 2 COMPETING EVENTS format 2 #####
set.seed(1996)
b0.pop<-15
b1.pop<- 0.3
alpha1<- 0.2
alpha2=0

omega_b0<-3
omega_b1<-0.1

sigma_a<-1

p1=0.02
g1 = 0.1

N=100

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
  #F1 = 1-exp(-p1*g1*exp(alpha1*r1)/(1-p1)*(exp((alpha1*r2-g1)*t)/(alpha1*r2-g1)*hyperg_2F1(1,1-alpha1*r2/g1,2-alpha1*r2/g1,-p1*exp(-g1*t)/(1-p1))-1/(alpha1*r2-g1)*hyperg_2F1(1,1-alpha1*r2/g1,2-alpha1*r2/g1,-p1/(1-p1))))
  F2 = (1-rev(F1)[1])*(1-exp(-t/10))
  # fait pareil que cumsum... 
  
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
longi$ytype=1

tte.data$ytype=2

alldata = rbind(longi,tte.data)

write.csv(alldata, "C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension/alexandra/joint_alex/datas/joint_comp_format2.csv",row.names = F)

## en nonlin 
set.seed(1996)
N=100

### erreur additive 

b0<-2
b1<- -0.15
b2 <- -0.2 
a<- 5

omega_b0<-sqrt(0.5)
omega_b1<-sqrt(0.05)
omega_b2<-sqrt(0.06)
omega_a<-sqrt(0.8)

alpha1 = 0.5
sigma_a<-0.9
p1 = 0.02
g1 = 0.1


r1 = rnorm(N,mean=b0 ,sd=omega_b0) # vector des random intercepts
r2 = rnorm(N,mean=b1 ,sd=omega_b1) 
r3 = rnorm(N,mean=b2, sd=omega_b2)
r4 = a*exp(rnorm(N,mean=0, sd=omega_a))
             
tt = 0:30	# design pour les observations longitudinales
longi = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt)))),time = rep(tt,N))
longi$obs = pmin(pmax(r1[longi$id] + r4[longi$id] *(exp(r2[longi$id]*(longi$time-r5[longi$id]))-exp(r3[longi$id]*(longi$time-r5[longi$id]))) + rnorm(nrow(longi),sd=sigma_a),0),24)
  
gettimes = function(r1,r2,r3,r4){ # Fonction qui renvoit un temps d'�v�nement pour un patient avec intercept r1 et slope r2
  t = seq(0,1000,length.out=10000)
  pas = t[2]-t[1]
  h = (p1*g1*exp(-g1*t)/(1-p1*(1-exp(-g1*t))))*exp(alpha1*(r1+r4*(exp(r2*(t))-exp(r3*(t)))))
  F1 = 1-exp(-cumsum(h)*pas)
  F2 = (1-rev(F1)[1])*(1-exp(-t/10))
  
  p1 = rev(F1)[1]
  p2 = 1-p1
    
  u = runif(1)
  event = ifelse(u<p1,1,2)
  
  if(event==1) z=(which(F1/p1>runif(1))[1])*pas 
  if(event==2) z=(which(F2/p2>runif(1))[1])*pas 
  c(z,event)
}
a=mapply(gettimes,r1,r2,r3,r4)

tte.data = data.frame(id=1:N,time = a[1,], obs = a[2,])
tte.data$obs[tte.data$time>30]=0	# Les temps > 30 seront censurées
tte.data$time[tte.data$time>30]=30	# Date de censure pour les temps >30
longi = longi[longi$time<=sapply(longi$id,function(x) tte.data$time[tte.data$id==x]),] # on supprime les observations postérieures au temps d'évènement
  
table(tte.data$obs)
summary(tte.data$time)
plot(survfit(Surv(time,obs)~1, data = tte.data))
  

# for monolix 

DATA_SURV<-data.frame(id=rep(tte.data$id[tte.data$obs!=0],each=2),time=0,obs=0,ytype=0)
for (i in tte.data$id[tte.data$obs!=0]){
  if (tte.data$time[tte.data$id==i]<=30 & tte.data$obs[tte.data$id==i]==1){
    DATA_SURV$ytype[DATA_SURV$id==i] = 2
    DATA_SURV$time[DATA_SURV$id==i][2]<-tte.data$time[tte.data$id==i]
    DATA_SURV$obs[DATA_SURV$id==i][2]<-1
  }
  if (tte.data$time[which(tte.data$id==i)]<=30 & tte.data$obs[which(tte.data$id==i)]==2){
    DATA_SURV$ytype[DATA_SURV$id==i] = 3
    DATA_SURV$time[DATA_SURV$id==i][2]<-tte.data$time[tte.data$id==i]
    DATA_SURV$obs[DATA_SURV$id==i][2]<-1
  }
}
DATA_SURV2<-data.frame(id=rep(tte.data$id[tte.data$obs==0],each=4),time=0,obs=0,ytype=rep(c(2,3),each=2,times=length(tte.data$id[tte.data$obs==0])))
for (i in tte.data$id[tte.data$obs==0]){
  DATA_SURV2$time[DATA_SURV2$id==i][2]<-30
  DATA_SURV2$time[DATA_SURV2$id==i][4]<-30
}

tte.data.final = rbind(DATA_SURV,DATA_SURV2)
alldata = rbind(longi,tte.data.final)

write.csv(alldata, "C:/Users/AlexandraLAVALLEY/Documents/Code_saemix/essais/joint_comp_format2.csv",row.names = F)


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


##### SIMU JOINT 3 mod longi + 2 COMPETING EVENTS #####
N=300

longi1 = data.frame(id=NA,obs=NA,time=NA)
b01 = 4.6
b11 = 0.05
omega_b01 = 2.3
omega_b11 = 0.3
sigma_a1 = 0.3
alpha1=0.6

longi2 = data.frame(id=NA,obs=NA,time=NA)
b02 = 7.4
b12 = 0.003
omega_b02 = 0.04
omega_b12 = 0.005
sigma_a2 = 0.05
alpha2 = -11

longi3 = data.frame(id=NA,obs=NA,time=NA)
b03 = 4.2
b13 = -0.16
omega_b03 = 0.9
omega_b13 = 0.15
sigma_a3 = 0.7
alpha3 = 0.14

p1 = 0.05
g1 = 0.1

rint1 = rnorm(N,mean=b01 ,sd=omega_b01) # vector des random intercepts
reff1 = rnorm(N,mean=b11 ,sd=omega_b11) # vector des random slopes
tt = 0:30	# design pour les observations longitudinales
longi1 = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt)))),time = rep(tt,N))
longi1$obs = rint1[longi1$id] + reff1[longi1$id]*longi1$time + rnorm(nrow(longi1),sd=sigma_b1) 

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
  F1 = 1-exp(-p1*g1*exp(alpha1*(r1-5.78)+alpha2*(r3-7.42)+alpha3*(r5-2.89))/(1-p1)*(exp((alpha1*r2+alpha2*r4+alpha3*r6-g1)*t)/(alpha1*r2+alpha2*r4+alpha3*r6-g1)*hyperg_2F1(1,1-(alpha1*r2+alpha2*r4+alpha3*r6)/g1,2-(alpha1*r2+alpha2*r4+alpha3*r6)/g1,-p1*exp(-g1*t)/(1-p1))-1/(alpha1*r2+alpha2*r4+alpha3*r6-g1)*hyperg_2F1(1,1-(alpha1*r2+alpha2*r4+alpha3*r6)/g1,2-(alpha1*r2+alpha2*r4+alpha3*r6)/g1,-p1/(1-p1))))
  F2 = (1-rev(F1)[1])*(1-exp(-t/10))
  # fait pareil que cumsum... 
  
  p1 = rev(F1)[1]
  p2 = 1-p1
  
  u = runif(1)
  event = ifelse(u<p1,1,2)
  
  if(event==1) z=(which(F1/p1>runif(1))[1])*pas 
  if(event==2) z=(which(F2/p2>runif(1))[1])*pas 
  c(z,event)
}

a=mapply(gettimes,rint1,reff1,rint2,reff2,rint3,reff3)

tte.data = data.frame(id=1:N,time = a[1,], obs = a[2,]) # dataset des données de survie
tte.data$obs[tte.data$time>30]=0	# Les temps > 30 seront censurées
tte.data$time[tte.data$time>30]=30	# Date de censure pour les temps >30
table(tte.data$obs)
summary(tte.data$time)
plot(cuminc(tte.data$time,tte.data$obs,cencode = 0))

longi1 = longi1[longi1$time<=sapply(longi1$id,function(x) tte.data$time[tte.data$id==x]),] # on supprime les observations postérieures au temps d'évènement
longi2 = longi2[longi2$time<=sapply(longi2$id,function(x) tte.data$time[tte.data$id==x]),] # on supprime les observations postérieures au temps d'évènement
longi3 = longi3[longi3$time<=sapply(longi3$id,function(x) tte.data$time[tte.data$id==x]),] # on supprime les observations postérieures au temps d'évènement

longi1$ytype=1
longi2$ytype=2
longi3$ytype=3
tte.data$ytype=4

alldata = rbind(longi1,longi2,longi3,tte.data)

write.csv(alldata,file="C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension/alexandra/joint_alex/datas/joint_multilongi_cmpr.csv",row.names = F)


##### SIMU BOUCLES #####

# comp risks 
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

for (it in 1:50){
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
  
  write.table(alldata, paste0("W:/saemix/dev/datas/JM_comprisks/jointCR",it,".txt"),row.names = F)
  
}

# joint lin+TTE 

set.seed(1996)
N=1000

b0=15
b1 = 0.3
omega_b0 = 0.5
omega_b1 = 0.1
sigma_a = 0.1
h0=0.01
alpha=0.10

for (w in 101:102){
  longi = data.frame(id=NA,obs=NA,time=NA)
  
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
  
  tte.data2<-data.frame(id=tte.data$id, time=0, ev=0)
  ttedata<-rbind(tte.data, tte.data2)
  ttedata <- ttedata[order(ttedata$id, ttedata$time),]
  
  names(ttedata)[3]="obs"
  
  longi$ytype=1
  ttedata$ytype=2
  
  alldata = rbind(longi,ttedata)
  
  write.table(alldata,file=paste0("W:/saemix/dev/comput_se/datas/data",w,".txt"),row.names = F)  
  #write.table(alldata,file=paste0("C:/Users/AlexandraLAVALLEY/Documents/Code_saemix/essais/comput_se/datas/data",w,".txt"),row.names = F) 
}


# joint lin + CR
set.seed(1996)
for (i in 5:100){
  b0.pop<-15
  b1.pop<- 0.1
  alpha1<- 0.15
  
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
    #h = (p1*g1*exp(-g1*t)/(1-p1*(1-exp(-g1*t))))*exp(alpha1*(r1+r2*t))
    #F1 = 1-exp(-cumsum(h)*pas)
    F1 = 1-exp(-p1*g1*exp(alpha1*r1)/(1-p1)*(exp((alpha1*r2-g1)*t)/(alpha1*r2-g1)*hyperg_2F1(1,1-alpha1*r2/g1,2-alpha1*r2/g1,-p1*exp(-g1*t)/(1-p1))-1/(alpha1*r2-g1)*hyperg_2F1(1,1-alpha1*r2/g1,2-alpha1*r2/g1,-p1/(1-p1))))
    F2 = (1-rev(F1)[1])*(1-exp(-t/10))
    # fait pareil que cumsum... 
    
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
  
  longi = longi[longi$time<=sapply(longi$id,function(x) tte.data$time[tte.data$id==x]),] # on supprime les observations postérieures au temps d'évènement
  longi$ytype=1
  
  tte.data$ytype=2
  
  alldata = rbind(longi,tte.data)
  
  write.table(alldata, paste0("W:/saemix/dev/comput_se/jointCR/datas/data",i,".txt"),row.names = F)
}


## joint lin 

library(survival)
## simu data

set.seed(1996)
N=100

### erreur additive 
for (i in 1:100){
  longi = data.frame(id=NA,obs=NA,time=NA)
  b0=15
  b1 = 0.3
  omega_b0 = 0.5
  omega_b1 = 0.1
  sigma_a = 1
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
  tte.data = data.frame(id=1:N,time = mapply(gettimes,rint,reff), obs = 1) # dataset des données de survie
  tte.data$obs[tte.data$time>30]=0	# Les temps > 30 seront censurées
  tte.data$time[tte.data$time>30]=30	# Date de censure pour les temps >30
  longi = longi[longi$time<=sapply(longi$id,function(x) tte.data$time[tte.data$id==x]),] # on supprime les observations postérieures au temps d'évènement
  
  table(tte.data$obs)
  summary(tte.data$time)
  
  longi$ytype=1
  tte.data$ytype=2
  
  # base saemix 
  alldata_saem = rbind(longi,tte.data)
  
  write.table(alldata_saem,file=paste0("C:/Users/AlexandraLAVALLEY/Documents/These/projet 3/calculs_se/joint_lin/datas/data",i,".txt"),row.names = F)
  
}


##### LIN 
set.seed(1996)
N=100

for (it in 1:100){
  longi = data.frame(id=NA,obs=NA,time=NA)
  b0=15
  b1 = 0.3
  omega_b0 = 3
  omega_b1 = 0.2
  sigma_a = 2
  
  rint = rnorm(N,mean=b0 ,sd=omega_b0) # vector des random intercepts
  reff = rnorm(N,mean=b1 ,sd=omega_b1) # vector des random slopes
  tt = 0:10	# design pour les observations longitudinales
  longi = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt)))),time = rep(tt,N))
  longi$obs = rint[longi$id] + reff[longi$id]*longi$time + rnorm(nrow(longi),sd=sigma_a) # on remplit avec les observations
  
  write.table(longi,file=paste0("C:/Users/AlexandraLAVALLEY/Documents/These/projet 3/calculs_se/lin/datas/data",it,".txt"),row.names = F)
  
}




#### SIMUS LASSO 1 #####
# on a 3 marqueurs, 1 seul associ� � l'ev, les 2 autres corr�l�s positivement avec le 1er

library(mvtnorm)

set.seed(1996)
N=100

longi1 = data.frame(id=NA,obs=NA,time=NA)
b01=15
b11 = 0.1
omega_b01 = 2
omega_b11 = 0.1
sigma_a1 = 1
alpha1=0.15

longi2 = data.frame(id=NA,obs=NA,time=NA)
b02 = 10
b12 = -0.1
omega_b02 = 4
omega_b12 = 0.1
sigma_a2 = 1

longi3 = data.frame(id=NA,obs=NA,time=NA)
b03 = 7
b13 = 0.8
omega_b03 = 2
omega_b13 = 0.4
sigma_a3 = 1

cov12 = 0.8 * sqrt(omega_b11) * sqrt(omega_b12)
cov13 = 0.8 * sqrt(omega_b11) * sqrt(omega_b13)

mcov = matrix(c(omega_b01,0,0,0,0,0,
                0,omega_b11,0,cov12,0,cov13,
                0,0,omega_b02,0,0,0,
                0,cov12,0,omega_b12,0,0,
                0,0,0,0,omega_b03,0,
                0,cov13,0,0,0,omega_b13), nrow = 6,byrow = T)
rmvnorm(10,rep(0,6),mcov)

h0 = 0.005
rvec = rmvnorm(N,mean=c(b01,b11,b02,b12,b03,b13), mcov)
rint1 = rvec[,1] # vector des random intercepts marq 1 
reff1 = rvec[,2] # vector des random slopes marq 2
tt = 0:30	# design pour les observations longitudinales
longi1 = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt)))),time = rep(tt,N))
longi1$obs = rint1[longi1$id] + reff1[longi1$id]*longi1$time + rnorm(nrow(longi1),sd=sigma_a1) 

rint2 = rvec[,3] # vector des random intercepts
reff2 = rvec[,4] # vector des random slopes
longi2 = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt)))),time = rep(tt,N))
longi2$obs = rint2[longi2$id] + reff2[longi2$id]*longi2$time + rnorm(nrow(longi2),sd=sigma_a2)

rint3 = rvec[,5] # vector des random intercepts
reff3 = rvec[,6] # vector des random slopes
longi3 = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt)))),time = rep(tt,N))
longi3$obs = rint3[longi3$id] + reff3[longi3$id]*longi3$time + rnorm(nrow(longi3),sd=sigma_a3)


gettimes = function(r1,r2){ # Fonction qui renvoit un temps d'�v�nement pour un patient avec intercept r1 et slope r2
  x = rexp(1)
  H = function(tt) h0*exp(alpha1*r1)/alpha1/r2*(exp(alpha1*r2*tt)-1) - x
  z = try(uniroot(H,interval=c(0,30))$root,silent=T)
  ifelse(class(z)=="try-error",Inf,z)
}
tte.data = data.frame(id=1:N,time = mapply(gettimes,rint1,reff1), obs = 1) # dataset des données de survie
tte.data$obs[tte.data$time>30]=0	# Les temps > 30 seront censurées
tte.data$time[tte.data$time>30]=30	# Date de censure pour les temps >30

longi1 = longi1[longi1$time<=sapply(longi1$id,function(x) tte.data$time[tte.data$id==x]),] # on supprime les observations postérieures au temps d'évènement
longi2 = longi2[longi2$time<=sapply(longi2$id,function(x) tte.data$time[tte.data$id==x]),] # on supprime les observations postérieures au temps d'évènement
longi3 = longi3[longi3$time<=sapply(longi3$id,function(x) tte.data$time[tte.data$id==x]),] # on supprime les observations postérieures au temps d'évènement

table(tte.data$obs)
summary(tte.data$time)
plot(survfit(Surv(time,obs)~1, data = tte.data))

longi1$ytype=1
longi2$ytype=2
longi3$ytype=3
tte.data$ytype=4


alldata = rbind(longi1,longi2,longi3,tte.data)

xyplot(obs ~ time,groups=id,data=longi1,
       type='l',xlim=c(0,30),xlab=("jours"))

write.csv(alldata,file="C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension/alexandra/joint_alex/datas/lasso3marq1ev.csv",row.names = F)





##### SIMU JOINT 1 mod longi iid 2 COMPETING EVENTS iid #####
set.seed(1996)
b0.pop<-15
b1.pop<- 0.1

omega_b0<-3
omega_b1<-0.1

sigma_a<-0.5

p=0.4
r = 0.07

N=500

rint = rnorm(N,mean=b0.pop ,sd=omega_b0) # vector des random intercepts
reff = rnorm(N,mean=b1.pop ,sd=omega_b1) # vector des random slopes
tt = 0:10	# design pour les observations longitudinales
longi = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt)))),time = rep(tt,N))
longi$obs = rint[longi$id] + reff[longi$id]*longi$time + rnorm(nrow(longi),sd=sigma_a) # on remplit avec les observations

  
gettimes=function(it){
  t = seq(0,1000,length.out=10000)
  pas = t[2]-t[1]
  F1 = p*(1-exp(-t*pas))
  F2 = (1-rev(F1)[1])*(1-exp(-t*r))
    
  p1 = rev(F1)[1]
  p2 = 1-p1
    
  u = runif(1)
  event = ifelse(u<p1,1,2)
    
  if(event==1) z=(which(F1/p1>runif(1))[1])*pas 
  if(event==2) z=(which(F2/p2>runif(1))[1])*pas 
  c(z,event)
}

a=mapply(gettimes,1:N)

tte.data = data.frame(id=1:N,time = a[1,], obs = a[2,]) # dataset des données de survie
tte.data$obs[tte.data$time>30]=0	# Les temps > 30 seront censurées
tte.data$time[tte.data$time>30]=30	# Date de censure pour les temps >30
table(tte.data$obs)
summary(tte.data$time)
plot(cuminc(tte.data$time,tte.data$obs,cencode = 0))

longi$ytype=1
tte.data$ytype=2

alldata = rbind(longi,tte.data)

write.csv(alldata, "C:/Users/AlexandraLAVALLEY/Documents/GitHub/saemixextension/alexandra/joint_alex/datas/comp_facile2.csv",row.names = F)


## format monolix

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

alldata = rbind(longi,DATA_SURV)
write.csv(alldata, "C:/Users/AlexandraLAVALLEY/Documents/Code_saemix/essais/CMPR_format1.csv",row.names = F)


# format monolix autre...
DATA_SURV<-data.frame(id=rep(tte.data$id[tte.data$obs!=0],each=2),time=0,obs=0,ytype=0)
for (i in tte.data$id[tte.data$obs!=0]){
  if (tte.data$time[tte.data$id==i]<=30 & tte.data$obs[tte.data$id==i]==1){
    DATA_SURV$ytype[DATA_SURV$id==i] = 2
    DATA_SURV$time[DATA_SURV$id==i][2]<-tte.data$time[tte.data$id==i]
    DATA_SURV$obs[DATA_SURV$id==i][2]<-1
  }
  if (tte.data$time[which(tte.data$id==i)]<=30 & tte.data$obs[which(tte.data$id==i)]==2){
    DATA_SURV$ytype[DATA_SURV$id==i] = 3
    DATA_SURV$time[DATA_SURV$id==i][2]<-tte.data$time[tte.data$id==i]
    DATA_SURV$obs[DATA_SURV$id==i][2]<-1
  }
}
DATA_SURV2<-data.frame(id=rep(tte.data$id[tte.data$obs==0],each=4),time=0,obs=0,ytype=rep(c(2,3),each=2,times=length(tte.data$id[tte.data$obs==0])))
for (i in tte.data$id[tte.data$obs==0]){
  DATA_SURV2$time[DATA_SURV2$id==i][2]<-30
  DATA_SURV2$time[DATA_SURV2$id==i][4]<-30
}

tte.data.final = rbind(DATA_SURV,DATA_SURV2)
alldata = rbind(longi,tte.data.final)

write.csv(alldata, "C:/Users/AlexandraLAVALLEY/Documents/Code_saemix/essais/CMPR_format2.csv",row.names = F)




####### SIMU JOINT 1 mod longi + 1 TTE INDEP (alpha=0) et NON INDEP ############
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

for (it in 1:100){
  rint = rnorm(N,mean=b0 ,sd=omega_b0) # vector des random intercepts
  reff = rnorm(N,mean=b1 ,sd=omega_b1) # vector des random slopes
  tt = 0:30	# design pour les observations longitudinales
  longi = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt)))),time = rep(tt,N))
  longi$obs = rint[longi$id] + reff[longi$id]*longi$time + rnorm(nrow(longi),sd=sigma_a) # on remplit avec les observations
  
  
  write.table(longi,file=paste0("W:/saemix/dev/comput_se/lin/datas/data",it,".txt"),row.names = F)
}





#### simu boucle 3mod longi + comp risks

library(mvtnorm)
library(stringr)

set.seed(1996)
N=300

for (it in 1:100){
  longi1 = data.frame(id=NA,obs=NA,time=NA)
  b01=7.4
  b11 = 0.003
  omega_b01 = 0.04**2
  omega_b11 = 0.005**2
  sigma_a1 = 0.05
  alpha1=-11
  
  longi2 = data.frame(id=NA,obs=NA,time=NA)
  b02 = 4.2
  b12 = -0.16
  omega_b02 = 0.9**2
  omega_b12 = 0.15**2
  sigma_a2 = 0.7
  alpha2 = 0.6
  
  longi3 = data.frame(id=NA,obs=NA,time=NA)
  b03 = 4.6
  b13 = -0.15
  b23 = -0.16
  a3 = 5.3
  omega_b03 = 2**2
  omega_b13 = 0.1**2
  omega_b23 = 0.07**2
  omega_a3 = 0.8**2
  sigma_a3 = 0.3
  alpha3 = 0.14
  
  longi4 = data.frame(id=NA,obs=NA,time=NA)
  b04 = 5.2
  b14 = 0.02
  omega_b04 = 1.7**2
  omega_b14 = 0.08**2
  sigma_a4 = 0.7
  alpha4 = 0
  
  longi5 = data.frame(id=NA,obs=NA,time=NA)
  b05 = 7
  b15 = -0.05
  omega_b05 = 0.8**2
  omega_b15 = 0.1**2
  sigma_b5 = 0.08
  alpha5 = 0
  
  longi6 = data.frame(id=NA,obs=NA,time=NA)
  b06 = 5.2
  b16 = -0.07
  omega_b06 = 0.7**2
  omega_b16 = 0.09**2
  sigma_b6 = 0.1
  alpha6 = 0
  
  longi7 = data.frame(id=NA,obs=NA,time=NA)
  b07 = 338
  b17 = -5
  omega_b07 = 107**2
  omega_b17 = 5**2
  sigma_b7 = 0.2
  alpha7 = 0
  
  # cov = cor * sigmax * sigmay
  # on fixe les cov à 0.8 
  # marq 4 corr au marq 1 sur la pente 
  cov14 = 0.8 * sqrt(omega_b11) * sqrt(omega_b14)
  
  # marq 5 corr au marq 2 sur la pente 
  cov25 = 0.8 * sqrt(omega_b12) * sqrt(omega_b15)
  
  mcov = matrix(c(omega_b01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,omega_b11,0,0,0,0,0,0,0,cov14,0,0,0,0,0,0,
                  0,0,omega_b02,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,omega_b12,0,0,0,0,0,0,0,cov25,0,0,0,0,
                  0,0,0,0,omega_b03,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,omega_b13,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,omega_b23,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,omega_a3,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,omega_b04,0,0,0,0,0,0,0,
                  0,cov14,0,0,0,0,0,0,0,omega_b14,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,omega_b05,0,0,0,0,0,
                  0,0,0,cov25,0,0,0,0,0,0,0,omega_b15,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,omega_b06,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,omega_b16,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,omega_b07,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,omega_b17),nrow = 16)
  
  colnames(mcov)  = c("omega_b01","omega_b11","omega_b02","omega_b12","omega_b03","omega_b13","omega_b23","omega_a3","omega_b04","omega_b14","omega_b05","omega_b15","omega_b06","omega_b16","omega_b07","omega_b17")
  rownames(mcov) = c("omega_b01","omega_b11","omega_b02","omega_b12","omega_b03","omega_b13","omega_b23","omega_a3","omega_b04","omega_b14","omega_b05","omega_b15","omega_b06","omega_b16","omega_b07","omega_b17")
  rvec = rmvnorm(N,mean=c(b01,b11,b02,b12,b03,b13,b23,a3,b04,b14,b05,b15,b06,b16,b07,b17), mcov)
  rint1 = rvec[,1] # vector des random intercepts marq 1 
  reff1 = rvec[,2] # vector des random slopes marq 1
  tt1 = seq(0,30,1.5)	# design pour les observations longitudinales
  longi1 = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt1)))),time = rep(tt1,N))
  longi1$obs = rint1[longi1$id] + reff1[longi1$id]*longi1$time + rnorm(nrow(longi1),sd=sigma_a1) 
  
  rint2 = rvec[,3] # vector des random intercepts
  reff2 = rvec[,4] # vector des random slopes
  tt2 = seq(0,30,2)
  longi2 = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt2)))),time = rep(tt2,N))
  longi2$obs = rint2[longi2$id] + reff2[longi2$id]*longi2$time + rnorm(nrow(longi2),sd=sigma_a2)
  
  rint3 = rvec[,5] # vector des random intercepts
  reff3 = rvec[,6] # vector des random slopes
  refff4 = rvec[,7]
  refff5 = a3*exp(rnorm(N,0,sqrt(omega_a3))) # attention param a lognormal !! on retire indpdemment comme il est pas lié aux autres pas grave 
  tt3 = seq(0,30,2)
  longi3 = data.frame(id = as.vector(sapply(1:N,function(x) rep(x,length(tt3)))),time = rep(tt3,N))
  longi3$obs = rint3[longi3$id] + refff5[longi3$id] *(exp(reff3[longi3$id]*longi3$time)-exp(refff4[longi3$id]*longi3$time)) + rnorm(nrow(longi3),sd=sigma_a3)
  
  
  med1 = median(longi1$obs)
  med2 = median(longi2$obs)
  med3 = median(longi3$obs)
  
  gettimes = function(r1,r2,r3,r4,r5,r6,r7,r8){ # Fonction qui renvoit un temps d'évènement pour un patient avec intercept r1 et slope r2
    
    t = seq(0,1000,length.out=10000)
    pas = t[2]-t[1]
    h = (p1*g1*exp(-g1*t)/(1-p1*(1-exp(-g1*t))))*exp(alpha1*(r1+r2*t-med1)+alpha2*(r3+r4*t-med2)+alpha3*(r5+r8*(exp(r6*t)-exp(r7*t))-med3))
    F1 = 1-exp(-cumsum(h)*pas)
    F2 = (1-rev(F1)[1])*(1-exp(-t/10))
    
    prob1 = rev(F1)[1]
    prob2 = 1-p1
    
    u = runif(1)
    event = ifelse(u<prob1,1,2)
    
    if(event==1) z=(which(F1/prob1>runif(1))[1])*pas 
    if(event==2) z=(which(F2/prob2>runif(1))[1])*pas 
    
    c(z,event)
    
  }
  
  p1 = 0.05
  g1 = 0.1
  
  a=mapply(gettimes,rint1,reff1,rint2,reff2,rint3,reff3,refff4,refff5)
  
  tte.data = data.frame(id=1:N,time = a[1,], obs = a[2,]) # dataset des donnÃ©es de survie
  tte.data$obs[tte.data$time>30 | is.na(tte.data$time)==T]=0	# Les temps > 30 seront censurÃ©es
  tte.data$time[tte.data$time>30 | is.na(tte.data$time)==T]=30 # Date de censure pour les temps >30
  tte.data$ytype=4
  #table(tte.data$obs)
  #summary(tte.data$time)
  #plot(cuminc(tte.data$time,tte.data$obs,cencode = 0))
  
  n_dcd = c(n_dcd,length(tte.data$id[tte.data$obs==1]))
  n_sor = c(n_sor,length(tte.data$id[tte.data$obs==2]))
  n_cen = c(n_cen,length(tte.data$id[tte.data$obs==0]))
  
  longi1 = longi1[longi1$time<=sapply(longi1$id,function(x) tte.data$time[tte.data$id==x]),] # on supprime les observations postÃ©rieures au temps d'Ã©vÃ¨nement
  longi2 = longi2[longi2$time<=sapply(longi2$id,function(x) tte.data$time[tte.data$id==x]),] # on supprime les observations postÃ©rieures au temps d'Ã©vÃ¨nement
  longi3 = longi3[longi3$time<=sapply(longi3$id,function(x) tte.data$time[tte.data$id==x]),] # on supprime les observations postÃ©rieures au temps d'Ã©vÃ¨nement
  
  longi1$ytype=1
  longi2$ytype=2
  longi3$ytype=3
  
  
  alldata = rbind(longi1,longi2,longi3,tte.data)
  
  #par(mfrow=c(1,2))
  #xyplot(obs ~ time,groups=id,data=longi2[longi2$id==3,],type='l',xlim=c(0,30),xlab=("jours"))
  #xyplot(obs ~ time,groups=id,data=longi5[longi5$id==3,],type='l',xlim=c(0,30),xlab=("jours"))
  
  write.table(alldata,file=paste0("W:/RISCOV/simu/back7/datas/data",it,".txt"),row.names = F)
  
}