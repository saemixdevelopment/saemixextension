# Run *after* saemix3_categoricalModel.R

library(xtable)

###### Exact FIM by AGQ (code by Sebastian Ueckert)
# Code Sebastian
source(file.path(dirAGQ,"default_settings.R"))
source(file.path(dirAGQ,"helper_functions.R"))
source(file.path(dirAGQ,"integration.R"))
source(file.path(dirAGQ,"model.R"))

saemix.fit <- ord.fit

# Setting up ordinal model
model <- Model$new(
  parameter_function = function(mu, b) list(alp1=mu[1]+b[1], alp2=mu[2], alp3=mu[3], alp4=mu[4], beta=mu[5] + b[2]),
  log_likelihood_function = function(y, design, alp1, alp2, alp3, alp4, beta) {
    logit1<-alp1 + beta*design$time
    logit2<-logit1+alp2
    logit3<-logit2+alp3
    logit4<-logit3+alp4
    pge1<-exp(logit1)/(1+exp(logit1))
    pge2<-exp(logit2)/(1+exp(logit2))
    pge3<-exp(logit3)/(1+exp(logit3))
    pge4<-exp(logit4)/(1+exp(logit4))
    pobs = (y==1)*pge1+(y==2)*(pge2 - pge1)+(y==3)*(pge3 - pge2)+(y==4)*(pge4 - pge3)+(y==5)*(1 - pge4)
    log(pobs)
  }, 
  simulation_function = function(design, alp1, alp2, alp3, alp4, beta) {
    logit1<-alp1 + beta*design$time
    logit2<-logit1+alp2
    logit3<-logit2+alp3
    logit4<-logit3+alp4
    pge1<-exp(logit1)/(1+exp(logit1))
    pge2<-exp(logit2)/(1+exp(logit2))
    pge3<-exp(logit3)/(1+exp(logit3))
    pge4<-exp(logit4)/(1+exp(logit4))
    x<-runif(length(time))
    ysim<-1+as.integer(x>pge1)+as.integer(x>pge2)+as.integer(x>pge3)+as.integer(x>pge4)
  },
  inverse_simulation_function = function(design, urand,alp1, alp2, alp3, alp4, beta) {
    if(is.null(urand)) return(seq_along(design$time))
    logit1<-alp1 + beta*design$time
    logit2<-logit1+alp2
    logit3<-logit2+alp3
    logit4<-logit3+alp4
    pge1<-exp(logit1)/(1+exp(logit1))
    pge2<-exp(logit2)/(1+exp(logit2))
    pge3<-exp(logit3)/(1+exp(logit3))
    pge4<-exp(logit4)/(1+exp(logit4))
    1+as.integer(urand>pge1)+as.integer(urand>pge2)+as.integer(urand>pge3)+as.integer(urand>pge4)
  },
  mu = saemix.fit@results@fixed.effects,
  omega = saemix.fit@results@omega[c(1,5),c(1,5)])


# define settings (agq with 3 grid points, quasi random monte-carlo and 500 samples)
settings <- defaults.agq(gq.quad_points = 3,  y_integration.method = "qrmc", y_integration.n_samples = 500, seed = 3257)

#### Design
# Checking whether everyone has the same visits - yes
time.patterns<-tapply(knee.saemix$time, knee.saemix$id, function(x) paste(x,collapse="-"))
unique(time.patterns)

# same 4 times for all subjects (0, 3, 7, 10)
design <- data.frame(time=sort(unique(knee.saemix$time)))
fim <- length(unique(knee.saemix$id)) * calc_fim(model, design, settings)
print(fim)
# calculate rse
rse <- calc_rse(model, fim)
print(rse)

est.se<-sqrt(diag(solve(fim)))
df <- data.frame(param=c(model$mu,diag(model$omega)),se=est.se)
df$rse <- abs(df$se/df$param*100)

print(df)

###### Comparing the SE with the different approaches
# Adding the exact FIM estimates to df2
l1<-paste0(format(par.estim, digits=2)," (",format(est.se,digits=2, trim=T),")")
ci.low <- par.estim - 1.96*est.se
ci.up <- par.estim + 1.96*est.se
l2<-paste0("[",format(ci.low, digits=2),", ",format(ci.up,digits=2, trim=T),"]")
df2<-cbind(df2, l1, l2)
colnames(df2)[7:8]<-paste0("FIM.",c("estimate","CI"))
print(df2)


