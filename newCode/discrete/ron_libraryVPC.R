library(vpc)

############ Categorical VPC

# simple function to simulate categorical data for single individual
sim_id <- function(id = 1) {
  n <- 10
  logit <- function(x) exp(x) / (1+exp(x))
  data.frame(id = id, time = seq(1, n, length.out = n),
             dv = round(logit((1:n) - n/2 + rnorm(n, 0, 1.5))) )
}
## simple function to simulate categorical data for a trial
sim_trial <- function(i = 1, n = 20) { # function to simulate categorical data for a trial
  data.frame(sim = i, do.call("rbind", lapply(1:n, sim_id)))
}

## simulate single trial for 20 individuals
obs <- sim_trial(n = 20)

## simulate 200 trials of 20 individuals
sim <- do.call("rbind", lapply(1:200, sim_trial, n = 20))

## Plot categorical VPC
vpc_cat(sim = sim, obs = obs)

############ Time-to-event VPC
data(rtte_obs_nm)
data(rtte_sim_nm)

# treat RTTE as TTE, no stratification
vpc_tte(sim = rtte_sim_nm[rtte_sim_nm$sim <= 20,],
        obs = rtte_obs_nm,
        rtte = FALSE,
        sim_cols=list(dv = "dv", idv = "t"), obs_cols=list(idv = "t"))

# Only first event
# idx<-which()


## Applied to tte.fit from saemix3_tteModel.Rmd notebook
n<-dim(ysim.tte@sim.data@datasim)[1]
simron<-cbind(ysim.tte@sim.data@datasim,dv=rep(c(0,1), n/2))
colnames(simron)<-c("id","sim","time","dv")

vpc_tte(sim = simron,
        obs = lung.saemix,
#        rtte = FALSE,
        sim_cols=list(id="id", dv = "dv", idv = "time"), obs_cols=list(id="id", dv="status", idv = "time"))

write.table(simron[simron$sim<21,],"/home/eco/work/saemix/hotline/ronvpc/simulatedTTE.csv", quote=F, row.names=F)

