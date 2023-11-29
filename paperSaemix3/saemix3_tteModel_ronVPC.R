# Run this after saemix3_tteModel.R to obtain VPC by Ron Keizer's package vpc

library(vpc)

n<-dim(simtte.fit@sim.data@datasim)[1]
simron<-cbind(simtte.fit@sim.data@datasim,dv=rep(c(0,1), n/2))
colnames(simron)<-c("id","sim","time","dv")

vpc_tte(sim = simron,
        obs = lung.saemix,
        #        rtte = FALSE,
        sim_cols=list(id="id", dv = "dv", idv = "time"), obs_cols=list(id="id", dv="status", idv = "time"))


if(FALSE) { # Setting late simulated times to the maximum event time (a censored event)
  simron$time[simron$time>max(lung.saemix$time)]<-max(lung.saemix$time)
  vpc_tte(sim = simron,
          obs = lung.saemix,
          #        rtte = FALSE,
          sim_cols=list(id="id", dv = "dv", idv = "time"), obs_cols=list(id="id", dv="status", idv = "time"))
}
