library(testthat)
saemixDir<-"/home/eco/work/saemix/saemixextension/newCode"
source(file.path(saemixDir,"errorModels.R"))

y1<-new(Class="SaemixDiscreteOutcome", name.outcome="success", type.outcome="categorical")
y2<-new(Class="SaemixDiscreteOutcome", name.outcome="relapse", type.outcome="event")
y3<-new(Class="SaemixContinuousOutcome", name.outcome="y1", error.model='proportional')
y4<-try(new(Class="SaemixContinuousOutcome", name.outcome="y2", error.model='user', error.function<-user.error1))
y5<-new(Class="SaemixContinuousOutcome", name.outcome="y2", error.model='user', error.npar=3, error.function<-user.error1)

print(y1)
print(y5)
show(y3)
showall(y3)

