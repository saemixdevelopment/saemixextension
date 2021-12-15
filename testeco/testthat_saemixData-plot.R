context("Plots for SaemixData")

test_that("Spaghetti plot", {
  x<-saemixData(name.data=file.path(datDir,"theo.saemix.tab"),header=T,na=".", name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),units=list(x="hr",y="mg/L"), name.X="Time",verbose=F)
  plot(x, type="b", col='blue', main="Spaghetti plot of theophylline data")
})


test_that("Individual plot", {
  x<-saemixData(name.data=file.path(datDir,"theo.saemix.tab"),header=T,na=".", name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),units=list(x="hr",y="mg/L"), name.X="Time",verbose=F)
  plot(x, individual=TRUE, ilist=1:6)
})

test_that("Individual plot, automatic reordering", {
  theo.saemix<-read.table(file.path(datDir,"theo.saemix.tab"),header=T,na=".")
  theo2<-theo.saemix[order(theo.saemix$Time),]
  x2<-saemixData(name.data=theo2,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L", covariates=c("kg","-")), name.X="Time")
  plot(x2, individual=TRUE, ilist=1:6)
})
