
saemix.model<-saemix.fit@model
idx.eps<-saemix.fit@model@indx.res
idx.iiv<-saemix.fit@model@indx.omega

### Fitting bootstrap - only population parameters and FIM for SE, no log-likelihood
saemix.bootOpt<-list(fix.seed=F,directory="current",displayProgress=F, save.graphs=F,map=F,ll.is=F)

res.boot<-data.frame()
for(iboot in 1:nboot) {
  data.bootCase <- dataGen.case(saemix.fit)
  fit.bootCase<-saemix(saemix.model,data.bootCase,saemix.bootOpt)
  res<-fit.bootCase@results
  l1<-c(iboot,res@fixed.effects, diag(res@omega)[idx.iiv],res@ll.lin)
  res.boot<-rbind(res.boot,l1) 
}
lnam<-c(saemix.model@name.fixed,saemix.model@name.random)
colnames(res.boot)<-c("Replicate",lnam,"LL")

