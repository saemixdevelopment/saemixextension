
namfich<-paste('param_',namsimdat,isim,".tab",sep="")
parsim<-read.table(file.path(datDir,namfich),header=T)
i<-1
par(mfrow=c(2,2))
for(i in 1:3) {
  plot(yfit@results@cond.mean.psi[,i],yfit@results@map.psi[,i],pch=20)
  abline(0,1)
}
par(mfrow=c(2,2))
for(i in 1:3) {
  plot(parsim[,(i+1)],yfit@results@map.psi[,i],pch=20)
  abline(0,1)
}
par(mfrow=c(2,2))
for(i in 1:3) {
  plot(yfit@results@cond.mean.psi[,i],yfit2@results@cond.mean.psi[,i],pch=20)
  abline(0,1)
}
par(mfrow=c(2,2))
for(i in 1:3) {
  plot(parsim[,(i+1)],yfit@results@cond.mean.psi[,i],pch=20)
  abline(0,1)
  x1<-mean((parsim[,(i+1)]-yfit@results@cond.mean.psi[,i])/parpop[i])
  x2<-sd((parsim[,(i+1)]-yfit@results@cond.mean.psi[,i])/parpop[i])
  cat("RBias=",x1,"   RMSE=",x2,"\n")
}
par(mfrow=c(2,2))
for(i in 1:3) {
  plot(parsim[,(i+1)],yfit2@results@cond.mean.psi[,i],pch=20)
  abline(0,1)
  x1<-mean((parsim[,(i+1)]-yfit2@results@cond.mean.psi[,i])/parpop[i])
  x2<-sd((parsim[,(i+1)]-yfit2@results@cond.mean.psi[,i])/parpop[i])
  cat("RBias=",x1,"   RMSE=",x2,"\n")
}

saemix.fit<-conddist.saemix(saemix.fit,nsamp = nsamp)
idx.eps<-saemix.fit@model@indx.res
idx.iiv<-saemix.fit@model@indx.omega
res<-saemix.fit@results
xcal<-calcul.FIM.recoded(saemix.fit)

l1<-c(isim,res@fixed.effects, diag(res@omega)[idx.iiv],res@respar[idx.eps], res@se.fixed, res@se.omega[idx.iiv], res@se.respar[idx.eps],res@ll.lin)
l2<-c(isim,res@fixed.effects, diag(res@omega)[idx.iiv],res@respar[idx.eps], sqrt(diag((solve(xcal$popFIM)))),res@ll.lin)
write(l1,namOrig,ncol=length((headersOrig)),append=T)
write(l2,namOrigSE,ncol=length((headersOrig)),append=T)

# Debugging FIM... again...
saemixObject<-yfit


