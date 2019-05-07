if(iscenar %in% 1:9) {
  c("Simulation",saemix.model@name.modpar,saemix.model@name.random,saemix.model@name.sigma[saemix.model@indx.res])
  
  namOrig<-file.path(resDir,paste("parpop_",namsimdat,".res",sep=""))
  namSE<-file.path(resDir,paste("sepop_",namsimdat,".res",sep=""))

  for (isim in 1:nsim) {
    namfich<-paste('data_',namsimdat,isim,".tab",sep="")
    saemix.data<-saemixData(file.path(datDir,namfich),header=T,name.group="id",name.predictors="dose",name.response="y",units=list(x="-",y="-"))
    
    yfit<-saemix(saemix.model,saemix.data,saemix.options)
    yfit2<-conddist.saemix(yfit)
    # Saving population parameters
    idx.eps<-yfit@model@indx.res
    idx.iiv<-yfit@model@indx.omega
    vec<-summary(yfit@results)
    l1<-unlist(c(isim,vec$fixed.effects[,2],vec$random.effects[,2],vec$logLik[1:2,2]))
    l2<-unlist(c(isim,vec$fixed.effects[,3],vec$random.effects[,3]))
    if(isim==1) {
      headersOrig<-c("Simulation",vec$fixed.effects[,1],as.character(vec$random.effects[,1]),"LL.lin","LL.IS")
      headersSE<-headersOrig[1:(length(headersOrig)-2)]
      headersSE[-c(1)]<-paste("SE.", headersSE[-c(1)],sep="")
      write(headersOrig,namOrig,ncol=length(headersOrig))
      write(headersSE,namSE,ncol=length(headersSE))
    }
    write(l1,namOrig,ncol=length((headersOrig)),append=T)
    write(l2,namSE,ncol=length((headersSE)),append=T)
    
    # Saving individual parameters
    npar<-dim(yfit@results@map.psi)[2]
    indpar<-data.frame(id=unique(yfit@data@data$id), yfit@results@map.psi, yfit@results@cond.mean.psi, yfit2@results@cond.mean.psi)
    colnames(indpar)[(npar+2):(3*npar+1)]<-c(paste(colnames(indpar)[2:(npar+1)],"cmean",sep="."),paste(colnames(indpar)[2:(npar+1)],"cmean2",sep="."))
    namfich<-paste('indpar_',namsimdat,isim,".tab",sep="")
    write.table(indpar,file.path(resDir,namfich),row.names=FALSE,quote=FALSE)
  }
}

rm(modfun)