########################################################################
# Debugging functions
########################################################################

# Predict for an saemixModel object

if(FALSE) {
  model1cpt<-function(psi,id,xidep) { 
    dose<-xidep[,1]
    tim<-xidep[,2]  
    ka<-psi[id,1]
    V<-psi[id,2]
    CL<-psi[id,3]
    k<-CL/V
    ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
    return(ypred)
  }
  object<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))))
  xidep<-data.frame(dose=100, tim=seq(0,24,2))
  xtim<-rep(xidep[,2],5)
  id<-rep(1:5, each=length(xidep[,2]))
  predictors<-do.call(rbind,rep(list(xidep),5))
  psi<-do.call(rbind,rep(list(c(2, 25, 0.5)),5))
}

########################################################################
# Plot function

library(ggplot2)

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  k<-CL/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}
xmod<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))))

xidep1<-data.frame(dose=100, tim=c(0,1,2,seq(4,48,4)))
xtim<-xidep1[,2]
id<-rep(1:12, each=length(xtim))
zedat<-do.call(rbind,rep(list(xidep1),12))
psi0<-c(2, 25, 1)
ypred<-c()
for(i in 1:12) {
  psi<-psi0*exp(rnorm(3,mean=0,sd=0.2))
  ypred<-c(ypred,
           xmod["model"](t(psi), rep(1,length(xtim)), xidep1)*(1+rnorm(length(xtim),sd=0.2)))
}
zedat<-cbind(id=id,zedat)
zedat$y<-ypred
  
ypred<-predict.saemixmodel(xmod, zedat[,c("dose","tim")], psi=psi0)
zedat$ppred<-ypred$predictions$pred

ggplot(data=zedat, aes(x=tim, y=y, group=id)) + geom_point() + geom_line(aes(x=tim, y=ppred)) + facet_wrap(.~id, ncol=3) + 
  labs(x="Time (h)", y="Concentrations (mg/L)") + theme_bw()

# Smoothing
ntim<-200
xtim2<-seq(min(zedat$tim), max(zedat$tim), length.out=ntim)
xidep2<-data.frame(dose=zedat$dose[1],tim=xtim2)
zedat2<-data.frame(id=rep(unique(zedat$id), each=ntim))
zedat2<-cbind(zedat2,
              do.call(rbind,rep(list(xidep2),length(unique(zedat$id)))))
zedat2$dose<-zedat$dose[match(zedat2$id, zedat$id)]

ypred<-predict.saemixmodel(xmod, zedat2[,c("dose","tim")], psi=psi0)
zedat2$ppred<-ypred$predictions$pred

ggplot(data=zedat, aes(x=tim, y=y, group=id)) + geom_point() + geom_line(data=zedat2,aes(x=tim, y=ppred)) + facet_wrap(.~id, ncol=3) + 
  labs(x="Time (h)", y="Concentrations (mg/L)") + theme_bw()

# Making a function... not that simple :-/ (lacking a lot of information only available when coupled with a dataset)

plot.saemixModel <- function(model, predictors, psi=NA, id=NA, name.x="x") {
  if(is.na(id) || length(id)!=dim(predictors)[1]) 
    id<-rep(1,dim(xidep)[1]) 
  
  ypred<-predict.saemixmodel(saemixModel, predictors, psi=psi, id=id)
  ypred<-ypred$predictions$pred
  
  
  xidep<-predictors
  idkeep<-id
  if(max(id)>length(unique(id))) { # indexes need to go from 1 to N
    id1<-1:length(unique(id))
    id2<-unique(id)
    id<-id1[match(id,id2)]
  }
  if(is.na(psi)) psi<-object["psi0"][1,,drop=FALSE]
  if(is.null(dim(psi))) psi<-as.data.frame(t(psi)) # psi given as a vector
  if(dim(psi)[2] != object@nb.parameters) {
    message(paste0("psi must have a number of columns equal to the number of parameters in the model (",object@nb.parameters,")\n")
    )
    return()
  }
  if(dim(psi)[1]==1 & length(unique(id))>1)
    psi<-do.call(rbind,rep(list(psi,length(unique(id)))))
  
}

# Making a function with an saemixData object to extract the proper predictors
plot.saemixModel <- function(smx.model, smx.data, psi=NA) {
  if(is.na(psi)) psi<-smx.model["psi0"][1,,drop=FALSE]
  if(is.null(dim(psi))) psi<-as.data.frame(t(psi)) # psi given as a vector
  if(dim(psi)[2] != smx.model@nb.parameters) {
    message(paste0("psi must have a number of columns equal to the number of parameters in the model (",smx.model@nb.parameters,")\n"))
    return()
  }
  # if(dim(psi)[1]>1 & dim(psi[1])!=smx.data@N) {
  #   message(paste0("psi must have a number of lines equal to the number of subjects (",smx.data@N,")\n"))
  #   return()
  # }
  if(dim(psi)[1]==1 || dim(psi)[1]<smx.data@N)
    psi<-do.call(rbind,rep(list(psi),length.out=smx.data@N))
  nvalues<-100
  xt<-seq(min(smx.data@data[,smx.data["name.X"]]), max(smx.data@data[,smx.data["name.X"]]), length.out=nvalues)
  xidep<-data.frame(x=xt)
  colnames(xidep)<-smx.data["name.X"]
  if(length(smx.data@name.predictors)>1) {
    id<-smx.data@data[,smx.data@name.group]
    otherpred<-smx.data@name.predictors[smx.data@name.predictors != smx.data["name.X"]]
    x1<-smx.data@data[match(unique(id), id), otherpred, drop=FALSE]
    dat1<-NULL
    for(i in 1:length(unique(id)))
      dat1<-rbind(dat1, 
                  do.call(rbind,rep(list(x1[i,,drop=FALSE]), nvalues)))
    xidep<-cbind(xidep, dat1)
    colnames(xidep[2:dim(xidep)[2]])<-otherpred
    xidep<-xidep[,smx.data["name.predictors"]] # Sort the predictors back in the correct order...
  }
  id<-rep(1:length(unique(id)), each=nvalues)
  y<-predict.saemixmodel(smx.model, predictors=xidep, psi=psi, id=id)
  gpred<-cbind(id=id,xidep,y=y$predictions$pred)
  colnames(gpred)[colnames(gpred)==smx.data@name.X]<-"x"

  gdat<-smx.data@data
  colnames(gdat)[colnames(gdat)==smx.data@name.X]<-"x"
  colnames(gdat)[colnames(gdat)==smx.data@name.response]<-"y"
  colnames(gdat)[colnames(gdat)==smx.data@name.group]<-"id"
  
  g1<-ggplot(data=gdat, aes(x=x, y=y, group=id)) + geom_point() + geom_line(data=gpred,aes(x=x, y=y)) + facet_wrap(.~id, nrow=3, ncol=4) + 
    labs(x=smx.data@name.X, y=smx.data@name.response) + theme_bw()
  return(g1)
} 

smx.data <- saemixData(name.data=file.path(datDir,"theo.saemix.tab"),header=T,na=".", name.group=c("Id"),
                       name.predictors=c("Dose","Time"),name.response=c("Concentration"),units=list(x="hr",y="mg/L"), name.X="Time",verbose=F)
model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  k<-CL/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}
smx.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))))

plot.saemixModel(smx.model, smx.data)
p1<-plot.saemixModel(smx.model, smx.data)

p2<-plot.saemixModel(smx.model, smx.data, psi=c(2, 30, 1.5))

p1
p2

########################################################################
