plotDebug <- function(x, y , range=c(0,1), psi=NULL, predictors=NULL, ...) {
            # If verbose=TRUE, print messages
            args1<-match.call(expand.dots=TRUE)
            list.args <- list(...)
            i1<-match("verbose",names(args1))
            if(is.na(i1)) verbose<-FALSE else verbose<-eval(args1[[i1]])
            # Set psi by default to the starting parameters given in the model (if not given as arguments)
            if(is.null(psi)) psi<-x@psi0[1,,drop=FALSE]
            if(is.null(dim(psi)[1])) psi<-matrix(psi,nrow=1) else psi<-psi[1,,drop=FALSE]
            npred<-length(x@name.predictors)
            if(npred==0 & is.null(predictors)) npred<-1 else {
              if(npred==0 & !missing(predictors)) {
                npred<-1+length(predictors)
              } else {
                if(npred>1 & (missing(predictors) || length(predictors)<(npred-1))) {
                  if(verbose) message("Please provide the value of the predictors other than X\n")
                  return("Missing predictors")
                }
              }
            }
            if(length(x@name.response)>1) {
              if(verbose) message("Currently the plot can only be obtained for single-response models.\n")
              return()
            }
            if(length(x@name.X)>0 & length(x@name.predictors)>0 && x@name.X != x@name.predictors[1]){
              if(verbose) message("Warning: X predictor supposed to be on the first axis, exiting without plot\n")
              return()
            }
            npts<-100
            id<-rep(1,npts+1)
            xval<-range[1]+(range[2]-range[1])*c(0:100)/100
            if(npred==1) {
              xdep<-matrix(xval,ncol=1)
            } else {
              xdep<-cbind(xval,matrix(rep(predictors[1:(npred-1)],(npts+1)), byrow=T,nrow=(npts+1)))
              if(length(x@name.X)>0) {
                colnames(xdep)<-c(x@name.X,x@name.predictors[x@name.predictors!=x@name.X])
                xdep<-xdep[,match(x@name.predictors,colnames(xdep))]
              } else colnames(xdep)<-paste("Predictor",1:npred)
            }
            ypred<-try(x@model(psi,id,xdep))
            if(!is.numeric(ypred) & verbose) {
              message("Problem when attempting to obtain predictions from the model.\n")
              message("Usage: plot(x, range=c(0,1), psi, predictors) \n")
              message("Possible solutions can be:\n")
              message("   1. provide suitable values for X (option range=c(<lower bound>, <upper bound>))\n")
              message("   2. provide values for additional predictors (option predictors=c(<value for predictor 1>, <value for predictor 2>, ...)).\n")
              message("   3. check values for the model parameters (defaults to component psi0[1,] of the model).\n")
              message("   4. the predictor used the X-axis is assumed to be in the first column; please check your model is written in a compatible way.\n")
            } else {
              if(length(x@name.X)==0 | length(x@name.predictors)==0) {
                if(verbose) message("Warning: X predictor supposed to be on the first axis\n")}
              if(verbose) message("Plot characteristics:\n")
              if(npred>1) {
                for(j in 1:dim(xdep)[2]) {
                  if(length(x@name.X)==0) {
                    if(j>1) message("   predictor:",colnames(xdep)[j],"=",xdep[1,j],"\n")
                  } else {
                    if(colnames(xdep)[j]!=x@name.X) message("    predictor:",colnames(xdep)[j],"=",xdep[1,j],"\n")
                  }
                }}
              if(verbose) message("   range for X-axis:",min(xval),"-",max(xval),"\n")
              if(verbose) message("   parameters used: ", paste(x@name.modpar,"=",psi[1,],collapse=", "),"\n")
              plot(xval,ypred,type="l",xlab=ifelse(length(x@name.X)==0, "X",x@name.X),ylab=ifelse(length(x@name.response)==0, "Response",x@name.response),...)
            }
}

model1cpt <-function(psi,id,xidep) { 
  tim<-xidep[,1]  
  dose<-xidep[,2]
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  k<-CL/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}
model1cpt.bis <-function(psi,id,xidep) { 
  tim<-xidep[,2]  
  dose<-xidep[,1]
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  k<-CL/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}
x<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", 
               psi0=matrix(c(1.5,30,1), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))))
x.bis<-saemixModel(model=model1cpt.bis,description="One-compartment model with first-order absorption", 
               psi0=matrix(c(1.5,30,1), ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))))
x.bis@name.X <- "tim"
x.bis@name.predictors <- c("dose","tim")

plotDebug(x, range=c(0,24), psi=c(1.5,20,2), predictors=350, verbose=TRUE)

plotDebug(x.bis, range=c(0,24), psi=c(1.5,20,2), predictors=350, verbose=TRUE)

