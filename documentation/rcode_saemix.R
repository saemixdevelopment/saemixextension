# Theophylline


#####################################################################################
# Paper submitted to JSS
#### Example 3.2 - Theophylline PK

library(saemix)

data(theo.saemix)

theo.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA,
   name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"), 
   name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), 
   name.X="Time")

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

###########################
# Model with a covariate effect (effect of Weight on CL)

theo.model<-saemixModel(model=model1cpt,
   description="One-compartment model with first-order absorption", 
   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, 
   c("ka","V","CL"))),transform.par=c(1,1,1), 
   covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE))

opt<-list(save=FALSE,save.graphs=FALSE)

theo.fit<-saemix(theo.model,theo.data,opt)
theo.fit<-llgq.saemix(theo.fit)

###########################
# Model without covariate

theo.model.base<-saemixModel(model=model1cpt,
   description="One-compartment model with first-order absorption", 
   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, 
   c("ka","V","CL"))),transform.par=c(1,1,1))

opt<-list(save=FALSE,save.graphs=FALSE)

theo.base<-saemix(theo.model.base,theo.data,opt)
theo.base<-llgq.saemix(theo.base)

ll1<-theo.base["results"]["ll.lin"]*(-2)
ll2<-theo.fit["results"]["ll.lin"]*(-2)
1-pchisq(ll1-ll2,1)

ll1<-theo.base["results"]["ll.is"]*(-2)
ll2<-theo.fit["results"]["ll.is"]*(-2)
1-pchisq(ll1-ll2,1)

ll1<-theo.base["results"]["ll.gq"]*(-2)
ll2<-theo.fit["results"]["ll.gq"]*(-2)
1-pchisq(ll1-ll2,1)

theo.fit["results"]["ll.is"]

######################################
# Diagnostic plots

# Plotting individual fits with selected options
par(mfrow=c(2,2))
plot(theo.fit,plot.type="individual.fit",new=FALSE,ilist=1:4,smooth=TRUE,ylog=T, pch=1, col="Blue",xlab="Time in hr",ylab="Theophylline concentrations (mg/L)")

# Plots of the observations versus predictions
plot(theo.fit, plot.type="observations.vs.predictions")

# Scatterplots and distribution of residuals
plot(theo.fit, plot.type="npde")

# VPC
plot(theo.fit, plot.type="vpc")

# Scatterplots and distribution of residuals
plot(theo.fit, plot.type="npde")

