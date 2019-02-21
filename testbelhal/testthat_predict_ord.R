
###################################################################################
# Using predict to return a vector of predictions

context("Testing predict.newdata for a likelihood model \n")

### FOR CAT DATA
psiM<-data.frame(id = seq(1,1000,1), alp1=seq(2.308,3,1),alp2 = seq(0.716,1,1),alp3 = seq(0.762,1,1))
logp<-saemixObject["model"]["model"](psiM[,2:4], saemix.fit@data@data[,c("ID")], saemix.fit@data@data[,c("TIME","Y")])
