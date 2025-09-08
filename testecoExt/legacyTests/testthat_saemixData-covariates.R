context("Environment issues with dataframes")

# From a data frame instead of a file on disk (doesn't work within test_that)
# theo.saemix<<-read.table(file.path(datDir,"theo.saemix.tab"),header=T,na=".")

test_that("Successful creation of a SaemixData object from dataframe object, full specification", {
  theo.saemix<-read.table(file.path(datDir,"theo.saemix.tab"),header=T,na=".")
  x<-saemixData(name.data=theo.saemix, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),units=list(x="hr",y="mg/L"), name.X="Time", verbose=FALSE)
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("Dose","Time"))
  expect_equal(x@name.group,c("Id"))
  expect_equal(x@name.response,"Concentration")
  expect_equal(x@header,TRUE)
  expect_equal(x@na,"NA")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"Time")
  expect_equal(x@N,12)
  expect_equal(x@ntot.obs,120)
})

context("Testing creation of SaemixData with covariates\n")
# From a data frame instead of a file on disk (doesn't work within test_that because of environment)
#theo.saemix<<-read.table(file.path(datDir,"theo.saemix.tab"),header=T,na=".")
#theo.saemix$Sex<<-ifelse(theo.saemix$Sex==1,"M","F")

test_that("Successful creation of a SaemixData object with covariates", {
  theo.saemix<-read.table(file.path(datDir,"theo.saemix.tab"),header=T,na=".")
  x<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L", covariates=c("kg","-")), name.X="Time")
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("Dose","Time"))
  expect_equal(x@name.group,c("Id"))
  expect_equal(x@name.response,"Concentration")
  expect_equal(x@name.covariates,c("Weight","Sex"))
  expect_equal(x@header,TRUE)
  expect_equal(x@sep," ")
  expect_equal(x@na,"NA")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"Time")
  expect_equal(x@N,12)
  expect_equal(x@ntot.obs,120)
  expect_equal(sort(unique(x@data$Sex)),c(0,1))
#  expect_equal(sort(unique(x@ocov$Sex)),c("F","M"))
})

#theo.saemix$Weight[theo.saemix$Id==3]<<-NA
#theo.saemix$Sex<<-ifelse(theo.saemix$Sex==1,"M","F")

test_that("Successful creation of a SaemixData object with covariates and one missing covariate", {
  theo.saemix<-read.table(file.path(datDir,"theo.saemix.tab"),header=T,na=".")
  theo.saemix$Weight[theo.saemix$Id==3]<-NA
  theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
  
  x<-saemixData(name.data=theo.saemix, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),
                name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
  
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("Dose","Time"))
  expect_equal(x@name.group,c("Id"))
  expect_equal(x@name.response,"Concentration")
  expect_equal(x@name.covariates,c("Weight","Sex"))
  expect_equal(x@header,TRUE)
  expect_equal(x@sep,"")
  expect_equal(x@na,"NA")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"Time")
  expect_equal(x@N,12)
  expect_equal(x@ntot.obs,120)
  expect_equal(sort(unique(x@data$Sex)),c(0,1))
  expect_equal(sort(unique(x@ocov$Sex)),c("F","M"))
  expect_equal(length(unique(x@data$Weight[!is.na(x@data$Weight)])),11)
})

context("Testing creation of SaemixData with censored or missing data\n")

test_that("Successful creation of a SaemixData object, presence of censored data", {
  PD1.saemix<-read.table(file.path(datDir,"PD1.saemix.tab"),header=T,na=".")
  tab1<-cbind(PD1.saemix,censored=0)
  tab1$censored[tab1$response<5]<-1
  tab1$response[tab1$response<5]<-5
  x<-saemixData(name.data=tab1,header=T, name.group=c("subject"),name.predictors=c("dose"),name.response=c("response"), name.covariates=c("gender"),name.cens="censored",units=list(x="mg",y="-",covariates="-"),verbose=F)
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("dose"))
  expect_equal(x@name.group,c("subject"))
  expect_equal(x@name.response,"response")
  expect_equal(x@name.covariates,"gender")
  expect_equal(x@header,TRUE)
  expect_equal(x@sep,"")
  expect_equal(x@na,"NA")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"dose")
  expect_equal(x@N,100)
  expect_equal(x@ntot.obs,300)
  expect_equal(sort(unique(x@data$gender)),c(0,1))
  expect_equal(sum(x@data$cens),2)
})

# Reading PD1.saemix and PD2.saemix in global environment (not needed ?)

test_that("Successful creation of a SaemixData object, two responses", {
  PD1.saemix<-read.table(file.path(datDir,"PD1.saemix.tab"),header=T,na=".")
  PD2.saemix<-read.table(file.path(datDir,"PD1.saemix.tab"),header=T,na=".")
  tab3<-rbind(cbind(PD1.saemix,type.rep=1),cbind(PD2.saemix,type.rep=2))
  tab3<-tab3[order(tab3$subject,tab3$dose),]
  
  x<-saemixData(name.data=tab3,header=T, name.group=c("subject"),name.predictors=c("dose"),name.response=c("response"), name.covariates=c("gender"),name.ytype="type.rep",units=list(x="mg",y="-",covariates="-"),verbose=F)
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("dose"))
  expect_equal(x@name.group,c("subject"))
  expect_equal(x@name.response,"response")
  expect_equal(x@name.covariates,"gender")
  expect_equal(x@header,TRUE)
  expect_equal(x@sep,"")
  expect_equal(x@na,"NA")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"dose")
  expect_equal(x@N,100)
  expect_equal(x@ntot.obs,600)
  expect_equal(sort(unique(x@data$gender)),c(0,1))
  expect_equal(sort(unique(x@data$ytype)),c(1,2))
})



test_that("Successful creation of a SaemixData object, missing data", {
  PD1.saemix<-read.table(file.path(datDir,"PD1.saemix.tab"),header=T,na=".")
  PD2.saemix<-read.table(file.path(datDir,"PD1.saemix.tab"),header=T,na=".")
  tab1<-PD1.saemix
  tab2<-PD2.saemix
  tab1<-cbind(tab1,censored=0)
  tab1$censored[tab1$response<5]<-1
  tab1$response[tab1$response<5]<-5
  tab3<-rbind(cbind(tab1,type.rep=1,missing=0),cbind(tab2,censored=0,type.rep=2,missing=0))
  tab3<-tab3[order(tab3$subject,tab3$dose),]
  tab3$missing[c(1,4,7)]<-1
  x<-saemixData(name.data=tab3,name.ytype="type.rep",name.cens="censored",name.mdv="missing", name.covariates=c("gender"), verbose=F)
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("dose"))
  expect_equal(x@name.group,c("subject"))
  expect_equal(x@name.response,"response")
  expect_equal(x@name.covariates,"gender")
  expect_equal(x@header,TRUE)
  expect_equal(x@sep,"")
  expect_equal(x@na,"NA")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"dose")
  expect_equal(x@N,100)
  expect_equal(x@ntot.obs,600)
  expect_equal(sort(unique(x@data$gender)),c(0,1))
  expect_equal(sum(x@data$mdv),3)
  expect_equal(sort(unique(x@data$ytype)),c(1,2))
})


test_that("Successful creation of a SaemixData object, several occasions", {
  PD1.saemix<-read.table(file.path(datDir,"PD1.saemix.tab"),header=T,na=".")
  PD2.saemix<-read.table(file.path(datDir,"PD1.saemix.tab"),header=T,na=".")
  tab1<-PD1.saemix
  tab2<-PD2.saemix
  tab1<-cbind(tab1,occas=1,censored=0)
  tab1$censored[tab1$response<5]<-1
  tab1$response[tab1$response<5]<-5
  tab1$occas[tab1$dose==90]<-2
  tab3<-rbind(cbind(tab1,type.rep=1,missing=0),cbind(tab2,occas=1,censored=0,type.rep=2,missing=0))
  tab3<-tab3[order(tab3$subject,tab3$dose),]
  tab3$missing[c(1,4,7)]<-1
  tab3<-tab3
  x<-saemixData(name.data=tab3,name.ytype="type.rep",name.cens="censored",name.mdv="missing",name.occ="occas", name.covariates=c("gender"), verbose=F)
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("dose"))
  expect_equal(x@name.group,c("subject"))
  expect_equal(x@name.response,"response")
  expect_equal(x@name.covariates,"gender")
  expect_equal(x@header,TRUE)
  expect_equal(x@sep,"")
  expect_equal(x@na,"NA")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"dose")
  expect_equal(x@N,100)
  expect_equal(x@ntot.obs,600)
  expect_equal(sort(unique(x@data$gender)),c(0,1))
  expect_equal(sum(x@data$mdv),3)
  expect_equal(sort(unique(x@data$ytype)),c(1,2))
  expect_equal(sort(unique(x@data$occ)),c(1,2))
})


test_that("Successful creation of a SaemixData object with covariates and genetic covariates", {
  pkpddat<-read.table(file.path(saemixDir,"testeco","pkpd_withcov_full.csv"),header=T,sep='\t')
  x<-saemixData(name.data=pkpddat,verbose=F,name.ytype="ytype", name.covariates = c("wt","age", "crcl", "trt"), name.genetic.covariates = paste("gen",1:9,sep=""))
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("time","dose"))
  expect_equal(x@name.group,c("id"))
  expect_equal(x@name.response,"y")
  expect_equal(x@name.covariates,c("wt","age", "crcl", "trt",paste("gen",1:9,sep="")))
  expect_equal(x@header,TRUE)
  expect_equal(x@sep,"")
  expect_equal(x@na,"NA")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"time")
  expect_equal(x@N,45)
  expect_equal(x@ntot.obs,540)
  expect_equal(sort(unique(as.character(x@data$trt))),c("C","trtA","trtB"))
  expect_equal(sort(unique(as.character(x@data$gen1))),c("AA","AT","TA","TT"))
})

test_that("Successful creation of a SaemixData object, missing values in predictor", {
  PD1.saemix<-read.table(file.path(datDir,"PD1.saemix.tab"),header=T,na=".")
  tab1<-PD1.saemix
  tab1[15,2]<-tab1[41,2]<-tab1[55,2]<-NA
  x<-saemixData(name.data=tab1,header=T, name.group="subject",name.predictors="dose", name.response="response",name.covariates="gender", units=list(x="mg",y="-",covariates="-"),verbose=F)
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("dose"))
  expect_equal(x@name.group,c("subject"))
  expect_equal(x@name.response,"response")
  expect_equal(x@name.covariates,"gender")
  expect_equal(x@header,TRUE)
  expect_equal(x@sep,"")
  expect_equal(x@na,"NA")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"dose")
  expect_equal(x@N,100)
  expect_equal(x@ntot.obs,297)
})

test_that("Successful creation of a SaemixData object, 2 subjects with all values missing", {
  PD1.saemix<-read.table(file.path(datDir,"PD1.saemix.tab"),header=T,na=".")
  tab1<-PD1.saemix
  tab1[7:12,3]<-NA
  x<-saemixData(name.data=tab1,header=T, name.group="subject",name.predictors="dose", name.response="response",name.covariates="gender", units=list(x="mg",y="-",covariates="-"),verbose=F)
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("dose"))
  expect_equal(x@name.group,c("subject"))
  expect_equal(x@name.response,"response")
  expect_equal(x@name.covariates,"gender")
  expect_equal(x@header,TRUE)
  expect_equal(x@sep,"")
  expect_equal(x@na,"NA")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"dose")
  expect_equal(x@N,98)
  expect_equal(x@ntot.obs,294)
})

test_that("Successful creation of a SaemixData object, missing values, should set mdv appropriately", {
  PD1.saemix<-read.table(file.path(datDir,"PD1.saemix.tab"),header=T,na=".")
  tab1<-PD1.saemix
  tab1[15,3]<-tab1[41,3]<-tab1[55,3]<-tab1[53,3]<-NA
  x<-saemixData(name.data=tab1,header=T, name.group="subject",name.predictors="dose", name.response="response",name.covariates="gender", units=list(x="mg",y="-",covariates="-"),verbose=F)
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("dose"))
  expect_equal(x@name.group,c("subject"))
  expect_equal(x@name.response,"response")
  expect_equal(x@name.covariates,"gender")
  expect_equal(x@header,TRUE)
  expect_equal(x@sep,"")
  expect_equal(x@na,"NA")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"dose")
  expect_equal(x@N,100)
  expect_equal(x@ntot.obs,300)
  expect_equal(sum(x@data$mdv),4)
})

test_that("Successful creation of a SaemixData object, missing values for covariate gender", {
  PD1.saemix<-read.table(file.path(datDir,"PD1.saemix.tab"),header=T,na=".")
  tab1<-PD1.saemix
  tab1[15,4]<-tab1[41,4]<-tab1[55,4]<-NA
  x<-saemixData(name.data=tab1,header=T, name.group="subject",name.predictors="dose", name.response="response",name.covariates="gender", units=list(x="mg",y="-",covariates="-"),verbose=F)
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(x@name.predictors,c("dose"))
  expect_equal(x@name.group,c("subject"))
  expect_equal(x@name.response,"response")
  expect_equal(x@name.covariates,"gender")
  expect_equal(x@header,TRUE)
  expect_equal(x@sep,"")
  expect_equal(x@na,"NA")
  expect_equal(x@name.mdv,"mdv")
  expect_equal(x@name.cens,"cens")
  expect_equal(x@name.occ,"occ")
  expect_equal(x@name.ytype,"ytype")
  expect_equal(x@name.X,"dose")
  expect_equal(x@N,100)
  expect_equal(x@ntot.obs,300)
  expect_equal(sort(unique(x@data$gender)),c(0,1))
  expect_equal(sum(is.na(x@data$gender)),3)
})

test_that("Successful creation of a SaemixData object, gender as a character string", {
  PD1.saemix<-read.table(file.path(datDir,"PD1.saemix.tab"),header=T,na=".")
  tab1<-PD1.saemix
  tab1$gender<-ifelse(tab1$gender==1,"W","M")
  tab1<-tab1
  x<-saemixData(name.data=tab1,header=T, name.group="subject",name.predictors="dose", name.response="response",name.covariates="gender", units=list(x="mg",y="-",covariates="-"),verbose=F)
  expect_is(x, "SaemixData") # tests for particular class
  expect_equal(sort(unique(x@data$gender)),c(0,1))
  expect_equal(sort(unique(x@ocov$gender)),c("M","W"))
})

