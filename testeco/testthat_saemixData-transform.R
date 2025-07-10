context("Testing covariate transformation - continuous \n")

test_that("Transforming variables - continuous numeric value", {
  expect_equal(transform.numeric(1:10,centering=0),1:10)
  expect_equal(transform.numeric(1:11,centering="mean"),seq(-5,5))
  expect_equal(transform.numeric(1:11,centering="median"),seq(-5,5))
  expect_equal(transform.numeric(1:10,centering=5),seq(-4,5))  
  expect_equal(transform.numeric(1:10,centering=0,transformation=log),log(1:10))
  expect_equal(transform.numeric(1:10,centering=(-5),transformation=log),log(6:15))
})

test_that("Transforming variables - continuous numeric values with NA", {
  vect<-c(1:10)
  vect[7]<-NA
  expect_equal(transform.numeric(vect,centering=0),vect)
  transform.numeric(vect,centering="mean")
})


test_that("Transforming datasets - continuous covariate", {
  cow.saemix<-read.table(file.path(datDir,"cow.saemix.tab"),header=T,na=".")
  tab1<-cow.saemix
  x<-saemixData(name.data=tab1,name.group="cow",name.predictors=c("time"), name.response="weight",name.covariates=c("birthyear","twin","birthrank"), units=list(x="d",y="kg",covariates=c("yr","-","-")),verbose=F)
  expect_is(x, "SaemixData") # tests for particular class
  x2<-transformContCov(x,birthyear,verbose=T)
  expect_equal(summary(x2@data$birthyear.mod),summary(x@data$birthyear-median(x@data$birthyear)))
  x2<-transformContCov(x,birthyear,centering="mean",verbose=T)
  expect_equal(summary(x2@data$birthyear.mod),summary(x@data$birthyear-mean(x@data$birthyear)))
  x2<-transformContCov(x,birthyear,centering=1998,verbose=T)
  expect_equal(summary(x2@data$birthyear.mod),summary(x@data$birthyear-1998))
  x2<-transformContCov(x,birthyear,centering=1900,transformation=function(x) log(x),verbose=T)
  expect_equal(summary(x2@data$birthyear.mod),summary(log(x@data$birthyear)-log(1900)))
})

test_that("Transforming datasets - continuous covariate renamed", {
  cow.saemix<-read.table(file.path(datDir,"cow.saemix.tab"),header=T,na=".")
  tab1<-cow.saemix
  x<-saemixData(name.data=tab1,name.group="cow",name.predictors=c("time"), name.response="weight",name.covariates=c("birthyear","twin","birthrank"), units=list(x="d",y="kg",covariates=c("yr","-","-")),verbose=F)
  expect_is(x, "SaemixData") # tests for particular class
  x2<-transformContCov(x,birthyear,centering=1900,transformation=log,verbose=F, newCovName = "logYear")
  expect_equal(summary(x2@data$logYear),summary(log(x@data$birthyear)-log(1900)))
  expect_equal(x2@trans.cov$logYear$transformation(c(1:10)),log(c(1:10))) # check log-transformation
  x2<-transformContCov(x,birthyear,centering="median", newCovName = "yr.cent",verbose=T)
  expect_equal(summary(x2@data$yr.cent),summary(x@data$birthyear-median(x@data$birthyear)))
  
})

test_that("Transforming datasets - continuous covariate with NA values", {
  cow.saemix<-read.table(file.path(datDir,"cow.saemix.tab"),header=T,na=".")
  tab1<-cow.saemix
  tab1$weight[tab1$weight<300]<-NA
  x<-saemixData(name.data=tab1,name.group="cow",name.predictors=c("time"), name.response="weight",name.covariates=c("birthyear","twin","birthrank","weight"), units=list(x="d",y="kg",covariates=c("yr","-","-")),verbose=F)
  expect_is(x, "SaemixData") # tests for particular class
  x2<-transformContCov(x,weight,centering="median",transformation=log,verbose=F, newCovName = "logWT")
  expect_equal(summary(x2@data$logWT),summary(log(x@data$weight)-log(median(x@data$weight,na.rm=T))))
  expect_equal(x2@trans.cov$logWT$transformation(c(1:10)),log(c(1:10))) # check log-transformation
})


test_that("Transforming variables - categorical covariates", {
  expect_equal(transform.numeric(1:10,centering=0),1:10)
  expect_equal(transform.numeric(1:11,centering="mean"),seq(-5,5))
  expect_equal(transform.numeric(1:11,centering="median"),seq(-5,5))
  expect_equal(transform.numeric(1:10,centering=5),seq(-4,5))  
  expect_equal(transform.numeric(1:10,centering=0,transformation=log),log(1:10))
  expect_equal(transform.numeric(1:10,centering=(-5),transformation=log),log(6:15))
})


test_that("Transforming datasets - categorical covariate", {
  cow.saemix<-read.table(file.path(datDir,"cow.saemix.tab"),header=T,na=".")
  tab1<-cow.saemix
  vec<-ifelse(tab1$twin==1,"Single","Twin")
  tab1$twin<-vec
  vec<-as.factor(tab1$birthrank)
  tab1$birthrank<-as.factor(vec)
  tab1<-tab1
  x<-saemixData(name.data=tab1,name.group="cow",name.predictors=c("time"), name.response="weight",name.covariates=c("birthyear","twin","birthrank"), units=list(x="d",y="kg",covariates=c("yr","-","-")),verbose=F)
  expect_is(x, "SaemixData") # tests for particular class
  x2<-transformCatCov(x, covariate=birthrank, newCat=c(1,2,2,3,3), verbose=TRUE)
  expect_equal(x2@name.covariates[grep("birthrank",x2@name.covariates)],c("birthrank.G2","birthrank.G3"))
  expect_equal(sum(x@data$birthrank==4)+sum(x@data$birthrank==5), sum(x2@data$birthrank.G2==1))
  expect_equal(sum(x@data$birthrank==6)+sum(x@data$birthrank==7), sum(x2@data$birthrank.G3==1))
  x2<-transformCatCov(x, covariate=birthrank, newCat=c(1,2,2,3,3), newCatName=c("ref","preg4-5","6-7"), verbose=TRUE)
  expect_equal(sum(x@data$birthrank==4)+sum(x@data$birthrank==5), sum(x2@data$'preg4-5'==1))
  expect_equal(sum(x@data$birthrank==6)+sum(x@data$birthrank==7), sum(x2@data$'6-7'==1))
  x2<-transformCatCov(x, covariate=birthrank, newCat=c(2,2,1,3,3), newCatName=c("preg3-4","ref5","preg6-7"), verbose=TRUE)
  expect_equal(sum(x@data$birthrank==4)+sum(x@data$birthrank==3), sum(x2@data$'preg3-4'==1))
  expect_equal(sum(x@data$birthrank==6)+sum(x@data$birthrank==7), sum(x2@data$'preg6-7'==1))
})

test_that("Transforming datasets - changing reference for binary covariate", {
  cow.saemix<-read.table(file.path(datDir,"cow.saemix.tab"),header=T,na=".")
  tab1<-cow.saemix
  vec<-ifelse(tab1$twin==1,"Single","Twin")
  tab1$twin<-vec
  vec<-as.factor(tab1$birthrank)
  tab1$birthrank<-as.factor(vec)
  tab1<-tab1
  x<-saemixData(name.data=tab1,name.group="cow",name.predictors=c("time"), name.response="weight",name.covariates=c("birthyear","twin","birthrank"), units=list(x="d",y="kg",covariates=c("yr","-","-")),verbose=F)
  expect_is(x, "SaemixData") # tests for particular class
  x2<-transformCatCov(x, covariate=twin, newCat=c(2,1), newCatName="Singleton", verbose=TRUE)
  expect_equal(sum(x2@data$Singleton),sum(x@ocov$twin=="Single"))
})


context("Testing subsetting\n")

