# Check class definitions and initialisation from SaemixCovariateModel.R

context("Creating default covariate models using the class")

test_that("Parent class", {
  x1<-new(Class="SaemixCovariateTransform", name="wt")
  expect_equal(x1@type, "continuous")
  expect_equal(x1@name, "wt")
  expect_equal(x1@name.orig, "wt")
})

test_that("Cont2cont covariate", {
  x1<-new(Class="covmodelCont2Cont", name="wt")
  expect_equal(x1@type, "continuous")
  expect_equal(x1@name, "wt")
  expect_equal(x1@name.orig, "wt")
})

test_that("Categorical covariate from continuous", {
  x1<-new(Class="covmodelCont2Cat", name="dwt", name.orig="wt")
  expect_equal(x1@name, "dwt")
  expect_equal(x1@name.orig, "wt")
  expect_equal(x1@ncat,2)
  expect_equal(length(x1@name.cat),0)
})

test_that("Categorical covariate from continuous", {
  x1<-new(Class="covmodelCont2Cat", name="dwt", name.orig="wt", ncat=4)
  expect_equal(x1@name, "dwt")
  expect_equal(x1@name.orig, "wt")
  expect_equal(x1@ncat,4)
  expect_equal(length(x1@name.cat),0)
  expect_equal(x1@reference,character())
})

test_that("Categorical covariate", {
  x1<-new(Class="covmodelCat2Cat", name="pgp")
  expect_equal(x1@name, "pgp")
  expect_equal(length(x1@ncat),0)
  expect_equal(length(x1@reference),0)
})

context("Creating complex covariate models using the class")

test_that("Continuous covariate", {
  x1<-new(Class="covmodelCont2Cont", name="wt")
  expect_equal(x1@centering.function(),1)
  x2<-new(Class="covmodelCont2Cont", name="lwt", name.orig="wt", transform.function=log)
  expect_equal(x2@centering.function(),1)
  x3<-new(Class="covmodelCont2Cont", name="lwt", name.orig="wt", transform.function=log, centering.value=60)
  expect_null(x3@centering.function())
})

test_that("Categorical covariate from a continuous covariate", {
  x1<-new(Class="covmodelCont2Cat", name="dwt", name.orig="wt", ncat=3)
  expect_equal(length(x1@breaks),0)
  expect_equal(x1@ncat,3)
  expect_equal(length(x1@name.cat),0)
})

test_that("Categorical covariate in 3 groups", {
  x1<-new(Class="covmodelCat2Cat", name="pgpMut", name.orig="pgp", groups=list(c("CC"), c("CT","TT")))
  expect_equal(length(x1@groups),2)
  expect_equal(x1@ncat,2)
  expect_equal(x1@name.cat,c("CC","CT-TT"))
})


test_that("Categorical covariate regrouped", {
  x1<-new(Class="covmodelCat2Cat", name="pgpMut", name.orig="pgp", groups=list(c("CC"), c("CT","TT")), name.cat=c("wild","mutant"))
  expect_equal(length(x1@groups),2)
  expect_equal(x1@ncat,2)
  expect_equal(x1@name.cat,c("wild","mutant"))
  x2<-new(Class="covmodelCat2Cat", name="pgpMut", name.orig="pgp", groups=list(c("CC"), c("CT","TT")))
  expect_equal(length(x2@groups),2)
  expect_equal(x2@ncat,2)
  expect_equal(x2@name.cat,c("CC","CT-TT"))
})


if(FALSE) { # Intialising to an empty function
  setClass(Class = "myClass",
           representation=representation(
             funtion = "function" # initialize
           ),
           validity=function(object){
             return(TRUE)
           }
  )
  
  x1<-new(Class="myClass")  
  print(x1)
}


