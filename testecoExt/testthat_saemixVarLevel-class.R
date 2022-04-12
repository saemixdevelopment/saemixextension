context("Creating SaemixVarLevel objects (representing a level of variability) using the class")

test_that("Variability levels", {
  x1<-new(Class="SaemixVarLevel")
  expect_equal(x1@name.level,"iiv")
  expect_equal(x1@variable,"id")
  x2<-new(Class="SaemixVarLevel", name.level="iov", variable="occ")
  expect_equal(x2@name.level,"iov")
  expect_equal(x2@variable,"occ")
  x1<-new(Class="SaemixVarLevel", variable="id", size=3, omega.model=diag(c(1,1,0)))
  expect_equal(x1@name.level,"iiv")
  expect_equal(x1@variable,"id")
  expect_equal(sum(x1@omega.model),2)
})

test_that("Diagonal IIV", {
  x1<-new(Class="SaemixVarLevel", variable="id", size=3)
  expect_equal(x1@name.level,"iiv")
  expect_equal(x1@variable,"id")
  expect_equal(sum(x1@omega.model),3)
  expect_equal(length(x1@omega.tri),6)
  expect_equal(length(x1@index.omega),3)
  expect_equal(length(x1@index.omega.fix),0)
})

test_that("Full var-cov", {
  omega1<-vec2mat(c(1,0.9,0,1,0,0.5))
  x1<-new(Class="SaemixVarLevel", variable="id", size=3, omega=omega1, omega.model=matrix(data=1, nrow=3, ncol=3))
  expect_equal(x1@name.level,"iiv")
  expect_equal(x1@variable,"id")
  expect_equal(sum(x1@omega.model),9)
  expect_equal(length(x1@omega.tri),6)
  expect_equal(length(x1@index.omega),6)
  expect_equal(length(x1@index.omega.fix),0)
})


test_that("Fixing some elements of the variance-covariance matrix", {
  x1<-new(Class="SaemixVarLevel", variable="id", size=3, omega.model.fix=diag(c(1,0,0)))
  expect_equal(x1@name.level,"iiv")
  expect_equal(x1@variable,"id")
  expect_equal(sum(x1@omega),3)
  expect_equal(sum(x1@omega.model),3)
  expect_equal(sum(x1@omega.model.fix),1)
  x1<-saemixVarNames(x1, c("E0","Emax","EC50"))
  showall(x1)
})

##################### Creating a list of variance levels through constructor function

context("Creating SaemixVarLevel objects using the constructor function")

test_that("One level, only size given", {
  x<-saemixVarModel(size=4)
  expect_is(x, "SaemixVarLevel")
  expect_equal(x@name.level,"iiv")
  expect_equal(x@variable,"id")
  expect_equal(x@omega, diag(4))
  expect_equal(x@omega.model, diag(4))
  expect_equal(sum(x@omega.model.fix),0)
})

test_that("One level, only omega given", {
  omega1<-vec2mat(c(1,0.9,0,1,0,0.5))
  x<-saemixVarModel(omega=omega1)
  expect_is(x, "SaemixVarLevel")
  expect_equal(x@name.level,"iiv")
  expect_equal(x@variable,"id")
  expect_equal(sum(x@omega.model),5)
  expect_equal(sum(x@omega.model.fix),0)
})

test_that("One level, size mismatch", {
  omega1<-vec2mat(c(1,0.9,0,1,0,0.5))
  omega.model1<-diag(c(1,1,0,1))
  x<-saemixVarModel(omega=omega1, omega.model=omega.model1, verbose=TRUE)
  expect_is(x, "SaemixVarLevel")
  expect_equal(x@name.level,"iiv")
  expect_equal(x@variable,"id")
  expect_identical(x@omega, omega1)
  expect_equal(sum(x@omega.model),5)
  expect_equal(sum(x@omega.model.fix),0)
})

test_that("One level, only omega given, 2 parameters with IIV", {
  omega1<-vec2mat(c(1,0.9,0,1,0,0))
  x<-saemixVarModel(omega=omega1)
  expect_is(x, "SaemixVarLevel")
  expect_equal(x@name.level,"iiv")
  expect_equal(x@variable,"id")
  expect_equal(sum(x@omega.model),4)
  expect_equal(sum(x@omega.model.fix),0)
})

test_that("Block diagonal var-cov, inferred from omega", {
  omega1<-vec2mat(c(1,0.9,0,1,0,0.5))
  x1<-saemixVarModel(omega=omega1)
  expect_equal(x1@name.level,"iiv")
  expect_equal(x1@variable,"id")
  expect_equal(sum(x1@omega.model),5)
  expect_equal(length(x1@omega.tri),6)
  expect_equal(length(x1@index.omega),4)
  expect_equal(x1@index.omega.var,c(1,4,6))
  expect_equal(x1@index.omega.covar,2)
  expect_equal(length(x1@index.omega.fix),0)
})

test_that("List of two levels, id and occ", {
  omega1<-vec2mat(c(1,0.9,0,1,0,0.5))
  omega2<-vec2mat(c(0.5,0.2,0,0.5,0,0))
  x<-list(saemixVarModel(omega=omega1, variable="id"), saemixVarModel(omega=omega2, variable="occ", verbose=T))
  expect_is(x[[1]], "SaemixVarLevel")
  expect_equal(x[[1]]@name.level,"iiv")
  expect_equal(x[[2]]@name.level,"iov")
  expect_equal(x[[1]]@variable,"id")
  expect_equal(x[[2]]@variable,"occ")
  expect_equal(sum(x[[1]]@omega.model),5)
  expect_equal(sum(x[[2]]@omega.model),4)
  expect_equal(sum(x[[1]]@omega.model.fix),0)
  expect_equal(sum(x[[2]]@omega.model.fix),0)
})

test_that("One level, only omega given, 2 parameters with IIV", {
  omega1<-vec2mat(c(1,0.9,0,1,0,0))
  x<-saemixVarModel(omega=omega1)
  x<-saemixVarNames(x, c("ka","V","CL"))
  expect_is(x, "SaemixVarLevel")
  expect_equal(x@name.level,"iiv")
  expect_equal(x@variable,"id")
  expect_equal(sum(x@omega.model),4)
  expect_equal(sum(x@omega.model.fix),0)
  expect_equal(length(x@omega.names[x@index.omega]),3)
})

test_that("2 parameters with IIV, parameter names given", {
  omega1<-vec2mat(c(0,0,0,1,0.9,1))
  x<-saemixVarModel(omega=omega1, parameter.names=c("ka","V","CL"))
  expect_is(x, "SaemixVarLevel")
  expect_equal(x@name.level,"iiv")
  expect_equal(x@variable,"id")
  expect_equal(sum(x@omega.model),4)
  expect_equal(sum(x@omega.model.fix),0)
  expect_equal(length(x@omega.names[x@index.omega]),3)
})

test_that("saemixVarModel with all defaults except size", {
  x<-saemixVarModel(size=3, omega=NULL, omega.model=NULL, omega.model.fix=NULL, parameter.names=NULL, verbose=FALSE)
  expect_is(x, "SaemixVarLevel")
  expect_equal(x@name.level,"iiv")
  expect_equal(x@variable,"id")
  expect_equal(sum(x@omega),3)
  expect_equal(sum(x@omega.model),3)
  expect_equal(sum(x@omega.model.fix),0)
})

##################### Testing for valid variance-covariance structures
context("Auxiliary functions for SaemixModel objects - resizing matrices \n")

test_that("Adjust matrices", {
  x1<-resizeMatrix(diag(c(1,1,1,0)), diag(c(1,1)))
  expect_equal(dim(x1)[1],4)
  x2<-resizeMatrix(diag(c(1,1,1,0)), diag(c(1,0,1,1,1)))
  expect_equal(dim(x2)[1],4)
  omega1<-matrix(1,ncol=4,nrow=4)
  x3<-resizeMatrix(omega1, diag(c(1,0,1,1,1)))
  expect_equal(sum(x3[4,]),1)
})

context("Auxiliary functions for SaemixModel objects - matrix to vector and back \n")

test_that("Switching from matrix to vector form", {
  omega1<-matrix(c(1:9), nrow=3)
  expect_equal(c(omega1),1:9)
  expect_equal(omega1[lower.tri(omega1, diag=TRUE)], c(1:3,5:6,9))
})

test_that("Switching from vector to matrix and back", {
  xvec<-c(1,0.9,0,1,0,0.5)
  omega1<-vec2mat(xvec)
  expect_equal(omega1[2,1],0.9)
  expect_equal(omega1[3,1],0)
  expect_equal(omega1[lower.tri(omega1,diag=TRUE)],xvec)
})

context("Auxiliary functions for SaemixModel objects - covariance model\n")

test_that("Covariance models - test matrix containing elements other than 0/1", {
  xmat<-matrix(c(1,0,0,2),ncol=2)
  expect_false(validate.covariance.model(xmat))
})

test_that("Covariance models - not a square matrix", {
  xmat<-matrix(c(1,0,0,2,3,4),ncol=2)
  expect_false(validate.covariance.model(xmat))
})

test_that("Covariance models - non symmetrical", {
  xmat<-matrix(c(1,0,1,1),ncol=2)
  expect_false(validate.covariance.model(xmat))
})

# Stuck at testing a matrix is block diagonal
# test_that("Covariance models - non block diagonal", {
#   xmat<-vec2mat(c(1,1,0,1,1,0,0,1,1,1))
#   expect_false(validate.covariance.model(xmat))
# })

test_that("Valid covariance models", {
  xmat<-diag(c(1,1,1,0))
  expect_true(validate.covariance.model(xmat))
  xmat<-vec2mat(c(1,1,0,1,0,1))
  expect_true(validate.covariance.model(xmat))
  xmat<-vec2mat(c(1,1,0,0,1,0,0,1,1,1))
  expect_true(validate.covariance.model(xmat))
})
