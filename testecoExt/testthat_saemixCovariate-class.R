# Check class definitions and initialisation from SaemixCovariate.R

context("Creating covariate model using Class")

test_that("Simple continuous covariate model", {
  x1<-new(Class="SaemixCovariate", name="wt")
  expect_equal(x1@covariate.transform@name, x1@name)
  expect_equal(x1@covariate.transform@name.orig, "wt")
  expect_equal(x1@covariate.transform@type.orig, "continuous")
  expect_equal(x1@type, "continuous")
  expect_equal(x1@covariate.transform@transform.function(1:5), 1:5)
  expect_equal(x1@covariate.transform@centering.function(1:5), 1)
})


test_that("Simple continuous covariate model", {
  covmodel<-new(Class="covmodelCont2Cont", name="lwt", name.orig="wt", transform.function=log, centering.value=60)
  x1<-new(Class="SaemixCovariate", name="lwt", unit="kg", beta=0.75, beta.fix=1, covariate.transform=covmodel)
  expect_equal(x1@covariate.transform@name, x1@name)
  expect_equal(x1@covariate.transform@name, "lwt")
  expect_equal(x1@covariate.transform@name.orig, "wt")
  expect_equal(x1@covariate.transform@type.orig, "continuous")
  expect_equal(x1@type, "continuous")
  expect_equal(x1@beta, 0.75)
  expect_equal(x1@beta.fix, 1)
})

test_that("Simple categorical covariate model", {
  x2<-new(Class="SaemixCovariate", name="pgp", type="categorical")
  expect_equal(x2@covariate.transform@name, x2@name)
  expect_equal(x2@covariate.transform@name.orig, "pgp")
  expect_equal(x2@covariate.transform@type.orig, "categorical")
  expect_equal(length(x2@covariate.transform@ncat), 0)
  expect_equal(x2@type, "categorical")
})

test_that("Categorical covariate model with 3 classes", {
  covmodel<-new(Class="covmodelCat2Cat", name="pgpMut", name.orig="pgp", groups=list(c("CC"), c("CT"),c("TT")))
  x2<-new(Class="SaemixCovariate", name="pgp2", type="categorical", beta=c(0.2, 0.5), covariate.transform=covmodel) 
  expect_equal(x2@covariate.transform@name, x2@name)
  expect_equal(x2@covariate.transform@name.orig, "pgp")
  expect_equal(x2@covariate.transform@type.orig, "categorical")
  expect_equal(x2@covariate.transform@ncat, 3)
  expect_equal(x2@type, "categorical")
  expect_equal(length(x2@beta), 2)
  expect_equal(sum(x2@beta.fix), 0)
  expect_equal(sum(x2@beta), 0.7)
})

context("Creating covariate models using constructor function - all defaults")

test_that("Creating a continuous covariate model - defaults", {
  lwt <- contCov(name="wt")
  expect_equal(lwt@name, "wt")
  expect_is(lwt, "SaemixCovariate")
  expect_equal(lwt@type, "continuous")
  expect_equal(lwt@unit, "")
  expect_equal(lwt@covariate.transform@transform.function(c(1,2)), c(1,2))
  expect_equal(lwt@covariate.transform@centering.function(1:10), 1)
  expect_equal(length(lwt@covariate.transform@centering.value), 0)
  expect_equal(lwt@beta, 0)
  expect_equal(lwt@beta.fix, 0)
})


test_that("Creating a categorical covariate model - defaults", {
  gender <- catCov(name="gender")
  expect_equal(gender@name, "gender")
  expect_is(gender, "SaemixCovariate")
  expect_equal(gender@type, "categorical")
  expect_equal(gender@unit, "")
  expect_equal(gender@beta, 0)
  expect_equal(gender@beta.fix, 0)
  expect_equal(length(gender@covariate.transform@ncat), 0)
  expect_equal(length(gender@covariate.transform@name.cat), 0)
  expect_equal(length(gender@covariate.transform@reference), 0)
})


context("Creating covariate models using constructor function")

# I want to set a full parameter model as eg
# cl=saemixPar(distribution="lognormal", mu.start=10, omega.start=0.5,
#              var.level=c("id","occ"),var.start=c(0.5,0.2), covariance=list(c("CL", "V"), c("CL")),
#              covariate=c(lwt=contCov(name="wt", transform=log, centering=mean, beta=0.75, beta.fix=1), 
#                          gender=catCov(name="sex", reference="F", beta=0.5), 
#                          mutPgp=catCov(name="pgp", groups=list(c("CC"), c("CT","TT")), beta=c(0.2, 0.5))
#              )
# )

# name="pgp2", type="categorical", beta=c(0.2, 0.5), covariate.transform=catCov(name="pgp", groups=list(c("CC"),c("CT"), c("TT"))))


test_that("Creating a continuous covariate model with initial CI", {
  lwt <- contCov(name="wt", transform=log, centering=mean, beta=0.75, beta.fix=1)
  expect_equal(lwt@name, "wt")
  expect_equal(lwt@covariate.transform@transform.function(c(1,2)), c(log(1), log(2)))
  expect_equal(lwt@covariate.transform@centering.function(1:10), mean(1:10))
  expect_equal(length(lwt@covariate.transform@centering.value), 0)
  expect_equal(lwt@beta, 0.75)
  expect_equal(lwt@beta.fix, 1)
})

test_that("Creating a binary covariate model, just name", {
  gender <- binCov(name="sex")
  expect_equal(gender@name, "sex")
  expect_equal(gender@covariate.transform@ncat,2)
  expect_equal(gender@beta, 0)
  expect_equal(gender@beta.fix, 0)
  expect_equal(length(gender@covariate.transform@reference),0)
})


test_that("Creating a binary covariate model with catCov, just name", {
  gender <- catCov(name="sex")
  expect_equal(gender@name, "sex")
  expect_equal(length(gender@covariate.transform@ncat),0)
  expect_equal(gender@beta, 0)
  expect_equal(gender@beta.fix, 0)
  expect_equal(length(gender@covariate.transform@reference),0)
})


test_that("Creating a binary covariate model with initial CI", {
  gender <- catCov(name="sex", reference="F", beta=0.5)
  expect_equal(gender@name, "sex")
  expect_equal(length(gender@covariate.transform@ncat),0)
  expect_equal(gender@beta, 0.5)
  expect_equal(gender@beta.fix, 0)
  expect_equal(gender@covariate.transform@reference,"F")
})


test_that("Creating a binary covariate model with initial CI", {
  gender <- binCov(name="sex", reference="F", beta=0.5)
  expect_equal(gender@name, "sex")
  expect_equal(gender@covariate.transform@ncat,2)
  expect_equal(gender@beta, 0.5)
  expect_equal(gender@beta.fix, 0)
  expect_equal(gender@covariate.transform@reference,"F")
})

test_that("Creating a covariate model with 3 categories", {
  pgp <- catCov(name="pgp", ncat=3)
  expect_equal(pgp@name, "pgp")
  expect_equal(pgp@covariate.transform@ncat,3)
})


test_that("Creating a covariate model with 3 categories", {
  pgp <- catCov(name="pgp", reference="CC", name.cat=c("CC","CT","TT"), beta=c(0.5,0.3))
  expect_equal(pgp@name, "pgp")
  expect_equal(pgp@covariate.transform@ncat,3)
  expect_equal(pgp@covariate.transform@name.cat,c("CC","CT","TT"))
  expect_equal(pgp@beta, c(0.5,0.3))
})


context("Creating a list of covariate models using the constructors")

test_that("List of covariates", {
  lcov<-list(lwt=contCov(name="wt", transform=log, centering=mean), 
             gender=catCov(name="sex", groups=list("F","M")), 
             mutPgp=catCov(name="pgp", groups=list(c("CC"), c("CT","TT"))))
  for(i in 1:length(lcov)) lcov[[i]]@name<-names(lcov)[i]
  expect_equal(length(lcov),3)
  expect_equal(lcov[[2]]@covariate.transform@type.orig, "categorical")
  expect_equal(length(lcov[[2]]@covariate.transform@ncat),1)
  expect_equal(lcov[[2]]@type, "categorical")
  expect_equal(lcov[[2]]@name, "gender")
  expect_equal(lcov[[2]]@covariate.transform@name.orig, "sex")
  expect_equal(lcov[[2]]@covariate.transform@reference, "F")
  expect_equal(lcov[[3]]@covariate.transform@reference, "CC")
})


