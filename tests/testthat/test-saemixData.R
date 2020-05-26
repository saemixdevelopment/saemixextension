test_that("saemixData can read from a file and from a data frame passed to it", {
  expect_output(
    {
      path <- system.file("testdata/theo.saemix.tab", package = "saemix")
      saemixData(path,
        header=TRUE,sep=" ",na=NA,
        name.group=c("Id"),name.predictors=c("Dose","Time"),
        name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
        units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
    },
    "Reading.*successfully")
  expect_output(
    {
      data(theo.saemix)
      saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA,
        name.group=c("Id"),name.predictors=c("Dose","Time"),
        name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
        units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
    },
    "successfully")

})
