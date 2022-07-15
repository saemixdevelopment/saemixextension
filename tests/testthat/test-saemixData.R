test_that("saemixData can read from a file and from a data frame passed to it", {
  data(theo.saemix)
  expect_output(
    {
      path <- tempfile(fileext = ".tab")
      write.table(theo.saemix, file = path,
        quote = FALSE, sep = " ", row.names = FALSE)
      saemixData(path,
        header=TRUE,sep=" ",na=NA,
        name.group=c("Id"),name.predictors=c("Dose","Time"),
        name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
        units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
    },
    "Reading.*successfully")
  expect_output(
    {
      saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA,
        name.group=c("Id"),name.predictors=c("Dose","Time"),
        name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
        units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
    },
    "successfully")

})
