context("Extraction and calculation of diagnostic messages")

test_that("Log-Likelihood, AIC and BIC can be extracted", {

  ll_is <- logLik(theo_fit_1)
  expect_equal(class(ll_is), "logLik")
  expect_equal(round(as.numeric(ll_is), 2), -172.81)
  expect_equal(attr(ll_is, "df"), 7)

  AIC_is <- AIC(theo_fit_1)
  expect_equal(round(AIC_is, 2), 359.63)

  BIC_is <- BIC(theo_fit_1)
  expect_equal(round(BIC_is, 2), 363.02)

  ll_lin <- logLik(theo_fit_1, method = "lin")
  expect_equal(round(as.numeric(ll_lin), 2), -171.98)

  AIC_lin <- AIC(theo_fit_1, method = "lin")
  expect_equal(round(AIC_lin, 2), 357.96)

  BIC_lin <- BIC(theo_fit_1, method = "lin")
  expect_equal(round(BIC_lin, 2), 361.35)

  expect_error(ll_gq <- logLik(theo_fit_1, method = "gq"), "not yet.*computed")
})
