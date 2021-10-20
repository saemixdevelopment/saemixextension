context("Model comparison")

test_that("Log-Likelihood, AIC and BIC can be extracted", {

  ll_is <- logLik(theo_fit)
  expect_equal(class(ll_is), "logLik")
  expect_equal(round(as.numeric(ll_is), 2), -172.81)
  expect_equal(attr(ll_is, "df"), 7)

  AIC_is <- AIC(theo_fit)
  expect_equal(round(AIC_is, 2), 359.63)

  BIC_is <- BIC(theo_fit)
  expect_equal(round(BIC_is, 2), 363.02)

  ll_lin <- logLik(theo_fit, method = "lin")
  expect_equal(round(as.numeric(ll_lin), 2), -171.98)

  AIC_lin <- AIC(theo_fit, method = "lin")
  expect_equal(round(AIC_lin, 2), 357.96)

  BIC_lin <- BIC(theo_fit, method = "lin")
  expect_equal(round(BIC_lin, 2), 361.35)

  expect_error(ll_gq <- logLik(theo_fit, method = "gq"), "not yet.*computed")
})

test_that("Model comparison works as expected", {

  # Importance sampling
  expect_message(
    comp_1 <- compare.saemix(theo_fit, theo_fit_1_cov, theo_fit_2_cov, theo_fit_3_cov),
    "by importance sampling")

  # AIC is lowest when using three covariates
  expect_true(all(comp_1[4, "AIC"] < comp_1[1:3, "AIC"]))

  # BIC and BIC.cov is lowest without covariates
  expect_true(all(comp_1[1, "BIC"] < comp_1[2:4, "BIC"]))
  expect_true(all(comp_1[1, "BIC.cov"] < comp_1[2:4, "BIC.cov"]))

  # Linearisation gives the same orders
  expect_message(
    comp_2 <- compare.saemix(theo_fit, theo_fit_1_cov, theo_fit_2_cov, theo_fit_3_cov,
      method = "lin"),
    "by linearisation")
  expect_true(all(comp_2[4, "AIC"] < comp_2[1:3, "AIC"]))
  expect_true(all(comp_2[1, "BIC"] < comp_2[2:4, "BIC"]))
  expect_true(all(comp_2[1, "BIC.cov"] < comp_2[2:4, "BIC.cov"]))

  # Gaussian quadrature also gives the same orders
  expect_error(compare.saemix(theo_fit, theo_fit_1_cov, method = "gq"),
    "not available for model 1")
  theo_fit_gq <- llgq.saemix(theo_fit)

  expect_error(compare.saemix(theo_fit_gq, theo_fit_1_cov, method = "gq"),
    "not available for model 2")

  theo_fit_1_cov_gq <- llgq.saemix(theo_fit_1_cov)

  expect_message(compare.saemix(theo_fit_gq, theo_fit_1_cov_gq, method = "gq"),
    "by Gaussian quadrature")
  theo_fit_2_cov_gq <- llgq.saemix(theo_fit_2_cov)
  theo_fit_3_cov_gq <- llgq.saemix(theo_fit_3_cov)

  expect_message(
    comp_3 <- compare.saemix(theo_fit_gq, theo_fit_1_cov_gq, theo_fit_2_cov_gq, theo_fit_3_cov_gq,
      method = "gq"),
    "by Gaussian quadrature")
  expect_true(all(comp_3[4, "AIC"] < comp_3[1:3, "AIC"]))
  expect_true(all(comp_3[1, "BIC"] < comp_3[2:4, "BIC"]))
  expect_true(all(comp_3[1, "BIC.cov"] < comp_3[2:4, "BIC.cov"]))

})
