context('Check stepwise-type methods via the function step.saemix')


test_that("Errors", {
  Theoph.4 <- subset(Theoph, Subject == 4)
  fm1 <- nls(conc ~ SSfol(Dose, Time, lKe, lKa, lCl),
             data = Theoph.4)  
  expect_error(step.saemix(fm1),"Invalid class for argument saemixObject.")
  expect_error(step.saemix(saemix.fit.nocov),"The saemixData object should contain covariates.")
  expect_error(step.saemix(saemix.fit,direction="stepwise"),"Invalid argument direction.")
  expect_warning(step.saemix(saemix.fit,direction="backward"),"Backward procedures are not recommended if initial complete model is too complex.")
})

