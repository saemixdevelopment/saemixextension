context('Check compare.saemix method')


test_that("Errors", {
  saemix.fit<-theo.fit1
  saemix.fit2<-theo.fit2
  expect_error(compare.saemix(saemix.fit),"'mod.list' must be a list.")
  expect_error(compare.saemix(list(saemix.fit)),"'compare.saemix' requires at least two models.")
  coplot(conc ~ Time | Subject, data = Theoph, show.given = FALSE)
  Theoph.4 <- subset(Theoph, Subject == 4)
  fm1 <- nls(conc ~ SSfol(Dose, Time, lKe, lKa, lCl),
             data = Theoph.4)  
  expect_error(compare.saemix(list(saemix.fit,fm1)),"All inputs should have class 'SaemixObject'.")
  expect_error(compare.saemix(list(saemix.fit,binary.fit)),"Compared models should be fitted on the same data.")
})

test_that("Appropriate use of BIC.covariate",{
  tab1 <- compare.saemix(list(saemix.fit, saemix.fit2,saemix.fit3),method="lin")
  expect_equal(dim(tab1)[2],2)
  tab2 <- compare.saemix(list(saemix.fit, saemix.fit2,saemix.fit2bis),method="is")
  expect_equal(dim(tab2)[2],3)
  tab3 <- compare.saemix(list(saemix.fit, saemix.fit2,saemix.fit2bis,saemix.fit5),method="is")
  expect_equal(dim(tab3)[2],2)
  tab4 <- compare.saemix(list(binary.fit, binary.fit2),method="is")
  expect_equal(dim(tab4)[2],3)
})