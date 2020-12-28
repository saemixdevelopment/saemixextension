library(saemix)

theo_data <- saemixData(name.data = theo.saemix,
  header = TRUE, sep = " ", na = NA,
  name.group = c("Id"), name.predictors = c("Dose", "Time"),
  name.response = c("Concentration"), name.covariates = c("Weight","Sex"),
  units = list(x = "hr", y = "mg/L", covariates = c("kg", "-")),
  name.X = "Time")

model1cpt<-function(psi, id, xidep) {
  dose <- xidep[, 1]
  tim  <- xidep[, 2]
  ka   <- psi[id, 1]
  V    <- psi[id, 2]
  CL   <- psi[id, 3]
  k    <- CL/V
  ypred <- dose * ka/(V * (ka - k)) * (exp(- k * tim) - exp(- ka * tim))
  return(ypred)
}

model_1cpt <- saemixModel(model = model1cpt,
  description = "One-compartment model with first-order absorption, no covariate",
  psi0 = matrix(c(1., 20, 0.5), ncol = 3, byrow = TRUE,
    dimnames = list(NULL, c("ka", "V", "CL"))),
  transform.par = c(1, 1, 1))

saemix.options <- list(seed = 123456, save = FALSE, save.graphs = FALSE)

theo_fit_1 <-saemix(model_1cpt, theo_data, saemix.options)
