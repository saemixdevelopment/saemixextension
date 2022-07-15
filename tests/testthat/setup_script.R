library(saemix)

data(theo.saemix)
theo_data <- saemixData(name.data = theo.saemix,
  header = TRUE, sep = " ", na = NA,
  name.group = c("Id"), name.predictors = c("Dose", "Time"),
  name.response = c("Concentration"), name.covariates = c("Weight","Sex"),
  units = list(x = "hr", y = "mg/L", covariates = c("kg", "-")),
  name.X = "Time",
  verbose = FALSE
)

# One-compartment model with first-order absorption
f_1cpt <- function(psi, id, xidep) {
  dose <- xidep[, 1]
  tim  <- xidep[, 2]
  ka   <- psi[id, 1]
  V    <- psi[id, 2]
  CL   <- psi[id, 3]
  k    <- CL/V
  ypred <- dose * ka/(V * (ka - k)) * (exp(- k * tim) - exp(- ka * tim))
  return(ypred)
}

m_1cpt <- saemixModel(model = f_1cpt,
  description = "One-compartment model",
  psi0 = matrix(c(1., 20, 0.5), ncol = 3, byrow = TRUE,
    dimnames = list(NULL, c("ka", "V", "CL"))),
  transform.par = c(1, 1, 1),
  verbose = FALSE
)

m_1cpt_1_cov <- saemixModel(model = f_1cpt,
  description = "One-compartment model, clearance dependent on weight",
  psi0 = matrix(c(1., 20, 0.5, 0.1, 0, -0.01), ncol = 3, byrow = TRUE,
    dimnames = list(NULL, c("ka", "V", "CL"))),
  transform.par = c(1, 1, 1),
  covariate.model = matrix(c(0, 0, 1, 0, 0, 0), ncol = 3, byrow = TRUE),
  verbose = FALSE
)

m_1cpt_2_cov <- saemixModel(model = f_1cpt,
  description = "One-compartment model, clearance dependent on weight and volume dependent on sex",
  psi0 = matrix(c(1., 20, 0.5, 0.1, 0, -0.01), ncol = 3, byrow = TRUE,
    dimnames = list(NULL, c("ka", "V", "CL"))),
  transform.par = c(1, 1, 1),
  covariate.model = matrix(c(0, 0, 1, 0, 1, 0), ncol = 3, byrow = TRUE),
  verbose = FALSE
)

m_1cpt_3_cov <- saemixModel(model = f_1cpt,
  description = "One-compartment model, clearance and absorption dependent on weight, volume dependent on sex",
  psi0 = matrix(c(1., 20, 0.5, 0.1, 0, -0.01), ncol = 3, byrow = TRUE,
    dimnames = list(NULL, c("ka", "V", "CL"))),
  transform.par = c(1, 1, 1),
  covariate.model = matrix(c(1, 0, 1, 0, 1, 0), ncol = 3, byrow = TRUE),
  verbose = FALSE
)

saemix.options <- list(seed = 123456, save = FALSE, save.graphs = FALSE,
  displayProgress = FALSE, print = FALSE)
theo_fit <- saemix(m_1cpt, theo_data, saemix.options)
theo_fit_1_cov <- saemix(m_1cpt_1_cov, theo_data, saemix.options)
theo_fit_2_cov <- saemix(m_1cpt_2_cov, theo_data, saemix.options)
theo_fit_3_cov <- saemix(m_1cpt_3_cov, theo_data, saemix.options)
