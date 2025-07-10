####################################################################################
####			SaemixModelPBPK class - User-level function			                      ####
####################################################################################

#' Function to create a SaemixModelPBPK object
#' 
#' 
#' @param model name of the function used to compute the structural model (default=modelBatch).
#' @param psi0 a matrix with a number of columns equal to the number d of
#' parameters in the model, and one (when no covariates are available) or two
#' (when covariates enter the model) giving the initial estimates for the fixed
#' effects. The column names of the matrix should be the names of the
#' parameters in the model, and will be used in the plots and the summaries.
#' When only the estimates of the mean parameters are given, psi0 may be a
#' named vector.
#' @param x.grid a list of d grids.
#' @param output a list with fields
#' \itemize{
#'   \item \code{name} path of the output
#'   \item \code{scaleFactor} ratio between the unit of data and the unit of model output (default=1)
#'   }
#' @param file file name with a list rgrid  that contains individual predictions previously computed).
#' @param simulations a list of simulations.
#' @param model.interpolate a boolean, whether interpolation should be used or not (default=T).
#' @param error.model type of residual error model. Valid types are constant,
#' proportional, combined and exponential (default=combined).
#' @param transform.par the distribution for each parameter (0=normal,
#' 1=log-normal, 2=probit, 3=logit). Defaults to a vector of 1s (all parameters
#' have a log-normal distribution)
#' @param fixed.estim whether parameters should be estimated (1) or fixed to
#' their initial estimate (0). Defaults to a vector of 1s
#' @param covariate.model a matrix giving the covariate model. Defaults to no
#' covariate in the model
#' @param covariance.model a square matrix of size equal to the number of
#' parameters in the model, giving the variance-covariance matrix of the model:
#' 1s correspond to estimated variances (in the diagonal) or covariances
#' (off-diagonal elements). Defaults to the identity matrix
#' @param omega.init a square matrix of size equal to the number of parameters
#' in the model, giving the initial estimate for the variance-covariance matrix
#' of the model.
#' @param error.init a vector of size 2 giving the initial value of a and b in
#' the error model. Defaults to 1 for each estimated parameter in the error
#' model
#' @return A SaemixModelPBPK object.
#' @author Marc Lavielle' 
#' @examples
#' 
#' output <- list(name = "Organism|PeripheralVenousBlood|Theophylline|Plasma (Peripheral Venous Blood)",
#'                scaleFactor = 180.17/1000)
#'
#' 
#' saemix.modeld <- saemixModelPBPKd(psi0=p0, 
#'                                   x.grid=x.grid05, 
#'                                   simulations=theophylline.simulations, 
#'                                   output=output)

#' 
#' @export saemixModelPBPK

saemixModelPBPK <-function(psi0=NULL, transform.par=NULL, x.grid=NULL, model=modelBatch,
                           simulations=NULL, output=NULL,
                           fixed.estim=NULL, omega.init=NULL, covariance.model=NULL, covariate.model=NULL,
                           error.model="combined", error.init=NULL, file=NULL, model.interpolate=TRUE,
                           modeltype=NULL , name.response = NULL) {
  
  
# 
# 
#   psi0=psi_ini
#   x.grid=x.grid05
#   simulations=unique(id1)
#   output=NULL
#   covariate.model = NULL
#   covariance.model=diag(1,5,5)
#   error.model = c("constant" , "constant")
#   modeltype = c("structural","structural")
# transform.par=NULL
# model=modelPKPD
#   fixed.estim=NULL
# omega.init=NULL
#  error.init=NULL
# file=NULL
# model.interpolate=TRUE
#  name.response = c("y1","y2")
# 

  
  # if (!is.null(output) && is.null(output[['scaleFactor']]))
  #   output$scaleFactor <- 1
  
  # ICI SERA A CHANGER DANS LE FUTUR IL AJOUTE NATURELEMENT UN SCALEFACTOR == 1
  
  #sim1 <- simulations[[1]]$simBatch
  # if (!"SimulationBatch" %in% class(sim1))
  #   stop("simulations[[i]][['simBatch']] must be a OSP SimulationBatch", call. = FALSE)
  if (is.null(names(psi0)))
    stop("names of psi0 are missing", call. = FALSE)
  if (is.matrix(psi0)) {
    d <- ncol(psi0)
  } else {
    d <- length(psi0)
    psi0 <- matrix(psi0, ncol=d, dimnames=list(NULL, names(psi0)))
  }
  # if (length(sim1$getVariableParameters())!=d)
  #   stop("length of psi0 does not match with the number of parameters required for the simulation ", call. = FALSE)
  if (!is.null(x.grid) && length(x.grid)!=d)
    stop("length of psi0 does not match with the dimension of the initial grid", call. = FALSE)
  if (!is.logical(model.interpolate))
    stop("model.interpolate should be logical", call. = FALSE)
  
  if (is.null(transform.par))
    transform.par <- rep(1, d)
  trans.parj <- transform.par
  if (is.null(fixed.estim))
    fixed.estim <- rep(1, d)
  if (is.null(error.init)) {
    if (any(grepl("combined", error.model))) #### MODIFICATION CAR NE PEUT PAS RECEVOIR PLUSIEURS MODELS D ERREUR 
      error.init <- NULL # error.init <- c(0.1, 0.2)
    else
      error.init <- NULL ### il faudra remettre à 0.3 je pense ??? 
  }
  if (is.null(omega.init))
    omega.init <- diag(rep(1, d))
  if (is.null(covariance.model))
    covariance.model <- diag(rep(1, d))
  
  
  
  N <- length(simulations)
  if (is.null(file) ) {
    if (model.interpolate) {
      xdim <- unlist(lapply(x.grid, length))
      if (d==1 & !is.list(x.grid))  x.grid <- list(x.grid)
      rgrid <- list(x=replicate(N, x.grid, FALSE), b=replicate(N, array(0, xdim), FALSE), f=replicate(N, NULL, FALSE), parameter=colnames(psi0))
    } else {
      rgrid <- FALSE
    }
  } else {
    if (!file.exists(file))
      stop(paste0("file ",file," does not exists"), call. = FALSE)
    load(file)
    if (!exists("rgrid") || length(rgrid) != 4)
      stop(paste0("file ",file," should contain a valid list rgrid obtained from a previous run"), call. = FALSE)
  }
  
  if (is.null(name.response)){
    if (length(output) > 1){
      name.response = vector(length = length(output))
      for (i in 1:length(output)){
        name.response[i] = output[[i]]$class
      }
    } else {
      name.response = ""
    }
  }
  
  assign("rgrid",rgrid, .GlobalEnv)
  assign("model.interpolate",model.interpolate, .GlobalEnv)
  assign("output",output, .GlobalEnv)
  
  
  model.approx <-function(psi, id, xidep) {
   # psi =  psiM
   # id = IdM
   # xidep = XM
    
    #start = Sys.time()

    rgrid <- get('rgrid',.GlobalEnv)
    model.interpolate <- get('model.interpolate',.GlobalEnv)
    output <- get('output',.GlobalEnv)
    
    if (model.interpolate){
      if (!exists("compute_simulation_env", envir = .GlobalEnv)) {
        compute_simulation_env <- new.env()
        compute_simulation_env$d <- length(rgrid$x[[1]])
        compute_simulation_env$xidep <- xidep
        compute_simulation_env$Liste_result_to_cpp <- vector("list", length(rgrid$x))
        
        assign("compute_simulation_env", compute_simulation_env, envir = .GlobalEnv)
      }
      
     # start_1 = Sys.time()
      
      sim_result = main_function(rgrid$b,rgrid$x,rgrid$f,psi,id,compute_simulation,trans.parj)
      
      #end_1 = Sys.time()
      #print(paste("temps C++ : " , end_1-start_1))
      
      rgrid$x = sim_result[[1]]
      rgrid$b = sim_result[[2]]
      rgrid$f = sim_result[[3]]
      
      assign("rgrid",rgrid, .GlobalEnv)
      
      #end = Sys.time()
      
      #print(paste("temps global : " , end - start))
      
      return(sim_result[[4]])
    } else {
      return(model(psi, id, xidep))
    }
  }
 
      
  saemix.model <-saemixModel(model=model.approx,
                             #modelPKSim=model,
                             psi0=psi0,
                             fixed.estim=fixed.estim,
                             transform.par=transform.par,
                             omega.init=omega.init,
                             covariance.model=covariance.model,
                             covariate.model=covariate.model,
                             error.model=error.model,  
                             error.init=error.init,
                             modeltype = modeltype,
                             name.response = name.response)
  
  return(saemix.model)
}




#' Function to run simulation batches with given parameter values
#' 
#' 
#' @param psi a vector or a matrix with a number of elements, resp. with the number of columns, equal to the number d of
#' parameters in the model
#' @param simulation a list with fields
#' \itemize{
#'   \item \code{simBatch} a SimulationBatch
#'   \item \code{indexTimeRemove} indexes of times to remove (optional)
#'   }
#' @param output a list with fields (optional)
#' \itemize{
#'   \item \code{name} path of the output
#'   \item \code{scaleFactor} ratio between the unit of data and the unit of model output (default=1)
#'   }
#' @return A vector of output values.
#' @author Marc Lavielle 
#' @export modelBatch
#' 
modelBatch <- function(psi, simulation, output=NULL) {
  simBatch <- simulation$simBatch
  indexTimeRemove <- simulation$indexTimeRemove
  if (!is.matrix(psi)) psi <- as.matrix(matrix(psi, nrow=1))
  #ypred <- NULL
  for (i in seq_len(nrow(psi))) {
    ids <- simBatch$addRunValues(parameterValues = psi[i,])
    simResults <- runSimulationBatches(simulationBatches = simBatch)
    simulatedData <- getOutputValues(simulationResults = simResults[[1]][[1]])$data
    if (is.null(output)){
      simulatedValues <- list(simulatedData[,3])
      if (length(indexTimeRemove)>0) simulatedValues =  simulatedValues[[1]][-indexTimeRemove]
      warning("By default, the first output column has been collected.")
    } else {
      if (length(output) < 2){
        output = output[[1]]
        simulatedValues <- simulatedData[[output$name]]*output$scaleFactor
        if (length(indexTimeRemove)>0) simulatedValues <- simulatedValues[-indexTimeRemove]
      } else {
      simulatedValues = list()
      for (output_element in output){
        simulatedValues[[output_element$class]] = simulatedData[[output_element$name]] * output_element$scaleFactor
        if (length(indexTimeRemove)>0) simulatedValues[[output_element$class]] <- simulatedValues[[output_element$class]][-indexTimeRemove]
      }
      }
    }
      ypred = as.data.frame(simulatedValues , check.names = FALSE) #potentiellement optionnel , à voir si plus simple avec une liste 
    #ypred <- c(ypred, simulatedValues)
  }
  return(ypred)
}


compute_simulation = function(psi,id){
  d = compute_simulation_env$d 
  xidep = compute_simulation_env$xidep
  Liste_result_to_cpp = compute_simulation_env$Liste_result_to_cpp 

  for (i in 1:length(psi)){
    if (length(psi[[i]]) > 0){
      longueur_sim = length(psi[[i]]) / d
      matrix_i = matrix(psi[[i]] , ncol = d , byrow = TRUE)
      id_i = rep(seq(longueur_sim),each = length(which(id == i)))
      xidep_i <- do.call(rbind, replicate(longueur_sim, xidep[id == i , ],       simplify = FALSE)) 
      Liste_result_to_cpp[[i]] = modelPKPD(matrix_i , id_i ,xidep_i)
    }
  }
  return(Liste_result_to_cpp)
}




