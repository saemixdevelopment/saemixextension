

####################################################################################
####			               createSimulations                                			####
####################################################################################

#' Function to create individuals from a dataset with the package \code{ospsuite}
#' 
#' 
#' @param data a dataframe or a dataset name containing information of the individuals. 
#' @param pkml.file link to the pkml file
#' @param parameter vector of parameters to estimate
#' @param individuals list of individuals. If missing the individual in the simulation file will be used for all simulations
#' @param timeFactor factor to modify the time units for unit conversion in minutes
#' @param individual.values specific parameter values that can be set up in the simulation
#' @param dose dose value
#' @param indexTimeRemove index of time to be removed (it can be set up to 1 to remove the time=0 simulated by default by PK-Sim if not present in the datset)
#' @param index vector of ID index that should be simulated. If NULL all IDs in the dataset will be simulated
#' @param id vector of ID that should be simulated. If NULL all IDs in the dataset will be simulated
#' @return A list of individuals for \code{ospsuite}.
#' @author Donato Teutonico.
#' 
#' @examples
#' 
#' myInd <- create_individuals(dataset = d, fromFile=FALSE)
#
#'
#' @export createSimulations                                
#' 
#' 

createSimulations <- function(data, pkml.file, parameter=NULL, individuals=NULL, timeFactor=1, pkmlFromFile=FALSE,
                              individual.values=NULL, dose=NULL, indexTimeRemove=NULL, index=NULL, id=NULL) {
  cat("\n // Create simulations //\n")
  
  names(data) <- tolower(names(data))
  data <- data %>% mutate(id=factor(id))
  
  # ------------------------------------
  # -----  Select id's -----------------
  list.id <- levels(data$id)
  list.index <- seq_len(length(list.id))
  if (!is.null(id))  list.index <- sort(list.index[match(id, list.id)])
  else if (!is.null(index)) list.index <- index
  
  list.id <- list.id[list.index]
  individuals <- individuals[list.index]
  
  # -----  iterate over the id's ----------------------
  list.sim <- list()
  pb <- progress::progress_bar$new(format = "   [:bar] :percent // id :what  ", total = length(list.id), width = 60)
  pb$tick(0)
  Sys.sleep(0.5)
  init.value <- NULL
  for (i in seq_along(list.id)) {
    pb$tick(tokens = list(what = i))
    
    # -----  extract individual data -----------------
    fi <- data %>% dplyr::filter(id==list.id[i]) %>% droplevels()
    
    # -----  create simulation with individual covariates -----------------
    # -----  If individuals= NULL the simulation vector is created cloning the input simulation 
    if (pkmlFromFile) pkml.file <- unique(fi$pkml.file)
    sim <- loadSimulation(pkml.file)
    indiv <- individuals[[i]]
    if(!is.null(indiv)){
    setParameterValuesByPath(parameterPaths = indiv$distributedParameters$paths,
                             values = indiv$distributedParameters$values,
                             simulation = sim)
    if (!is.null(individual.values))
      setParameterValuesByPath(parameterPaths = names(individual.values),
                               values = individual.values[i,],
                               simulation = sim)
    }
    # -----  use dose information  -----------------
    if (!is.null(dose)) {
      doseParam <- getParameter(dose$path, sim)
      myDoseBS<- toBaseUnit(quantity = doseParam, values = fi$amt[1], unit = dose$unit)
      setParameterValues(parameters = doseParam, values = myDoseBS)
    }
    
    # -----  use observation times  -----------------
    clearOutputIntervals(simulation = sim)
    time.i <- round(fi$time*timeFactor, 4)
    for (j in seq_along(time.i))
      addOutputInterval(simulation = sim, startTime = time.i[j], endTime = time.i[j], resolution = 0)
    
    # -----  create batch simulation  -----------------
    simBatch <- createSimulationBatch(sim, parametersOrPaths = parameter)
    
    # -----  find output times to remove  -----------------
    if (is.null(indexTimeRemove)) {
      if (is.null(init.value))
        init.value <- unlist(lapply(parameter, function(x) {getParameter(x, sim)$value}))
      ids <- simBatch$addRunValues(parameterValues = init.value)
      simResults <- runSimulationBatches(simulationBatches = simBatch)
      simulatedTimes <- getOutputValues(simulationResults = simResults[[1]][[1]])$data[["Time"]]
      removeTime <- rep(TRUE, length(simulatedTimes))
      for (j in seq_along(time.i))
        removeTime[which.min(abs(time.i[j]-simulatedTimes))] <- FALSE
      itr <- which(removeTime)
      #simulatedTimes <- round(getOutputValues(simulationResults = simResults[[1]][[1]])$data[["Time"]], 4)
      #itr <- which(is.na(match(simulatedTimes, time.i)))
    } else
      itr <- indexTimeRemove
    list.sim[[i]] <- list(simBatch=simBatch, indexTimeRemove=itr, time=time.i)
  }
  
  # ----  use id names
  names(list.sim) <- list.id
  return(list.sim)
}