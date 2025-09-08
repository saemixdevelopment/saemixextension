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
  
  # if (!is.null(output) && is.null(output[['scaleFactor']]))
  #   output$scaleFactor <- 1
  
  # ICI SERA A CHANGER DANS LE FUTUR IL AJOUTE NATURELEMENT UN SCALEFACTOR == 1
  
  sim1 <- simulations[[1]]$simBatch
  if (!"SimulationBatch" %in% class(sim1))
    stop("simulations[[i]][['simBatch']] must be a OSP SimulationBatch", call. = FALSE)
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
   
    
    rgrid <- get('rgrid',.GlobalEnv)
    model.interpolate <- get('model.interpolate',.GlobalEnv)
    output <- get('output',.GlobalEnv)
    
    index_batch_a_simuler = c()
    
    alg <- 0.5
    if (is.list(psi) & !is.data.frame(psi))   psi <- psi[[1]]
    if (!is.matrix(psi))   psi=matrix(psi, nrow=1)
    ypred <- NULL
    N <- nrow(psi)
    d <- ncol(psi)
    st <- 0
    
    counter = 0
    
    Liste_id_a_visiter = vector("list",N)
    Liste_batch_a_simuler <- c()
    Liste_poids_id = vector("list",N)
 
    for (k in seq_len(N)) {
      j <- ((unique(id)[k]-1) %% length(rgrid$x)) +1 #### il y aura probablement une erreur plus tard ici 
      psii <- psi[k,]
      #if (model.interpolate) { j'ai enlevé le model interpolate pour le moment
      Wi <- rgrid$x[[j]]
      lWi <- rgrid$b[[j]]
        #fi <- rgrid$f[[j]]
        
      ii <- pii <- vector(mode="list", length=d)
        
        ##### PARTIE RAPIDE SANS APPEL A MODEL BATCH 
        for (l in seq_len(d)) {
          Wil <- Wi[[l]]
          if (length(Wil)>1) {
            i2l <- detect_index(Wil, function(x) x>psii[l])
            
            if (i2l==0) {
              ni <- length(Wil)
              m <- 0
              while (psii[l] > Wil[ni]) {
                if (trans.parj[l]==0) 
                  Wil <- c( Wil, Wil[ni]+(Wil[ni] - Wil[ni-1]))
                else
                  Wil <- c(Wil, ((2*Wil[ni]^alg) - Wil[ni-1]^alg)^(1/alg))
                ni <- length(Wil)
                m <- m+1
              }
              di <- dim(lWi)
              di[l] <- m
              lWi <- abind::abind(lWi, array(0,dim=di), along=l)
              i2l <- length(Wil)
              Wi[[l]] <- Wil
            } else if (i2l==1) { # dans ce cas des lignes vont être ajoutées au début donc ça change toute l'organisation !!!! 
              m <- 0
              while (psii[l] < Wil[1]) {
                if (trans.parj[l]==0) 
                  Wil <- c(Wil[1]-(Wil[2] - Wil[1]), Wil)
                else
                  Wil <- c((Wil[1]^2)/Wil[2], Wil)
                m <- m+1
              }
              if (trans.parj[l]==1 && Wil[1]<=0)  
                Wil[1] <- psii[l]/5
              
              di <- dim(lWi)
              di[l] <- m
              lWi <- abind::abind(array(0,dim=di), lWi, along=l)
              i2l <- 2
              Wi[[l]] <- Wil
              
              if (k %/% length(rgrid$x) > 0){
                for (i in 0:(k%/%length(rgrid$x) - 1)){
                  vec = rep(0 , d)
                  vec[l] = m
                  if (d == 1){
                    vec_add_position = rep(vec,1)
                  } else {
                    vec_add_position = rep(vec,2^(d))
                  }
                  Liste_id_a_visiter[[j+length(rgrid$x)*i]] = Liste_id_a_visiter[[j+length(rgrid$x)*i]] + vec_add_position
                }
              }
            }
            ii[[l]] <- c(i2l-1, i2l)
            pil <- (Wil[i2l] - psii[l])/(Wil[i2l] - Wil[i2l-1])
            pii[[l]] <- c(pil, 1-pil)
          } else {
            ii[[l]] <- pii[[l]] <- 1
          }
        }
        
        #### FIN DE LA PARTIE RAPIDE SANS APPEL A MODEL BATCH 
        ig <- as.matrix(expand.grid(ii))
        colnames(ig) <- NULL
        whi <- apply(as.matrix(expand.grid(pii)), 1, prod)
        
        Liste_poids_id[[k]] = c(whi)
        
        new_value = FALSE
       
        
        for (l in 1:nrow(ig)) {
          igl <- ig[l,]
          iwl <- do.call(`[`, c(list(lWi), igl)) # = lWi[igl[1], igl[2], ...]
          if (iwl == 0) {
            m.psi <- Wi[[1]][igl[1]]
            if (d >= 2) {
              for (m in 2:d){
                m.psi <- c(m.psi, Wi[[m]][igl[m]])
                }
            }
            Liste_id_a_visiter[[k]] <- c(Liste_id_a_visiter[[k]] , igl)
            simulations[[j]]$simBatch$addRunValues(m.psi)
            counter = counter + 1
            new_value = TRUE
            lWi = do.call(`[<-`, c(list(lWi), igl , -1)) 
           
          } else {
            Liste_id_a_visiter[[k]] = c(Liste_id_a_visiter[[k]] , igl)
          }
        }
        
        if (new_value){
           rgrid$b[[j]] = lWi 
           rgrid$x[[j]] = Wi 
           index_batch_a_simuler = c(index_batch_a_simuler , j)
           
        }
    }

    index_batch_a_simuler = unique(index_batch_a_simuler)
    simulationBatches = list()
    numberOfCores = 10
    #counter = 0 
   # all_simResults = list()
    
    simulationRunOptions <- SimulationRunOptions$new(
      numberOfCores,
      checkForNegativeValues = NULL,
      showProgress = TRUE
    )
    
    for (index in index_batch_a_simuler){
        simBatch = simulations[[index]]$simBatch
        simulationBatches =  append(simulationBatches, simBatch)
    }
   
    
    if (!is.null(simulationBatches)){
    
    # Passer l'objet instancié à la fonction runSimulationBatches
    simResults <- runSimulationBatches(
      simulationBatches = simulationBatches,
      simulationRunOptions = simulationRunOptions  # Passer l'instance ici
    )
    }
    
    ypredi = c()
    hash_table <- new.env(hash = TRUE)

    for (k in seq_len(N)) {
      j <- ((unique(id)[k]-1) %% length(rgrid$x)) +1 

      if (!exists(as.character(j),envir = hash_table)){
        hash_table[[as.character(j)]] = 1
        iterator_j = hash_table[[as.character(j)]]
      } else {
        iterator_j = hash_table[[as.character(j)]]
      }
      
      
      col_pour_prediction = c()
      matrice_position_grille <- rgrid$b[[j]]
      copy_matrice =  rgrid$b[[j]]
      predict_grille <- rgrid$f[[j]]
      
      indexTimeRemove <- simulations[[j]]$indexTimeRemove
      
      iterator_i = which(index_batch_a_simuler == j)
      
      for (i  in 0:as.numeric((length(Liste_id_a_visiter[[k]])/d)-1)){
        if (do.call(`[`, c(list(matrice_position_grille), Liste_id_a_visiter[[k]][(i*d+1):((i+1)*d)] )) > 0 ){
            col <- do.call(`[`, c(list(matrice_position_grille), Liste_id_a_visiter[[k]][(i*d+1):((i+1)*d)] ))
        } else {
          check_null = FALSE
          if (iterator_j > length(simResults[[iterator_i]])) {
            cat("Erreur détectée : i =", i, ", j =", j, "\n")
            browser()  # Arrêt uniquement quand l'erreur va se produire
          }
          if (is.null(simResults[[iterator_i]][[iterator_j]])){
            longueur = sum(id == k)
            result = rep(Inf,longueur)
            check_null = TRUE
          }
          
          if (!check_null){
            if (is.null(simResults[[iterator_i]][[iterator_j]])){
              browser()
            }
          simulatedData <- getOutputValues(simulationResults=simResults[[iterator_i]][[iterator_j]])$data
          
          if (is.null(output)){
            simulatedValues <- list(simulatedData[,3])
            if (length(indexTimeRemove)>0) simulatedValues =  simulatedValues[[1]][-indexTimeRemove]
            warning("By default, the first output column has been collected.")
          } else {
            if (length(output) < 2){
              result <- simulatedData[[output[[1]]$name]]*output[[1]]$scaleFactor
              if (length(indexTimeRemove)>0) result <- result[-indexTimeRemove]
            } else {
            simulatedValues = list()
            for (output_element in output){
              simulatedValues[[output_element$class]] = simulatedData[[output_element$name]] * output_element$scaleFactor
              if (length(indexTimeRemove)>0) simulatedValues[[output_element$class]] <- simulatedValues[[output_element$class]][-indexTimeRemove]
            }
            y_pred = as.data.frame(simulatedValues)
            
            result <- as.character(xidep[id == k , "ytype"])
            for(element in unique(xidep$ytype)){
              result[element == result] = y_pred[[element]]
            }
            }
          }
          }
      
          f.m <- as.numeric(result)
          predict_grille <- cbind(predict_grille, f.m)
          col <- ncol(predict_grille)
          matrice_position_grille <- do.call(`[<-`, c(list(matrice_position_grille), Liste_id_a_visiter[[k]][(i*d+1):((i+1)*d)] , col)) 
          iterator_j = iterator_j + 1
        }
        col_pour_prediction = c(col_pour_prediction , col )
      
      }
      
      
      if (iterator_j > hash_table[[as.character(j)]] ){
        rgrid$f[[j]] <- predict_grille
        rgrid$b[[j]] <- matrice_position_grille
        hash_table[[as.character(j)]] = iterator_j
      }
      
      ypredi <- c(ypredi,predict_grille[,col_pour_prediction]%*% Liste_poids_id[[k]])
    }
    
    assign("rgrid",rgrid, .GlobalEnv)
    
    return(ypredi)
    
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





