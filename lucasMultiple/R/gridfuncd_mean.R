
####################################################################################
####			               gridfuncd                                        			####
####################################################################################

#' Function to create the grid used by saemixPBPK.
#' 
#' 
#' @param simulations a list of simulations.
#' @param data a data frame with the data.
#' @param output a list with fields
#' \itemize{
#'   \item \code{name} path of the output
#'   \item \code{scaleFactor} ratio between the unit of data and the unit of model output (default=1)
#'   }
#' @param psi0 a vector of initial values.
#' @param psi.min a vector of minimum parameter values (default = 0.05*psi0).
#' @param psi.max a vector of maximum parameter values (default = 10*psi0).
#' @param e.obj the maximum relative error between interpolations and model predictions (default=0.05).
#' @param h.obj the minimum RMSE between model predictions and observations (default=0.25). Used to eliminate unrealistic values in the grid.
#' @param rh a coefficient between 0 and 1 (default = 0.9).
#' @param nw0 number of points in the initial grid (default = 11).
#' @param alpha power used to define the initial grid (default = 0.5).
#' @param N0 number of id's used to define the initial grid (default = 6).
#' @param id a vector of selected id's (default = NULL)
#' @param index a vector of selected indexes (default = NULL).
#' @param model model function (default = modelBatch).
#' @param seed seed used to randomly select N0 id's for the initial grid (default = 12345).
#' @return A list of d grids.
#' @author Marc Lavielle.
#' 
#' @examples
#' 
#' x.grid05 <- gridfuncd(simulations=theod1.simulations, data=theo.data, psi0=c(Cl=0.005))
#' x.grid02 <- gridfuncd(simulations=theod2.simulations, data=theo.data, psi0=c(Cl=0.005, lipo=1), nw0=7, N.select=5, e.obj=0.02)
#
#'
#' @export gridfuncd
#' 
#' 

# 
# 
# 
# 
# simulations=itra.simulations
# output=info.output
# data=itra.data_4
# psi0=itra.psi_0
# e.obj = 0.05
# class = "type"
# psi.min=NULL
# psi.max=NULL
# e.obj=0.05
# h.obj=0.25
# rh=0.9
# nw0=11
# alpha=0.5
# N0=6
# model=modelBatch
# seed=12345
# index=NULL
# id=NULL


gridfuncd <- function(simulations=NULL, data=NULL, output=NULL, psi0=NULL, psi.min=NULL, psi.max=NULL, e.obj=0.05, 
                      h.obj=0.25, rh=0.9, nw0=11, alpha=0.5, N0=6, model=modelBatch, seed=12345,
                      index=NULL, id=NULL , class=NULL) {
  
  # func_pred <- function(model , p , sim , output ){
  #   pred_df = model(p,sim,output)
  #   pred_df[,2]  = pred_df[,2] * 10000
  #  return(pred_df) 
  # }
  # 
  # simulations=theophylline.simulations
  # output=info.output
  # data=theo.data
  # psi0=theo.psi0
  # e.obj = 0.05
  # class = NULL
  # psi.min=NULL
  # psi.max=NULL
  # h.obj=0.25
  # rh=0.9
  # nw0=11
  # alpha=0.5
  # N0=6
  # seed=12345
  # model=modelBatch
  # index=NULL
  # id=NULL


  model <<- model 
  
  
  # browser()
  # 
  # if (identical(model,modelBatch)) {
  #   
  # }
  
  names(data) <- tolower(names(data))
  if (!("id" %in% names(data)) )
    stop("'id' column (identifiers) is missing in the data file", call. = FALSE)
  if (!("y" %in% names(data)) )
    stop("'y' column (observations) is missing in the data file", call. = FALSE)
  if (!is.null(psi.min)) {
    if (length(psi.min) != length(psi0))
      stop("psi.min and psi0 should have the same length", call. = FALSE)
    if (min(psi0-psi.min) <= 0)
      stop("psi0 should be greater than psi.min", call. = FALSE)
  }
  if (!is.null(psi.max)) {
    if (length(psi.max) != length(psi0))
      stop("psi.max and psi0 should have the same length", call. = FALSE)
    if (min(psi.max-psi0) <= 0)
      stop("psi.max should be greater than psi0", call. = FALSE)
  }
  if (e.obj<=0 | e.obj>=1) 
    stop("e.obj should be strictly between 0 and 1", call. = FALSE)
  if (h.obj<=0 | h.obj>=1) 
    stop("h.obj should be strictly between 0 and 1", call. = FALSE)
  if (rh<=0 | rh>=1) 
    stop("rh should be strictly between 0 and 1", call. = FALSE)
  if (alpha<=0) 
    stop("alpha should be strictly positive", call. = FALSE)
  
  cat(" // Compute grid //\n")
  set.seed(seed)
  N <- length(simulations)
  if (N0>N) 
    stop(paste0("N0 should be less than or equal to ", N), call. = FALSE)
  
  if (length(output) < 2){
    data <- data %>% select(id, y) %>% mutate(id=factor(id))
  } else {
    data <- data %>% select(id, y , class) %>% mutate(id=factor(id))
  }
 
  
  list.id <- levels(data$id)
  
  if (!is.null(id)) {
    index <- match(id, list.id)
  if (anyNA(index))
    stop("input list of id's does not match with the id's of the data file", call. = FALSE)
    index <- sort(index)
  }
  if (!is.null(index)) {
    if (min(index)<1 | max(index)>N)
      stop(paste0("input list of indexes should be between 0 and ", N), call. = FALSE)
    select.index <- index
  } else
    select.index <- sort(sample(N, N0))
  select.sim <- simulations[select.index]
  
  sim1 <- select.sim[[1]]$simBatch
  if (!"SimulationBatch" %in% class(sim1))
    stop("simulation must be a OSP SimulationBatch", call. = FALSE)
  # if (length(sim1$getVariableParameters())!=length(psi0))
  #   stop("length of psi0 does not match with the number of parameters required for the simulation ", call. = FALSE)
  
  select.data <- data %>% dplyr::filter(id %in% list.id[select.index]) %>% droplevels()
  if (length(output) < 2){
    y <- select.data$y 
  } else {
    y <- select.data[, c("y" , class)] 
  }
 
  d <- length(psi0)
  psi0.names <- names(psi0)
  names(psi0) <- NULL
  
  assign("psi0.names" , psi0.names , .GlobalEnv)
  
  if (is.null(psi.min))
    psi.min <- 0.05*psi0 
  if (is.null(psi.max))
    psi.max <- 10*psi0 
  
  for (jp in seq_len(d)) {
    w.min <- w.max <- psi0
    w.min[jp] <- psi.min[jp]
    list_w.min = list(w.min)
    for (k in 1:50){
      if (k == 50){
        stop("Unable to find a reasonable fit... try another parameter value", call. = FALSE)
      }
      for (i in 1:5){
        position.min = list_w.min[[i]]
        position.min[jp] <- (1-rh)*psi0[jp] + rh*list_w.min[[i]][jp]
        list_w.min = append(list_w.min , list(position.min))
      }
      h.min <- try(pred.mse.mean(y, select.sim, list_w.min, output), silent = T) ##### premier calcul 
      if (detect_index(h.min, function(x) x>h.obj) > 0){
        index = detect_index(h.min, function(x) x > h.obj)
        psi.min[jp] <- list_w.min[[index]][jp]
        break
      } else {
        list_w.min[[length(list_w.min)]][jp] = (1-rh)*psi0[jp] + rh*list_w.min[[length(list_w.min)]][jp]
        list_w.min = list(list_w.min[[length(list_w.min)]])
      }
    }
    
    # PARTIE POUR TROUVER PSI MAX 
    w.max[jp] <- psi.max[jp]
    list_w.max = list(w.max)
    for (k in 1:50){
      if (k == 50){
        stop("Unable to find a reasonable fit... try another parameter value", call. = FALSE)
      }
      for (i in 1:5){
        position.max = list_w.max[[i]]
        position.max[jp] <-  (1-rh)*psi0[jp] + rh*list_w.max[[i]][jp]
        list_w.max = append(list_w.max , list(position.max))
      }
      h.max <- try(pred.mse.mean(y, select.sim, list_w.max, output), silent = T) ##### premier calcul 
      
      if (detect_index(h.max, function(x) x >h.obj) > 0){
        index = detect_index(h.max, function(x) x > h.obj)
        psi.max[jp] <- list_w.max[[index]][jp]
        break
      } else {
        list_w.max[[length(list_w.max)]][jp] =  (1-rh)*psi0[jp] + rh*list_w.max[[length(list_w.max)]][jp]
        list_w.max = list(list_w.max[[length(list_w.max)]])
      }
    }
  }
  
  xg <- list()
  list_all_positions = vector("list" , length(d))
  list_all_rmse = vector("list" , length(d))

  
  for (jp in seq_len(d)) {
    
    lwa <- seq(psi.min[jp]^alpha, psi.max[jp]^alpha, length=nw0)^(1/alpha)
    qwa <- rep(TRUE, nw0)
    lw1 <- lw0 <- NULL
    test.w <- F
    list_positions = vector("list" , length(lwa))
    list_rmse = vector("list" , length(lwa))
    iterator_k = 1
    
    while (TRUE) {
      print(iterator_k)
      list_position_a_explorer = c()
    for (i in 1:length(lwa)) {
      if (qwa[i] == TRUE){
        wd <- wm <- psi0
        wd[jp] <- lwa[i]
        wb <- wd
        w.max <- (2^(iterator_k))*lwa[i]
        wb[jp] <- w.max
        wm[jp] <- (lwa[i] + w.max)/2
        
        list_position_a_explorer = c(list_position_a_explorer , TRUE)
        
        for (sim in select.sim) {
            sim$simBatch$addRunValues(wd)
            sim$simBatch$addRunValues(wb)
            sim$simBatch$addRunValues(wm)
        }
        list_positions[[i]] = append(list_positions[[i]] , list(wb) )
      } else {
        list_position_a_explorer = c(list_position_a_explorer , FALSE)
        next}
    }
      
      simulationBatches = list()
    
      for (sim in select.sim) {
        simulationBatches = append(simulationBatches , sim$simBatch)
      }
      
      simResults = sim_para(simulationBatches)
      iterator_i = 0
      
      for (position in 1:length(list_position_a_explorer)) { # de taille de lwa
        if (list_position_a_explorer[position] == TRUE){ # si c'est égale TRUE , il y a eu un calcul à cette position 
      
          pred0_df = data.frame()
          predm_df = data.frame() # initialisation de mes dataframes. 
          predb_df = data.frame()
          prediction_list = list(pred0_df,predb_df,predm_df) #je les mets dans une liste
          
          for (index_result in 1:length(simResults)){ #pour chaque simulation dans simResults (nb sim = nb d'individus )
             for (j in (3*iterator_i+1):(3*iterator_i+3)){
      
               simulatedValues = get_sim_result(index_result,j,select.sim,simResults,output)
               
               if (j %% 3 == 1){
                 prediction_list[[1]] = rbind(prediction_list[[1]] , as.data.frame(simulatedValues,check.names = FALSE))
               }
               if (j %% 3 == 2){
                 prediction_list[[2]] = rbind(prediction_list[[2]] , as.data.frame(simulatedValues,check.names = FALSE))
               }
               if (j %% 3 == 0){
                 prediction_list[[3]] = rbind(prediction_list[[3]] , as.data.frame(simulatedValues,check.names = FALSE))
               }
             }
          }
          
        pred0 = prediction_list[[1]]
        predb = prediction_list[[2]]
        predm = prediction_list[[3]]
        
        predm[predm==0] <- 0.000001
        preda <- (pred0 + predb)/2
        e <- max(abs((predm - preda)/predm))
        
        if (e > e.obj){
          list_rmse[[position]] = append(list_rmse[[position]] , e)
          qwa[position] = FALSE
        }
        iterator_i = iterator_i + 1
        }
      }
      
     iterator_k = iterator_k + 1
     if (!any(qwa)){
       break
     }
     
    }
    list_all_positions[[jp]] = list_positions
    list_all_rmse[[jp]] = list_rmse
  }
  
  
  
  optimization_positions = vector("list" , d)
  optimization_rmse = vector("list" , d)
  lw0  = vector("list" , d)
  
  for (jp in 1:d){
    lwa <- seq(psi.min[jp]^alpha, psi.max[jp]^alpha, length=nw0)^(1/alpha)
    lw0[[jp]] = append(lw0[[jp]] , lwa[1])
    p_length = length(list_all_positions[[jp]][[1]])
    optimization_positions[[jp]] = append(optimization_positions[[jp]] , list(c(lwa[1] , list_all_positions[[jp]][[1]][[p_length]][[jp]])))
    optimization_rmse[[jp]] = append(optimization_rmse[[jp]] , list(list_all_rmse[[jp]][[1]]))
    
    for (i in 3:(length(lwa)+1)) {
      if (i < (length(lwa)-1) ){
      p_length = length(list_all_positions[[jp]][[i-2]])
      if (p_length > 1){
        if (list_all_positions[[jp]][[i-2]][[p_length-1]][jp] > lwa[i]){
          print(i)
          next
        }
      }
      }
      lw0[[jp]] = append(lw0[[jp]] , lwa[i-1])
      p_length = length(list_all_positions[[jp]][[i-1]])
      optimization_positions[[jp]] = append(optimization_positions[[jp]] , list( c(lwa[i-1] , list_all_positions[[jp]][[i-1]][[p_length]][[jp]])))
      optimization_rmse[[jp]] = append(optimization_rmse[[jp]] , list(list_all_rmse[[jp]][[i-1]]))
    }
  }
  
  
     # JE VEUX A PARTIR DE CETTE FONCTION OBTENIR LE VECTEUR LW1 PUIS FAIRE L'interpolation après 
     lw1 <- my.optimize(select.sim, d , optimization_positions , e.obj, output , lw0 , 
                        optimization_rmse , psi0)
     lw1 = lapply(lw1,function(sous_liste) {unlist(sous_liste)})
 
    x.grid <- vector("list" , d)
    
    for (jp in 1:d){
      
      lw1_j = lw1[[jp]]
      lw0_j = lw0[[jp]]
      w = lw0_j[1]
      
      while (!is.na(w)) {
        x.grid[[jp]] <- c(x.grid[[jp]], w)
        w <- approx(lw0_j, lw1_j, w)$y
      }
    }
      return(x.grid)
    }
      
    
    


test.grid <- function(w0, w1, select.sim, output) {
  K <- 11
  wx <- seq(w0, w1, length=11)
  pred <- NULL
  for (w in wx) {
    predw <- NULL
    for (sim in select.sim) 
      predw <- c(predw, model(w, sim, output))
    pred <- cbind(pred, predw)
  }
  app <- NULL
  for (k in 1:K) {
    app <- cbind(app, ((K-k)*pred[,1]+(k-1)*pred[,K])/(K-1) )
  }
  e <- apply((abs((pred-app)/app)), MARGIN = 2, max)
  return(data.frame(x=wx, e=e))
}



my.optimize <- function(select.sim, d , optimization_positions , e.obj, output , lw0 , 
                        optimization_rmse , psi0) {
  e.tol <- e.obj/10
  de <- 0.3
  list_wk = vector("list" , d)
  
  for (para_j in 1:d){
    for (para_i in 1:length(optimization_positions[[para_j]])){
      list_wk[[para_j]] = append(list_wk[[para_j]],list(optimization_positions[[para_j]][[para_i]][2]))
    }
  }
  
  for (jp in 1:d){

    lw0_j = lw0[[jp]]
    opti_pos = optimization_positions[[jp]]
    ek = optimization_rmse[[jp]]
    list_wk[[jp]]
    qwa = rep(TRUE,length(lw0_j))
    ea = lapply(1:length(lw0_j) , function(x) c(0))
    eb = vector("list" , length(lw0_j))
    iterator_k = 1
    
    while(TRUE){
      
      print(iterator_k)
      list_position_a_explorer = c()
      print(qwa)
      
      for (i in 1:length(lw0_j)) {
        if (qwa[i] == TRUE){
          
          if (ek[[i]] > e.obj){
            opti_pos[[i]][2] = list_wk[[jp]][[i]]
            eb[[i]] = ek[[i]]
          } else {
            opti_pos[[i]][1] = list_wk[[jp]][[i]]
            ea[[i]] = ek[[i]]
          }
          
          pa = ((eb[[i]] - e.obj)^(de)) / (((e.obj-ea[[i]])^(de)) + ((eb[[i]]-e.obj)^(de)))
          w.max = pa * opti_pos[[i]][1] + (1-pa)*opti_pos[[i]][2]
          
          list_wk[[jp]][[i]]= w.max
          
          matrix_0 <- matrix_max <- matrix_m <- psi0
          matrix_0[jp] = lw0_j[[i]]
          matrix_max[jp] = w.max
          matrix_m[jp] = (lw0_j[[i]] + w.max)/ 2
          
          list_position_a_explorer = c(list_position_a_explorer , TRUE)
          
          for (sim in select.sim) {
            sim$simBatch$addRunValues(matrix_0)
            sim$simBatch$addRunValues(matrix_max)
            sim$simBatch$addRunValues(matrix_m)
          }
        } else {
          list_position_a_explorer = c(list_position_a_explorer , FALSE)
          next}
      } 
      
      simulationBatches = list()
      
      for (sim in select.sim) {
        simulationBatches = append(simulationBatches , sim$simBatch)
      }
      
      simResults = sim_para(simulationBatches)
    iterator_i = 0
    
    for (position in 1:length(list_position_a_explorer)) { # de taille de lwa
      if (list_position_a_explorer[position] == TRUE){ # si c'est égale TRUE , il y a eu un calcul à cette position 
        
        pred0_df = data.frame()
        predm_df = data.frame() # initialisation de mes dataframes. 
        predb_df = data.frame()
        prediction_list = list(pred0_df,predb_df,predm_df) #je les mets dans une liste
        
        for (index_result in 1:length(simResults)){ #pour chaque simulation dans simResults (nb sim = nb d'individus )
          for (j in (3*iterator_i+1):(3*iterator_i+3)){
            
            simulatedValues = get_sim_result(index_result,j,select.sim,simResults,output)
           
            if (j %% 3 == 1){
              prediction_list[[1]] = rbind(prediction_list[[1]] , as.data.frame(simulatedValues,check.names = FALSE))
            }
            if (j %% 3 == 2){
              prediction_list[[2]] = rbind(prediction_list[[2]] , as.data.frame(simulatedValues,check.names = FALSE))
            }
            if (j %% 3 == 0){
              prediction_list[[3]] = rbind(prediction_list[[3]] , as.data.frame(simulatedValues,check.names = FALSE))
            }
          }
        }
        
        pred0 = prediction_list[[1]]
        predb = prediction_list[[2]]
        predm = prediction_list[[3]]
        
        predm[predm==0] <- 0.000001
        preda <- (pred0 + predb)/2
        e <- max(abs((predm - preda)/predm))
        
        ek[[position]] = e
        
        if (abs(e - e.obj) < e.tol){
          qwa[position] = FALSE
        }
        iterator_i = iterator_i + 1
      }
    }
    iterator_k = iterator_k + 1
    if (!any(qwa)){
      break
    }
}
  }
  return (list_wk)
}



sim_para <- function(simulationBatches){
  numberOfCores = 10
  
  simulationRunOptions <- SimulationRunOptions$new(
    numberOfCores,
    checkForNegativeValues = NULL,
    showProgress = NULL
  )
  
  simResults <- runSimulationBatches(
    simulationBatches = simulationBatches,
    simulationRunOptions = simulationRunOptions 
  )
  
  return(simResults)
  
}

get_sim_result <- function(index_indi,index_simu,list.sim,simResults,output){
  
  indexTimeRemove = list.sim[[index_indi]]$indexTimeRemove
  
  simulatedData <- getOutputValues(simulationResults=simResults[[index_indi]][[index_simu]])$data
  
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
  
  return(simulatedValues)
}


pred.mse.mean <- function(y, list.sim, p, output) {

  
  simulationBatches = list()
  
  for (sim in list.sim) {
    for (position in p){
      sim$simBatch$addRunValues(position)
    }
    simulationBatches = append(simulationBatches , sim$simBatch)
  }
  
  simResults = sim_para(simulationBatches)
  
  list_result = list()
  
  for (i in 1:(length(p))){
    pred_df = data.frame()  
    for (index_result in 1:length(simResults)){
      
      simulatedValues = get_sim_result(index_result,i,list.sim,simResults,output)
      
      pred_df = rbind(pred_df , as.data.frame(simulatedValues,check.names = FALSE))
    }
    list_result = append(list_result , list(pred_df))
  }
  vec_h = c()
  
  for (col in list_result){
    pred_df = col 

    if (length(col) < 2){
      sum_square_normalized = mean(sapply(pred_df , function(x){
        mean((y - x)^(2)) / mean(y^(2))
      }))
    } else {
      sum_square_normalized = mean(sapply(colnames(pred_df) , function(col_pred){
        mean((y[y$type == col_pred,"y"] - pred_df[,col_pred])^(2)) / mean(y[y$type == col_pred, "y"]^(2))
      }))
    }
    h <- 1 - sum_square_normalized
    vec_h = c(vec_h , h)
  }
    
  return(vec_h)
}

  