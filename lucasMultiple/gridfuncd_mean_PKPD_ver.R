
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

gridfuncd <- function(simulations=NULL, data=NULL, output=NULL, psi0=NULL, psi.min=NULL, psi.max=NULL, e.obj=0.05, 
                      h.obj=0.25, rh=0.9, nw0=11, alpha=0.5, N0=6, model=NULL, seed=12345,
                      index=NULL, id=NULL , class=NULL) {
# # 
  # simulations=NULL
  # data= data_sim1_cpt
  # output=NULL
  # psi0= psi_ini
  # psi.min=NULL
  # psi.max=NULL
  # e.obj=0.05
  # h.obj=0.25
  # rh=0.9
  # nw0=11
  # alpha=0.5
  # N0=6
  # model=modelPKPD.cpt
  # seed=12345
  # index=NULL
  # id=NULL
  # class= "ytype"

  model <<- model 
  
  
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
  
  N <- length(unique(data$id))
  if (N0>N) 
    stop(paste0("N0 should be less than or equal to ", N), call. = FALSE)

  if (!is.factor(data$id)){
    data$id = as.factor(data$id)
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
  
  
  #select.sim <- simulations[select.index]
  
  # sim1 <- select.sim[[1]]$simBatch
  # if (!"SimulationBatch" %in% class(sim1))
  #   stop("simulation must be a OSP SimulationBatch", call. = FALSE)
  # if (length(sim1$getVariableParameters())!=length(psi0))
  #   stop("length of psi0 does not match with the number of parameters required for the simulation ", call. = FALSE)
  
  select.data <- data %>% dplyr::filter(id %in% list.id[select.index]) %>% droplevels()
  
  if (length(unique(data[[class]])) < 2){
    y <- select.data %>% select(id, y) %>% mutate(id=factor(id))
  } else {
    y <- select.data %>% select(id, y , class) %>% mutate(id=factor(id))
  }
  
  id = c()
  for (i in 1:length(select.index)){
    id = c(id , rep(i,length(select.data$id[select.data$id == select.index[i]])))
  }
  
  # id = select.data$id si on fait ça les id ne seront pas consécutives 
  xidep = select.data[,c("dose","tim","ytype")]

  d <- length(psi0)
  psi0.names <- names(psi0)
  names(psi0) <- NULL
  
 # assign("psi0.names" , psi0.names , .GlobalEnv)
  
  if (is.null(psi.min))
    psi.min <- 0.05*psi0 
  if (is.null(psi.max))
    psi.max <- 10*psi0
  
  
  for (jp in seq_len(d)) {
    print(jp)
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
      
      h.min <- try(pred.mse.mean(y, list_w.min, id , xidep), silent = T) ##### premier calcul 
      print("min")
      print(h.min)
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
      h.max <- try(pred.mse.mean(y, list_w.max, id , xidep), silent = T) ##### premier calcul 
      print("max")
      print(h.max)
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
  
  
  ##### FIN DE PARTIE POUR TROUVER LES BORNES PSI MIN ET PSI MAX
  
  list_all_positions = vector("list" , d )
  list_all_rmse = vector("list" , d )
  


  
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
    for (i in 1:length(lwa)) {
      if (qwa[i] == TRUE){
        wd <- wm <- psi0
        wd[jp] <- lwa[i]
        wb <- wd
        w.max <- (2^(iterator_k))*lwa[i]
        wb[jp] <- w.max
        wm[jp] <- (lwa[i] + w.max)/2
        
        list_positions[[i]] = append(list_positions[[i]] , list(wb) )
        
        list_position = list(wd,wb,wm)
        prediction_list = list()
        
        for (k in 1:length(list_position)){
          psi =  matrix(rep(list_position[[k]] , length(unique(id))) , ncol = length(list_position[[k]]) , nrow = length(unique(id)) , byrow = TRUE)
          prediction_list = append(prediction_list , list(model(psi,id,xidep)) )
        }
        
        pred0 = prediction_list[[1]]
        predb = prediction_list[[2]]
        predm = prediction_list[[3]]
        
        predm[predm==0] <- 0.000001
        preda <- (pred0 + predb)/2
        preda[preda == 0] = 0.000001 
        e <- max(abs((predm - preda)/predm))
        if (e > e.obj){
          list_rmse[[i]] = append(list_rmse[[i]] , e)
          qwa[i] = FALSE
        }
      } else {
        next
      }
    }
    
    if (!any(qwa)){
      break
    }
  
    index = which(qwa)
    if (list_positions[[index[1]]][[iterator_k]][jp] > psi.max[jp]) {
       list_rmse[[index[1]]] = 0 
       break
    }
    iterator_k = iterator_k + 1
  }
    list_all_positions[[jp]] = list_positions
    list_all_rmse[[jp]] = list_rmse
  }
  
  #FIN DE LA RECHERCHE DE POSITIONS QUI DEPASSENT LES 5% 
  
  optimization_positions = vector("list" , d)
  optimization_rmse = vector("list" , d)
  lw0  = vector("list" , d)
  
  for (jp in 1:d){
    
    lwa <- seq(psi.min[jp]^alpha, psi.max[jp]^alpha, length=nw0)^(1/alpha)
    iterator_i = 1
    
    while(iterator_i < length(lwa)){
      
      p_length = length(list_all_positions[[jp]][[iterator_i]])
      lw0[[jp]] = append(lw0[[jp]] , lwa[iterator_i])
      optimization_positions[[jp]] = append(optimization_positions[[jp]] , list(c(lwa[iterator_i] , list_all_positions[[jp]][[iterator_i]][[p_length]][[jp]])))
      optimization_rmse[[jp]] = append(optimization_rmse[[jp]] , list(list_all_rmse[[jp]][[iterator_i]]))
      
      if (list_all_rmse[[jp]][[iterator_i]] == 0){
        break
      }

      if (p_length > 1){
        index = which(list_all_positions[[jp]][[iterator_i]][[p_length-1]][jp] > lwa)
        pos_inf = which(index > (iterator_i))
        if (length(pos_inf) > 1){
          pos_inf = pos_inf[-length(pos_inf)]
          iterator_i = max(pos_inf)
          next
        }
      }
      iterator_i = iterator_i + 1
    }
  }
  
     # JE VEUX A PARTIR DE CETTE FONCTION OBTENIR LE VECTEUR LW1 PUIS FAIRE L'interpolation après 
 
     lw1 <- my.optimize(d , optimization_positions , lw0 , 
                        optimization_rmse , psi0, xidep , y , id ,e.obj )
     
     
     lw1 = lapply(lw1,function(sous_liste) {unlist(sous_liste)})
 
    x.grid <- vector("list" , d)
    
    for (jp in 1:d){
      lw1_j = lw1[[jp]]
      lw0_j = lw0[[jp]]
      w = lw0_j[1]

      if (length(lw0_j) == 1){
        x.grid[[jp]] = c(lw0_j[1],psi0[jp],lw1_j[1])
        next
      }
      
      while (!is.na(w)) {
        x.grid[[jp]] <- c(x.grid[[jp]], w)
        w <- approx(lw0_j, lw1_j, w)$y
      }
    }
      return(x.grid)
    }
      

my.optimize <- function(d , optimization_positions , lw0 , 
                        optimization_rmse , psi0, xidep , y , id , e.obj ) {
  
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
    ek = lapply(optimization_rmse[[jp]],  function(x) x[1])
   # list_wk[[jp]]
    qwa = rep(TRUE,length(lw0_j))
    ea = lapply(1:length(lw0_j) , function(x) c(0))
    eb = vector("list" , length(lw0_j))
    iterator_k = 1
    
   
    
    while(TRUE){
      print(iterator_k)
      print(qwa)
      for (i in 1:length(lw0_j)) {
        if(ek[[i]] == 0){
          qwa[i] = FALSE
        }
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
          
          list_position = list(matrix_0,matrix_max,matrix_m)
          prediction_list = list()
          
          for (k in 1:length(list_position)){
            psi =  matrix(rep(list_position[[k]] , length(unique(id))) , ncol = length(list_position[[k]]) , nrow = length(unique(id)) , byrow = TRUE)
            prediction_list = append(prediction_list , list(model(psi,id,xidep)) )
          }
          
          pred0 = prediction_list[[1]]
          predb = prediction_list[[2]]
          predm = prediction_list[[3]]
          
          predm[predm==0] <- 0.000001
          preda <- (pred0 + predb)/2
          preda[preda==0] = 0.000001
          e <- max(abs((predm - preda)/predm))
          ek[[i]] = e
          
          if (abs(e - e.obj) < e.tol){
            qwa[i] = FALSE
          }
          
          if ((e < e.obj) & (iterator_k > 30) ){
            qwa[i] = FALSE
          }
        } else {
          next
        }
      }
      iterator_k = iterator_k + 1
      if (!any(qwa)){
        break
      }
      # fin while 
    }
  }
    return(list_wk)
  }
          

    
pred.mse.mean <- function(y, p , id , xidep) {
  vec_h = c()
  for (i in 1:(length(p))){
    psi =  matrix(rep(p[[i]] , length(unique(id))) , ncol = length(p[[i]]) , nrow = length(unique(id)) , byrow = TRUE)
    pred = model(psi,id,xidep)
    pred[is.nan(pred)] = Inf
    sum_square_normalized = mean(sapply(unique(y$ytype) , function(col_pred){
      mean((y[y$ytype == col_pred,"y"] - pred[y$ytype == col_pred])^(2)) / mean(y[y$ytype == col_pred, "y"]^(2))
    }))
    h <- 1 - sum_square_normalized
    vec_h = c(vec_h , h)
  }
  return(vec_h)
}

  