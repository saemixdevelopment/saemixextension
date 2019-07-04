default_settings <- list(
  agq = list(
    derivative_step = 1E-6,
    adaptive = TRUE,
    adaptation_mode = 'mode',
    scaled_likelihood = TRUE,
    gq.quad_points = 5,
    y_integration.n_samples = 1000,
    y_integration.method = 'mc',
    b_integration.method = 'gq',
    jit = FALSE,
    seed = 12345
  ),
  laplace = list(
    derivative_step = 1E-6,
    adaptive = TRUE,
    adaptation_mode = 'mode',
    scaled_likelihood = TRUE,
    gq.quad_points = 1,
    y_integration.n_samples = 1000,
    y_integration.method = 'mc',
    b_integration.method = 'gq',
    jit = FALSE,
    seed = 12345
  )
)

for(group_name in names(default_settings)){
  settings <- default_settings[[group_name]]
  f <- function(){
    for(name in names(as.list(formals()))) settings[[name]] <- get(name)
    return(settings)
  }
  formals(f) <- as.pairlist(default_settings[[group_name]]) 
  assign(paste0("defaults.", group_name), f)
}

filter_settings <- function(settings, prefix){
  filtered <- settings[grepl(paste0(prefix,"."), names(settings))]
  names(filtered) <- gsub(paste0(prefix,"."),"", names(filtered))
  filtered
}