
set_default_options <- function(options){
  if(!exists(".current.options")) .current.options <<- list()
  .current.options<<-add_list_elements(.current.options, options)
  return(NULL)
}

set_options  <- function(options){
  if(!exists(".current.options")) .current.options <<- list()
  .current.options<<-add_list_elements(.current.options, options, replace=T)
  return(NULL)
}

add_list_elements <- function(existing, additional, replace=F){
  for(n in names(additional)){
    if(! n %in% names(existing)) {
      existing[[n]] <- additional[[n]]
    }else{
      if(!is.list(additional[[n]]) & replace) existing[[n]] <- additional[[n]]
      if(is.list(additional[[n]])) {
        existing[[n]] <- add_list_elements(existing[[n]], additional[[n]], replace)
      }
    }
  }
  return(existing)
}

# trace
tr <- function(m) sum(diag(m))