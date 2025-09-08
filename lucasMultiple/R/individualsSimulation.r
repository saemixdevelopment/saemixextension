

####################################################################################
####			               create_individuals                                			####
####################################################################################

#' Function to create individuals from a dataset with the package \code{ospsuite}
#' 
#' 
#' @param dataset a dataframe or a dataset name containing information of the individuals. The dataset should contain information for:
#' \itemize{
#'   \item species (If missing it will be assumed to be \code{Human})
#'   \item population (required if species is \code{Human})
#'   \item gender (required if species is \code{Human})
#'   \item weight
#'   \item age (required if species is \code{Human})
#'   \item height (required if species is \code{Human})
#'   }
#' 
#' @param sep the field separator string. Values within each row of \code{dataset} are separated by this string. This argument is passed to the function \code{read.table}
#' @param weight_unit the unit of the weights contained in \code{dataset}. Default is "kg".
#' @param height_unit the unit of the heights contained in \code{dataset}. Default is "cm".
#' @param age_unit the unit of the ages contained in \code{"dataset"}. Default is "year(s)".
#' @param progress Logical. If \code{TRUE} an output on the console will be provided showing the progress of the simulation.
#' @param fromFile Logical. If \code{TRUE} the function will try to read a dataset file. If \code{FALSE} (default) the function will use the dataframe \code{dataset} 
#' @param seed seed which will be set up for the simulation process
#' @return A list of individuals for \code{ospsuite}.
#' @author Donato Teutonico.
#' 
#' @examples
#' 
#' myInd <- create_individuals(dataset = d, fromFile=FALSE)
#
#'
#' @export create_individuals
#' 
#' 


create_individuals <- function(dataset,
                               sep = ",",
                               weight_unit = "kg",
                               height_unit = "cm",
                               age_unit = "year(s)",
                               progress = TRUE,
                               fromFile = FALSE,
                               seed = NULL) {
  
  if(!is.null(seed)) set.seed(seed)

  if(fromFile) data <- read.table(dataset, header = TRUE, sep = sep)
  else data <- dataset
  header <- colnames(data)
  header <- toupper(header)
  header <- trimws(header)
  
  weight_unit <- tolower(weight_unit)
  height_unit <- tolower(height_unit)
  age_unit <- tolower(age_unit)
  
  # find columns
  id_col <- grepl("^ID$", header)
  species_col <- grepl("^SPECIES$", header)
  weight_col <- grepl("^WEIGHT", header)
  
  isThereAnyHuman <- any(toupper(unique(data[species_col]))=="HUMAN")
  if (sum(species_col)==0) {
    isThereAnyHuman <- TRUE
    warning("Species information was not identified in the dataset, it will assumed to be Human", call. = FALSE, immediate. = TRUE, noBreaks. = FALSE)
  }
  if(isThereAnyHuman){
  age_col <- grepl("^AGE", header)
  sex_col <- grepl("^SEX|^GENDER", header)
  height_col <- grepl("^HEIGHT", header)
  pop_col <- grepl("^POPULATION$", header)
 
  # error handling
  if (sum(id_col) != 1)
    stop("Error identifying ID column", call. = FALSE)
  
  if (sum(age_col) != 1)
    stop("Error identifying AGE column", call. = FALSE)
  
  if (sum(sex_col) != 1)
    stop("Error identifying SEX/GENDER column", call. = FALSE)
  
  if (sum(weight_col) != 1)
    stop("Error identifying WEIGHT column", call. = FALSE)
  
  if (sum(height_col) != 1)
    stop("Error identifying HEIGHT column", call. = FALSE)
  
  if (sum(pop_col) != 1)
    stop("Error identifying POPULATION column", call. = FALSE)
}

  id_col <- which(id_col)
  species_col <- which(species_col)
  weight_col <- which(weight_col)
  
  if(isThereAnyHuman){
  age_col <- which(age_col)
  sex_col <- which(sex_col)
  height_col <- which(height_col)
  pop_col <- which(pop_col)
}
  
  # create individuals
  data <- data %>% 
    group_by_at(id_col) %>% 
    slice(1) %>% 
    ungroup() %>%
    arrange_at(id_col) 
  data <- data %>% mutate(seed=seed)
  
  sex_matches <- toupper(unname(unlist((ospsuite::Gender))))
  pop_matches <- toupper(unname(unlist((ospsuite::HumanPopulation))))
  species_matches <- toupper(unname(unlist(ospsuite::Species)))
  
  error_ids <- c()
  good_ids <- c()
  
  
  .create_ind <- function(x) {
    
    species <- toupper( trimws( x[species_col]) )
    
    if (length(species)==0) species <- "HUMAN" # if species information is missing in the dataset assumes "Human"
    species_id <- which(grepl(species, species_matches))
    weight <- as.numeric( x[weight_col] )
    
    if(species=="HUMAN"){
    age <- as.numeric( x[age_col] )
    id <- as.numeric( x[id_col] )
    height <- as.numeric( x[height_col] )
    sex <- toupper( trimws( x[sex_col] ) )
    pop <- toupper( trimws( x[pop_col]) )
    
    sex_id <- which(grepl(paste0("^", sex, "$"), sex_matches))
    pop_id <- which(grepl(pop, pop_matches))
    
    sex <- ospsuite::Gender[[sex_id]]
    pop <- ospsuite::HumanPopulation[[pop_id]]
    }
    species <- ospsuite::Species[[species_id]]
    
    res <- tryCatch(
      create_individual(species=species,
                        population = pop,
                        gender = sex,
                        age = age,
                        age_unit = age_unit,
                        weight = weight,
                        weight_unit = weight_unit,
                        height = height,
                        height_unit = height_unit,
                        seed = seed),
      error = function(e) { 
        error_ids <- c(error_ids, id)
        return(NULL) 
      }
    )
    
    if (progress)
      pb$tick()
    
    return(res)
  }
  
  if (progress)
    pb <- progress::progress_bar$new(total = nrow(data), 
                                     format = "  Individuals [:bar] :percent",
                                     width = 80)
  pb$tick(0)
  Sys.sleep(0.5)
  
  individuals <- apply(data, 1, .create_ind)
  
  ids <- data[[id_col]]
  if (length(error_ids) > 0) {
    warning(paste("Could not create individuals with ids:", paste(error_ids, collapse = ", ")))
    ids <- setdiff(ids, error_ids)
  }

  names(individuals) <- paste("ID", ids)
  
  return(individuals)
}





create_individual <- function(species = "Human",
                              population = "European_ICRP_2002",
                              gender = "MALE",
                              weight = 73,
                              age = 30,
                              height = 176,
                              gestational_age = 40,
                              molecule_ontogenies = NULL, 
                              weight_unit = "kg",
                              height_unit = "cm",
                              age_unit = "year(s)",
                              gestational_age_unit = "week(s)",
                              seed = NULL) {
  
  ##  Check imputs validity
  if(!species %in% unlist(ospsuite::Species))
    stop("Species input should be among the options in ospsuite::Species", call. = FALSE)
  
  if (species=="Human"){  
    if(!population %in% unlist(ospsuite::HumanPopulation))
      stop("Population input should be among the options in ospsuite::HumanPopulation", call. = FALSE)
    
    if(!gender %in% unlist(ospsuite::Gender))
      stop("Gender input should be among the options in ospsuite::Gender", call. = FALSE)
  }
  
  
  
  if (species=="Human"){
    myInd <- createIndividualCharacteristics(species = species,
                                             population = population,
                                             gender = gender,
                                             weight = weight,
                                             height = height,
                                             age = age,
                                             gestationalAge = gestational_age,
                                             weightUnit  = weight_unit,
                                             heightUnit  = height_unit,
                                             ageUnit = age_unit,
                                             gestationalAgeUnit = gestational_age_unit,
                                             moleculeOntogenies = molecule_ontogenies, 
                                             seed = seed)
  } else {
    myInd <- createIndividualCharacteristics(species = species,
                                             weight = weight,
                                             weightUnit  = weight_unit,
                                             seed = seed)
    
  }
  individual <- createIndividual(individualCharacteristics = myInd)
  return(individual)
}




















