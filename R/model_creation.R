#' Create model function pointer
#'
#' Creaters a function pointer for \code{\link{model_func_isolated}} or \code{\link{model_func_additive}} to save passing multiple vectors constantly. This function will return the same thing, but requires only a vector of (unnamed) parameters and a vector of times providing the parameter vector is the same order as in the given parTab argument.
#' @param parTab the full parameter table - see example csv file
#' @param form string to indicate if this uses the isolated or competitive version of the model. \code{\link{model_func_isolated}}, \code{\link{model_func_competitive}}
#' @return a function pointer for \code{\link{model_func_isolated}} or \code{\link{model_func_competitive}}
#' @export
#' @useDynLib antibodyKinetics
create_model_func <- function(parTab, form="isolated"){
  strains <- unique(parTab$strain)
  strains <- strains[complete.cases(strains)]
  
  exposures <- parTab[parTab$names =="t_i",]
  exposure_indices <- which(parTab$names =="t_i")
  
  cr_table <- parTab[parTab$names == "x",]
  cr_indices <- which(parTab$names == "x")
  
  order_tab <- parTab[parTab$names == "mod",]
  order_indices <- which(parTab$names == "mod")
  
  parTab1 <- parTab[!(parTab$names %in% c("t_i","x","mod")),]
  parTab_indices <- which(!(parTab$names %in% c("t_i","x","mod")))

  if(form == "isolated") ver <- 0
  if(form == "competitive") ver <- 1

  f <- function(pars,times){
      parTab1$values <- pars[parTab_indices]
      cr_table$values <- pars[cr_indices]
      order_tab$values <- pars[order_indices]
      exposures$values <- pars[exposure_indices]
      y <- model_func(parTab1,cr_table,order_tab,exposures,strains,times,ver)
      return(y)
  }
  return(f)  
}

#' Create model function pointer groups
#'
#' Creaters a function pointer for \code{\link{model_func_groups}} to save passing multiple vectors constantly. This function will return the same thing, but requires only a vector of (unnamed) parameters and a vector of times providing the parameter vector is the same order as in the given parTab argument.
#' @param parTab the full parameter table - see example csv file
#' @param form a string to indicate the form of the model ("isolated" or "competitive")
#' @return a function pointer for \code{\link{model_func_isolated}} or \code{\link{model_func_competitive}}
#' @export
#' @useDynLib antibodyKinetics
#' @seealso \code{\link{create_model_func}}
create_model_group_func <- function(parTab, form = "isolated"){
  strains <- unique(parTab$strain)
  strains <- strains[complete.cases(strains)]
  
  exposures <- parTab[parTab$names =="t_i",]
  exposure_indices <- which(parTab$names =="t_i")
  
  cr_table <- parTab[parTab$names == "x",]
  cr_indices <- which(parTab$names == "x")
  
  order_tab <- parTab[parTab$names == "mod",]
  order_indices <- which(parTab$names == "mod")
  
  parTab1 <- parTab[!(parTab$names %in% c("t_i","x","mod")),]
  parTab_indices <- which(!(parTab$names %in% c("t_i","x","mod")))

  f_model <- NULL
  if(form=="isolated") f_model <- model_func_isolated
  if(form=="competitive") f_model <- model_func_competitive

  f <- function(pars,times){
      parTab1$values <- pars[parTab_indices]
      cr_table$values <- pars[cr_indices]
      order_tab$values <- pars[order_indices]
      exposures$values <- pars[exposure_indices]
      y <- model_func_groups(parTab1,cr_table,order_tab,exposures,strains,times, f_model)
      return(y)
  }
  return(f)  
}

#' Create model function pointer cpp implementation
#'
#' Creaters a function pointer for \code{\link{model_func_group_cpp}} to save passing multiple vectors constantly. This function will return the same thing, but requires only a vector of (unnamed) parameters and a vector of times providing the parameter vector is the same order as in the given parTab argument. This function can also be used to create a pointer to the same model function, but solving a likelihood function. This currently uses the "isolated" form of the model. Support could be added to allow either the competitive or isolated form.
#' @param parTab the full parameter table - see example csv file
#' @param dat if posterior function, need the matrix of data. First row is model times, and subsequent rows are trajectories (each row is trajectory of antibodies for one strain, grouped by exposure group)
#' @param PRIOR_FUNC optional pointer to prior calculating function that takes current parameter vector
#' @param version string of either "model" (for pure model function) or "posterior" (for posterior calculation)
#' @param convert_types optionally, a named vector converting strings of infection types to integers. The Cpp funciton needs these as integers, but the default arguments should be fine
#' @param convert_strains as for convert_types, but relating to the infection strain names
#' @param convert_groups if the groups are named, used to convert names to integers
#' @param individuals vector indicating how many individuals are in each group ie. relating to rows in the data matrix
#' @param form string of "isolated" or "competitive" indicating whether we're using the isolation boosting or competitive boosting version of the model
#' @return a function pointer for \code{\link{model_func_group_cpp}} or \code{\link{posterior_func_group_cpp}}
#' @export
#' @useDynLib antibodyKinetics
create_model_group_func_cpp <- function(parTab, dat=NULL, PRIOR_FUNC = NULL,
                                        version="model",
                                        convert_types = c("all"=0,"infection"=1,"vacc"=2,"adj"=3,"mod"=4,"NA"=5),
                                        convert_strains = c("A"=1,"B"=2,"C"=3,"D"=4,"E"=5),
                                        convert_groups = c("1"=1,"2"=2,"3"=3,"4"=4,"5"=5),
                                        individuals = c(1,1,1,1,1),
                                        form = "isolated"
                                        ){
##########################################################
    ## Firstly, we need to isolate group specific exposures
##########################################################
    ## Get unique groups
    groups <- unique(parTab$group)
    groups <- groups[groups != "all"]
    exposure_indices <- NULL
    exposure_i_lengths <- NULL
    
    ## For each group, isolate the parameter table indices that correspond to exposures for
    ## that group. Save these indices in a contiguous vector, and also store the indices
    ## of THIS vector that each group corresponds to eg. first 10 members correspond to
    ## group 1 would give the first entry to the exposure_i_lengths vector of 10
    for(group in groups){
        tmp <- which(parTab$names == "t_i" & parTab$group == group)
        exposure_indices <- c(exposure_indices, tmp)
        exposure_i_lengths <- c(exposure_i_lengths, length(tmp))
    }
    
    ## The length of this vector is the number of groups plus 1
    ## Index starts at 0
    exposure_i_lengths <- c(0,cumsum(exposure_i_lengths))
    
#########################################################
    ## Order modifier parameters
#########################################################
    ## Which indices of parTab correspond to order modifiers?
    order_indices <- which(parTab$names == "mod")
    
#########################################################
    ## Cross reactivity parameters
#########################################################
    strains <- unique(parTab$strain)
    strains <- strains[!is.na(strains)]
    cr_inds <- NULL
    
    ## In blocks of length(strains),
    cr_lengths <- rep(length(strains),length(strains))
    cr_lengths <- c(0,cumsum(cr_lengths))

    ## This is just creating a vectorised matrix, where the size of each "block"
    ## corresponds to the number of strains. Thus, to find the cross reactivity between
    ## say the second and third strain, we can look at the 3rd element of the 2nd block
    for(strain1 in strains){
        for(strain2 in strains){
            ## The parameter table should only have one row entry for each cross reactivity because I'm assuming
            ## symmetric distance
            tmpStrains <- sort(c(strain1,strain2))
            cr_inds <- c(cr_inds,which(parTab$exposure == tmpStrains[1] & parTab$strain == tmpStrains[2] & parTab$names == "x"))
        }
    }
    
#########################################################
    ## Model parameters
#########################################################
    ## Indices of model parameters. 1 = infection, 2 = vaccine, 3 = adjuvant
    param_indices <- which(parTab$index == "parameter")

    ## How many indices correspond to each exposure type?
    par_lengths <- c(
        length(param_indices[which(parTab$type %in% c("all","infection") & parTab$index == "parameter")]),
        length(param_indices[which(parTab$type %in% c("all","vacc") & parTab$index == "parameter")]),
        length(param_indices[which(parTab$type %in% c("all","adj") & parTab$index == "parameter")]))

    ## Contiguous vector of parameters corresponding to infection, vaccine and adjuvanted vaccine
    par_type_ind <- c(param_indices[which(parTab$type %in% c("all","infection") & parTab$index == "parameter")],
                      param_indices[which(parTab$type %in% c("all","vacc") & parTab$index == "parameter")],
                      param_indices[which(parTab$type %in% c("all","adj") & parTab$index == "parameter")])

    ## Add 0 for first index to par_lengths
    par_lengths <- c(0,cumsum(par_lengths))
    
    
#########################################################
    ## Model parameters
#########################################################
    ## Now we can loop through each group, each strain and
    ## each exposure and get the correct model parameters
    pars <- parTab$values

    ## What are the exposure types, strains, order and priming?
    exposure_types <- parTab$type
    exposure_strains <- parTab$exposure
    measured_strains <- parTab$strain
    exposure_orders <- parTab$order
    exposure_primes <- parTab$primed

    convert_types_back <- names(convert_types)
    convert_strains_back <- names(convert_strains)
    

    ## Convert any named elements to integers
    exposure_types <- convert_types[exposure_types]
    exposure_strains <- convert_strains[exposure_strains]
    measured_strains <- convert_strains[measured_strains]
    
    strains <- convert_strains[strains]
    groups <- convert_groups[groups]

    ## Because we are dealing with Cpp where everything is 0 indexed, need to subtract 1 from
    ## the indices calculated in R
    exposure_indices <- exposure_indices - 1
    cr_inds <- cr_inds - 1
    par_type_ind <- par_type_ind - 1
    order_indices <- order_indices - 1
    exposure_i_lengths <- exposure_i_lengths
    par_lengths <- par_lengths
    cr_lengths <- cr_lengths

    f <- NULL
    ## Return the posterior calculation or just model calculation depending on what
    ## is asked for

    ver <- 0
    if(form == "isolated") ver <- 0
    if(form == "competitive") ver <- 1

    
    if(version == "posterior"){
        times <- dat[1,]
        dat <- dat[2:nrow(dat),]
        f <- function(pars){
            ln <- posterior_func_group_cpp(pars, times, groups, individuals, strains,
                                           exposure_types, exposure_strains, measured_strains,
                                           exposure_orders, exposure_primes, 
                                           exposure_indices, cr_inds, par_type_ind, order_indices,
                                           exposure_i_lengths,  par_lengths, cr_lengths,
                                           ver, dat)

            if(!is.null(PRIOR_FUNC)) ln <- ln + PRIOR_FUNC(pars)
            ln
        }
    } else {
        f <- function(pars, times){
            model_func_group_cpp(pars, times, groups, strains,
                                 exposure_types, exposure_strains, measured_strains,
                                 exposure_orders, exposure_primes, 
                                 exposure_indices, cr_inds, par_type_ind, order_indices,
                                 exposure_i_lengths,  par_lengths, cr_lengths,
                                 ver)
        }
    }
    return(f)
}
