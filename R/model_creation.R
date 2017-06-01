#' Create model function pointer
#'
#' Creaters a function pointer for \code{\link{model_func_isolated}} or \code{\link{model_func_additive}} to save passing multiple vectors constantly. This function will return the same thing, but requires only a vector of (unnamed) parameters and a vector of times providing the parameter vector is the same order as in the given parTab argument.
#' @param parTab the full parameter table - see example csv file
#' @param form string to indicate if this uses the isolated or competitive version of the model. \code{\link{model_func_isolated}}, \code{\link{model_func_competitive}}
#' @return a function pointer for \code{\link{model_func_isolated}} or \code{\link{model_func_competitive}}
#' @export
#' @useDynLib antibodyKinetics
create_model_func <- function(parTab, exposures,form="isolated",cross_reactivity=FALSE,typing=FALSE){
  strains <- unique(parTab$strain)
  strains <- strains[complete.cases(strains)]
  strains <- strains[strains != "all"]
  
 # exposures <- parTab[parTab$names =="t_i",]
#  exposure_indices <- which(parTab$names =="t_i")
  
  cr_table <- parTab[parTab$names == "x",]
  cr_indices <- which(parTab$names == "x")
  
  order_tab <- parTab[parTab$names == "mod",]
  order_indices <- which(parTab$names == "mod")
  
  parTab1 <- parTab[!(parTab$names %in% c("t_i","x","mod")),]
  parTab_indices <- which(!(parTab$names %in% c("t_i","x","mod")))

  if(form == "isolated") ver <- 0
  if(form == "competitive") ver <- 1

  f <- function(pars,times){
      parTab$values <- pars
      cr_table$values <- pars[cr_indices]
      order_tab$values <- pars[order_indices]
      
      y <- model_func(parTab,cr_table,order_tab,exposures,strains,
                      times,ver,cross_reactivity=FALSE,typing=FALSE)
      return(y)
  }
  return(f)  
}

#' Create model function pointer groups
#'
#' Creaters a function pointer for \code{\link{model_func_groups}} to save passing multiple vectors constantly. This function will return the same thing, but requires only a vector of (unnamed) parameters and a vector of times providing the parameter vector is the same order as in the given parTab argument.
#' @param parTab the full parameter table - see example csv file
#' @param form a string to indicate the form of the model ("isolated" or "competitive")
#' @param cross_reactivity if TRUE, then uses cross-reactive boosting rather than ID based boosting
#' @param typing if TRUE, then uses the type specific parameters rather than universal parameters
#' @return a function pointer for \code{\link{model_func_isolated}} or \code{\link{model_func_competitive}}
#' @export
#' @useDynLib antibodyKinetics
#' @seealso \code{\link{create_model_func}}
create_model_group_func <- function(parTab, exposures, form = "isolated",
                                    cross_reactivity=FALSE,typing=FALSE){
  strains <- unique(parTab$strain)
  strains <- strains[complete.cases(strains)]
  strains <- strains[strains != "all"]
  
 #exposures <- parTab[parTab$names =="t_i",]
#  exposure_indices <- which(parTab$names =="t_i")
  
  cr_table <- parTab[parTab$names == "x",]
  cr_indices <- which(parTab$names == "x")
  
  order_tab <- parTab[parTab$names == "mod",]
  order_indices <- which(parTab$names == "mod")
  
  parTab1 <- parTab[!(parTab$names %in% c("t_i","x","mod")),]
  parTab_indices <- which(!(parTab$names %in% c("t_i","x","mod")))

  f_model <- NULL
  if(form=="isolated") ver <- 0
  if(form=="competitive") ver <- 1

  f <- function(pars,times){
      parTab$values <- pars
      cr_table$values <- pars[cr_indices]
      order_tab$values <- pars[order_indices]
      
      y <- model_func_groups(parTab, cr_table, order_tab, exposures,
                             strains, times, ver, cross_reactivity, typing)
      return(y)
  }
  return(f)  
}

#' Create model function pointer cpp implementation
#'
#' Creaters a function pointer for \code{\link{model_func_group_cpp}} to save passing multiple vectors constantly. This function will return the same thing, but requires only a vector of (unnamed) parameters and a vector of times providing the parameter vector is the same order as in the given parTab argument. This function can also be used to create a pointer to the same model function, but solving a likelihood function. This currently uses the "isolated" form of the model. Support could be added to allow either the competitive or isolated form.
#' @param parTab the full parameter table - see example csv file
#' @param exposureTab table of exposure times etcs - see example csv file
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
create_model_group_func_cpp <- function(parTab, exposureTab, 
                                        dat=NULL, PRIOR_FUNC = NULL,
                                        version="model",
                                        convert_types = c("all"=0,"infection"=1,"vacc"=2,"adj"=3,"mod"=4,"NA"=5),
                                        convert_strains = c("A"=1,"B"=2,"C"=3,"D"=4,"E"=5),
                                        convert_groups = c("1"=1,"2"=2,"3"=3,"4"=4,"5"=5),
                                        individuals = c(1,1,1,1,1),
                                        form = "isolated",
                                        typing=FALSE,
                                        cross_reactivity=FALSE
                                        ){
##########################################################
    ## Firstly, we need to isolate group specific exposures
##########################################################
    ## Get unique groups
    groups <- unique(exposureTab$group)
    groups <- groups[groups != "all"]
    groups <- groups[!is.na(groups)]
    
    ## Get unique strains
    strains <- unique(exposureTab$strain)
    strains <- strains[strains != "all"]
    strains <- strains[!is.na(strains)]
    
    exposure_indices <- NULL
    exposure_i_lengths <- NULL
    
    strain_indices <- NULL
    strain_i_lengths <- NULL
    
    
    ## For each group, isolate the exposure table indices that correspond to exposures for
    ## that group. Save these indices in a contiguous vector, and also store the indices
    ## of THIS vector that each group corresponds to eg. first 10 members correspond to
    ## group 1 would give the first entry to the exposure_i_lengths vector of 10
    for(group in groups){
        tmp <- which(exposureTab$group == group)
        exposure_indices <- c(exposure_indices, tmp)
        exposure_i_lengths <- c(exposure_i_lengths, length(tmp))
        ## For each strain
        tmpExposures <- exposureTab[tmp,]
        strain_i_lengths_tmp <- NULL
        
        for(strain in strains){
            tmp <- which(tmpExposures$strain == strain)
            strain_indices <- c(strain_indices, tmp)
            strain_i_lengths_tmp <- c(strain_i_lengths_tmp, length(tmp))
        }
        strain_i_lengths <- c(strain_i_lengths, c(0,cumsum(strain_i_lengths_tmp)))
    }

    ## The length of this vector is the number of groups plus 1
    ## Index starts at 0
    exposure_i_lengths <- c(1,exposure_i_lengths)
    exposure_i_lengths <- cumsum(exposure_i_lengths)
    
#########################################################
    ## Order modifier parameters
#########################################################
    ## Which indices of parTab correspond to order modifiers?
    order_inds <- which(parTab$names == "mod")
    
#########################################################
    ## Cross reactivity parameters
#########################################################
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
    ## For each exposure, we need to get the correct model parameters
    ## If using typing, get the parameter indices that correspond to each type in turn
    ## Convert types to indices as in the convert_types argument
    param_indices <- 1:nrow(parTab)
    
    par_inds <- NULL
    par_lengths <- NULL
    for(i in 1:nrow(exposureTab)){
        tmpType <- exposureTab[i,"type"]
        tmpExp <- exposureTab[i,"exposure"]
        tmpStrain <- exposureTab[i,"strain"]
        tmpID <- exposureTab[i, "id"]
        if(typing){
            if(!cross_reactivity){
                tmp <- param_indices[which(parTab$type %in% c("all",tmpType) & parTab$exposure %in% c(NA,tmpExp) & 
                                           parTab$strain %in% c(NA,tmpStrain,"all") & !(parTab$names %in% c("mod","x")))]
            } else {
                tmp <- param_indices[which(parTab$type %in% c("all",tmpType) & 
                                           parTab$strain %in% c(NA,tmpStrain,"all") & !(parTab$names %in% c("mod","x")))]
                
            }
        } else {
            tmp <- param_indices[which(parTab$id %in% c("all",tmpID) & parTab$strain %in% c("all",tmpStrain,NA) & !(parTab$names %in% c("mod","x")))]
        }
        par_inds <- c(par_inds, tmp)
        par_lengths <- c(par_lengths, length(tmp))
        
    }
    par_lengths <- c(0, cumsum(par_lengths))
    
    
#########################################################
    ## Model parameters
#########################################################
    ## Now we can loop through each group, each strain and
    ## each exposure and get the correct model parameters
    pars <- parTab$values

    ## What are the exposure types, strains, order and priming?
    exposure_types <- exposureTab$type
    exposure_times <- exposureTab$values
    exposure_strains <- exposureTab$exposure
    exposure_measured <- exposureTab$strain
    exposure_orders <- exposureTab$order
    exposure_primes <- exposureTab$primed
    exposure_next <- exposureTab$next_t

    convert_types_back <- names(convert_types)
    convert_strains_back <- names(convert_strains)
    

    ## Convert any named elements to integers
    exposure_types <- convert_types[exposure_types]
    exposure_strains <- convert_strains[exposure_strains]
    exposure_measured <- convert_strains[exposure_measured]
    
    strains <- convert_strains[strains]
    groups <- convert_groups[groups]

    ## Because we are dealing with Cpp where everything is 0 indexed, need to subtract 1 from
    ## the indices calculated in R
    exposure_indices <- exposure_indices - 1
    cr_inds <- cr_inds - 1
    par_inds <- par_inds - 1
    order_inds <- order_inds - 1
    strain_indices <- strain_indices - 1
    exposure_i_lengths <- exposure_i_lengths - 1
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
            ln <- posterior_func_group_cpp(pars, times, groups, strains,
                                           exposure_indices, exposure_i_lengths, strain_indices, strain_i_lengths,
                                           exposure_times, exposure_strains, exposure_next, exposure_measured,
                                           exposure_orders, exposure_primes, cr_inds, par_inds,
                                           order_inds, par_lengths, cr_lengths, ver, individuals, dat)
            if(!is.null(PRIOR_FUNC)) ln <- ln + PRIOR_FUNC(pars)
            ln
        }
    } else {
        f <- function(pars, times){
            model_func_group_cpp(pars, times, groups, strains,
                                 exposure_indices, exposure_i_lengths, strain_indices, strain_i_lengths,
                                 exposure_times, exposure_strains, exposure_next, exposure_measured,
                                 exposure_orders, exposure_primes, cr_inds, par_inds,
                                 order_inds, par_lengths, cr_lengths, ver)
        }
    }
    return(f)
}
