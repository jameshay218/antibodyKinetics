#' Read in ferret data
#'
#' Reads in the csv of ferret antibody titre data from a specified file
#' @param data_file the full file path of the data file
#' @return the melted antibody titre data to be used in the model fitting
#' @export
#' @useDynLib antibodyKinetics
read_ferret_data <- function(data_file){
    raw_data <- read.csv(data_file, stringsAsFactors=FALSE,na.strings=c("NA","-","?"),header=0)
    colnames(raw_data) <- raw_data[1,]
    raw_data <- raw_data[2:nrow(raw_data),]

    melted.data <- reshape2::melt(raw_data)
    melted.data$variable <- as.integer(as.character(melted.data$variable))
    return(melted.data)

}

#' Add noise to a titre observation
#'
#' Adds observation noise to the provided vector of titres
#' @param pars vector of parameters with S, EA and MAX_TITRE
#' @param y the titre to add observation error
#' @param normal indicates if data is from the discretised normal or not (uses a shouldered observation error function if not)
#' @return a single, modified titre value
#' @export
#' @useDynLib
add_noise <- function(pars, y, normal=FALSE){
    MAX_TITRE <- pars["MAX_TITRE"]
    S <- pars["S"]
    EA <- pars["EA"]

    if(y < 0) y <- 0
    if(y > MAX_TITRE) y <- MAX_TITRE

    if(!normal){
        probs <- numeric(MAX_TITRE+1)
        
        probs[] <- (1.0/(MAX_TITRE-2.0))*(1.0-S-EA)
        
        if(y == MAX_TITRE){
            probs[y+1] <- S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA)
            probs[y] <- EA/2.0
        } else if (y == 0){
            probs[y+1] <- S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA)
            probs[y+2] <- EA/2.0        
        } else {
            probs[y+1] <- S
            probs[y] <- EA/2.0
            probs[y+2] <- EA/2.0
        }
    } else {
        probs <- numeric(MAX_TITRE + 1)
        for(i in 1:length(probs)){
            probs[i] <- norm_error(y, i-1, S, MAX_TITRE)
        }
    }
    probs <- cumsum(probs)

    tmp <- runif(1,0,1)

    i <- 1
    while(probs[i] < tmp){
        i <- i + 1
    }
    return(i-1)
}


#' Modify parameter table
#'
#' Modifies the parTab data frame depending on the model settings
#' @param parTab the data frame of parTab as read in by a typical parTab file
#' @param options list with the following elements:
#' \enumerate{
#'     \item typing: boolean to indicate if "typing" is used in the model
#'     \item cr: boolean to indicate if cross reactivity is included
#'     \item priming: boolean for priming
#'     \item monophasic_waning: boolean for monophasic waning. FALSE if biphasic waning
#'     \item y0_mod: boolean if titre dependent boosting is used
#'     \item antigenic_seniority: boolean if antigenic seniority is included
#'     \item form: string argument for which model form is used
#' }
#' Set the contents of these elements as described above (mostly booleans) to set the parameter table to include or exclude each model feature
#' @param fixed_S boolean, if TRUE then the standard deviation of the observation error function is estimated
#' @return the correct parTab data frame for the given model structure
#' @export
parTab_modification <- function(parTab, options,fixed_S=FALSE){
    form=options$form
    antigenic_seniority=options$antigenic_seniority
    cr=options$cr
    typed_cr=options$typed_cr
    priming=options$priming
    types=options$types
    monophasic_waning=options$monophasic_waning
    y0_mod=options$y0_mod
    
    ## Modify if monophasic waning
    if(monophasic_waning){
        parTab[parTab$names %in% c("ts","dp"),"fixed"] <- 1
        parTab[parTab$names %in% c("ts","dp"),"values"] <- 0
    } else {
        parTab[parTab$names %in% c("ts","dp"),"fixed"] <- 0
    }

    ## Modify if titre dependent boosting
    if(y0_mod){
        parTab[parTab$names =="y0_mod","fixed"] <- 0
        parTab[parTab$names == "y0_mod","upper_bound"] <- 1
        parTab[parTab$names == "y0_mod","lower_bound"] <- -1
        
        parTab[parTab$names =="boost_limit","fixed"] <- 0
        parTab[parTab$names =="boost_limit","upper_bound"] <- 12
        parTab[parTab$names =="boost_limit","lower_bound"] <- 0
    } else {
        parTab[parTab$names =="y0_mod","fixed"] <- 1
        parTab[parTab$names =="boost_limit","fixed"] <- 1
        parTab[parTab$names =="y0_mod","values"] <- -1000
        parTab[parTab$names =="boost_limit","values"] <- -20

        parTab[parTab$names == "y0_mod","upper_bound"] <- 10000
        parTab[parTab$names == "y0_mod","lower_bound"] <- -10000
        parTab[parTab$names =="boost_limit","upper_bound"] <- 1000
        parTab[parTab$names =="boost_limit","lower_bound"] <- -1000
    }

    ## Modify if antigenic seniority
    if(antigenic_seniority){
        parTab[parTab$names == "mod","fixed"] <- 0
        parTab[parTab$names == "mod","fixed"][1] <- 1
    } else {
        parTab[parTab$names == "mod","fixed"] <- 1
        parTab[parTab$names == "mod","values"] <- 1   
    }

    ## Modify if estimating error variance
    if(!fixed_S){
        parTab[parTab$names == "S","fixed"] <- 0
    } else {
        parTab[parTab$names == "S","fixed"] <- 1
    }

    ## Modify if priming
    if(priming){
        parTab[parTab$names %in% c("beta","c"),"fixed"] <- 0
    } else {
        parTab[parTab$names %in% c("beta","c"),"fixed"] <- 1
    }
    parTab[parTab$names == "m","upper_bound"] <- 12
    return(parTab)
}

#' Convert runName string
#'
#' Given a string describing the run name (eg. "CYAY3BN"), generates a list of model options to be fed into \code{\link{parTab_modification}}
#' @param runName a string with 7 characters, each corresponding to: form (C or I); antigenic seniority (Y or N); cross reactivity (A or T); priming (Y or N); number of exposure types (3 or 6); waning (B or M); titre dependent boosting (Y or N)
#' @return a list to be fed into the options argument of \code{\link{parTab_modification}}
#' @export
convert_runName_to_options <- function(runName){
    options <- substring(runName, seq(1,nchar(runName),1), seq(1,nchar(runName),1))
    names(options) <- c("form","as","cr","priming","types","wane","y0")

    form <- "isolated"
    if(options["form"] == "C") form <- "competitive"

    antigenic_seniority <- FALSE
    if(options["as"] == "Y") antigenic_seniority <- TRUE
    
    cr <- FALSE
    typed_cr <- FALSE
    if(options["cr"] == "A"){
        cr <- TRUE
    } else if(options["cr"] == "T"){
        cr <- TRUE
        typed_cr <- TRUE
    }

    priming <- FALSE
    if(options["priming"] == "Y") priming <- TRUE

    types <- 6
    if(options["types"] == 3) types <- 3
    if(options["types"] == 0) types <- 0

    monophasic_waning <- FALSE
    if(options["wane"] == "M") monophasic_waning <- TRUE

    y0_mod <- FALSE
    if(options["y0"] == "Y") y0_mod <- TRUE

    return(list(form=form,antigenic_seniority=antigenic_seniority,cr=cr,typed_cr=typed_cr,
                priming=priming,types=types,monophasic_waning=monophasic_waning,
                y0_mod=y0_mod))
}

#' Convert options to runName
#'
#' The inverse function of \code{\link{convert_runName_to_options}}
#' @param options list with the following elements:
#' 1. typing: boolean to indicate if "typing" is used in the model
#' 2. cr: boolean to indicate if cross reactivity is included
#' 3. priming: boolean for priming
#' 4. monophasic_waning: boolean for monophasic waning. FALSE if biphasic waning
#' 5. y0_mod: boolean if titre dependent boosting is used
#' 6. antigenic_seniority: boolean if antigenic seniority is included
#' 7. form: string argument for which model form is used
#' Set the contents of these elements as described above (mostly booleans) to set the parameter table to include or exclude each model feature
#' @return a 7 character long string corresponding to a run name identifier
#' @export
convert_options_to_runName <- function(options){

    form <- "C"
    if(options$form == "isolated") form <- "I"

    as <- "Y"
    if(options$antigenic_seniority == FALSE) as <- "N"
    
    cr <- "N"
    if(options$cr == TRUE){
        cr <- "T"
        if(options$typed_cr == FALSE) cr <- "A"
    }

    priming <- "Y"
    if(options$priming == FALSE) priming <- "N"

    types <- as.character(options$types)

    wane <- "B"
    if(options$monophasic_waning == TRUE) wane <- "M"

    y0_mod <- "Y"
    if(options$y0_mod == FALSE) y0_mod <- "N"

    runName <- paste(c(form,as,cr,priming,types,wane,y0_mod),sep="",collapse="")
    return(runName)
}

#' Describe model parameters
#'
#' Prints out descriptions of the model parameters corresponding to \code{\link{model_trajectory}} and \code{\link{model_trajectory_cpp}}
#' @return nothing
#' @export
parameter_descriptions <- function(){
    message("============================")
    message("== Parameter descriptions ==")
    message("============================")
    message(" -- Note: parameters are described in the order they should be passed to model_trajectory_cpp")
    message(cat("Name","Description", sep="\t"))
    message(cat("lower_bound", "assumed true lower titre bound on a log scale (a value of 0 is represents no haemagglutination at 1:1 dilution, and therefore not necessarily no antibodies", sep="\t"))
    message(cat("S", "standard deviation of the discritised normally distributed error function OR see `obs_error` function for alternative measurement error function", sep="\t"))
    message(cat("EA", "deprecated  - see `obs_error` function", sep="\t"))
    message(cat("MAX_TITRE", "maximum observable antibody titre on log scale", sep="\t"))
    message(cat("mu", "maximum homologous boost after exposure", sep="\t"))
    message(cat("tp", "time from exposure to peak boost", sep="\t"))
    message(cat("dp", "proportion of initial boost lost following initial waning phase", sep="\t"))
    message(cat("ts", "duration in days of initial waning phase", sep="\t"))
    message(cat("m", "gradient of long term waning phase in log units lost per day", sep="\t"))
    message(cat("beta", "cross reactivity gradient of additional boost from priming", sep="\t"))
    message(cat("c", "additional homologous boosting from primed exposure", sep="\t"))
    message(cat("sigma", "cross reactivity gradient", sep="\t"))
    message(cat("y0_mod", "gradient of titre-dependent boosting term. -10000 turns this feature off; -1 is maximum titre dependent suppression; 1 is maximum titre dependent enhancement", sep="\t"))
    message(cat("boost_limit", "maximum titre below which titre-dependent boosting takes place. If y0_mod is -10000, then this is ignored", sep="\t"))
    message(cat("primed", "1 if this is a primed exposure, 0 otherwise", sep="\t"))
    message(cat("mod", "corresponding to rho in model description - scaling parameter for antigenic seniority between 0 and 1", sep="\t"))
    message(cat("x", "antigenic distance between the exposure strain and the strain being measured here (0 for homologous", sep="\t"))
    message(cat("t_i", "time of exposure in days", sep="\t"))
    message(cat("y0", "initial titre at time of exposure.", sep="\t"))
    message(cat("eff_y0", "initial titre for the purpose of calculating titre dependent boosting. If competitive version, this is the actual titre at the time of exposure. If isolated, this should be 0", sep="\t"))
    NULL
}
