#' Read in ferret data
#'
#' Reads in the csv of ferret antibody titre data from a specified file
#' @param data_file the full file path of the data file
#' @return the melted antibody titre data
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
#' @param pars vector of parameters with S, EA and MAX_TITRE
#' @param y the titre to add observation error
#' @param normal indicates if data is from the discretised normal or not
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
#' @param typing boolean to indicate if "typing" is used in the model
#' @param cr boolean to indicate if cross reactivity is included
#' @param priming boolean for priming
#' @param monophasic_waning boolean for monophasic waning. FALSE if biphasic waning
#' @param y0_mod boolean if titre dependent boosting is used
#' @param antigenic_seniority boolean if antigenic seniority is included
#' @param form string argument for which model form is used
#' @return the correct parTab data frame
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
        parTab[parTab$names %in% c("ts","dp"),"fixed"] <- 1
    }

    ## Modify if titre dependent boosting
    if(y0_mod){
        parTab[parTab$names =="y0_mod","fixed"] <- 0
    } else {
        parTab[parTab$names =="y0_mod","fixed"] <- 1
        parTab[parTab$names =="y0_mod","values"] <- -20
    }

    ## Modify if antigenic seniority
    if(antigenic_seniority){
        parTab[parTab$names == "mod","fixed"] <- 0      
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
    return(parTab)
}

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
