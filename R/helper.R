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
