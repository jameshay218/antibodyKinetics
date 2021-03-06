#' Example Exposure Table
#' Exposure input table used to generate trajectories for a model with the following properties:
#' \itemize{
#'     \item 5 experimental groups
#'     \item 5 influenza strains being tracked (A-E)
#'     \item Competitive interactions of multiple exposures
#'     \item Antigenic seniority present
#'     \item Type-specific cross reactivity parameters
#'     \item Impact of priming infection subsequent vaccine response
#'     \item 6 distinct exposure types
#'     \item Bi-phasic antibody waning
#'     \item No titre-dependent boosting
#' }
#' 
#' Each row corresponds to an observable antibody titre trajectory in the model. Columns are as follows:
#' 
#' \describe{
#'     \item{id}{Unique identifier describing which group (G), exposure number (E) and exposure strain (S) this row corresponds to}
#'     \item{values}{Time in days of this exposure}
#'     \item{type}{If using type-specific parameters, which exposure type does this parameter correspond to? either infection, vacc or adj if 3 types, or infection1, infection2, vacc1, vacc2, adj1, adj2 if using 6 types.}
#'     \item{exposure}{f using exposure-strain specific interactions, which strain is in this exposure?}
#'     \item{strain}{If using exposure-strain specific interactions, which observed strain are we currently describing?}
#'     \item{order}{Which number exposure is this, starting from 1 (ie. the nth exposure)}
#'     \item{primed}{1 if this is a vaccination following priming infection, 0 otherwise}
#'     \item{group}{Which experimental group does this exposure correspond to?}
#'     \item{end}{At what time in days does the experimental protocol end?}
#'     \item{next_t}{Time at which the next exposure in this experimental group occurs}
#' }
#' @docType data
#' 
#' @usage data(exampleExposureTab)
#' 
#' @format An object of class \code{"data.frame"}
"exampleExposureTab"