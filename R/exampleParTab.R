#' Example Parameter Table
#' 
#' Parameter input table used to generate trajectories for a model with the following properties:
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
#' Each row corresponds to a parameter in the model. Columns of the parameter table are as follows:
#' 
#' \describe{
#'     \item{names}{name of the parameter, see parameter_descriptions()}
#'     \item{id}{Which exposure ID does this parameter correspond to? This will be "all" for most model variants}
#'     \item{values}{Value of this parameter}
#'     \item{type}{If using type-specific parameters, which exposure type does this parameter correspond to? either infection, vacc or adj if 3 types, or infection1, infection2, vacc1, vacc2, adj1, adj2 if using 6 types.}
#'     \item{exposure}{f using exposure-strain specific interactions, which strain is in this exposure?}
#'     \item{strain}{If using exposure-strain specific interactions, which observed strain are we currently describing?}
#'     \item{order}{only corresponds to the "mod" parameters, with 1 being first exposure and 4 being the 4th exposure}
#'     \item{fixed}{1 if this is a fixed parameter, 0 if to be estimated in the MCMC procedure}
#'     \item{steps}{initial step size for the univariate MCMC sampler}
#'     \item{lower_bound}{lower allowable bound for this parameter}
#'     \item{upper_bound}{upper allowable bound for this parameter}
#' }
#' @docType data
#' 
#' @usage data(exampleParTab)
#' 
#' @format An object of class \code{"data.frame"}
"exampleParTab"