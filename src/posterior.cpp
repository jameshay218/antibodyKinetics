#include "model.h"
using namespace Rcpp;


//' Observation error matrix solver
//' 
//' Calculates a single likelihood for an observed titre given the true titre and observation error parameters
//' @param actual integer for the believed true titre
//' @param obs integer for the observed titre
//' @param S probability of observing the true titre
//' @param EA probability of a + or - 1 observation
//' @param MAX_TITRE integer for the maximum observable titre
//' @return a single probability value
//' @export
//' @useDynLib antibodyKinetics
//[[Rcpp::export]]
double obs_error(int actual, int obs, double S, double EA, int MAX_TITRE){
  if(actual == (MAX_TITRE) && obs == (MAX_TITRE)) return(S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA));
  else if(actual==0 && obs==0) return(S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA));
  else if(actual==obs) return(S);
  else if(actual == (obs + 1) || actual==(obs-1)) return(EA/2.0);
  return((1.0/(MAX_TITRE-2.0))*(1.0-S-EA));
}


//' Discretised normal error
//'
//' Gives the probability of a titre observation given a true titre.
//' @param actual the assumed true titre
//' @param obs the observed titre
//' @param sd standard deviation of the observation error function
//' @param MAX_TITRE the maximum observable titre
//' @export
//[[Rcpp::export]]
double norm_error(double actual, int obs, double sd, int MAX_TITRE){
  double lik = 0;
  
  if(obs >= MAX_TITRE){
    lik = R::pnorm(MAX_TITRE, actual, sd, 0, 0);
  } 
  if(obs < 1.0){
    lik = R::pnorm(1, actual, sd, 1, 0);    
  }
  else {
    lik = R::pnorm(obs+1, actual, sd, 1, 0) - R::pnorm(obs, actual, sd, 1, 0);
  }
  /*
    double lik = 0;
  
    if(obs > MAX_TITRE || obs < 0){
    return 0;
    }
  
    lik = R::pnorm(obs+1, actual, sd, 1, 0) - R::pnorm(obs, actual, sd, 1, 0);
  */
  
  //lik = lik/(R::pnorm(MAX_TITRE,actual,sd,1,0) - R::pnorm(0, actual, sd, 1, 0));
  return lik;  
}

double norm_error_james(double actual, int obs, double sd, int MAX_TITRE){
  double lik = 0;
  
  if(obs > MAX_TITRE || obs < 0){
    return 0;
  }

  if(actual > MAX_TITRE) actual = MAX_TITRE;
  if(actual < 0) actual = 0;
  
  lik = R::pnorm(obs+1, actual, sd, 1, 0) - R::pnorm(obs, actual, sd, 1, 0);
  
  // Normalise by observable range
  lik = lik/(R::pnorm(MAX_TITRE,actual,sd,1,0) - R::pnorm(0, actual, sd, 1, 0));
  return lik;  
}

double norm_error_adam(double actual, int obs, double sd, int MAX_TITRE){
  double lik = 0;
  
  if(obs >= MAX_TITRE){
    lik = R::pnorm(MAX_TITRE, actual, sd, 0, 0);
  } 
  if(obs < 1.0){
    lik = R::pnorm(1, actual, sd, 1, 0);    
  }
  else {
    lik = R::pnorm(obs+1, actual, sd, 1, 0) - R::pnorm(obs, actual, sd, 1, 0);
  }
  
  return lik;  
}


//' Observation error function
//'
//' Given a vector of believed true titres and a vector of observed data, 
//' calculates a likelihood based on a given observation error matrix
//' @param y NumericVector of believed true titres
//' @param data NumericVector of observed data, matching y
//' @param params observation error matrix paramters in order S, EA, MAX_TITRE
//' @export
//[[Rcpp::export]]
double obs_likelihood(NumericVector y, NumericVector data, NumericVector params){
  double ln = 0;
  double tmp = 0;
  int MAX_TITRE = params(3);
  for(int i = 0; i < y.length();++i){
    //if(y(i) < 0) y(i) = 0;
    //if(y(i) >= MAX_TITRE) y(i) = MAX_TITRE;
    //ln += R::dnorm(data(i),y(i),params(1),1);
    if(!NumericVector::is_na(data(i))){
      tmp = norm_error(y(i),data(i),params(1),MAX_TITRE);
    }
    if(tmp == 0){
      ln += -10000;
    } else {
      ln += log(tmp);
    }
    //Rcpp::Rcout << "Data: " << data(i) << "; Model: " << y(i) << "; Lik: " << log(norm_error(y(i),data(i),1.2,MAX_TITRE)) << std::endl;
    //ln += log(obs_error(floor(y(i)), floor(data(i)),params(1),params(2),MAX_TITRE));
    //Rcpp::Rcout << data(i) << " " << y(i) << ": " << R::dpois(data(i),floor(y(i)),1) << std::endl;
    
    //ln += R::dpois(y(i),data(i),1);
  }
  return ln;
}

//' Posterior calculation cpp implementation
//'
//' Solves the antibody kinetics model for given parameters, and then calculates a likelihood
//' for the given data set. This is a complex function, so should only really be called through
//' \code{\link{create_model_group_func_cpp}}. Lookx at this code to really understand what's
//' going on! The key confusing thing is that the length of the vectors has to match the number
//' of rows from the overall parameter table
//' @param pars the vector of model parameters as in parTab
//' @param times the vector of times to solve the model over
//' @param groups IntegerVector of the exposure groups (starting at group 1)
//' @param strains IntegerVector of strains involved in exposures (ie. observed and exposed), starting at 1
//' @param exposure_indices IntegerVector of indices from the exposure table
//' @param exposure_i_lengths IntegerVector of lengths describing the size of blocks in the exposure_indices vector that relate to each exposure group
//' @param strain_indices IntegerVector of indices relating to the exposure table, with contiguous indices for each group and then each strain
//' @param strain_i_lengths IntegerVector of lengths describing the size of blocks in the strain_indices vector that relate to each strain and group
//' @param exposure_times NumericVector infection time of each exposure
//' @param exposure_strains IntegerVector of exposure strains for each exposure
//' @param exposure_next NumericVector specifying the time of the exposure after the current one (for subsetting times)
//' @param exposure_measured IntegerVector of measured strain for each exposure
//' @param exposure_orders IntegerVector of order of exposures
//' @param exposure_primes IntegerVector of whether each exposure was primed or not (for priming modifier)
//' @param cr_inds IntegerVector of indices describing which parTab rows relate to cross reactivity parameters
//' @param par_inds IntegerVector of indices describing which parTab rows relate to model parameters
//' @param par_names CharacterVector names from parTab
//' @param order_inds IntegerVector of indices describing which parTab rows relate to order modifer parameters
//' @param exposure_i_lengths IntegerVector of lengths describing the size of blocks in the exposure_indices vector that relate to each exposure group
//' @param par_lengths IntegerVector of lengths describing the size of blocks in the par_type_ind vector that relate to each exposure type
//' @param cr_lengths IntegerVector of lengths describing the size of blocks in the cr_ind vector that relate to each strain
//' @param version integer indicating which version of the model will be solve. 0 solves the isolated form, and 1 solves the competitive exposure form.
//' @param individuals IntegerVector indicating how many individuals there are in each group
//' @param data NumericMatrix of antibody titre data. Each row should be complete observations of titres against a given strain for a given group. If we have 5 strains measured and 5 groups, rows 1:5 should be titres in group 1, rows 6:10 titres in group 3 etc. If more than 1 individual in each group, multiply these criteria by number of inidividuals in that group (ie., rows 1:15 for group 1)
//' @return a single likelihood of observing the data given the model parameters
//' @export
//[[Rcpp::export]]
double posterior_func_group_cpp(
				NumericVector pars, NumericVector times, 
				IntegerVector groups, IntegerVector strains,
				IntegerVector exposure_indices, IntegerVector exposure_i_lengths,
				IntegerVector strain_indices, IntegerVector strain_i_lengths,
				NumericVector exposure_times, IntegerVector exposure_strains,
				NumericVector exposure_next, IntegerVector exposure_measured,
				IntegerVector exposure_orders, IntegerVector exposure_primes, 
				IntegerVector cr_inds, IntegerVector par_inds,
				CharacterVector par_names,
				IntegerVector order_inds, IntegerVector par_lengths,
				IntegerVector cr_lengths, int version, IntegerVector individuals,
				NumericMatrix data){

  NumericMatrix result = model_func_group_cpp(pars, times, groups, strains, exposure_indices,
					      exposure_i_lengths,strain_indices, strain_i_lengths,
					      exposure_times, exposure_strains, exposure_next,
					      exposure_measured, exposure_orders, exposure_primes,
					      cr_inds, par_inds, par_names, order_inds, par_lengths,
					      cr_lengths, version);
  double ln = 0;
  int index_data = 0;
  int index_model = 0;

  for(int i = 0; i < groups.size(); ++i){
    for(int j = 0; j < strains.size(); ++j){
      for(int k = 0; k < individuals[i]; ++k){
	ln += obs_likelihood(result(index_model,_), data(index_data,_), pars);
	index_data++;
      }
      index_model++;
    }
  }
  return ln;
}






//[[Rcpp::export]]
NumericMatrix posterior_func_group_cpp_matrix(
					      NumericVector pars, NumericVector times, 
					      IntegerVector groups, IntegerVector strains,
					      IntegerVector exposure_indices, IntegerVector exposure_i_lengths,
					      IntegerVector strain_indices, IntegerVector strain_i_lengths,
					      NumericVector exposure_times, IntegerVector exposure_strains,
					      NumericVector exposure_next, IntegerVector exposure_measured,
					      IntegerVector exposure_orders, IntegerVector exposure_primes, 
					      IntegerVector cr_inds, IntegerVector par_inds,
					      CharacterVector par_names,
					      IntegerVector order_inds, IntegerVector par_lengths,
					      IntegerVector cr_lengths, int version, IntegerVector individuals,
					      NumericMatrix data){
  
  NumericMatrix result = model_func_group_cpp(pars, times, groups, strains, exposure_indices,
                                              exposure_i_lengths,strain_indices, strain_i_lengths,
                                              exposure_times, exposure_strains, exposure_next,
                                              exposure_measured, exposure_orders, exposure_primes,
                                              cr_inds, par_inds, par_names, order_inds, par_lengths,
                                              cr_lengths, version);
  int n_col = data.ncol();
  int n_row = data.nrow();
  NumericMatrix likelihoods(n_row,n_col);
  int index_data = 0;
  int index_model = 0;
  int MAX_TITRE = pars(3);
  double sd = pars(1);
  for(int i = 0; i < groups.size(); ++i){
    for(int j = 0; j < strains.size(); ++j){
      for(int k = 0; k < individuals[i]; ++k){
        for(int x = 0; x < times.size(); ++x){
          likelihoods(index_data,x) = log(norm_error(result(index_model,x), data(index_data,x), sd,MAX_TITRE));
        }
        index_data++;
      }
      index_model++;
    }
  }
  return likelihoods;
}
