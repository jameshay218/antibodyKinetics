#include <Rcpp.h>
using namespace Rcpp;

//' Model trajectory calc cpp
//'
//' Calculates the ferret model trajectory for a single infection event.
//' Uses a Cpp implementation for speed. Need to obey order of parameters as in the
//' example
//' @param pars the vector of model parameters
//' @param times the vector of time in days to solve over
//' @return a vector of antibody titres
//' @export
//' @family model functions
//' @useDynLib antibodyKinetics
//' @examples
//' pars <- c("mu"=8,"tp"=12,"dp"=0.5,"ts"=10,"m"=0.003,
//'           "sigma"=0.01, "beta"=0.02,"c"=4,"y0_mod"=0,
//'            "primed"=0,"mod"=1,
//'           "x"=0,"t_i"=10,"y0"=0,"eff_y0"=0)
//' times <- seq(0,100,by=10)
//' y <- model_trajectory_cpp(pars,times)
//' @export
//[[Rcpp::export]]
NumericVector model_trajectory_cpp(NumericVector pars, NumericVector times, bool logSigma){
  double lower_bound = pars[0];
  double mu = pars[4];
  double tp = pars[5];
  double dp = pars[6];
  double ts = pars[7];
  double m = pars[8];
  double c = pars[11];
  double primed = pars[13];
  double mod = pars[14];
  double x = pars[15];
  double t_i = pars[16];

  double sigma, beta, y0_mod;
  if(logSigma){
    sigma = exp(pars[9]);
    beta = exp(pars[10]);
    y0_mod = exp(pars[12]);
  } else {
    sigma = pars[9];
    beta = pars[10];
    y0_mod = pars[12];

}

  // We have y0 twice. In the non-additive version (one antibody producing process),
  // eff_y0 is zero. Otherwise, it's the titre at the time of exposure.
  double y0 = pars[17];
  double eff_y0 = pars[18];

  double t = 0;
  double tmp = 0;

  double cr = exp(-sigma*x);
  double prime_cr = c*exp(-beta*x)*primed;
  double mod_boost = exp(-y0_mod*y0);

  mu = mu*cr*mod*mod_boost + prime_cr;

  NumericVector y(times.size());
  
  for(int i = 0; i < times.size(); ++i){
    t = times[i];
    if(t <= t_i) tmp = 0;
    else if(t > t_i && t <= (t_i + tp)) tmp = (mu/tp)*(t-t_i);
    else if(t > (tp+t_i) && t <=(ts + t_i+tp)) tmp = ((-(dp*mu)/ts)*(t) + ((mu*dp)/ts)*(t_i+tp) + mu);
    else tmp = (-m*(t)+m*(t_i+tp+ts)+(1-dp)*mu);
    tmp += eff_y0;
    if(tmp < lower_bound) tmp = lower_bound;
    y[i] = tmp;
  }
  return y; 
}




//' Model calculation cpp implementation
//'
//' Solves the antibody kinetics model for given parameters. This is a complex function, 
//' so should only really be called through \code{\link{create_model_group_func_cpp}}. 
//' Look at this code to really understand what's going on! The key confusing thing is 
//' that the length of the vectors has to match the number of rows from the overall 
//' parameter table
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
//' @param order_inds IntegerVector of indices describing which parTab rows relate to order modifer parameters
//' @param exposure_i_lengths IntegerVector of lengths describing the size of blocks in the exposure_indices vector that relate to each exposure group
//' @param par_lengths IntegerVector of lengths describing the size of blocks in the par_type_ind vector that relate to each exposure type
//' @param cr_lengths IntegerVector of lengths describing the size of blocks in the cr_ind vector that relate to each strain
//' @param version integer indicating which version of the model will be solve. 0 solves the isolated form, and 1 solves the competitive exposure form.
//' @return a matrix of antibody kinetic trajectories, with rows for group and then strain
//' @export
//' @useDynLib antibodyKinetics
//[[Rcpp::export]]
NumericMatrix model_func_group_cpp(NumericVector pars, NumericVector times, 
				   IntegerVector groups, IntegerVector strains,
				   IntegerVector exposure_indices, IntegerVector exposure_i_lengths,
				   IntegerVector strain_indices, IntegerVector strain_i_lengths,
				   NumericVector exposure_times, IntegerVector exposure_strains,
				   NumericVector exposure_next, IntegerVector exposure_measured,
				   IntegerVector exposure_orders, IntegerVector exposure_primes, 
				   IntegerVector cr_inds, IntegerVector par_inds, 
				   IntegerVector order_inds, IntegerVector par_lengths,
				   IntegerVector cr_lengths, int version){
  int group; // Index of group
  int exposure_strain, measured_strain; // Index of strain
  int A, B; // Used to get vector subsets within a range
  int type, order; // Used for indexing exposure properties
  double t_i, mod, cr,isPrimed, next_t, y0; // Temporary model parameters

  IntegerVector tmp_exposures_group, tmp_exposures_strain; // For storing subset of exposures
  IntegerVector tmp_strains, tmp_strain_lengths;
  IntegerVector tmp_exposures;
  IntegerVector time_indices = seq_len(times.size()); // Sequence along the time vector for subsetting other stuff
  IntegerVector tmp_time_indices = time_indices;  // For storing subset of time vector indices

  NumericVector fullPars; // For storing subset of model parameters
  NumericMatrix results(strains.size()*groups.size(), times.size()); // Storing results
  NumericVector y; // Storing a single trajectory
  NumericVector tmpTimes = times; // Storing a subset of times to solve over

  tmpTimes = times;
  tmp_time_indices = time_indices;

  int index = 0;
  int index_dat = 0;

  
  // For each group
  for(int i = 0; i < groups.size(); ++i){
    group = groups[i];

    // Get exposures for this group
    A = exposure_i_lengths[i];
    B = exposure_i_lengths[i+1] - 1;
    tmp_exposures_group = exposure_indices[Range(A,B)];

    // For each strain, a subset of these exposures apply
    tmp_strains = strain_indices[Range(A,B)];
    tmp_strain_lengths = strain_i_lengths[Range((i*5),(i*5)+5)];
    
    // For each strain in this group
    for(int j = 0; j < tmp_strain_lengths.size(); ++j){
      A = tmp_strain_lengths[j];
      B = tmp_strain_lengths[j+1] - 1;
      tmp_exposures = tmp_exposures_group[tmp_strains[Range(A,B)]];
      
      y0 = 0;
      
      // For each exposure for this strain
      for(int k = 0; k < tmp_exposures.size(); ++k){
	tmpTimes = times;
	/* The par_lengths vector should be the same size as
	   the number of exposures (+1 for the first index of 0)
	   Use this to get subset of parameters.
	*/
	A = par_lengths[tmp_exposures[k]];
	B = par_lengths[tmp_exposures[k+1]] - 1;
	fullPars = pars[par_inds[Range(A,B)]];
	
	index = tmp_exposures[k];
	t_i = exposure_times[index];
	exposure_strain = exposure_strains[index]-1;
	measured_strain = exposure_measured[index]-1;
	order = exposure_orders[index];
	mod = pars[order_inds[order]];
	isPrimed = exposure_primes[index];

	cr = pars[cr_inds[cr_lengths[measured_strain] + exposure_strain]];
	
	fullPars.push_back(isPrimed);
	fullPars.push_back(mod);
	fullPars.push_back(cr);
	fullPars.push_back(t_i);
	fullPars.push_back(y0);
	
	// Depending on model version, use effective y0 of zero or 
	// actual y0 at start of exposure
	if(version == 0){
	  fullPars.push_back(0);
	} else {
	  fullPars.push_back(y0);
	}

	/* Which times we solve over depends on the version. If the isolated boosting version (0), 
	   solve over all times. If it's the competitive boosting version (1), we need to subset the
	   times vector to those times between the current and next infection. */
	if(version == 1){
	  if(next_t == times[times.size()-1]){
	    tmpTimes = times[times >= t_i & times <= next_t];
	    tmp_time_indices = time_indices[times >= t_i & times <= next_t];
	  } else {
	    tmpTimes = times[times >= t_i & times < next_t];
	    tmp_time_indices = time_indices[times >= t_i & times < next_t];
	  }
	}
	tmpTimes.push_back(next_t);

	
	// Solve the model
	y = model_trajectory_cpp(fullPars, tmpTimes, TRUE);

	/* Add model solutions to the correct row in the results matrix. If version 0, this will
	   be all time points. If version 1, this will correspond to the times relating to the
	   current exposure */
	for(int ii = 0; ii < tmp_time_indices.size(); ++ii){
	  results(index_dat, tmp_time_indices[ii]-1) += y[ii];
	}
	// Get the last element from the y vector, as this will be y0 for the next exposure
	y0 = y[y.size()-1];
      }
      index_dat++;
    }
  }
  return results;
}

 
