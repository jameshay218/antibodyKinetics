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
NumericVector model_trajectory_cpp(NumericVector pars, NumericVector times){
  double lower_bound = pars[0];
  double mu = pars[4];
  double tp = pars[5];
  double dp = pars[6];
  double ts = pars[7];
  double m = pars[8];
  double sigma = exp(pars[9]);
  double beta = exp(pars[10]);
  double c = pars[11];
  double y0_mod = exp(pars[12]);
  double primed = pars[13];
  double mod = pars[14];
  double x = pars[15];
  double t_i = pars[16];

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
    if(tmp < lower_bound) tmp = lower_bound;
    tmp += eff_y0;
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
//' @param exposure_types IntegerVector of exposure types matching the exposure table. "all"=0, "infection"=1,"vacc"=2,"adj"=3,"mod"=4,"NA"=5
//' @param exposure_strains IntegerVector of exposure strains for each exposure. Note that this goes from 1 to 5.
//' @param measured_strains IntegerVector of measured strains for each exposure (ie. second index in cross reactivity calculation)
//' @param exposure_orders IntegerVector of order of exposures
//' @param exposure_primes IntegerVector of whether each exposure was primed or not (for priming modifier)
//' @param exposure_indices IntegerVector of indices describing which parTab rows relate to exposure parameters
//' @param cr_inds IntegerVector of indices describing which parTab rows relate to cross reactivity parameters
//' @param par_type_ind IntegerVector of indices describing which parTab rows relate to model parameters
//' @param order_indices IntegerVector of indices describing which parTab rows relate to order modifer parameters
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
				   IntegerVector exposure_types, IntegerVector exposure_strains, 
				   IntegerVector measured_strains, IntegerVector exposure_orders, 
				   IntegerVector exposure_primes, IntegerVector exposure_indices, 
				   IntegerVector cr_inds, IntegerVector par_type_ind, 
				   IntegerVector order_indices, IntegerVector exposure_i_lengths, 
				   IntegerVector par_lengths, IntegerVector cr_lengths,
				   int version){
  int group; // Index of group
  int strain; // Index of strain
  int A, B; // Used to get vector subsets within a range
  int type, order, exposure_strain, old_strain, ii; // Used for indexing exposure properties
  double old_ti, t_i, mod, cr,isPrimed, old_cr, next_t, y0; // Temporary model parameters

  IntegerVector tmp_exposures; // For storing subset of exposures
  IntegerVector time_indices = seq_len(times.size()); // Sequence along the time vector for subsetting other stuff
  IntegerVector tmp_time_indices = time_indices;  // For storing subset of time vector indices

  NumericVector fullPars; // For storing subset of model parameters
  NumericMatrix results(strains.size()*groups.size(), times.size()); // Storing results
  NumericVector y; // Storing a single trajectory
  NumericVector tmpTimes = times; // Storing a subset of times to solve over

  tmpTimes = times;
  tmp_time_indices = time_indices;

  int index = 0;

  // For each group
  for(int i = 0; i < groups.size(); ++i){
    group = groups[i];

    // Get exposures for this group
    A = exposure_i_lengths[i];
    B = exposure_i_lengths[i+1] - 1;
    tmp_exposures = exposure_indices[Range(A,B)];

    // For each strain
    for(int j = 0; j < strains.size(); ++j){
      strain = strains[j] -1;
      y0 = 0.0;
      
      // For each exposure
      for(int k = 0; k < tmp_exposures.size(); ++k){
	t_i = pars[tmp_exposures[k]];
	order = exposure_orders[tmp_exposures[k]] - 1;
	exposure_strain = exposure_strains[tmp_exposures[k]] - 1;
	mod = pars[order_indices[order]];
	isPrimed = exposure_primes[tmp_exposures[k]];

	cr = pars[cr_inds[cr_lengths[strain] + exposure_strain]];

	// **********************************************************************************
	/* Trying to work out best way to deal with multiple exposures at the same time.
	   I think the best bet is to see if the new t_i is the same as the old one,
	   and to replace the trajectory if the new exposure has lower antigenic distance. 
	   This might make inference more difficult as we get switching of dominant processes,
	   but should work out...
	*/
	// If this isn't the first infection, check the previous one
	if(k > 0){
	  old_ti = pars[tmp_exposures[k-1]];
	  // If this new exposure was at the same time as the previous one
	  if(old_ti == t_i){
	    // Get the previous exposure strain
	    old_strain = exposure_strains[tmp_exposures[k-1]] - 1;
	    // Get cross reactivity to old strain
	    old_cr = pars[cr_inds[cr_lengths[strain] + old_strain]];
	    // If old strain was more closely related, then skip this exposure
	    if(old_cr < cr) continue;
	    else {
	      // Otherwise, remove the previously recorded trajectory
	      //Rcpp::Rcout << "Removing old... " << std::endl;
	      for(int ii = 0; ii < tmp_time_indices.size(); ++ii){
		results(index, tmp_time_indices[ii]) -= y[ii];
	      }
	    }
	  }
	}
	
	// **********************************************************************************
	
	// **********************************************************************************
	/* If this isn't the last infection, we need to find the time of the next infection
	   and solve the model for this time, as this will be used as y0 for the next exposure.
	   This extra time is appended to the vector of times to solve over, and we pull this out 
	   and use it to get y0 for the next exposure. 
	*/
	ii = k; // Save current exposure
	// If this is the last exposure, the "next infection time" is simply the end of the times vector
	if((k+1) == tmp_exposures.size()){
	  next_t = times[times.size()-1];
	  // Otherwise
	} else {
	  // Next exposure time might be next exposure in vector
	  next_t = pars[tmp_exposures[ii + 1]];
	  // However, next exposure time might be the same as this one. If this is the case.
	  // need to keep looking until we find a later exposure or the end of the times vector
	  while(next_t == t_i && (ii+1) < tmp_exposures.size()){
	    next_t = pars[tmp_exposures[ii + 1]];
	    ii++;
	  }
	  // If we got to the last exposure and it's still the same, use the last time
	  if((ii+1) == tmp_exposures.size()){
	    next_t = times[times.size()-1];
	  }
	}
	
	/* Which times we solve over depends on the version. If the isolated boosting version (0), 
	   solve over all times. If it's the competitive boosting version (1), we need to subset the
	   times vector to those times between the current and next infection. */
	if(version == 1){
	  if(next_t == times[times.size()-1]){
	    tmpTimes = times[times >= t_i & times <= next_t];
	    tmp_time_indices = time_indices[times >= t_i & times <= next_t];
	    tmpTimes.push_back(next_t);
	  } else {
	    tmpTimes = times[times >= t_i & times < next_t];
	    tmp_time_indices = time_indices[times >= t_i & times < next_t];
	  }
	}
	tmpTimes.push_back(next_t);

	// **********************************************************************************
	
	type = exposure_types[tmp_exposures[k]] - 1;
	A = par_lengths[type];
	B = par_lengths[type+1] - 1;

	fullPars = pars[par_type_ind[Range(A,B)]];
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


	// Solve the model
	y = model_trajectory_cpp(fullPars, tmpTimes);

	/* Add model solutions to the correct row in the results matrix. If version 0, this will
	   be all time points. If version 1, this will correspond to the times relating to the
	   current exposure */
	for(int ii = 0; ii < tmp_time_indices.size(); ++ii){
	  results(index, tmp_time_indices[ii]-1) += y[ii];
	}
	// Get the last element from the y vector, as this will be y0 for the next exposure
	y0 = y[y.size()-1];
      }
      index++;
    }
  }
  return results;
}
