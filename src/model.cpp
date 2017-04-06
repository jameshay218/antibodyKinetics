#include <Rcpp.h>
using namespace Rcpp;


//' @export
//' @useDynLib antibodyKinetics
//[[Rcpp::export]]
NumericMatrix predict_titres(NumericVector params, NumericVector times){
  double y0 = 0;
  double final_t;
  double t_i;
  double mu, dRF, tp, ts, m;
  double tmp;
  double tmp2;
  NumericVector::iterator t = times.begin();
  NumericVector Y(times.size());
  


  int j = 0;
  double trunc_lower = params[0];
  double first_infection = params[1];
  double max_t = times[times.size()-1];
  int no_infections=params.size()/6;
  NumericMatrix out(times.size(),2);

  for(int i=1;i <= no_infections; ++i){
    if(i==no_infections){
      final_t = max_t;
    }
    else {
      final_t = params[6*(i)+1];
    }

    t_i = params[6*(i-1)+1];
    mu = params[6*(i-1)+2];
    dRF = params[6*(i-1)+3];
    tp = params[6*(i-1)+4];
    ts = params[6*(i-1)+5];
    m = params[6*(i-1)+6];

    while(t != times.end() && *t <= final_t){
      tmp = 0;
      if(*t <= t_i) tmp = 0;
      else if(*t > t_i && *t <= (t_i+tp)) tmp = (mu/tp)*(*t-t_i);
      else if(*t > (tp+t_i) && *t <=(ts + t_i+tp)) tmp = ((-(dRF*mu)/ts)*(*t) + ((mu*dRF)/ts)*(t_i+tp) + mu);
      else tmp = (-m*(*t)+m*(t_i+tp+ts)+(1-dRF)*mu);
      tmp += y0;
      
      Y[j] = tmp;
      ++t;
      ++j;
    }
    
    if((6*i + 1) < params.size()){
      tmp2 = params[6*i + 1];
      if(tmp2 <= t_i) tmp = 0;
      else if(tmp2 > t_i && tmp2 <= (t_i+tp)) tmp = (mu/tp)*(tmp2-t_i);
      else if(tmp2 > (tp+t_i) && tmp2 <=(ts + t_i+tp)) tmp = ((-(dRF*mu)/ts)*(tmp2) + ((mu*dRF)/ts)*(t_i+tp) + mu);
      else tmp = (-m*(tmp2)+m*(t_i+tp+ts)+(1-dRF)*mu);
      y0 = y0 + tmp;
    }
    
    if(y0 < trunc_lower){
      y0 = trunc_lower;
    }
  }
  
  out(_,0) = times;
  out(_,1)=Y;
  return out;
}

//' Single trajectory function
//' @export
//[[Rcpp::export]]
NumericVector simple_model(NumericVector pars, NumericVector times){
  double mu = pars(0)*pars(1)*pars(2)*(pars(3)*pars(4) + 1*!pars(4));
  double dp = pars(5);
  double tp = pars(6);
  double ts = pars(7);
  double m = pars(8);
  double y0 = pars(9);
  double t_i = pars(10);
  double t;
  double tmp;
  
  NumericVector y(times.size());

  for(int i = 0; i < times.size(); ++i){
    t = times(i);
    tmp = 0;
    if(t <= t_i) tmp = y0;
    else if(t > t_i && t <= (t_i + tp)) tmp = (mu/tp)*(t-t_i);
    else if(t > (tp+t_i) && t <=(ts + t_i+tp)) tmp = ((-(dp*mu)/ts)*(t) + ((mu*dp)/ts)*(t_i+tp) + mu);
    else tmp = (-m*(t)+m*(t_i+tp+ts)+(1-dp)*mu);
    y(i) = tmp;
  }
  return(y);
}



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
//'           "sigma"=0.01, "beta"=0.02,"c"=4, "primed"=0,"mod"=1,
//'           "x"=0,"t_i"=10)
//' times <- seq(0,100,by=10)
//' y <- model_trajectory_cpp(pars,times)
//' @export
//[[Rcpp::export]]
NumericVector model_trajectory_cpp(NumericVector pars, NumericVector times){
  double mu = pars[3];
  double tp = pars[4];
  double dp = pars[5];
  double ts = pars[6];
  double m = pars[7];
  double sigma = exp(pars[8]);
  double beta = exp(pars[9]);
  double c = pars[10];
  double primed = pars[11];
  double mod = pars[12];
  double x = pars[13];
  double t_i = pars[14];

  double y0 = 0;
  double t = 0;
  double tmp = 0;

  double cr = exp(-sigma*x);
  double prime_cr = c*exp(-beta*x)*primed;

  mu = mu*cr*mod + prime_cr;

  NumericVector y(times.size());
  
  for(int i = 0; i < times.size(); ++i){
    t = times[i];
    if(t <= t_i) tmp = 0;
    else if(t > t_i && t <= (t_i + tp)) tmp = (mu/tp)*(t-t_i);
    else if(t > (tp+t_i) && t <=(ts + t_i+tp)) tmp = ((-(dp*mu)/ts)*(t) + ((mu*dp)/ts)*(t_i+tp) + mu);
    else tmp = (-m*(t)+m*(t_i+tp+ts)+(1-dp)*mu);
    tmp += y0;
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
				   IntegerVector par_lengths, IntegerVector cr_lengths){
  int group;
  int strain;
  int A, B, type, order, exposure_strain;
  double t_i, mod, cr,isPrimed;
  
  IntegerVector tmp_exposures;
  NumericVector fullPars;
  NumericMatrix results(strains.size()*groups.size(), times.size());
  NumericVector y;
  int index = 0;
  for(int i = 0; i < groups.size(); ++i){
    group = groups[i];
    A = exposure_i_lengths[i];
    B = exposure_i_lengths[i+1] - 1;
    tmp_exposures = exposure_indices[Range(A,B)];
    for(int j = 0; j < strains.size(); ++j){
      strain = strains[j] -1;
      for(int k = 0; k < tmp_exposures.size(); ++k){
	t_i = pars[tmp_exposures[k]];
	order = exposure_orders[tmp_exposures[k]] - 1;
	exposure_strain = exposure_strains[tmp_exposures[k]] - 1;
	mod = pars[order_indices[order]];
	isPrimed = exposure_primes[tmp_exposures[k]];

	cr = pars[cr_inds[cr_lengths[strain] + exposure_strain]];
	type = exposure_types[tmp_exposures[k]] - 1;
	A = par_lengths[type];
	B = par_lengths[type+1] - 1;

	fullPars = pars[par_type_ind[Range(A,B)]];
	fullPars.push_back(isPrimed);
	fullPars.push_back(mod);
	fullPars.push_back(cr);
	fullPars.push_back(t_i);

	y = model_trajectory_cpp(fullPars, times);

	for(int ii = 0; ii < y.size(); ++ii){
	  results(index, ii) += y[ii];
	}
      }
      index++;
    }
  }
  return results;
}
