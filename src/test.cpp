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


//' Converts to unit scale
//'
//' @param x the double to be converted
//' @param min the minimum value on the linear scale
//' @param max the maximum value on the linear scale
//' @return the value converted to a unit scale
//' @export
//[[Rcpp::export]]
double toUnitScale(double x, double min, double max){
  return((x-min)/(max-min));
}

//' Converts to linear scale
//'
//' @param x the double to be converted back to linear scale
//' @param min the minimum value on the linear scale
//' @param max the maximum value on the linear scale
//' @return the value converted to a linear scale
//' @export
//[[Rcpp::export]]
double fromUnitScale(double x, double min, double max){
  return(min + (max-min)*x);
}


//' @export
//[[Rcpp::export]]
double obs_error(int actual, int obs, double S, double EA){
  int MAX_TITRE = 13;
  if(actual == (MAX_TITRE) && obs == (MAX_TITRE)) return(S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA));
  else if(actual==0 && obs==0) return(S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA));
  else if(actual==obs) return(S);
  else if(actual == (obs + 1) || actual==(obs-1)) return(EA/2.0);
  return((1.0/(MAX_TITRE-2.0))*(1.0-S-EA));
}


//' @export
//[[Rcpp::export]]
double obs_likelihood(NumericVector y, NumericVector data, NumericVector params){
  double ln = 0;
  int MAX_TITRE = params(2);
  for(int i = 0; i < y.length();++i){
    if(y(i) < 0) y(i) = 0;
    if(y(i) >= MAX_TITRE) y(i) = MAX_TITRE;
    ln += log(obs_error(floor(y(i)), floor(data(i)),params(0),params(1)));
  }
  return ln;
}



//' @export
//[[Rcpp::export]]
NumericVector model_trajectory_cpp(NumericVector pars, NumericVector times){
  double mu = pars[3];
  double tp = pars[4];
  double dp = pars[5];
  double ts = pars[6];
  double m = pars[7];
  double sigma = pars[8];
  double beta = pars[9];
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
  //Rcpp::Rcout << "primed: " << primed << std::endl;
  //Rcpp::Rcout << "beta: " << beta << std::endl;
  //Rcpp::Rcout << "c: " << c << std::endl;
  //Rcpp::Rcout << "x: " << x << std::endl;
  //Rcpp::Rcout << "mu: " << mu << std::endl;
  //Rcpp::Rcout << "cr: " << prime_cr << std::endl;
  
  mu = mu*cr*mod + prime_cr;
  //Rcpp::Rcout << "Combined mu: " << mu << std::endl;


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


//' @export
//[[Rcpp::export]]
double posterior_func_group_cpp(NumericVector pars, NumericVector times, IntegerVector groups, IntegerVector strains,
				  IntegerVector exposure_types, IntegerVector exposure_strains, IntegerVector measured_strains, IntegerVector exposure_orders, IntegerVector exposure_primes,
				  IntegerVector exposure_indices, IntegerVector cr_inds, IntegerVector par_type_ind, IntegerVector order_indices,
				       IntegerVector exposure_i_lengths, IntegerVector par_lengths, IntegerVector cr_lengths, NumericMatrix data){
  int group;
  int strain;
  int A, B, type, order, exposure_strain;
  double t_i, mod, cr,isPrimed;
  double ln = 0;
  IntegerVector tmp_exposures;
  NumericVector fullPars;
  NumericMatrix results(strains.size()*groups.size(), times.size());
  NumericVector y;
  int index = 0;
  for(int i = 0; i < groups.size(); ++i){
    group = groups[i];
    //Rcpp::Rcout << "Group: " << group << std::endl;
    A = exposure_i_lengths[i];
    B = exposure_i_lengths[i+1] - 1;
    tmp_exposures = exposure_indices[Range(A,B)];
    for(int j = 0; j < strains.size(); ++j){
      strain = strains[j] -1;
      //Rcpp::Rcout << "Strain: " << strain << std::endl;
      for(int k = 0; k < tmp_exposures.size(); ++k){
	t_i = pars[tmp_exposures[k]];
	//Rcpp::Rcout << "Exposure: " << t_i << std::endl;
	order = exposure_orders[tmp_exposures[k]] - 1;
	exposure_strain = exposure_strains[tmp_exposures[k]] - 1;
	mod = pars[order_indices[order]];
	isPrimed = exposure_primes[tmp_exposures[k]];

	cr = pars[cr_inds[cr_lengths[strain] + exposure_strain]];
	//Rcpp::Rcout << "Exposure strain: " << exposure_strain << std::endl;
	//Rcpp::Rcout << "cr: " << cr << std::endl;

	type = exposure_types[tmp_exposures[k]] - 1;
	//Rcpp::Rcout << par_lengths[type] << " " << par_lengths[type+1] << std::endl;
	A = par_lengths[type];
	B = par_lengths[type+1] - 1;

	fullPars = pars[par_type_ind[Range(A,B)]];
	fullPars.push_back(isPrimed);
	fullPars.push_back(mod);
	fullPars.push_back(cr);
	fullPars.push_back(t_i);

	//Rcpp::Rcout << "Pars: " << fullPars << std::endl;

	y = model_trajectory_cpp(fullPars, times);
	//Rcpp::Rcout << y << std::endl;
	//Rcpp::Rcout << "Y: " << y << std::endl;
	for(int ii = 0; ii < y.size(); ++ii){
	  results(index, ii) += y[ii];
	}
      }
      ln += obs_likelihood(results(index,_), data(index,_), fullPars);
      index++;
    }
  }
  return ln;
}


//' @export
//[[Rcpp::export]]
NumericMatrix model_func_group_cpp(NumericVector pars, NumericVector times, IntegerVector groups, IntegerVector strains,
				  IntegerVector exposure_types, IntegerVector exposure_strains, IntegerVector measured_strains, IntegerVector exposure_orders, IntegerVector exposure_primes,
				  IntegerVector exposure_indices, IntegerVector cr_inds, IntegerVector par_type_ind, IntegerVector order_indices,
				  IntegerVector exposure_i_lengths, IntegerVector par_lengths, IntegerVector cr_lengths){
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
    //Rcpp::Rcout << "Group: " << group << std::endl;
    A = exposure_i_lengths[i];
    B = exposure_i_lengths[i+1] - 1;
    tmp_exposures = exposure_indices[Range(A,B)];
    for(int j = 0; j < strains.size(); ++j){
      strain = strains[j] -1;
      //Rcpp::Rcout << "Strain: " << strain << std::endl;
      for(int k = 0; k < tmp_exposures.size(); ++k){
	t_i = pars[tmp_exposures[k]];
	//Rcpp::Rcout << "Exposure: " << t_i << std::endl;
	order = exposure_orders[tmp_exposures[k]] - 1;
	exposure_strain = exposure_strains[tmp_exposures[k]] - 1;
	mod = pars[order_indices[order]];
	isPrimed = exposure_primes[tmp_exposures[k]];

	cr = pars[cr_inds[cr_lengths[strain] + exposure_strain]];
	//Rcpp::Rcout << "Exposure strain: " << exposure_strain << std::endl;
	//Rcpp::Rcout << "cr: " << cr << std::endl;

	type = exposure_types[tmp_exposures[k]] - 1;
	//Rcpp::Rcout << par_lengths[type] << " " << par_lengths[type+1] << std::endl;
	A = par_lengths[type];
	B = par_lengths[type+1] - 1;

	fullPars = pars[par_type_ind[Range(A,B)]];
	fullPars.push_back(isPrimed);
	fullPars.push_back(mod);
	fullPars.push_back(cr);
	fullPars.push_back(t_i);

	//Rcpp::Rcout << "Pars: " << fullPars << std::endl;

	y = model_trajectory_cpp(fullPars, times);
	//Rcpp::Rcout << y << std::endl;
	//Rcpp::Rcout << "Y: " << y << std::endl;
	for(int ii = 0; ii < y.size(); ++ii){
	  results(index, ii) += y[ii];
	}
      }
      index++;
    }
  }
  return results;
}
