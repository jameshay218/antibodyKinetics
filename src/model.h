#include <Rcpp.h>
using namespace Rcpp;

#ifndef MODEL_TRAJECTORY_CPP
#define MODEL_TRAJECTORY_CPP
NumericVector model_trajectory_cpp(NumericVector pars, NumericVector times);
#endif

#ifndef MODEL_FUNC_GROUP_CPP
#define MODEL_FUNC_GROUP_CPP
NumericMatrix model_func_group_cpp(NumericVector pars, NumericVector times, 
				   IntegerVector groups, IntegerVector strains,
				   IntegerVector exposure_indices, IntegerVector exposure_i_lengths,
				   IntegerVector strain_indices, IntegerVector strain_i_lengths,
				   NumericVector exposure_times, IntegerVector exposure_strains,
				   NumericVector exposure_next, IntegerVector exposure_measured,
				   IntegerVector exposure_orders, IntegerVector exposure_primes, 
				   IntegerVector cr_inds, IntegerVector par_inds,
				   CharacterVector par_names,
				   IntegerVector order_inds, IntegerVector par_lengths,
				   IntegerVector cr_lengths, int version);
#endif
