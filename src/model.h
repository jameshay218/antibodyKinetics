#include <Rcpp.h>
using namespace Rcpp;

NumericVector model_trajectory_cpp(NumericVector pars, NumericVector times);
NumericMatrix model_func_group_cpp(NumericVector pars, NumericVector times, 
				   IntegerVector groups, IntegerVector strains,
				   IntegerVector exposure_types, IntegerVector exposure_strains, 
				   IntegerVector measured_strains, IntegerVector exposure_orders, 
				   IntegerVector exposure_primes, IntegerVector exposure_indices, 
				   IntegerVector cr_inds, IntegerVector par_type_ind, 
				   IntegerVector order_indices, IntegerVector exposure_i_lengths, 
				   IntegerVector par_lengths, IntegerVector cr_lengths, 
				   int version);
