// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// toUnitScale
double toUnitScale(double x, double min, double max);
RcppExport SEXP _antibodyKinetics_toUnitScale(SEXP xSEXP, SEXP minSEXP, SEXP maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type min(minSEXP);
    Rcpp::traits::input_parameter< double >::type max(maxSEXP);
    rcpp_result_gen = Rcpp::wrap(toUnitScale(x, min, max));
    return rcpp_result_gen;
END_RCPP
}
// fromUnitScale
double fromUnitScale(double x, double min, double max);
RcppExport SEXP _antibodyKinetics_fromUnitScale(SEXP xSEXP, SEXP minSEXP, SEXP maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type min(minSEXP);
    Rcpp::traits::input_parameter< double >::type max(maxSEXP);
    rcpp_result_gen = Rcpp::wrap(fromUnitScale(x, min, max));
    return rcpp_result_gen;
END_RCPP
}
// model_trajectory_cpp
NumericVector model_trajectory_cpp(NumericVector pars, NumericVector times);
RcppExport SEXP _antibodyKinetics_model_trajectory_cpp(SEXP parsSEXP, SEXP timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type times(timesSEXP);
    rcpp_result_gen = Rcpp::wrap(model_trajectory_cpp(pars, times));
    return rcpp_result_gen;
END_RCPP
}
// model_func_group_cpp
NumericMatrix model_func_group_cpp(NumericVector pars, NumericVector times, IntegerVector groups, IntegerVector strains, IntegerVector exposure_indices, IntegerVector exposure_i_lengths, IntegerVector strain_indices, IntegerVector strain_i_lengths, NumericVector exposure_times, IntegerVector exposure_strains, NumericVector exposure_next, IntegerVector exposure_measured, IntegerVector exposure_orders, IntegerVector exposure_primes, IntegerVector cr_inds, IntegerVector par_inds, IntegerVector order_inds, IntegerVector par_lengths, IntegerVector cr_lengths, int version);
RcppExport SEXP _antibodyKinetics_model_func_group_cpp(SEXP parsSEXP, SEXP timesSEXP, SEXP groupsSEXP, SEXP strainsSEXP, SEXP exposure_indicesSEXP, SEXP exposure_i_lengthsSEXP, SEXP strain_indicesSEXP, SEXP strain_i_lengthsSEXP, SEXP exposure_timesSEXP, SEXP exposure_strainsSEXP, SEXP exposure_nextSEXP, SEXP exposure_measuredSEXP, SEXP exposure_ordersSEXP, SEXP exposure_primesSEXP, SEXP cr_indsSEXP, SEXP par_indsSEXP, SEXP order_indsSEXP, SEXP par_lengthsSEXP, SEXP cr_lengthsSEXP, SEXP versionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type times(timesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strains(strainsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type exposure_indices(exposure_indicesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type exposure_i_lengths(exposure_i_lengthsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strain_indices(strain_indicesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strain_i_lengths(strain_i_lengthsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type exposure_times(exposure_timesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type exposure_strains(exposure_strainsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type exposure_next(exposure_nextSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type exposure_measured(exposure_measuredSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type exposure_orders(exposure_ordersSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type exposure_primes(exposure_primesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cr_inds(cr_indsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type par_inds(par_indsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type order_inds(order_indsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type par_lengths(par_lengthsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cr_lengths(cr_lengthsSEXP);
    Rcpp::traits::input_parameter< int >::type version(versionSEXP);
    rcpp_result_gen = Rcpp::wrap(model_func_group_cpp(pars, times, groups, strains, exposure_indices, exposure_i_lengths, strain_indices, strain_i_lengths, exposure_times, exposure_strains, exposure_next, exposure_measured, exposure_orders, exposure_primes, cr_inds, par_inds, order_inds, par_lengths, cr_lengths, version));
    return rcpp_result_gen;
END_RCPP
}
// obs_error
double obs_error(int actual, int obs, double S, double EA, int MAX_TITRE);
RcppExport SEXP _antibodyKinetics_obs_error(SEXP actualSEXP, SEXP obsSEXP, SEXP SSEXP, SEXP EASEXP, SEXP MAX_TITRESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type actual(actualSEXP);
    Rcpp::traits::input_parameter< int >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< double >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type EA(EASEXP);
    Rcpp::traits::input_parameter< int >::type MAX_TITRE(MAX_TITRESEXP);
    rcpp_result_gen = Rcpp::wrap(obs_error(actual, obs, S, EA, MAX_TITRE));
    return rcpp_result_gen;
END_RCPP
}
// norm_error
double norm_error(double actual, int obs, double sd, int MAX_TITRE);
RcppExport SEXP _antibodyKinetics_norm_error(SEXP actualSEXP, SEXP obsSEXP, SEXP sdSEXP, SEXP MAX_TITRESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type actual(actualSEXP);
    Rcpp::traits::input_parameter< int >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< int >::type MAX_TITRE(MAX_TITRESEXP);
    rcpp_result_gen = Rcpp::wrap(norm_error(actual, obs, sd, MAX_TITRE));
    return rcpp_result_gen;
END_RCPP
}
// obs_likelihood
double obs_likelihood(NumericVector y, NumericVector data, NumericVector params);
RcppExport SEXP _antibodyKinetics_obs_likelihood(SEXP ySEXP, SEXP dataSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(obs_likelihood(y, data, params));
    return rcpp_result_gen;
END_RCPP
}
// posterior_func_group_cpp
double posterior_func_group_cpp(NumericVector pars, NumericVector times, IntegerVector groups, IntegerVector strains, IntegerVector exposure_indices, IntegerVector exposure_i_lengths, IntegerVector strain_indices, IntegerVector strain_i_lengths, NumericVector exposure_times, IntegerVector exposure_strains, NumericVector exposure_next, IntegerVector exposure_measured, IntegerVector exposure_orders, IntegerVector exposure_primes, IntegerVector cr_inds, IntegerVector par_inds, IntegerVector order_inds, IntegerVector par_lengths, IntegerVector cr_lengths, int version, IntegerVector individuals, NumericMatrix data);
RcppExport SEXP _antibodyKinetics_posterior_func_group_cpp(SEXP parsSEXP, SEXP timesSEXP, SEXP groupsSEXP, SEXP strainsSEXP, SEXP exposure_indicesSEXP, SEXP exposure_i_lengthsSEXP, SEXP strain_indicesSEXP, SEXP strain_i_lengthsSEXP, SEXP exposure_timesSEXP, SEXP exposure_strainsSEXP, SEXP exposure_nextSEXP, SEXP exposure_measuredSEXP, SEXP exposure_ordersSEXP, SEXP exposure_primesSEXP, SEXP cr_indsSEXP, SEXP par_indsSEXP, SEXP order_indsSEXP, SEXP par_lengthsSEXP, SEXP cr_lengthsSEXP, SEXP versionSEXP, SEXP individualsSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type times(timesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strains(strainsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type exposure_indices(exposure_indicesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type exposure_i_lengths(exposure_i_lengthsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strain_indices(strain_indicesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strain_i_lengths(strain_i_lengthsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type exposure_times(exposure_timesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type exposure_strains(exposure_strainsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type exposure_next(exposure_nextSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type exposure_measured(exposure_measuredSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type exposure_orders(exposure_ordersSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type exposure_primes(exposure_primesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cr_inds(cr_indsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type par_inds(par_indsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type order_inds(order_indsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type par_lengths(par_lengthsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cr_lengths(cr_lengthsSEXP);
    Rcpp::traits::input_parameter< int >::type version(versionSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type individuals(individualsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(posterior_func_group_cpp(pars, times, groups, strains, exposure_indices, exposure_i_lengths, strain_indices, strain_i_lengths, exposure_times, exposure_strains, exposure_next, exposure_measured, exposure_orders, exposure_primes, cr_inds, par_inds, order_inds, par_lengths, cr_lengths, version, individuals, data));
    return rcpp_result_gen;
END_RCPP
}
// posterior_func_group_cpp_matrix
NumericMatrix posterior_func_group_cpp_matrix(NumericVector pars, NumericVector times, IntegerVector groups, IntegerVector strains, IntegerVector exposure_indices, IntegerVector exposure_i_lengths, IntegerVector strain_indices, IntegerVector strain_i_lengths, NumericVector exposure_times, IntegerVector exposure_strains, NumericVector exposure_next, IntegerVector exposure_measured, IntegerVector exposure_orders, IntegerVector exposure_primes, IntegerVector cr_inds, IntegerVector par_inds, IntegerVector order_inds, IntegerVector par_lengths, IntegerVector cr_lengths, int version, IntegerVector individuals, NumericMatrix data);
RcppExport SEXP _antibodyKinetics_posterior_func_group_cpp_matrix(SEXP parsSEXP, SEXP timesSEXP, SEXP groupsSEXP, SEXP strainsSEXP, SEXP exposure_indicesSEXP, SEXP exposure_i_lengthsSEXP, SEXP strain_indicesSEXP, SEXP strain_i_lengthsSEXP, SEXP exposure_timesSEXP, SEXP exposure_strainsSEXP, SEXP exposure_nextSEXP, SEXP exposure_measuredSEXP, SEXP exposure_ordersSEXP, SEXP exposure_primesSEXP, SEXP cr_indsSEXP, SEXP par_indsSEXP, SEXP order_indsSEXP, SEXP par_lengthsSEXP, SEXP cr_lengthsSEXP, SEXP versionSEXP, SEXP individualsSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type times(timesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strains(strainsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type exposure_indices(exposure_indicesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type exposure_i_lengths(exposure_i_lengthsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strain_indices(strain_indicesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strain_i_lengths(strain_i_lengthsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type exposure_times(exposure_timesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type exposure_strains(exposure_strainsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type exposure_next(exposure_nextSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type exposure_measured(exposure_measuredSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type exposure_orders(exposure_ordersSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type exposure_primes(exposure_primesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cr_inds(cr_indsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type par_inds(par_indsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type order_inds(order_indsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type par_lengths(par_lengthsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cr_lengths(cr_lengthsSEXP);
    Rcpp::traits::input_parameter< int >::type version(versionSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type individuals(individualsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(posterior_func_group_cpp_matrix(pars, times, groups, strains, exposure_indices, exposure_i_lengths, strain_indices, strain_i_lengths, exposure_times, exposure_strains, exposure_next, exposure_measured, exposure_orders, exposure_primes, cr_inds, par_inds, order_inds, par_lengths, cr_lengths, version, individuals, data));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_antibodyKinetics_toUnitScale", (DL_FUNC) &_antibodyKinetics_toUnitScale, 3},
    {"_antibodyKinetics_fromUnitScale", (DL_FUNC) &_antibodyKinetics_fromUnitScale, 3},
    {"_antibodyKinetics_model_trajectory_cpp", (DL_FUNC) &_antibodyKinetics_model_trajectory_cpp, 2},
    {"_antibodyKinetics_model_func_group_cpp", (DL_FUNC) &_antibodyKinetics_model_func_group_cpp, 20},
    {"_antibodyKinetics_obs_error", (DL_FUNC) &_antibodyKinetics_obs_error, 5},
    {"_antibodyKinetics_norm_error", (DL_FUNC) &_antibodyKinetics_norm_error, 4},
    {"_antibodyKinetics_obs_likelihood", (DL_FUNC) &_antibodyKinetics_obs_likelihood, 3},
    {"_antibodyKinetics_posterior_func_group_cpp", (DL_FUNC) &_antibodyKinetics_posterior_func_group_cpp, 22},
    {"_antibodyKinetics_posterior_func_group_cpp_matrix", (DL_FUNC) &_antibodyKinetics_posterior_func_group_cpp_matrix, 22},
    {NULL, NULL, 0}
};

RcppExport void R_init_antibodyKinetics(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
