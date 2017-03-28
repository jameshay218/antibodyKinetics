// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// predict_titres
NumericMatrix predict_titres(NumericVector params, NumericVector times);
RcppExport SEXP antibodyKinetics_predict_titres(SEXP paramsSEXP, SEXP timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type times(timesSEXP);
    rcpp_result_gen = Rcpp::wrap(predict_titres(params, times));
    return rcpp_result_gen;
END_RCPP
}
// toUnitScale
double toUnitScale(double x, double min, double max);
RcppExport SEXP antibodyKinetics_toUnitScale(SEXP xSEXP, SEXP minSEXP, SEXP maxSEXP) {
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
RcppExport SEXP antibodyKinetics_fromUnitScale(SEXP xSEXP, SEXP minSEXP, SEXP maxSEXP) {
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
// obs_error
double obs_error(int actual, int obs, double S, double EA);
RcppExport SEXP antibodyKinetics_obs_error(SEXP actualSEXP, SEXP obsSEXP, SEXP SSEXP, SEXP EASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type actual(actualSEXP);
    Rcpp::traits::input_parameter< int >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< double >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type EA(EASEXP);
    rcpp_result_gen = Rcpp::wrap(obs_error(actual, obs, S, EA));
    return rcpp_result_gen;
END_RCPP
}
// posterior
double posterior(NumericVector y, NumericVector data, NumericVector params);
RcppExport SEXP antibodyKinetics_posterior(SEXP ySEXP, SEXP dataSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(posterior(y, data, params));
    return rcpp_result_gen;
END_RCPP
}
