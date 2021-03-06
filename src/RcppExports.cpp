// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// boot_ardl_c
Rcpp::List boot_ardl_c(arma::mat r_in, arma::mat GAMMAX, arma::mat A, arma::mat start_z, arma::vec omegat, arma::vec interc, arma::vec trend);
RcppExport SEXP _bootCT_boot_ardl_c(SEXP r_inSEXP, SEXP GAMMAXSEXP, SEXP ASEXP, SEXP start_zSEXP, SEXP omegatSEXP, SEXP intercSEXP, SEXP trendSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type r_in(r_inSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type GAMMAX(GAMMAXSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type start_z(start_zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omegat(omegatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type interc(intercSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trend(trendSEXP);
    rcpp_result_gen = Rcpp::wrap(boot_ardl_c(r_in, GAMMAX, A, start_z, omegat, interc, trend));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bootCT_boot_ardl_c", (DL_FUNC) &_bootCT_boot_ardl_c, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_bootCT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
