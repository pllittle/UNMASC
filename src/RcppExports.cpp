// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Rcpp_LL
arma::vec Rcpp_LL(const arma::mat& RD, const arma::vec& DP, const arma::vec& LBC, const double& pi, const double& psi);
RcppExport SEXP _UNMASC_Rcpp_LL(SEXP RDSEXP, SEXP DPSEXP, SEXP LBCSEXP, SEXP piSEXP, SEXP psiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type RD(RDSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type DP(DPSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type LBC(LBCSEXP);
    Rcpp::traits::input_parameter< const double& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const double& >::type psi(psiSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_LL(RD, DP, LBC, pi, psi));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_BAF_clust
Rcpp::List Rcpp_BAF_clust(const arma::mat& RD, const arma::vec& DP, const arma::vec& log_DP, const arma::vec& LBC, const arma::vec& props0, const double& mm0, const double& psi0, const arma::uword& max_iter, const double& eps, const bool& show);
RcppExport SEXP _UNMASC_Rcpp_BAF_clust(SEXP RDSEXP, SEXP DPSEXP, SEXP log_DPSEXP, SEXP LBCSEXP, SEXP props0SEXP, SEXP mm0SEXP, SEXP psi0SEXP, SEXP max_iterSEXP, SEXP epsSEXP, SEXP showSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type RD(RDSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type DP(DPSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type log_DP(log_DPSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type LBC(LBCSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type props0(props0SEXP);
    Rcpp::traits::input_parameter< const double& >::type mm0(mm0SEXP);
    Rcpp::traits::input_parameter< const double& >::type psi0(psi0SEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type show(showSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_BAF_clust(RD, DP, log_DP, LBC, props0, mm0, psi0, max_iter, eps, show));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_logSumExp
double Rcpp_logSumExp(const arma::vec& log_x);
RcppExport SEXP _UNMASC_Rcpp_logSumExp(SEXP log_xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type log_x(log_xSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_logSumExp(log_x));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_get_target
Rcpp::StringVector Rcpp_get_target(const arma::mat& intervals, const arma::vec& positions, const double& pad);
RcppExport SEXP _UNMASC_Rcpp_get_target(SEXP intervalsSEXP, SEXP positionsSEXP, SEXP padSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type intervals(intervalsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type positions(positionsSEXP);
    Rcpp::traits::input_parameter< const double& >::type pad(padSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_get_target(intervals, positions, pad));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_RDQ
arma::mat Rcpp_RDQ(const arma::mat& RDQ, const double& qq);
RcppExport SEXP _UNMASC_Rcpp_RDQ(SEXP RDQSEXP, SEXP qqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type RDQ(RDQSEXP);
    Rcpp::traits::input_parameter< const double& >::type qq(qqSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_RDQ(RDQ, qq));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_UNMASC_Rcpp_LL", (DL_FUNC) &_UNMASC_Rcpp_LL, 5},
    {"_UNMASC_Rcpp_BAF_clust", (DL_FUNC) &_UNMASC_Rcpp_BAF_clust, 10},
    {"_UNMASC_Rcpp_logSumExp", (DL_FUNC) &_UNMASC_Rcpp_logSumExp, 1},
    {"_UNMASC_Rcpp_get_target", (DL_FUNC) &_UNMASC_Rcpp_get_target, 3},
    {"_UNMASC_Rcpp_RDQ", (DL_FUNC) &_UNMASC_Rcpp_RDQ, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_UNMASC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
