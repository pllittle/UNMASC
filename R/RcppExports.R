# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

Rcpp_LL <- function(RD, DP, LBC, pi, psi) {
    .Call('_UNMASC_Rcpp_LL', PACKAGE = 'UNMASC', RD, DP, LBC, pi, psi)
}

Rcpp_BAF_clust <- function(RD, DP, log_DP, LBC, props0, mm0, psi0, max_iter = 4e3L, eps = 1e-5, show = TRUE) {
    .Call('_UNMASC_Rcpp_BAF_clust', PACKAGE = 'UNMASC', RD, DP, log_DP, LBC, props0, mm0, psi0, max_iter, eps, show)
}

Rcpp_logSumExp <- function(log_x) {
    .Call('_UNMASC_Rcpp_logSumExp', PACKAGE = 'UNMASC', log_x)
}

Rcpp_get_target <- function(intervals, positions, pad = 0) {
    .Call('_UNMASC_Rcpp_get_target', PACKAGE = 'UNMASC', intervals, positions, pad)
}

Rcpp_RDQ <- function(RDQ, qq = 0.25) {
    .Call('_UNMASC_Rcpp_RDQ', PACKAGE = 'UNMASC', RDQ, qq)
}

