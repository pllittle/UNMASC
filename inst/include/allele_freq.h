#include <RcppArmadillo.h>
#include "minor.h"

#ifndef _AF_
#define _AF_

arma::vec Rcpp_LL(const arma::mat& RD,
	const arma::vec& DP,const arma::vec& LBC,
	const double& pi,const double& psi);
Rcpp::List Rcpp_BAF_clust(const arma::mat& RD,
	const arma::vec& DP,const arma::vec& log_DP,
	const arma::vec& LBC,const arma::vec& props0,
	const double& mm0,const double& psi0,
	const arma::uword& max_iter,
	const double& eps,const bool& show);

#endif
