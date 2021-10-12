#include <RcppArmadillo.h>

#ifndef _MINOR_
#define _MINOR_

template<typename T>
void printR_obj(const T& obj);

double Rcpp_norm(const arma::vec&);

double Rcpp_logSumExp(const arma::vec&);

Rcpp::StringVector Rcpp_get_target(const arma::mat& intervals,
	const arma::vec& positions,const double& pad);

#endif