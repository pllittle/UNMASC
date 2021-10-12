#include "../inst/include/minor.h"

// [[Rcpp::depends("RcppArmadillo")]]

// ----------
// Minor functions
// ----------

template<typename T>
void printR_obj(const T& obj){
	Rcpp::Rcout << obj << std::endl;	
}

double Rcpp_norm(const arma::vec& a){
	return arma::norm(a);
}

// [[Rcpp::export]]
double Rcpp_logSumExp(const arma::vec& log_x){
	if( log_x.n_elem == 1 ){
		return log_x.at(0);
	} else {
		double max_val = arma::max(log_x);
		arma::vec log_x_2 = log_x - max_val;
		return std::log(arma::sum(arma::exp(log_x_2))) + max_val;
	}
}

// [[Rcpp::export]]
Rcpp::StringVector Rcpp_get_target(const arma::mat& intervals,
	const arma::vec& positions,const double& pad = 0){
	
	arma::uword num_elem = positions.n_elem, ii;
	Rcpp::StringVector xx(num_elem);
	arma::uvec tmp_index = arma::zeros<arma::uvec>(1);
	
	for(ii = 0; ii < num_elem; ii++){
		tmp_index = arma::find( intervals.col(0) - pad <= positions.at(ii) 
			&& positions.at(ii) <= intervals.col(1) + pad );
		if(tmp_index.n_elem > 0){
			xx(ii) = "YES";
		} else {
			xx(ii) = "NO";
		}
	}
	
	return xx;
}


// ----------
// RDQ/Quantile functions
// ----------

double my_quantile(const arma::vec& xx,const double& pp,
	const arma::uword& type = 7){
	
	double mm, gg, gamma;
	arma::uword jj, nn = xx.n_elem;
	arma::vec s_xx = arma::sort(xx);
	
	// Calculate mm
	if(type == 1 || type == 2){
		mm = 0;
	} else if(type == 3){
		mm = -0.5;
	} else if(type == 4){
		mm = 0;
	} else if(type == 5){
		mm = 0.5;
	} else if(type == 6){
		mm = pp;
	} else if(type == 7){
		mm = 1.0 - pp;
	} else if(type == 8){
		mm = (pp + 1.0) / 3.0;
	} else {
		mm = 0.25 * pp + 3.0/8.0;
	}
	
	// Calculate jj and gg
	jj = std::floor(nn * pp + mm);
	gg = nn * pp + mm - jj;
	
	// Calculate gamma
	gamma = 1.0;
	if(type == 1){
		if(gg == 0.0) gamma = 0;
	} else if(type == 2){
		if(gg == 0.0) gamma = 0.5;
	} else if(type == 3){
		if(gg == 0.0 && jj % 2 == 0) gamma = 0;
	} else {
		gamma = gg;
	}
	
	return (1.0 - gamma) * s_xx.at(jj - 1) + gamma * s_xx.at(jj);
}

// [[Rcpp::export]]
arma::mat Rcpp_RDQ(const arma::mat& RDQ,const double& qq = 0.25){
	// matrix RDQ columns = c(mutID,tAD,tDP,nDP,Q)
	// qq = quantile
	
	// Get unique mutID
	arma::vec uniq_mutID = arma::unique(RDQ.col(0));

	// Create output matrix
	arma::mat out = arma::zeros<arma::mat>(uniq_mutID.n_elem,5);
	arma::uword ii;
	for(ii = 0; ii < uniq_mutID.n_elem; ii++){
		arma::uvec tmp_index = arma::find(RDQ.col(0) == uniq_mutID.at(ii));
		arma::mat tmp_mat = RDQ.rows(tmp_index);
		out.at(ii,0) = uniq_mutID.at(ii);
		out.at(ii,1) = tmp_mat.at(0,1);
		out.at(ii,2) = tmp_mat.at(0,2);
		out.at(ii,3) = my_quantile(tmp_mat.col(3),qq);
		out.at(ii,4) = my_quantile(tmp_mat.col(4),qq);
	}
	
	return out;
}


// ----------
// Unused functions
// ----------

double Rcpp_LL_binom(const double& xx,const double& nn,
	const double& lbc,const double& pi){
	
	return lbc + xx * std::log(pi) + 
		(nn - xx) * std::log(1.0 - pi);
}

double Rcpp_LL_betaBinom(const double& xx,const double& nn,
	const double& lbc,const double& aa,const double& bb){
	
	// note: pi = aa / (aa + bb), psi = 1 / (aa + bb)
	// note: aa = pi / psi, bb = (1 - pi) / psi
	
	return lbc + std::lgamma(xx + aa) +
		std::lgamma(nn - xx + bb) +
		std::lgamma(aa + bb) -
		std::lgamma(nn + aa + bb) -
		std::lgamma(aa) - std::lgamma(bb);
}

arma::vec Rcpp_LL_betaBinom2(const arma::mat& RD,const arma::vec& DP,
	const arma::vec& LBC,const double& aa,const double& bb){
	
	// note: pi = aa / (aa + bb), psi = 1 / (aa + bb)
	// note: aa = pi / psi, bb = (1 - pi) / psi
	
	return LBC + arma::lgamma(RD.col(0) + aa) +
		arma::lgamma(RD.col(1) + bb) +
		std::lgamma(aa + bb) -
		arma::lgamma(DP + aa + bb) -
		std::lgamma(aa) - std::lgamma(bb);
}

double ZOISB_LLZ(const arma::mat& RD,
	const arma::vec& DP,const arma::vec& log_DP,
	const arma::vec& LBC,const arma::mat& log_IND,
	arma::mat& Z,arma::vec& params){

	// E-step
	Z.col(0) = std::log(params.at(0)) - log_DP;
	Z.col(1) = std::log(params.at(1)) + LBC + DP * std::log(0.5);
	Z.col(2) = std::log(params.at(2)) + log_IND.col(0);
	Z.col(4) = std::log(params.at(4)) + log_IND.col(1);
	if( params.at(6) == 0.0 ){
		Z.col(3).fill(-1.0 * arma::datum::inf);
		Z.col(5).fill(-1.0 * arma::datum::inf);
	} else {
		Z.col(3) = std::log(params.at(3)) + LBC + 
			RD.col(0) * std::log(params.at(6)) +
			RD.col(1) * std::log(1.0 - params.at(6));
		Z.col(5) = std::log(params.at(5)) + LBC + 
			RD.col(0) * std::log(1.0 - params.at(6)) +
			RD.col(1) * std::log(params.at(6));		
	}
	
	double LL = 0.0, tmp_num;
	arma::uword bb, B = DP.n_elem;
	for(bb = 0; bb < B; bb++){
		tmp_num = Rcpp_logSumExp(Z.row(bb).t());
		LL += tmp_num;
		Z.row(bb) = arma::exp(Z.row(bb).t() - tmp_num).t();
	}
	
	// M-step
	// update alpha
	params.subvec(0,5) = arma::mean(Z,0).t();
	// update pp
	double pp_den = arma::sum( (Z.col(3) + Z.col(5)) % DP );
	if( pp_den > 0 ){
		params.at(6) = arma::sum(Z.col(3) % RD.col(0) + Z.col(5) % RD.col(1)) /
			pp_den;
	} else {
		params.at(6) = 0.0;
	}
	
	return LL;
}

arma::vec ZOISB_alpha2props(const arma::vec& alpha){
	arma::vec props = arma::zeros<arma::vec>(5);
	
	props.at(0) = alpha.at(0);
	if( alpha.at(0) < 1 ){
		props.at(1) = alpha.at(1) / 
			(1.0 - alpha.at(0));
	}
	if( alpha.at(0) + alpha.at(1) < 1 ){
		props.at(2) = ( alpha.at(2) + alpha.at(3) ) / 
			(1.0 - (alpha.at(0) + alpha.at(1)));
	}
	if( alpha.at(2) + alpha.at(3) > 0 ){
		props.at(3) = alpha.at(2) / 
			(alpha.at(2) + alpha.at(3));
	}
	if( alpha.at(4) + alpha.at(5) > 0 ){
		props.at(4) = alpha.at(4) /
			(alpha.at(4) + alpha.at(5));
	}
	
	return props;
}

Rcpp::List ZOISB_EM(const arma::mat& RD, 
	const arma::vec& DP,const arma::vec& log_DP,
	const arma::vec& LBC, const arma::mat& log_IND,
	const arma::vec& params0, 
	const arma::uword& max_iter = 4e3,
	const double& conv = 1e-7,const bool& show = false){

	if( arma::min(DP) <= 0 ){
		Rcpp::stop("min(DP) = 0");
	}
	
	arma::uword B = DP.n_elem, cc, ms,
		iter = 0;
	arma::mat Z = arma::zeros<arma::mat>(B,6);
	arma::vec params = params0;
	double old_LL = 0.0, BIC, curr_LL = 0.0;
	arma::vec curr_params = arma::zeros<arma::vec>(7);
	arma::uvec classify = arma::zeros<arma::uvec>(B);
	
	while(iter < max_iter){
		// EM
		old_LL = ZOISB_LLZ(RD, DP, log_DP, LBC, log_IND, Z, params);

		if( iter > 0 ){
			// Adjust
			classify = arma::index_max(Z,1);
			for(cc = 0; cc < 6; cc++){
				if( arma::sum(classify == cc) == 0 ){
					params.at(cc) = 0.0;
				}
			}
			params.subvec(0,5) = params.subvec(0,5) /
				arma::sum(params.subvec(0,5));
			
			// Check convergence
			if( std::abs(old_LL - curr_LL) < conv && 
				Rcpp_norm(params - curr_params) < conv ){
				break;
			}
		}
		
		if(show){
			Rcpp::Rcout << "Iter=" << iter+1 << ": LL=" << old_LL << "\n";
			Rcpp::Rcout << "\tParams=" << params.t() << "\n";
		}
		
		curr_LL = old_LL;
		curr_params = params;
		iter += 1;
	}
	
	std::string converge;
	BIC = 2*old_LL - std::log(B)*ms;
	if(iter < max_iter){
		converge = "yes";
		old_LL = ZOISB_LLZ(RD, DP, log_DP, LBC, log_IND, Z, params);
		ms = arma::sum( params.subvec(0,5) > 0 ) - 1 +
			1*( params.at(6) > 0 && params.at(6) < 1 
			&& params.at(3) + params.at(5) > 0 );
		if( ms == 0 ) ms = 1;
		BIC = 2*old_LL - std::log(B)*ms;
		classify = arma::index_max(Z,1) + 1;
		arma::vec props = ZOISB_alpha2props(params.subvec(0,5));
		return Rcpp::List::create(
			Rcpp::Named("iter",iter),
			Rcpp::Named("converge",converge),
			Rcpp::Named("params0",Rcpp::NumericVector(params0.begin(),params0.end())),
			Rcpp::Named("LL",old_LL),
			Rcpp::Named("ms",ms),
			Rcpp::Named("BIC",BIC),
			Rcpp::Named("Z",Z),
			Rcpp::Named("class",Rcpp::NumericVector(classify.begin(),classify.end())),
			Rcpp::Named("params",Rcpp::NumericVector(params.begin(),params.end())),
			Rcpp::Named("props",Rcpp::NumericVector(props.begin(),props.end())),
			Rcpp::Named("maf",params.at(6))
			);
	} else {
		converge = "no";
		return Rcpp::List::create(
			Rcpp::Named("BIC",BIC),
			Rcpp::Named("converge",converge));
	}

}

// ---
