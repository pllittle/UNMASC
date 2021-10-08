#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <omp.h>

// [[Rcpp::depends("RcppArmadillo")]]

template<typename T>
void printR_obj(const T& obj){
	Rcpp::Rcout << obj << std::endl;	
}


// --------------------
// Intermediate Functions
// --------------------

// [[Rcpp::export]]
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


// --------------------
// Binomial vs Beta-Binomial functions and optimizations

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

// [[Rcpp::export]]
arma::vec Rcpp_LL(const arma::mat& RD,const arma::vec& DP,
	const arma::vec& LBC,const double& pi,const double& psi){
	
	// Binomial or Beta-Binomial
	if( psi == 0.0 ){
		return LBC + RD.col(0) * std::log(pi) + RD.col(1) * std::log(1.0 - pi);
	} else {
		double aa = pi / psi, bb = (1.0 - pi) / psi;
		return LBC + arma::lgamma(RD.col(0) + aa) +
			arma::lgamma(RD.col(1) + bb) +
			std::lgamma(aa + bb) - arma::lgamma(DP + aa + bb) -
			std::lgamma(aa) - std::lgamma(bb);
	}
}

double Rcpp_ecLL(const arma::mat& RD,const arma::vec& DP,
	const arma::vec& LBC,const arma::mat& Z,const arma::vec& PARS){
	
	arma::uword cc;
	double psi = std::exp(PARS.at(1)),mm = std::exp(PARS.at(0));
	mm /= 1.0 + mm;
	double LL = 0.0;
	
	for(cc = 2; cc < 5; cc++){
		if(cc == 2){
			LL += arma::sum(Z.col(cc) % Rcpp_LL(RD,DP,LBC,mm,psi));
		} else if(cc == 3){
			LL += arma::sum(Z.col(cc) % Rcpp_LL(RD,DP,LBC,0.5,psi));
		} else {
			LL += arma::sum(Z.col(cc) % Rcpp_LL(RD,DP,LBC,1.0 - mm,psi));
		}
	}
	
	return LL;
}

arma::vec Rcpp_ecGRAD(const arma::mat& RD,const arma::vec& DP,
	const arma::mat& Z,const arma::vec& PARS){
	
	arma::uword ii, cc, NN = RD.n_rows;
	double psi = std::exp(PARS.at(1)),mm = std::exp(PARS.at(0));
	mm /= 1.0 + mm;
	double ipsi = 1.0 / psi, part_ipsi = -1.0 * ipsi * ipsi,
		aa = mm * ipsi,bb = (1.0 - mm) * ipsi,
		dd = 0.5 * ipsi,part_aa_mm = ipsi,part_bb_mm = -1.0 * part_aa_mm,
		dg_aa = R::digamma(aa),dg_bb = R::digamma(bb),dg_dd = R::digamma(dd),
		dg_ipsi = R::digamma(ipsi),
		part_aa_psi = part_ipsi * mm,part_bb_psi = part_ipsi * (1.0 - mm),
		part_dd_psi = 0.5 * part_ipsi,
		dg_xx_aa,dg_nx_bb,dg_xx_bb,dg_nx_aa,dg_xx_dd,dg_nx_dd,dg_nn_ipsi;
	arma::vec GRAD = arma::zeros<arma::vec>(2); // mm,psi
	
	for(ii = 0; ii < NN; ii++){
		dg_xx_aa = R::digamma(RD.at(ii,0) + aa);
		dg_nx_bb = R::digamma(RD.at(ii,1) + bb);
		dg_xx_bb = R::digamma(RD.at(ii,0) + bb);
		dg_nx_aa = R::digamma(RD.at(ii,1) + aa);
		dg_xx_dd = R::digamma(RD.at(ii,0) + dd);
		dg_nx_dd = R::digamma(RD.at(ii,1) + dd);
		dg_nn_ipsi = R::digamma(DP.at(ii) + ipsi);
	for(cc = 2; cc < 5; cc++){
		// part unc_mm
		if(cc == 2){
			GRAD.at(0) += Z.at(ii,cc) * ( dg_xx_aa * part_aa_mm + 
				dg_nx_bb * part_bb_mm - dg_aa * part_aa_mm - dg_bb * part_bb_mm );
		}
		if(cc == 4){
			GRAD.at(0) += Z.at(ii,cc) * ( dg_xx_bb * part_bb_mm +
				dg_nx_aa * part_aa_mm - dg_aa * part_aa_mm - dg_bb * part_bb_mm );
		}
		
		// part log(psi)
		if(cc == 2){
			GRAD.at(1) += Z.at(ii,cc) * ( dg_xx_aa * part_aa_psi +
				dg_nx_bb * part_bb_psi - dg_nn_ipsi * part_ipsi +
				dg_ipsi * part_ipsi - dg_aa * part_aa_psi -
				dg_bb * part_bb_psi );
		} else if(cc == 3){
			GRAD.at(1) += Z.at(ii,cc) * ( dg_xx_dd * part_dd_psi +
				dg_nx_dd * part_dd_psi - dg_nn_ipsi * part_ipsi +
				dg_ipsi * part_ipsi - 2.0 * dg_dd * part_dd_psi );
		} else {
			GRAD.at(1) += Z.at(ii,cc) * ( dg_xx_bb * part_bb_psi +
				dg_nx_aa * part_aa_psi - dg_nn_ipsi * part_ipsi +
				dg_ipsi * part_ipsi - dg_aa * part_aa_psi -
				dg_bb * part_bb_psi );
		}
		
	}}
	
	GRAD.at(0) *= mm * (1.0 - mm);
	GRAD.at(1) *= psi;
	return GRAD;
}

void Rcpp_ecBFGS(const arma::mat& RD,const arma::vec& DP,
	const arma::vec& LBC,const arma::mat& Z,arma::vec& PARS,
	arma::uword& converge,const arma::mat& I_np,
	const arma::uword& max_iter = 4e3,const double& eps = 1e-10,
	const double& gr_eps = 1e-2,const bool& show = true){
	
	converge = 0;
	arma::uword iter = 0, jj, uu, reset_Bk = 0,np = 2;
	arma::mat inv_Bk = I_np, ISYT = inv_Bk;
	arma::vec xk = PARS, curr_xk = arma::zeros<arma::vec>(np),
		new_xk = curr_xk,gr_k = curr_xk,p_k = curr_xk,
		s_k = curr_xk,y_k = curr_xk;
	double old_LL,new_LL,inv_norm_p_k,tmp_step,ys,fnscale = -1.0,
		curr_LL = 0.0,mm,psi,iBG;
	mm = std::exp(PARS.at(0));
	mm /= 1.0 + mm;
	double orig_LL = Rcpp_ecLL(RD,DP,LBC,Z,xk);
	
	while(iter < max_iter){
		// Calculate Direction p_k
		gr_k = fnscale * Rcpp_ecGRAD(RD,DP,Z,xk);
		
		p_k = -1.0 * inv_Bk * gr_k;
		inv_norm_p_k = 1.0 / std::max(1.0,Rcpp_norm(p_k));

		// Line search for new xk
		uu = 0;
		old_LL = fnscale * Rcpp_ecLL(RD,DP,LBC,Z,xk);
		
		// Less than 0 b/c fnscale = -1
		if( false && old_LL < 0.0 ){
			break; 
		}
		
		for(jj = 0; jj <= 30; jj++){
			tmp_step 	= inv_norm_p_k / std::pow(2,jj);
			new_xk 		= xk + tmp_step * p_k;
			new_LL = fnscale * Rcpp_ecLL(RD,DP,LBC,Z,new_xk);
			if( new_LL < old_LL ){ // minimizing
				s_k = tmp_step * p_k;
				y_k = fnscale * Rcpp_ecGRAD(RD,DP,Z,new_xk) - gr_k;
				if( y_k.has_nan() ) continue;
				ys = arma::dot(y_k,s_k);
				if( ys > 0.0 ){
					ISYT = I_np - (s_k * y_k.t()) / ys;
					inv_Bk = ISYT * inv_Bk * ISYT.t() + s_k * s_k.t() / ys;
				}
				
				xk = new_xk;
				old_LL = new_LL;
				uu = 1;
				break;
			}
		}
		
		if( uu == 0 ) { // aka no update
			if( Rcpp_norm(gr_k) > 1.0 ){
				if( show ) printR_obj("\tReset inv_Bk");
				inv_Bk = I_np;
				reset_Bk++;
			} else {
				if( show ) printR_obj("\tFailed line search");
				break;
			}
		}
		
		if( reset_Bk > 5 ) break;
		
		// Check Convergence
		if( iter > 0 ){
			if( std::abs(curr_LL - old_LL) < eps &&
				Rcpp_norm(curr_xk - xk) < eps ){
				gr_k = Rcpp_ecGRAD(RD,DP,Z,xk);
				if( Rcpp_norm(gr_k) < eps ){
					break;
				}
			}
		}
		
		curr_xk = xk;
		curr_LL = old_LL;
		iter++;
	
	}
	
	// Update parameters
	PARS = xk;
	mm = std::exp(PARS.at(0));
	mm /= 1.0 + mm;
	psi = std::exp(PARS.at(1));
	
	// Calculate LL, GRAD
	old_LL = Rcpp_ecLL(RD,DP,LBC,Z,PARS);
	gr_k = Rcpp_ecGRAD(RD,DP,Z,PARS);
	p_k = -1.0 * inv_Bk * gr_k;
	iBG = Rcpp_norm(p_k);
	if( iBG < gr_eps ){
		converge = 1;
	}
	
	if( show ){
		Rcpp::Rcout << "\tLL_old = " << orig_LL << "; Iter = " << iter 
			<< "; LL_new = " << old_LL << "\n";
		double nG = Rcpp_norm(gr_k);
		Rcpp::Rcout << "\tConvergence Indicators: \n"
			<< "\t   NormGrad = " << nG << ";\n"
			<< "\t   NormiBkGr = " << iBG << ";\n";
		Rcpp::Rcout << "\t   MAF = " << mm << "\n";
		Rcpp::Rcout << "\t   PSI = " << psi << "\n";
		Rcpp::Rcout << "\t   Convergence Status = ";
		if( converge == 1 ){
			Rcpp::Rcout << "YES\n";
		} else {
			Rcpp::Rcout << "NO\n";
		}
		Rcpp::Rcout << "\n";
	}
	
}

double Rcpp_LLZ(const arma::mat& RD,const arma::vec& DP,
	const arma::vec& log_DP,const arma::vec& LBC,
	const arma::mat& log_IND,const arma::mat& I_np,
	arma::mat& Z,arma::vec& props,
	double& mm,double& psi,const bool& show = true){
	
	arma::uword ii, NN = RD.n_rows;
	
	// E-step
	if( show ) Rcpp::Rcout << "  Run E step ...\n";
	Z.col(0) = std::log(props.at(0)) - log_DP;
	Z.col(1) = std::log(props.at(1)) + log_IND.col(0);
	Z.col(2) = std::log(props.at(2)) + Rcpp_LL(RD,DP,LBC,mm,psi);
	Z.col(3) = std::log(props.at(3)) + Rcpp_LL(RD,DP,LBC,0.5,psi);
	Z.col(4) = std::log(props.at(4)) + Rcpp_LL(RD,DP,LBC,1.0 - mm,psi);
	Z.col(5) = std::log(props.at(5)) + log_IND.col(1);
	
	// Calculate observed LL
	double LL = 0.0, tmp_num;
	for(ii = 0; ii < NN; ii++){
		tmp_num = Rcpp_logSumExp(Z.row(ii).t());
		LL += tmp_num;
		Z.row(ii) = arma::exp(Z.row(ii).t() - tmp_num).t();
	}
	
	// M-steps
	if( show ) Rcpp::Rcout << "  Run M step ...\n";
	
	// update proportions
	props = arma::mean(Z,0).t();
	
	// If psi = 0, binomial mixtures, use code below
	if( psi == 0.0 ){
		// update mm
		double pp_den = arma::sum( (Z.col(2) + Z.col(4)) % DP );
		if( pp_den > 0 ){
			mm = arma::sum(Z.col(2) % RD.col(0) + Z.col(4) % RD.col(1)) /
				pp_den;
		} else {
			mm = 0.0;
		}
	} else {
		arma::vec PARS = arma::zeros<arma::vec>(2);
		PARS.at(0) = std::log(mm / (1.0 - mm));
		PARS.at(1) = std::log(psi);
		arma::uword converge = 0;
		Rcpp_ecBFGS(RD,DP,LBC,Z,PARS,converge,I_np,4e3,1e-5,1e-2,show);
		mm = std::exp(PARS.at(0));
		mm /= 1.0 + mm;
		psi = std::exp(PARS.at(1));
	}
	
	return LL;
}

Rcpp::List Rcpp_ZOISB_EM(const arma::mat& RD,const arma::vec& DP,
	const arma::vec& log_DP,const arma::vec& LBC,const arma::mat& log_IND,
	const arma::vec& props0,const double& mm0,const double& psi0,
	const arma::uword& max_iter = 4e3,const double& conv = 1e-7,
	const bool& show = false){

	if( arma::min(DP) <= 0 ){
		Rcpp::stop("min(DP) = 0");
	}
	
	if( show ) Rcpp::Rcout << "Set constants ...\n";
	double old_LL = 0.0,BIC,curr_LL = 0.0,mm = mm0,psi = psi0;
	arma::uword NN = DP.n_elem,cc,ms = 1,iter = 0;
	arma::mat Z = arma::zeros<arma::mat>(NN,6);
	arma::vec props = props0,curr_props = props0,
		PARS = arma::zeros<arma::vec>(2),curr_PARS = PARS;
	arma::uvec classify = arma::zeros<arma::uvec>(NN);
	PARS.at(0) = std::log(mm0 / (1.0 - mm0));
	if( psi0 > 0.0 ) PARS.at(1) = std::log(psi0);
	arma::mat I_np = arma::eye<arma::mat>(2,2);
	
	if( show ) Rcpp::Rcout << "Begin EM ...\n";
	while(iter < max_iter){
		// EM
		old_LL = Rcpp_LLZ(RD,DP,log_DP,LBC,log_IND,I_np,Z,props,mm,psi,show);
		PARS.at(0) = std::log(mm / (1.0 - mm));
		if( psi > 0.0 ) PARS.at(1) = std::log(psi);
		
		if( iter > 0 ){
			// Adjust
			classify = arma::index_max(Z,1);
			for(cc = 0; cc < 6; cc++){
				if( arma::sum(classify == cc) == 0 ){
					props.at(cc) = 0.0;
				}
			}
			props /= arma::sum(props);
			
			// Check convergence
			if( std::abs(old_LL - curr_LL) < conv && 
				Rcpp_norm(props - curr_props) < conv && 
				Rcpp_norm(PARS - curr_PARS) < conv ){
				break;
			}
		}
		
		if( show ){
			Rcpp::Rcout << "Iter = " << iter+1 << ": LL = " << old_LL << "; ";
			if(iter > 0) Rcpp::Rcout << "diff_LL = " << old_LL - curr_LL << ";\n";
			Rcpp::Rcout << "\tProps = " << props.t();
			Rcpp::Rcout << "\tPARS = " << PARS.t() << "\n";
		}
		
		curr_LL = old_LL;
		curr_props = props;
		curr_PARS = PARS;
		iter += 1;
	}
	
	std::string converge;
	BIC = 2*old_LL - std::log(NN)*ms;
	if(iter < max_iter){
		converge = "yes";
		old_LL = Rcpp_LLZ(RD,DP,log_DP,LBC,log_IND,I_np,Z,props,mm,psi,show);
		
		// Calculate model size
		ms = arma::sum( props > 0.0 ) - 1; 							// number of mix props > 0
		ms += 1 * ( props.at(0) > 0.0 );								// account for noise distribution
		ms += 1 * ( props.at(3) > 0.0 ); 								// some RD classify to MAF = 0.5
		ms += 1 * ( props.at(2) + props.at(4) > 0.0 ); 	// estimate mm
		ms += 1 * ( props.at(0) < 1.0 && psi > 0.0 ); 	// estimate psi, overdispersion
		
		BIC = 2*old_LL - std::log(NN)*ms;
		classify = arma::index_max(Z,1) + 1;
		return Rcpp::List::create(
			Rcpp::Named("iter",iter),Rcpp::Named("converge",converge),
			Rcpp::Named("props0",Rcpp::NumericVector(props0.begin(),props0.end())),
			Rcpp::Named("props",Rcpp::NumericVector(props.begin(),props.end())),
			Rcpp::Named("LL",old_LL),Rcpp::Named("ms",ms),Rcpp::Named("BIC",BIC),
			Rcpp::Named("class",Rcpp::NumericVector(classify.begin(),classify.end())),
			Rcpp::Named("maf",mm),Rcpp::Named("psi",psi)
			);
	} else {
		converge = "no";
		return Rcpp::List::create(
			Rcpp::Named("BIC",BIC),
			Rcpp::Named("converge",converge));
	}

}



// --------------------
// Get TARGET indicator
// --------------------

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


// --------------------
// Model VAF general mixture of zero/one inflation, symmetry, balance (ZOISB)
// --------------------

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


// R quantiles
// Reference: https://www.rdocumentation.org/packages/stats/versions/3.5.2/topics/quantile

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

arma::uvec Rcpp_quantile(const arma::vec& xx,const arma::vec& pp){
	arma::uvec uvec_qt = arma::zeros<arma::uvec>(pp.n_elem);
	arma::uword ii;
	
	for(ii = 0; ii < pp.n_elem; ii++){
			uvec_qt.at(ii) = std::floor(my_quantile(xx,pp.at(ii),7));
	}
	
	return uvec_qt;
}

// ----------
// Code to increase runtime for RDQ
// ----------

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
// BAF Clustering, modeling with binomial or beta-binomial density
// ----------

double Rcpp_BAF_binom_LL_Z(const arma::mat& RD,const arma::vec& DP,
	const arma::vec& log_DP,const arma::vec& LBC,arma::mat& Z, 
	arma::vec& pi,double& pp,arma::uvec& counts){
	
	arma::uword B = DP.n_elem;
	arma::uword bb,cc;
	for(cc = 0; cc < 4; cc++){
		if(cc == 0){
			Z.col(cc) = std::log(pi.at(cc)) - log_DP;
		} else if(cc == 1){
			Z.col(cc) = std::log(pi.at(cc)) + LBC + DP * log(0.5);
		} else if(cc == 2){
			Z.col(cc) = std::log(pi.at(cc)) + LBC + 
				RD.col(0) * std::log(pp) + RD.col(1) * std::log(1.0 - pp);
		} else {
			Z.col(cc) = std::log(pi.at(cc)) + LBC + 
				RD.col(0) * std::log(1.0 - pp) + RD.col(1) * std::log(pp);
		}
	}
	
	double LL = 0.0;
	double tmp_num,tmp_sum,tmp_pi_p;
	counts.zeros();
	for(bb = 0; bb < B; bb++){
		tmp_num = Rcpp_logSumExp(Z.row(bb).t());
		LL += tmp_num;
		Z.row(bb) = arma::exp(Z.row(bb).t() - tmp_num).t();
		counts.at(arma::index_max(Z.row(bb).t())) += 1;
	}
	
	// M-steps
	pi = arma::sum(Z,0).t();
	for(cc = 0; cc < 4; cc++){
		if( counts.at(cc) == 0 ) pi.at(cc) = 0.0;
	}
	pi /= arma::sum(pi);
	tmp_sum = arma::sum( arma::sum(Z.cols(2,3),1) % DP );
	if( tmp_sum > 0.0 ){
		pp = arma::sum( (Z.col(2) % RD.col(0) + Z.col(3) % RD.col(1)) ) / tmp_sum;
		tmp_pi_p = pi.at(2) / arma::sum(pi.subvec(2,3));
		if( tmp_pi_p < 0.5 ){
			pp = 1.0 - pp;
			pi.at(2) = pi.at(3);
			pi.at(3) = 1.0 - arma::sum(pi.subvec(0,2));
		}
	} else {
		pp = 0.5;
	}
	
	return LL;
}

Rcpp::List Rcpp_BAF_binom(const arma::mat& RD,
	const arma::vec& DP,const arma::vec& log_DP,
	const arma::vec& LBC,const arma::vec& pi0,
	const double& pp0,const arma::uword& max_iter = 4e3,
	const double& eps = 1e-5,const bool& show = true){
	
	arma::uword B = RD.n_rows;
	arma::vec pi = pi0;
	double pp = pp0;
	arma::mat Z = arma::zeros<arma::mat>(B,4);
	arma::uword iter = 0; 
	double old_LL,BIC;
	double curr_LL = 0.0;
	double curr_pp = 0.0;
	arma::vec curr_pi = arma::zeros<arma::vec>(4);
	arma::uvec counts = arma::zeros<arma::uvec>(4);
	std::string converge;
	
	// Optimize
	while(iter < max_iter){
		// EM-steps
		old_LL = Rcpp_BAF_binom_LL_Z(RD,DP,log_DP,LBC,Z,pi,pp,counts);
		
		if(show){
			Rcpp::Rcout << "LL=" << old_LL << "; pp=" << pp << "\n";
			Rcpp::Rcout << pi.t();
		}
		
		// Check convergence
		if(iter > 0){
			if( std::abs(old_LL - curr_LL) < eps ){
				if( Rcpp_norm(curr_pi - pi) < eps
					&& std::abs(curr_pp - pp) < eps ){
					break;
				}
			}
		}
		
		curr_LL = old_LL;
		curr_pi = pi;
		curr_pp = pp;
		iter++;
	}
	
	arma::uword ms = arma::sum(pi > 0.0) - 1 + 
		1*(pi.at(1) > 0) +
		1*(arma::sum(pi.subvec(2,3)) > 0.0);
	if( ms == 0 ) ms = 1;
	/*
	arma::uword ms = 1*(params.at(0) > 0.0 && params.at(0) < 1.0) +
		1*(params.at(0) < 1.0 && params.at(1) > 0.0 && params.at(1) < 1.0) +
		1*(params.at(0) < 1.0 && params.at(1) < 1.0 && params.at(2) > 0.0 && params.at(2) < 1.0) + 
		1*(params.at(0) < 1.0 && params.at(1) < 1.0 && params.at(3) > 0.0 && params.at(3) < 1.0);
	*/
	
	old_LL = Rcpp_BAF_binom_LL_Z(RD,DP,log_DP,LBC,Z,pi,pp,counts);
	BIC = 2*old_LL - log(B)*ms;
	arma::uvec classify = arma::index_max(Z,1) + 1;
	if(iter < max_iter){
		converge = "yes";
	} else {
		converge = "no";
	}
	arma::vec params = arma::zeros<arma::vec>(4);
	params.at(0) = pi.at(0);
	if( pi.at(0) < 1.0 ) params.at(1) = pi.at(1) / (1.0 - pi.at(0));
	params.at(2) = pp;
	if( arma::sum(pi.subvec(2,3)) > 0.0 ){
		params.at(3) = pi.at(2) / arma::sum(pi.subvec(2,3));
	} else {
		params.at(3) = 0.5;
	}
	
	return Rcpp::List::create(
		Rcpp::Named("iter",iter),
		Rcpp::Named("converge",converge),
		Rcpp::Named("LL",old_LL),
		Rcpp::Named("ms",ms),
		Rcpp::Named("BIC",BIC),
		Rcpp::Named("Z",Z),
		Rcpp::Named("class",Rcpp::NumericVector(classify.begin(),classify.end())),
		Rcpp::Named("pi",Rcpp::NumericVector(pi.begin(),pi.end())),
		Rcpp::Named("pp",pp),
		Rcpp::Named("params",Rcpp::NumericVector(params.begin(),params.end())));
}

// ----------
// NEW BAF Functions
// ----------

double Rcpp_BAF_ecLL(const arma::mat& RD,const arma::vec& DP,
	const arma::vec& LBC,const arma::mat& Z,const arma::vec& props,
	const arma::vec& PARS){
	
	double psi = std::exp(PARS.at(1)),mm = std::exp(PARS.at(0));
	mm /= 1.0 + mm;
	double LL = 0.0;
	
	if( props.at(1) > 0.0 ) LL += arma::sum(Z.col(1) % Rcpp_LL(RD,DP,LBC,0.5,psi));
	if( props.at(2) > 0.0 ) LL += arma::sum(Z.col(2) % Rcpp_LL(RD,DP,LBC,mm,psi));
	if( props.at(3) > 0.0 ) LL += arma::sum(Z.col(3) % Rcpp_LL(RD,DP,LBC,1.0 - mm,psi));
	
	return LL;
}

arma::vec Rcpp_BAF_ecGRAD(const arma::mat& RD,const arma::vec& DP,
	const arma::mat& Z,const arma::vec& props,const arma::vec& PARS){
	
	arma::uword ii, cc, NN = RD.n_rows;
	double psi = std::exp(PARS.at(1)),mm = std::exp(PARS.at(0));
	mm /= 1.0 + mm;
	double ipsi = 1.0 / psi, part_ipsi = -1.0 * ipsi * ipsi,
		aa = mm * ipsi,bb = (1.0 - mm) * ipsi,
		dd = 0.5 * ipsi,part_aa_mm = ipsi,part_bb_mm = -1.0 * part_aa_mm,
		dg_aa = R::digamma(aa),dg_bb = R::digamma(bb),dg_dd = R::digamma(dd),
		dg_ipsi = R::digamma(ipsi),
		part_aa_psi = part_ipsi * mm,part_bb_psi = part_ipsi * (1.0 - mm),
		part_dd_psi = 0.5 * part_ipsi,
		dg_xx_aa,dg_nx_bb,dg_xx_bb,dg_nx_aa,dg_xx_dd,dg_nx_dd,dg_nn_ipsi;
	arma::vec GRAD = arma::zeros<arma::vec>(2); // mm,psi
	
	for(ii = 0; ii < NN; ii++){
		dg_xx_aa = R::digamma(RD.at(ii,0) + aa);
		dg_nx_bb = R::digamma(RD.at(ii,1) + bb);
		dg_xx_bb = R::digamma(RD.at(ii,0) + bb);
		dg_nx_aa = R::digamma(RD.at(ii,1) + aa);
		dg_nn_ipsi = R::digamma(DP.at(ii) + ipsi);
	for(cc = 1; cc < 4; cc++){
		// part unc_mm
		if( props.at(cc) > 0.0 && cc == 2 ){
			GRAD.at(0) += Z.at(ii,cc) * ( dg_xx_aa * part_aa_mm + 
				dg_nx_bb * part_bb_mm - dg_aa * part_aa_mm - dg_bb * part_bb_mm );
		}
		if( props.at(cc) > 0.0 && cc == 3 ){
			GRAD.at(0) += Z.at(ii,cc) * ( dg_xx_bb * part_bb_mm +
				dg_nx_aa * part_aa_mm - dg_aa * part_aa_mm - dg_bb * part_bb_mm );
		}
		
		// part log(psi)
		if( props.at(cc) > 0.0 && cc == 1){
			dg_xx_dd = R::digamma(RD.at(ii,0) + dd);
			dg_nx_dd = R::digamma(RD.at(ii,1) + dd);
			GRAD.at(1) += Z.at(ii,cc) * ( dg_xx_dd * part_dd_psi +
				dg_nx_dd * part_dd_psi - dg_nn_ipsi * part_ipsi +
				dg_ipsi * part_ipsi - 2.0 * dg_dd * part_dd_psi );
		} else if( props.at(cc) > 0.0 && cc == 2){
			GRAD.at(1) += Z.at(ii,cc) * ( dg_xx_aa * part_aa_psi +
				dg_nx_bb * part_bb_psi - dg_nn_ipsi * part_ipsi +
				dg_ipsi * part_ipsi - dg_aa * part_aa_psi -
				dg_bb * part_bb_psi );
		} else if( props.at(cc) > 0.0 && cc == 3){
			GRAD.at(1) += Z.at(ii,cc) * ( dg_xx_bb * part_bb_psi +
				dg_nx_aa * part_aa_psi - dg_nn_ipsi * part_ipsi +
				dg_ipsi * part_ipsi - dg_aa * part_aa_psi -
				dg_bb * part_bb_psi );
		}
		
	}}
	
	GRAD.at(0) *= mm * (1.0 - mm);
	GRAD.at(1) *= psi;
	return GRAD;
}

void Rcpp_BAF_ecBFGS(const arma::mat& RD,const arma::vec& DP,
	const arma::vec& LBC,const arma::mat& Z,const arma::vec& props,
	arma::vec& PARS,arma::uword& converge,const arma::mat& I_np,
	const arma::uword& max_iter = 4e3,const double& eps = 1e-10,
	const double& gr_eps = 1e-2,const bool& show = true){
	
	converge = 0;
	arma::uword iter = 0, jj, uu, reset_Bk = 0,np = 2;
	arma::mat inv_Bk = I_np, ISYT = inv_Bk;
	arma::vec xk = PARS, curr_xk = arma::zeros<arma::vec>(np),
		new_xk = curr_xk,gr_k = curr_xk,p_k = curr_xk,
		s_k = curr_xk,y_k = curr_xk;
	double old_LL,new_LL,inv_norm_p_k,tmp_step,ys,fnscale = -1.0,
		curr_LL = 0.0,mm,psi,iBG;
	mm = std::exp(PARS.at(0));
	mm /= 1.0 + mm;
	double orig_LL = Rcpp_BAF_ecLL(RD,DP,LBC,Z,props,xk);
	
	while(iter < max_iter){
		// Calculate Direction p_k
		gr_k = fnscale * Rcpp_BAF_ecGRAD(RD,DP,Z,props,xk);
		
		p_k = -1.0 * inv_Bk * gr_k;
		inv_norm_p_k = 1.0 / std::max(1.0,Rcpp_norm(p_k));

		// Line search for new xk
		uu = 0;
		old_LL = fnscale * Rcpp_BAF_ecLL(RD,DP,LBC,Z,props,xk);
		
		// Less than 0 b/c fnscale = -1
		if( false && old_LL < 0.0 ){
			break; 
		}
		
		for(jj = 0; jj <= 30; jj++){
			tmp_step 	= inv_norm_p_k / std::pow(2,jj);
			new_xk 		= xk + tmp_step * p_k;
			new_LL = fnscale * Rcpp_BAF_ecLL(RD,DP,LBC,Z,props,new_xk);
			if( new_LL < old_LL ){ // minimizing
				s_k = tmp_step * p_k;
				y_k = fnscale * Rcpp_BAF_ecGRAD(RD,DP,Z,props,new_xk) - gr_k;
				if( y_k.has_nan() ) continue;
				ys = arma::dot(y_k,s_k);
				if( ys > 0.0 ){
					ISYT = I_np - (s_k * y_k.t()) / ys;
					inv_Bk = ISYT * inv_Bk * ISYT.t() + s_k * s_k.t() / ys;
				}
				
				xk = new_xk;
				old_LL = new_LL;
				uu = 1;
				break;
			}
		}
		
		if( uu == 0 ) { // aka no update
			if( Rcpp_norm(gr_k) > 1.0 ){
				if( show ) printR_obj("\tReset inv_Bk");
				inv_Bk = I_np;
				reset_Bk++;
			} else {
				if( show ) printR_obj("\tFailed line search");
				break;
			}
		}
		
		if( reset_Bk > 5 ) break;
		
		// Check Convergence
		if( iter > 0 ){
			if( std::abs(curr_LL - old_LL) < eps &&
				Rcpp_norm(curr_xk - xk) < eps ){
				gr_k = Rcpp_BAF_ecGRAD(RD,DP,Z,props,xk);
				if( Rcpp_norm(gr_k) < eps ){
					break;
				}
			}
		}
		
		curr_xk = xk;
		curr_LL = old_LL;
		iter++;
	
	}
	
	// Update parameters
	PARS = xk;
	mm = std::exp(PARS.at(0));
	mm /= 1.0 + mm;
	psi = std::exp(PARS.at(1));
	
	// Calculate LL, GRAD
	old_LL = Rcpp_BAF_ecLL(RD,DP,LBC,Z,props,PARS);
	gr_k = Rcpp_BAF_ecGRAD(RD,DP,Z,props,PARS);
	p_k = -1.0 * inv_Bk * gr_k;
	iBG = Rcpp_norm(p_k);
	if( iBG < gr_eps ){
		converge = 1;
	}
	
	if( show ){
		Rcpp::Rcout << "\tLL_old = " << orig_LL << "; Iter = " << iter 
			<< "; LL_new = " << old_LL << "\n";
		double nG = Rcpp_norm(gr_k);
		Rcpp::Rcout << "\tConvergence Indicators: \n"
			<< "\t   NormGrad = " << nG << ";\n"
			<< "\t   NormiBkGr = " << iBG << ";\n";
		Rcpp::Rcout << "\t   MAF = " << mm << "\n";
		Rcpp::Rcout << "\t   PSI = " << psi << "\n";
		Rcpp::Rcout << "\t   Convergence Status = ";
		if( converge == 1 ){
			Rcpp::Rcout << "YES\n";
		} else {
			Rcpp::Rcout << "NO\n";
		}
		Rcpp::Rcout << "\n";
	}
	
}

double Rcpp_BAF_LL_Z(const arma::mat& RD,const arma::vec& DP,
	const arma::vec& log_DP,const arma::vec& LBC,const arma::uword& BB,
	const arma::mat& I_np,arma::mat& Z,arma::vec& props,double& mm,
	double& psi,const bool& show = true){
	
	arma::uword ii;
	
	// E-step
	Z.fill(-1.0 * arma::datum::inf);
	if( props.at(0) > 0.0 ) Z.col(0) = std::log(props.at(0)) - log_DP;
	if( props.at(1) > 0.0 ) Z.col(1) = std::log(props.at(1)) + Rcpp_LL(RD,DP,LBC,0.5,psi);
	if( props.at(2) > 0.0 ) Z.col(2) = std::log(props.at(2)) + Rcpp_LL(RD,DP,LBC,mm,psi);
	if( props.at(3) > 0.0 ) Z.col(3) = std::log(props.at(3)) + Rcpp_LL(RD,DP,LBC,1.0 - mm,psi);
	
	// Calculate observed LL
	double LL = 0.0, tmp_num;
	for(ii = 0; ii < BB; ii++){
		tmp_num = Rcpp_logSumExp(Z.row(ii).t());
		LL += tmp_num;
		Z.row(ii) = arma::exp(Z.row(ii).t() - tmp_num).t();
	}
	
	// M-steps
	
	// update proportions
	props = arma::mean(Z,0).t();
	
	// If psi = 0, binomial mixtures, use code below
	if( psi == 0.0 ){
		// update mm
		double pp_den = arma::sum( (Z.col(2) + Z.col(3)) % DP );
		if( pp_den > 0 ){
			mm = arma::sum(Z.col(2) % RD.col(0) + Z.col(3) % RD.col(1)) / pp_den;
		} else {
			mm = 0.5;
		}
	} else {
		arma::vec PARS = arma::zeros<arma::vec>(2);
		if( arma::sum( props.subvec(2,3) ) == 0.0 ){
			mm = 0.5;
		}
		PARS.at(0) = std::log(mm / (1.0 - mm));
		PARS.at(1) = std::log(psi);
		arma::uword converge = 0;
		Rcpp_BAF_ecBFGS(RD,DP,LBC,Z,props,PARS,
			converge,I_np,4e3,1e-5,1e-2,show);
		mm = std::exp(PARS.at(0));
		mm /= 1.0 + mm;
		psi = std::exp(PARS.at(1));
	}
	
	return LL;
}

// [[Rcpp::export]]
Rcpp::List Rcpp_BAF_clust(const arma::mat& RD,const arma::vec& DP,
	const arma::vec& log_DP,const arma::vec& LBC,const arma::vec& props0,
	const double& mm0,const double& psi0,const arma::uword& max_iter = 4e3,
	const double& eps = 1e-5,const bool& show = true){
	
	arma::uword BB = RD.n_rows,iter = 0,cc,ms = 1;
	double mm = mm0,old_LL,BIC,curr_LL = 0,psi = psi0;
	arma::vec props = props0,PARS = arma::zeros<arma::vec>(2),
		curr_PARS = PARS,curr_props = props;
	arma::mat Z = arma::zeros<arma::mat>(BB,4),
		I_np = arma::eye<arma::mat>(2,2);
	arma::uvec classify = arma::zeros<arma::uvec>(BB);
	
	// Optimize
	while(iter < max_iter){
		// EM-steps
		old_LL = Rcpp_BAF_LL_Z(RD,DP,log_DP,LBC,BB,I_np,Z,props,mm,psi,show);
		PARS.at(0) = std::log(mm / (1.0 - mm));
		if( psi > 0.0 ) PARS.at(1) = std::log(psi);
		if( false && psi < 1e-5 ){
			psi = 0.0;
			PARS.at(1) = -1.0 * arma::datum::inf;
		}
		
		if(iter > 0){
			// Adjust
			classify = arma::index_max(Z,1);
			for(cc = 0; cc < 6; cc++){
				if( arma::sum(classify == cc) == 0 ){
					props.at(cc) = 0.0;
				}
			}
			props /= arma::sum(props);
			if( arma::sum(props.subvec(2,3)) == 0.0 ){
				PARS.at(0) = 0.0;
				mm = 0.5;
			}
			
			// Check convergence
			if( std::abs(old_LL - curr_LL) < eps ){
				if( Rcpp_norm(curr_props - props) < eps
					&& Rcpp_norm(curr_PARS - PARS) < eps ){
					break;
				}
			}
		}
		
		if( show ){
			Rcpp::Rcout << "Iter = " << iter+1 << ": LL = " << old_LL << "; ";
			if(iter > 0) Rcpp::Rcout << "diff_LL = " << old_LL - curr_LL << ";";
			Rcpp::Rcout << "\n";
			Rcpp::Rcout << "\tProps = " << props.t();
			Rcpp::Rcout << "\tPARS = " << PARS.t() << "\n";
		}
		
		curr_LL = old_LL;
		curr_props = props;
		curr_PARS = PARS;
		iter++;
	}
	
	std::string converge;
	BIC = 2 * old_LL - std::log(BB) * ms;
	if(iter < max_iter){
		converge = "yes";
		old_LL = Rcpp_BAF_LL_Z(RD,DP,log_DP,LBC,BB,I_np,Z,props,mm,psi,show);
		if( props.at(0) == 1.0 ) psi = 0.0;
		
		// Calculate model size
		ms = arma::sum( props > 0.0 ) - 1; 							// number of mix props > 0
		ms += 1 * ( props.at(0) > 0.0 );								// account for noise distribution
		ms += 1 * ( props.at(1) > 0.0 ); 								// some RD classify to MAF = 0.5
		ms += 1 * ( props.at(2) + props.at(3) > 0.0 ); 	// estimate mm
		ms += 1 * ( props.at(0) < 1.0 && psi > 0.0 ); 	// estimate psi, overdispersion
		
		BIC = 2 * old_LL - std::log(BB) * ms;
		classify = arma::index_max(Z,1) + 1;
		
		// Defining vector of output parameters: params = (prop_noise,prop_half,MAF,prop_MAF,psi)
		arma::vec params = arma::zeros<arma::vec>(5);
		params.at(0) = props.at(0);
		if( props.at(0) < 1.0 ) params.at(1) = props.at(1) / (1.0 - props.at(0));
		params.at(2) = mm;
		if( arma::sum(props.subvec(2,3)) > 0.0 ){
			params.at(3) = props.at(2) / arma::sum(props.subvec(2,3));
		} else {
			params.at(3) = 0.5;
		}
		params.at(4) = psi;
		
		return Rcpp::List::create(Rcpp::Named("iter",iter),Rcpp::Named("converge",converge),
			Rcpp::Named("props",Rcpp::NumericVector(props.begin(),props.end())),
			Rcpp::Named("LL",old_LL),Rcpp::Named("ms",ms),Rcpp::Named("BIC",BIC),
			Rcpp::Named("class",Rcpp::NumericVector(classify.begin(),classify.end())),
			Rcpp::Named("maf",mm),Rcpp::Named("psi",psi),
			Rcpp::Named("params",Rcpp::NumericVector(params.begin(),params.end())));
	} else {
		converge = "no";
		return Rcpp::List::create(Rcpp::Named("BIC",BIC),Rcpp::Named("converge",converge));
	}
	
}




