#include "../inst/include/allele_freq.h"

// [[Rcpp::depends("RcppArmadillo")]]

// ----------
// BAF Clustering, modeling with binomial or beta-binomial density
// ----------

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
				if( show ) Rcpp::Rcout << "\tReset inv_Bk\n";
				inv_Bk = I_np;
				reset_Bk++;
			} else {
				if( show ) Rcpp::Rcout << "\tFailed line search\n";
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

