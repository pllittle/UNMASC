# Smart Functions
smart_uniq_df = function(input_df,vars){
	input_df[!duplicated(input_df[,vars]),]
}
smart_sprintf = function(...){
	orig = sprintf(...)
	orig = gsub("\n","",orig)
	orig = gsub("\t","",orig)
	orig
}
smart_header = function(OBJ,req_COLS){
	for(COL in req_COLS){
		if( !(COL %in% colnames(OBJ)) )
			stop(sprintf("%s header missing",COL))
	}
}
smart_RefAlt = function(vcs){
	smart_header(OBJ = vcs,req_COLS = c("Ref","Alt"))
	
	vcs$INDEL = ifelse(nchar(vcs$Ref) == 1 & nchar(vcs$Alt) == 1,FALSE,TRUE)
	idx_bs = which(vcs$INDEL == FALSE)
	vcs$RefAlt = "INDEL"
	vcs$RefAlt[idx_bs] = paste0(vcs$Ref[idx_bs],">",vcs$Alt[idx_bs])
	rm(idx_bs)
	# smart_table(vcs$RefAlt)
	vcs
}

# Picard hsMetric normal sample clustering
gauss_unif = function(DATA,VARS,ID,show_plot=FALSE){
	if(FALSE){
		DATA = hsMetrics
		# VARS = interest_vars[c(1)]
		# VARS = interest_vars[c(1,2,3)]
		VARS = interest_vars
		show_plot = FALSE
	}
	# DATA: columns should include the continuous variable and studynumber

	# Identify 1 dimensional(marginal) outliers (1 Gaussian, 1 noise cluster)
	calc_obsLL_Z = function(XX,range_XX,mu,sigma2,epsilon){
		if(FALSE){
			mu = init_mu
			sigma2 = init_sigma2
			epsilon = init_eps
		}
		
		N = nrow(XX)
		D = ncol(XX)
		Z_mat = matrix(NA,N,2) # either gaussian or uniform noise
		run_E_step = TRUE

		if(D == 1){
			inv_sigma2 = 1/sigma2
			log_abs_det_sigma2 = log(sigma2)
		} else {
			tmp_rcond = rcond(sigma2)
			if( is.na(tmp_rcond) || tmp_rcond == 0){
				run_E_step = FALSE
			} else {
				inv_sigma2 = solve(sigma2,tol=0.1*rcond(sigma2))
				log_abs_det_sigma2 = log(abs(det(sigma2)))
			}
		}
		
		if( run_E_step ){
		
			for(dd in seq(2)){
				if(dd == 1){
					Z_mat[,dd] = log(epsilon) + 
						sum(apply(range_XX,2,function(x) log(1/(x[2]-x[1]))))
				} else {
					Z_mat[,dd] = log(1-epsilon) + (-D/2) * 
						log(2*pi) - 0.5 * log_abs_det_sigma2 - 
						0.5 * apply(XX,1,function(x) t(x - mu) 
						%*% inv_sigma2 %*% (x-mu))
				}
			}
			
			LL = 0
			for(bb in seq(N)){
				tmp_num = Rcpp_logSumExp(Z_mat[bb,])
				Z_mat[bb,] = exp(Z_mat[bb,] - tmp_num)
				LL = LL + tmp_num
			}
			
			list(LL=LL,Z=Z_mat)
			
		} else {
			list(LL=NA,Z=NA)
		}
	}

	# Get Data
	XX = as.matrix(DATA[,VARS])
	N = nrow(XX)
	D = ncol(XX)
	range_XX = apply(XX,2,range)

	# Init params
	init_mu = as.numeric(apply(range_XX,2,function(x) runif(1,x[1],x[2])))
	if(D == 1){
		init_sigma2 = var(XX)
	} else {
		init_sigma2 = diag(as.numeric(apply(XX,2,var)))
		#	keep_init_sigma2 = init_sigma2
	}
	init_eps = runif(1,0,0.1)
	keep_init_eps = init_eps

	iter = 0; max_iter = 1e3; conv_eps = 1e-6; curr_LL = NA
	while(TRUE){
		# E steps
		old_stuff = calc_obsLL_Z(XX,range_XX,init_mu,init_sigma2,init_eps)
		old_LL = old_stuff$LL; # print(old_LL)
		old_Z = old_stuff$Z
		
		if( is.na(old_LL) ){
			init_eps = 0.9
			break
		}
		
		if(show_plot && iter %% 5 == 0){
			out_class = apply(old_Z,1,which.max)
			plot(DATA[,VARS],col=ifelse(out_class == 1,"red","blue"))
		}
		
		# M steps
		init_eps = mean(old_Z[,1])
		# init_mu = sum(old_Z[,2] * XX) / sum(old_Z[,2])
		tmp_mu = rep(0,D)
		for(ii in seq(N)){
			tmp_mu = tmp_mu + old_Z[ii,2] * XX[ii,]
		}
		init_mu = tmp_mu / sum(old_Z[,2])
		# init_sigma2 = sum(old_Z[,2] * (XX - init_mu)^2) / sum(old_Z[,2])
		tmp_sigma2 = matrix(0,D,D)
		for(ii in seq(N)){
			tmp_sigma2 = tmp_sigma2 + old_Z[ii,2] * (XX[ii,] - init_mu) %*% t((XX[ii,] - init_mu))
		}
		init_sigma2 = tmp_sigma2 / sum(old_Z[,2])

		if(iter > 0){ # Check convergence
			if( old_LL >= curr_LL && abs(old_LL - curr_LL) < conv_eps ){
				# print("Converged")
				break
			}
		}

		curr_LL = old_LL
		iter = iter + 1
	}

	# Recursive part if optimization runs into a bad local optimum
	if( init_eps > 0.5 ){
		# print("Re-initializing ...")
		gauss_unif(DATA = DATA,VARS = VARS,ID = ID,show_plot = show_plot)
	} else {
		# print("Converged")
		rownames(init_sigma2) = colnames(init_sigma2)
		output = list(mean=init_mu,covariance=init_sigma2,eps=init_eps)
		output = sapply(output,function(x) round(x,7))
		out_class = apply(old_Z,1,which.max)
		if( show_plot ) plot(DATA[,VARS],col=ifelse(out_class == 1,"red","blue"))

		# output outlier samples
		output$outliers = DATA[[ID]][out_class == 1]
		output$Z = old_Z
		output$LL = old_LL
		output
	}
	
}
get_max_LL_result = function(DATA,VARS,ID,TRIALS){
	if(FALSE){
		DATA = hsMetrics; VARS = interest_vars; TRIALS = 10
	}
	
	all_result = list()
	for(ii in seq(TRIALS)){
		if( ii %% 2 == 0 ) message(".",appendLF = FALSE)
		if( ii %% 20 == 0 || ii == TRIALS ) message(sprintf("%s out of %s\n",ii,TRIALS),appendLF = FALSE)
		all_result[[ii]] = gauss_unif(DATA = DATA,
			VARS = VARS,ID = ID)
	}
	
	vec_LL = sapply(all_result,function(x) x$LL)
	max_index = which.max(vec_LL)
	# print(vec_LL)
	all_result[[max_index]]
}

#' @title run_hsMetric_clustRank
#' @description Clusters sample metrics to identify outlier and 
#'	reliable samples to establish a normal control cohort
#' @param DATA A data.frame of sample ID and metrics to cluster on
#' @param VARS A character vector of at least two metrics to cluster on
#'	within \code{DATA}
#' @param ID A character string specifying the column name in \code{DATA}
#'	that corresponds to sample ID.
#' @param TRIALS An integer number of times to re-run the EM algorithm,
#'	more trials help to avoid using a sub-optimal solution.
#' @param note A character string to label the plot title
#' @return R data.frame of cluster results and sample labeling.
#'
#' @export
run_hsMetric_clustRank = function(DATA,VARS,ID,TRIALS=50,note=""){
	if( !is(object = DATA,class2 = "data.frame") ) stop("DATA should be a data.frame")
	if( length(VARS) < 2 ) stop("Add more variables")
	if( any(!(VARS %in% colnames(DATA))) ) stop("Missing columns")
	
	# Plot
	# plot(DATA[,VARS])
	clust_out = get_max_LL_result(DATA = DATA,
		VARS = VARS,ID = ID,TRIALS = TRIALS)
	my_cols = hue_pal()(2)
	DATA$samp_col = ifelse(DATA[[ID]] %in% clust_out$outliers,
		my_cols[1],my_cols[2])
	print(clust_out[c("mean","covariance","eps","outliers")]) 
	par(cex.main=1.2)
	my_main = "Samples colored by cluster class"
	if( note != "" ) my_main = paste0(my_main,": ",note)
	plot(DATA[,VARS],main=my_main,col = DATA$samp_col)
	par(cex.main=1.2)

	# Label/rank samples
	DATA$class = "good_sample"
	DATA$class[which(DATA[[ID]] %in% clust_out$outliers)] = "outlier"
	prob_df = smart_df(clust_out$Z)
	names(prob_df) = paste0("prob_",c("outlier","Gaussian"))
	DATA = cbind(DATA,prob_df)
	DATA = DATA[order(DATA$prob_outlier,-DATA$prob_Gaussian),]
	rownames(DATA) = NULL
	DATA$rank = seq(nrow(DATA))
	out = DATA[,c(ID,"samp_col","class","prob_outlier","prob_Gaussian","rank")]
	
	# Preview
	message("Top ranked samples ...\n",appendLF = FALSE)
	print(head(out))
	message("Low ranked samples ...\n",appendLF = FALSE)
	print(tail(out))
	
	# Output
	out
}

# Artifact Functions
artifact_filter = function(xx){
	if(FALSE){
		xx = fPN_segANNO$FFPE[1]
		xx = "noise(11,64.71%);signal_like(6,35.29%)"
		xx = "noise(6,35.29%);signal_like(11,64.71%)"
	}
	
	# Search for signal, if nothing return TRUE, else get percent
	src_xx = grepl("signal",xx)
	if(src_xx){
		xx2 = strsplit(xx,";")[[1]]
		xx2 = xx2[grepl("signal",xx2)]
		xx2 = as.numeric(strsplit(xx2,"[,%]")[[1]][2]) 
		ifelse(xx2 >= 50,FALSE,TRUE)
	} else {
		TRUE
	}
}

# OXOG/FFPE Functions
BINOM_NOISE = function(mat_RD,binom = TRUE){
	if(FALSE){
		mat_RD = as.matrix(ALL_VCS[tmp_index,c("tAD","tRD")])
		mat_RD = as.matrix(nn2[,c("tAD","tRD")])
	}
	
	DATA = smart_df(mat_RD)
	names(DATA) = c("tAD","tRD")
	DATA$tDP = rowSums(DATA[,c("tAD","tRD")])
	DATA$log_DP = log(DATA$tDP)
	DATA$LBC = lchoose(DATA$tDP,DATA$tAD)
	grid_maf0 = c(seq(2,8,2)/100,seq(1,3)/10)
	iPSI_vec = c(0,0.005)
	if( binom == TRUE ) iPSI_vec = 0
	
	final_clust = NULL
	
	# noise
	for(iPSI in iPSI_vec){
		tmp_clust = Rcpp_BAF_clust(RD = as.matrix(DATA[,c("tAD","tRD")]),
			DP = DATA$tDP,log_DP = DATA$log_DP,LBC = DATA$LBC,props0 = c(1,0,0,0),
			mm0 = 0.04,psi0 = iPSI,max_iter = 4e3,eps = 1e-5,show = FALSE)
		if( tmp_clust$converge == "yes" ){
			if( is.null(final_clust) || tmp_clust$BIC > final_clust$BIC ){
				final_clust = tmp_clust
			}
		}
		rm(tmp_clust)
	}
	# final_clust[c("iter","converge","BIC","props","maf","psi")]
	
	props0_mat = matrix(c(
		# 1,1,0,0,
		1,0,1,1,
		# 0,1,1,1,
		1,1,1,1),ncol=4,byrow=TRUE)
	props0_mat = props0_mat / rowSums(props0_mat)
	# props0_mat
	
	for(iPSI in iPSI_vec){
	for(maf0 in grid_maf0){
	for(row_props0 in seq(nrow(props0_mat))){
		# iPSI = iPSI_vec[2]; iPSI
		# maf0 = grid_maf0[5]; maf0
		props0 = props0_mat[row_props0,]
		tmp_clust = Rcpp_BAF_clust(RD = as.matrix(DATA[,c("tAD","tRD")]),
			DP = DATA$tDP,log_DP = DATA$log_DP,LBC = DATA$LBC,props0 = props0,
			mm0 = maf0,psi0 = iPSI,max_iter = 1e3,eps = 1e-4,show = FALSE)
		# tmp_clust[c("iter","converge","BIC","props","maf","psi")]
		if( tmp_clust$converge == "yes" ){
			if( is.null(final_clust) || tmp_clust$BIC > final_clust$BIC ){
				final_clust = tmp_clust
			}
		}
		rm(tmp_clust)
	}}}
	# final_clust[c("iter","converge","BIC","props","maf","psi")]
	
	if( binom == FALSE && sum(final_clust$props[1:2]) == 0 ){
		return(BINOM_NOISE(mat_RD = mat_RD,binom = TRUE))
	}
	
	infer = as.character(final_clust$class)
	infer[infer %in% c("1","2")] = "noise"
	infer[infer %in% c("3","4")] = "signal"
	# smart_table(infer)
	list(CLUST = final_clust,INFER = infer)
}
BINOM_INFER = function(vec_RD,params){

	if(params[1] > 1) params[1] = 1
	if(params[2] > 1) params[2] = 1
	vec_alpha = c(params[1],
		(1-params[1])*params[2],
		(1-params[1])*(1-params[2])*params[4],
		(1-params[1])*(1-params[2])*(1-params[4]))
	VAF = vec_RD[1] / sum(vec_RD)
	# Infer segment status
	stat = get_ALLELE_STAT(params)
	# Infer VC status
	DP = sum(vec_RD)
	LBC = lchoose(DP,vec_RD[1])
	log_vec = rep(NA,4)
	log_vec[1] = log(vec_alpha[1]) - log(DP)
	log_vec[2] = log(vec_alpha[2]) + Rcpp_LL(RD = t(vec_RD),DP = DP,LBC = LBC,pi = 0.5,psi = params[5])
	log_vec[3] = log(vec_alpha[3]) + Rcpp_LL(RD = t(vec_RD),DP = DP,LBC = LBC,pi = params[3],psi = params[5])
	log_vec[4] = log(vec_alpha[4]) + Rcpp_LL(RD = t(vec_RD),DP = DP,LBC = LBC,pi = 1-params[3],psi = params[5])
	log_vec = as.numeric(exp(log_vec - Rcpp_logSumExp(log_vec)))
	log_vec
}
ARTI_find = function(vcs,strand,nSEG,qscore_thres,ad_thres,rd_thres,
	binom = TRUE,max_tVAF = 0.45,ID = "STUDYNUMBER",count_thres = 20,
	ARTI,eps_thres=0.5,psi_thres=0.01){
	
	if(FALSE){
		vcs = vcs # [which(vcs$STUDYNUMBER == SN),]
		strand = rds$strand
		nSEG = nSEG
		
		vcs = vcf; nSEG = nSEG
		
		binom = TRUE
		max_tVAF = 0.45; ID = "STUDYNUMBER"; count_thres = 20
		ARTI = c("OXOG","FFPE")[1]
		eps_thres = 0.5; psi_thres = 0.01
		
	}
	
	req_COLS = c("RefAlt","tAD","tDP","tVAF",
		"mutID","Chr","Position","ThsdG_AF","EXAC_AF")
	smart_header(OBJ = vcs,req_COLS = req_COLS)
	
	if( ARTI == "OXOG" ){
		message(sprintf("%s: Finding oxoG artifacts ...\n",date()),appendLF = FALSE)
	} else if( ARTI == "FFPE" ){
		message(sprintf("%s: Finding FFPE artifacts ...\n",date()),appendLF = FALSE)
	} else if( ARTI == "ARTI" ){
		message(sprintf("%s: Finding ARTI artifacts ...\n",date()),appendLF = FALSE)
	}
	vcs = vcs[,names(vcs) != ARTI]
	vcs$INFER = NA
	
	# Append strand bias column
	message(sprintf("%s: Merge strand info ...\n",date()),appendLF = FALSE)
	vcs = smart_rmcols(vcs,"Strand_bias_P")
	vcs = smart_merge(vcs,strand[,c("mutID","Chr",
		"Position","Ref","Alt","Strand_bias_P")],all.x=TRUE)
	
	# Append H2M info
	vcs = anno_nSEG(DAT = vcs,nSEG = nSEG,
		eps_thres = eps_thres,psi_thres = psi_thres)
	
	# Filters
	if( ARTI == "OXOG" ){
		art_filter = vcs$RefAlt %in% c("C>A","G>T")
	} else if( ARTI == "FFPE" ){
		art_filter = vcs$RefAlt %in% c("C>T","G>A","INDEL")
	} else if( ARTI == "ARTI" ){
		art_filter = vcs$RefAlt %in% c("A>C","T>G")
	}
	dp_filter 			= vcs$tDP >= rd_thres & vcs$tAD >= ad_thres # ad_thres eliminates artifact cluster near 0
	vaf_filter 			= vcs$tVAF > 0 & vcs$tVAF < min(c(1,max_tVAF))
	germ_filter 		= vcs$ThsdG_AF == 0 & vcs$EXAC_AF == 0
	sb_filter				= is.na(vcs$Strand_bias_P) | vcs$Strand_bias_P > 0.05
	qscore_filter		= vcs$Qscore >= qscore_thres
	# qscore_filter		= ( vcs$Qscore >= qscore_thres & vcs$Qscore < 1000 )
	# qscore_filter		= vcs$Qscore >= 1e3 # seems to work best
	# h2m_filter			= is.na(vcs$nANNO) | vcs$nANNO == "mappable"
	h2m_filter 				= rep(TRUE,nrow(vcs))
	
	advst_filter		= ( dp_filter & vaf_filter
										& sb_filter & h2m_filter)
	advsta_filter		= ( advst_filter & art_filter )
	glob_filter			= ( advsta_filter & germ_filter
										& qscore_filter )
	
	IDs = sort(unique(vcs[[ID]]))
	clust = NULL
	cnt = 1
	for(id in IDs){
		# id = IDs[1]
		message(paste0(id," "),appendLF = FALSE)
		if( cnt %% 3 == 0 ) message("\n",appendLF = FALSE)
		sn_filter = vcs[[ID]] == id
		idx = which(sn_filter & glob_filter)
		length(idx)
			
		if(FALSE){ # Check distribution
			idx = which(sn_filter 		# control-specific
				& dp_filter							# depth filter
				& vaf_filter						# vaf > 0 & vaf < 1
				& sb_filter							# strand bias filter
				& h2m_filter
				& art_filter						# artifact RefAlt filter
				& germ_filter						# germ filter
				& qscore_filter					# qscore filter
				)
			# idx = which(glob_filter)
			length(idx)
			smart_hist(vcs$tVAF[idx],breaks = 50,xlim = c(0,1))
			plot(vcs[idx,c("x_coor","tVAF","tAD","tDP","Qscore")])
			# plot(vcs[idx,c("x_coor","tVAF")],col=ifelse(vcs$Qscore[idx] > 1000,"red","black"))
			
			mat_RD = as.matrix(vcs[idx,c("tAD","tRD")])
			b_out = BINOM_NOISE(mat_RD = mat_RD,binom = !binom)
			# names(b_out)
			b_out$CLUST[c("iter","converge","BIC","props","maf","psi")]
			smart_table(b_out$CLUST$class)
			plot(vcs[idx,c("x_coor","tVAF","Qscore")[1:2]],ylim = c(0,1),
				# col = factor(b_out$CLUST$class),
				col = factor(b_out$INFER),pch = 16,
				cex = vcs$DP_cex[idx])
			
			tmp_vcs = vcs[idx,]
			plot(tmp_vcs[,c("tVAF","tDP","Qscore")])
			idx2 = which(tmp_vcs$Qscore > 1e3)
			tmp_vcs[idx2,]
			smart_table(tmp_vcs$EXONIC,tmp_vcs$tVAF > 0.4)
			
		}
		
		if( length(idx) >= count_thres ){
			mat_RD = as.matrix(vcs[idx,c("tAD","tRD")])
			bin_noise_out = BINOM_NOISE(mat_RD = mat_RD,binom = binom)
			clust[[id]] = bin_noise_out$CLU
			if(FALSE){
				names(bin_noise_out$CLUST)
				bin_noise_out$CLUST[c("iter","converge","props","maf","psi")]
				table(bin_noise_out$INFER)
			}
			
			if( sum(bin_noise_out$CLUST$props[1:2]) == 0 ){
				next
			}
			
			vcs$INFER[idx] = bin_noise_out$INFER
			
			# Label all other qualified variants meeting RD, Qscore, tVAF thresholds
			tmp_index2 = which(sn_filter
				& advsta_filter
				& is.na(vcs$INFER)
				)
			length(tmp_index2)
			if(FALSE){
				plot(vcs[tmp_index2,c("x_coor","tVAF","Qscore")[1:2]],ylim = c(0,1),
					# col = factor(b_out$CLUST$class),
					# col = factor(b_out$INFER),pch = 16,
					cex = vcs$tDP_cex[tmp_index2])
			}
			ss = t(apply(vcs[tmp_index2,c("tAD","tRD")],1,function(xx) 
				BINOM_INFER(xx,bin_noise_out$CLUST$params)))
			ss = apply(ss,1,function(xx) ifelse(which.max(xx) %in% c(1,2),"noise","signal"))
			ss = as.character(ss)
			# table(ss)
			vcs$INFER[tmp_index2] = ss
			
			# Label all other RD, tVAF threshold variants as oxog_like or noise
			tmp_index3 = which(sn_filter
				& advst_filter
				& is.na(vcs$INFER)
				)
			ss = t(apply(vcs[tmp_index3,c("tAD","tRD")],1,function(xx) 
				BINOM_INFER(xx,bin_noise_out$CLUST$params)))
			ss = apply(ss,1,function(xx) ifelse(which.max(xx) %in% c(1,2),"noise","signal_like"))
			ss = as.character(ss)
			vcs$INFER[tmp_index3] = ss
		}
		
		cnt = cnt + 1
	}
	message("\n",appendLF = FALSE)
	names(vcs)[names(vcs) == "INFER"] = ARTI
	vcs = vcs[,names(vcs) != "Strand_bias_P"]
	
	list(vcs = vcs,clust = clust)
}

FFPE_find = function(vcs,strand,ad_thres,rd_thres,binom = TRUE,
	max_tVAF = 0.45,ID = "STUDYNUMBER",count_thres = 20){
	
	if(FALSE){
		vcs = vcf; strand = strand; ad_thres = 5; rd_thres = 10; binom = FALSE
		max_tVAF = 0.45; ID = "STUDYNUMBER"
	}
	
	req_COLS = c("RefAlt","tAD","tDP","tVAF",
		"mutID","Chr","Position","ThsdG_AF","EXAC_AF")
	smart_header(OBJ = vcs,req_COLS = req_COLS)
	
	message(sprintf("%s: Finding FFPE artifacts ...\n",date()),appendLF = FALSE)
	vcs = vcs[,names(vcs) != "FFPE"]
	vcs$INFER = NA
	
	# Append strand bias column
	message(sprintf("%s: Merge strand info ...\n",date()),appendLF = FALSE)
	vcs = vcs[,names(vcs) != "Strand_bias_P"]
	vcs = smart_merge(vcs,strand[,c("mutID","Chr",
		"Position","Ref","Alt","Strand_bias_P")],all.x=TRUE)
	
	# Filters
	art_filter 			= vcs$RefAlt %in% c("C>T","G>A","INDEL")
	ad_filter 			= vcs$tAD > 0
	dp_filter 			= vcs$tDP >= rd_thres
	vaf_filter 			= vcs$tVAF > 0 & vcs$tVAF < max_tVAF
	germ_filter 		= vcs$ThsdG_AF == 0 & vcs$EXAC_AF == 0
	sb_filter				= is.na(vcs$Strand_bias_P) | vcs$Strand_bias_P > 0.05
	target_filter 	= vcs$TARGET == "YES"
	
	IDs = sort(unique(vcs[[ID]]))
	cnt = 1
	for(id in IDs){
		# id = IDs[1]
		message(paste0(id," "),appendLF = FALSE)
		if(cnt %% 5 == 0) message("\n",appendLF = FALSE)
		sn_filter = vcs[[ID]] == id
		tmp_index = which(sn_filter
			& ad_filter & dp_filter
			& vaf_filter & sb_filter
			& germ_filter
			& target_filter
			& art_filter
			)
			
		if( length(tmp_index) >= count_thres ){
			bin_noise_out = BINOM_NOISE(mat_RD = as.matrix(vcs[tmp_index,c("tAD","tRD")]),binom = binom)
			vcs$INFER[tmp_index] = bin_noise_out$INFER
			
			# Label all other variants C>T/G>A variants meeting RD, Qscore, tVAF thresholds
			tmp_index2 = which(sn_filter
				& ad_filter & dp_filter
				& vaf_filter & sb_filter
				& target_filter
				& art_filter
				& is.na(vcs$INFER)
				)
			ss = t(apply(vcs[tmp_index2,c("tAD","tRD")],1,function(xx) 
				BINOM_INFER(xx,bin_noise_out$CLUST$params)))
			ss = apply(ss,1,function(xx) ifelse(which.max(xx) == 1,"noise","signal"))
			ss = as.character(ss)
			# smart_table(ss)
			vcs$INFER[tmp_index2] = ss
			
			# Label all other RD, Qscore, tVAF threshold variants as oxog_like or noise
			tmp_index3 = which(sn_filter
				& ad_filter & dp_filter
				& vaf_filter & sb_filter
				& target_filter
				& is.na(vcs$INFER)
				)
			ss = t(apply(vcs[tmp_index3,c("tAD","tRD")],1,function(xx) 
				BINOM_INFER(xx,bin_noise_out$CLUST$params)))
			ss = apply(ss,1,function(xx) ifelse(which.max(xx) == 1,"noise","signal_like"))
			ss = as.character(ss)
			vcs$INFER[tmp_index3] = ss
		}
		
		cnt = cnt + 1
	}
	message("\n",appendLF = FALSE)
	names(vcs)[names(vcs) == "INFER"] = "FFPE"
	vcs = vcs[,names(vcs) != "Strand_bias_P"]
	
	vcs
}

# Parse/Prep Tumor Only MultiCaller VCFs
filter_E_snpEff = function(in_df){
	# More robust function
	# in_df = one_npn_RDQ
	srch_set = c(
		"CHROMOSOME_LARGE DELETION","CHROMOSOME_LARGE_DELETION",
		"CHROMOSOME_LARGE_DUPLICATION","CHROMOSOME_LARGE_INVERSION",
		
		"CODON_CHANGE","CODON_DELETION","CODON_INSERTION",
		"CODON_CHANGE_PLUS_CODON_DELETION",
		"CODON_CHANGE_PLUS_CODON_INSERTION",
		
		"EXON","EXON_DELETED","EXON_DELETED_PARTIAL",
		"EXON_DUPLICATION","EXON_DUPLICATION_PARTIAL",
		"EXON_INVERSION","EXON_INVERSION_PARTIAL",
		
		"FRAME_SHIFT",
		
		"GENE_DELETED","GENE_DUPLICATION","GENE_FUSION","GENE_FUSION_HALF",
		"GENE_FUSION_REVERESE","GENE_REARRANGEMENT","GENE_FUSION_REVERSE",
		
		"NON_SYNONYMOUS_CODING","NON_SYNONYMOUS_START","NON_SYNONYMOUS_STOP",
		
		"PROTEIN_PROTEIN_INTERACTION_LOCUS","PROTEIN_STRUCTURAL_INTERACTION_LOCUS",
		"RARE_AMINO_ACID",
		
		"SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","SPLICE_SITE_REGION",
		"SPLICE_SITE_BRANCH","SPLICE_SITE_BRANCH_U12",
		
		"STOP_LOST","START_GAINED","START_LOST","STOP_GAINED",
		
		"SYNONYMOUS_CODING","SYNONYMOUS_START","SYNONYMOUS_STOP",
		
		"TRANSCRIPT_DELETED","UTR_3_DELETED","UTR_5_DELETED","NEXT_PROT"
		)
	tmp_index = sort(unique(unlist(sapply(srch_set,function(x) 
		grep(x,in_df$snpEff),USE.NAMES=FALSE))))
	in_df[tmp_index,]
}
effect_12 = function(eff){
	# transform input onco_effect into snpEff
	# there's an issue with onco's intron effect, doesn't always match snpEff
	if(is.na(eff)){
		eff1 = ""; eff2 = eff1
	} else if(eff == "3_prime_UTR_variant"){
		eff1 = "MODIFIER"; eff2 = "UTR_3_PRIME"
	} else if(eff == "5_prime_UTR_variant"){
		eff1 = "MODIFIER"; eff2 = "UTR_5_PRIME"
	} else if(eff == "frameshift_variant"){
		eff1 = "HIGH"; eff2 = "FRAME_SHIFT"
	} else if(eff == "inframe_deletion"){
		eff1 = "MODERATE"; eff2 = "CODON_DELETION"
	} else if(eff == "inframe_insertion"){
		eff1 = "MODERATE"; eff2 = "CODON_INSERTION"
	} else if(eff == "initiator_codon_variant"){
		eff1 = "LOW"; eff2 = "NON_SYNONYMOUS_START"
	} else if(eff == "intergenic_variant"){
		eff1 = "MODIFIER"; eff2 = "INTERGENIC"
	} else if(eff == "intragenic_variant"){
		eff1 = "MODIFIER"; eff2 = "INTRAGENIC"
	} else if(eff == "intron_variant"){
		eff1 = "MODIFIER"; eff2 = "INTRON"
	} else if(eff %in% c("missense","missense_variant")){
		eff1 = "MODERATE"; eff2 = "NON_SYNONYMOUS_CODING"
	} else if(eff == "splice_region_variant"){
		eff1 = "LOW"; eff2 = "SPLICE_SITE_REGION"
	} else if(eff == "stop_gained"){
		eff1 = "HIGH"; eff2 = "STOP_GAINED"
	} else if(eff == "stop_lost"){
		eff1 = "HIGH"; eff2 = "START_LOST"
	} else if(eff == "stop_retained_variant"){
		eff1 = "LOW"; eff2 = "SYNONYMOUS_STOP"
	} else if(eff == "synonymous_variant"){
		eff1 = "LOW"; eff2 = "SYNONYMOUS_CODING"
	} else if(eff == "upstream_gene_variant"){
		eff1 = "MODIFIER"; eff2 = "UPSTREAM"
	}
	c(eff1,eff2)
}
parse_EFF_string = function(effect_string){
	# Documentation of the EFF string from:  http://snpeff.sourceforge.net/SnpEff_manual.html#input
	# ##INFO=<ID=EFF,Number=.,Type=String,Description="Predicted effects for this variant.Format: 'Effect ( Effect_Impact |  Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Exon [ | ERRORS | WARNINGS ] )' \">
	# Effect_Impact can be {High, Moderate, Low, Modifer} and we need to return the HIGHEST of those, in that order
	if(FALSE){ # Testing
		effect_string = "INTRON(MODIFIER||||ARID1A|nonsense_mediated_decay|CODING|ENST00000466382|4),INTRON(MODIFIER||||ARID1A|protein_coding|CODING|ENST0000045
			7599|17),INTRON(MODIFIER||||ARID1A|protein_coding|CODING|ENST00000540690|4),STOP_GAINED(HIGH|NONSENSE|Cag/Tag|p.Gln1136*|ARID1A|protein_coding|CODING|ENST00000
			374152|17),STOP_GAINED(HIGH|NONSENSE|Cag/Tag|p.Gln1519*|ARID1A|protein_coding|CODING|ENST00000324856|18),UPSTREAM(MODIFIER||||ARID1A|nonsense_mediated_decay|CO
			DING|ENST00000532781|)"
		effect_string = "EXON(MODIFIER||||TUBB8P11|unprocessed_pseudogene|NON_CODING|ENST00000415481|1),INTRON(MODIFIER||||FAM41C|lincRNA|NON_CODING|ENST00000427857|2),INTRON(MODIFIER||||FAM41C|lincRNA|NON_CODING|ENST00000446136|2),UPSTREAM(MODIFIER||||FAM41C|lincRNA|NON_CODING|ENST00000432963|)"
		effect_string = check_second
	}
	
	if(FALSE){ # Older code from Stuart and Alan, this will return a list of one element
		effects = strsplit(effect_string,',')
			# INTRON(MODIFIER||||ARID1A|nonsense_mediated_decay|CODING|ENST00000466382|4)
		effStart = sapply(effects, "sub", pattern="\\|.*", replacement="")

		effList = strsplit (effStart, "\\(")

		if (length(grep ("HIGH",effList) )>0) {
			effValue = effList[grep("HIGH",effList)[1]]
		} else if (length (grep ("MODERATE",effList)) >0) {
			effValue = effList[grep("MODERATE",effList)[1]]
		} else if (length (grep ("LOW",effList)) >0) {
			effValue = effList[grep("LOW",effList)[1]]
		} else if (length (grep ("MODIFIER",effList)) >0) {
			effValue = effList[grep("MODIFIER",effList)[1]] 
		} else {
			effValue = c("UNKNOWN","UNKNOWN")
		}

		effValue
	}
	
	# Newer Code: Returns a two-element vector with Effect and Effect2
	effects = strsplit(effect_string,',')[[1]]
	# effMat = t(sapply(effects,function(x) strsplit(strsplit(x,"[|]")[[1]][1],"\\(")[[1]],USE.NAMES=FALSE))
	effMat = t(sapply(effects,function(x) strsplit(x,"[|()]")[[1]][c(1,2,5,9,10)],USE.NAMES=FALSE))
	effMat = unique(effMat)
	
	srch_HIGH = grep("HIGH",effMat[,2])
	srch_MODERATE = grep("MODERATE",effMat[,2])
	srch_LOW = grep("LOW",effMat[,2])
	srch_MODIFIER = grep("MODIFIER",effMat[,2])
	srch_index = NA
	effValue = c("UNKNOWN","UNKNOWN")
	
	if( length(srch_HIGH) > 0 ) {
		srch_index = srch_HIGH
		effValue[2] = "HIGH"
	} else if( length(srch_MODERATE) > 0 ) {
		srch_index = srch_MODERATE
		effValue[2] = "MODERATE"
	} else if( length(srch_LOW) > 0 ) {
		srch_index = srch_LOW
		effValue[2] = "LOW"
	} else if( length(srch_MODIFIER) > 0 ) {
		srch_index = srch_MODIFIER
		effValue[2] = "MODIFIER"
	}
	
	# effMat = as.matrix(effMat[,-2])
	if( length(srch_index) == 1 ){
		effValue[1] = paste(effMat[srch_index,-2],collapse="|")
	} else {
		effValue[1] = paste(apply(effMat[srch_index,-2],1,function(x) paste(x,collapse="|")),collapse=";")
	}
	effValue
	
}
parse_onco = function(one_info){
	if(FALSE){
		one_info = onco$INFO[1]
	}
	
	part_info0 = strsplit(one_info,";")[[1]]
	desired_annos = c("^gene=",
		"^COSMIC_n_overlapping_mutations",
		"^1000gp3_AF=",
		"^ExAC_AF=",
		"^variant_classification=",
		"^Ensembl_so_term=",
		"^EFF=",
		"^gc_content=",
		"^ref_context=",
		"^dbNSFP_MutationTaster_pred=",
		"^CGC_Cancer_Somatic_Mut=",
		"^CGC_Cancer_Germline_Mut=",
		"^ExAC_ESP_AF_GLOBAL=",
		#"^ExAC_KG_AF_GLOBAL=",
		"^HGNC_Ensembl_Gene_ID=",
		"^blarg=")
	part_info = part_info0[unlist(sapply(desired_annos,
		function(x) grep(x,part_info0),USE.NAMES=FALSE))]
	
	# Gene
	check_gene = grep("^gene=",part_info)
	tmp_gene = NA
	if( length(check_gene) > 0 ){
		tmp_gene = paste0("m_",part_info[check_gene])
		tmp_gene = strsplit(tmp_gene[grep("m_gene",tmp_gene)],"=")[[1]][2]
		if(tmp_gene == "Unknown") tmp_gene = NA
	}
	
	# Cosmic
	check_cosm = grep("^COSMIC_n_overlapping_mutations",part_info)
	tmp_cosm = 0
	if(length(check_cosm) > 0){
		tmp_cosm = strsplit(part_info[check_cosm],"=")[[1]][2]
	}
	
	# 1000G_AF
	check_1000G = grep("^1000gp3_AF=",part_info)
	tmp_1k = NA
	if(length(check_1000G) > 0){
		tmp_1k = strsplit(part_info[check_1000G],"=")[[1]][2]
		tmp_1k = min(strsplit(tmp_1k,"[|]")[[1]])
	}
	
	# Exac_AF
	check_exac = grep("^ExAC_AF=",part_info)
	tmp_exac = NA
	if(length(check_exac) > 0){
		tmp_exac = strsplit(part_info[check_exac],"=")[[1]][2]
	}
	
	# Variant classification
	check_var_class = grep("^variant_classification=",part_info)
	tmp_varclass = NA
	if(length(check_var_class) > 0){
		tmp_varclass = strsplit(part_info[check_var_class],"=")[[1]][2]
	}
	
	# Effect2
	check_eff = grep("^Ensembl_so_term=",part_info)
	tmp_eff = NA
	if(length(check_eff) > 0){
		tmp_eff = strsplit(part_info[check_eff],"=")[[1]][2]
	}
	
	# SnpEff
	check_snpEff = grep("^EFF=",part_info)
	tmp_snpEff = c("NA","NA-impact")
	if(length(check_snpEff) > 0){
		check_second = strsplit(part_info[check_snpEff],"=")[[1]][2]
		if( !is.na(check_second) ){
			tmp_snpEff = parse_EFF_string(check_second)
		}
	}
	
	# gc_content
	check_gc = grep("^gc_content=",part_info)
	tmp_gc = NA
	if(length(check_gc) > 0){
		tmp_gc = strsplit(part_info[check_gc],"=")[[1]][2]
	}
	
	# ref_context
	check_ref = grep("^ref_context=",part_info)
	tmp_ref = NA
	if(length(check_ref) > 0){
		tmp_ref = strsplit(part_info[check_ref],"=")[[1]][2]
	}
	
	# dbNSFP_MutationTaster_pred
	check_mutpred = grep("^dbNSFP_MutationTaster_pred=",part_info)
	tmp_mutpred = NA
	if(length(check_mutpred) > 0){
		tmp_mutpred = strsplit(part_info[check_mutpred],"=")[[1]][2]
	}
	
	# CGC Somatic
	check_cgcsom = grep("^CGC_Cancer_Somatic_Mut=",part_info)
	tmp_cgcsom = NA
	if(length(check_cgcsom) > 0){
		tmp_cgcsom = strsplit(part_info[check_cgcsom],"=")[[1]][2]
	}
	
	# CGC Germline
	check_cgcgerm = grep("^CGC_Cancer_Germline_Mut=",part_info)
	tmp_cgcgerm = NA
	if(length(check_cgcgerm) > 0){
		tmp_cgcgerm = strsplit(part_info[check_cgcgerm],"=")[[1]][2]
	}
	
	# ESP AF
	check_espAF = grep("^ExAC_ESP_AF_GLOBAL=",part_info)
	tmp_espAF = NA
	if(length(check_espAF) > 0){
		tmp_espAF = strsplit(part_info[check_espAF],"=")[[1]][2]
		tmp_espAF = ifelse(tmp_espAF == "NA",NA,tmp_espAF)
	}
	
	# ENSG
	check_ensg = grep("^HGNC_Ensembl_Gene_ID=",part_info)
	tmp_ensg = NA
	if(length(check_ensg) > 0){
		tmp_ensg = strsplit(part_info[check_ensg],"=")[[1]][2]
	}
	
	# Output
	c(tmp_gene,tmp_cosm,tmp_1k,tmp_exac,
		tmp_varclass,tmp_eff,tmp_snpEff,
		tmp_gc,tmp_ref,tmp_mutpred,
		tmp_cgcsom,tmp_cgcgerm,tmp_ensg,
		tmp_espAF)

}
get_counts = function(in_counts,vcf_format="GT:AD:DP"){
	# in_counts = onco[1,"NORMAL"]
	# in_counts = onco$NORMAL[1]
	# vcf_format = "GT:DP:AD"
	if(vcf_format == "GT:AD:DP"){
		spl_counts = strsplit(strsplit(in_counts,":")[[1]][2],"[|]")[[1]]
	} else if(vcf_format == "GT:DP:AD"){
		spl_counts = strsplit(strsplit(in_counts,":")[[1]][3],"[|]")[[1]]
	}
	as.numeric(spl_counts)
}
create_target = function(target_bed,input_df,pad,verbose=TRUE){
	if(FALSE){
		target_bed = aa$tumor_target_bed
		input_df = tmp_df
		pad = 0
	}
	
	# Newer Code
	if(verbose) message("Getting target info = ",appendLF = FALSE)
	input_df$TARGET = NA
	for(one_contig in sort(unique(input_df$Chr))){
		# one_contig = unique(input_df$Chr)[1]
		if(verbose) message(paste0(gsub("chr","",one_contig)," "),appendLF = FALSE)
		sub_target_bed = target_bed[which(target_bed$chrom == one_contig),]
		one_index = which(input_df$Chr == one_contig)
		input_df$TARGET[one_index] = Rcpp_get_target(intervals = as.matrix(sub_target_bed[,c("start","end")]),
			positions = input_df$Position[one_index],pad = pad)
	}
	if(verbose) message("\n",appendLF = FALSE)
	
	input_df
}
import_bed = function(bed_fn){
	message(sprintf("Import bed file: %s\n",bed_fn),appendLF = FALSE)
	bed_df = fread(file = bed_fn,
		sep = '\t',header = TRUE,
		data.table = FALSE,
		showProgress = FALSE)
	bed_df = bed_df[,c("chrom","start","end")]
	names(bed_df) = c("chrom","start","end")
	bed_df
}

# Annotation functions
get_INDEL_RefAlt_TRINUC = function(ALL_VCS){
	ALL_VCS$INDEL = ifelse(nchar(ALL_VCS$Ref) == 1 & nchar(ALL_VCS$Alt) == 1,FALSE,TRUE)
	index_bs = which(ALL_VCS$INDEL == FALSE)
	ALL_VCS$RefAlt = "INDEL"
	ALL_VCS$RefAlt[index_bs] = paste0(ALL_VCS$Ref[index_bs],">",ALL_VCS$Alt[index_bs])
	
	get_tri = function(alt,refcon){
		if(FALSE){
			alt = ALL_VCS$Alt[1]
			refcon = ALL_VCS$ref_context[1]
		}
		
		refcon = toupper(refcon)
		rr = strsplit(refcon,"")[[1]][10:12]
		paste0(rr[1],"[",rr[2],">",alt,"]",rr[3])
	}
	ALL_VCS$TRINUC = "INDEL"
	ALL_VCS$TRINUC[index_bs] = apply(ALL_VCS[index_bs,c("Alt","ref_context")],1,
		function(xx) get_tri(xx[1],xx[2]))
	ALL_VCS
}

#' @title new_STRAND
#' @description Extracts strand-specific read counts from a 
#'	bam file.
#' @param BAM_fn A character string for the bam's file path.
#' @param vcs A R data.frame of variant calls containing 
#'	columns \code{Chr}, \code{Position}, \code{Ref}, and 
#'	\code{Alt} corresponding to the chromosome,
#'	genomic position, reference base(s), and alternate base(s),
#'	respectively.
#' @param minBQ An integer for minimum base quality.
#' @param minMQ An integer for minimum mapping quality.
#' @param NT A positive integer for the number of threads.
#' 
new_STRAND = function(BAM_fn,vcs,minBQ = 13,minMQ = 40,NT = 1){
	if(FALSE){
		BAM_fn = nBAM_fn; vcs = svi[1:20,]; NT = smart_ncores()
	}
	
	vcs = unique(vcs[,c("mutID","Chr","Position","Ref","Alt")])
	message(sprintf("%s: Total number of loci for strand = %s\n",date(),nrow(vcs)),appendLF = FALSE)
	bf = BamFile(BAM_fn)
	p_param = PileupParam(max_depth = 1e4,
		min_base_quality = minBQ,min_mapq = minMQ,
    min_nucleotide_depth = 1,min_minor_allele_depth = 0,
    distinguish_strands = TRUE,distinguish_nucleotides = TRUE,
    ignore_query_Ns = TRUE,include_deletions = TRUE,
		include_insertions = TRUE,
    left_bins = NULL,query_bins = NULL,cycle_bins = NULL)
	cl = makeCluster(NT) # NT = num threads
	registerDoParallel(cl)
	ii = NULL
	xx = foreach(ii=seq(nrow(vcs)),.combine="rbind") %dopar% {
		# ii = 2
		chr = vcs$Chr[ii]; pos = vcs$Position[ii]
		ref = vcs$Ref[ii]; alt = vcs$Alt[ii]
		param = ScanBamParam(which=GRanges(chr,IRanges(start = pos,end = pos)))
		pp = pileup(bf,scanBamParam = param,pileupParam = p_param)
		pp = pp[which(pp$nucleotide %in% c(ref,alt)),]
		
		rev_ref = pp$count[which(pp$strand == "-" & pp$nucleotide == ref)]
		rev_alt = pp$count[which(pp$strand == "-" & pp$nucleotide == alt)]
		for_ref = pp$count[which(pp$strand == "+" & pp$nucleotide == ref)]
		for_alt = pp$count[which(pp$strand == "+" & pp$nucleotide == alt)]
		
		rev_ref = ifelse(length(rev_ref) > 0,rev_ref,0)
		rev_alt = ifelse(length(rev_alt) > 0,rev_alt,0)
		for_ref = ifelse(length(for_ref) > 0,for_ref,0)
		for_alt = ifelse(length(for_alt) > 0,for_alt,0)
		
		mm = matrix(c(rev_ref,rev_alt,for_ref,for_alt),ncol=2,byrow=TRUE)
		colnames(mm) = c("REF","ALT"); rownames(mm) = c("REV","FOR")
		fish_out = fisher.test(mm)
		c(mm[2,1],mm[2,2],mm[1,1],mm[1,2],fish_out$p.value)
	}
	stopCluster(cl)
	
	rownames(xx) = NULL
	xx = smart_df(xx)
	names(xx) = c("For_Ref","For_Alt","Rev_Ref",
		"Rev_Alt","Strand_bias_P")
	vcs = smart_df(vcs,xx)
	vcs
}

get_EXONIC = function(npn){
	if(FALSE){
		npn = all_vcs
	}
	
	uniq_npn = unique(npn[,c("mutID","snpEff")])
	uniq_npn$EXONIC = "NO"
	uniq_npn$EXONIC[which(uniq_npn$mutID %in% filter_E_snpEff(uniq_npn)$mutID)] = "YES"
	uniq_npn = uniq_npn[,c("mutID","EXONIC")]
	npn = smart_merge(npn[,which(names(npn) != "EXONIC")],uniq_npn)
	npn
}
get_INFER = function(vec_RD,params){
	if(FALSE){
		vec_RD = c(530,443); vec_RD; sum(vec_RD)
			vec_RD = round(vec_RD / 2); vec_RD; sum(vec_RD)
		params = c(0.15403,0,0.34672,0.50387)
	}
	
	if( is.na(params[1]) ){
		return(NA)
	}
	
	# sometimes approx 0 values are like -7e-17 or 1.00000 > 1
	params[which(params > 1)] = 1
	params[which(params < 0)] = 0
	
	if( params[1] == 0 ){
		params[1] = 1e-3 
		# Since oxoG variants were pre-filtered, this might have been our "absent noise"
	}
	
	vec_alpha = c(params[1],
		(1-params[1])*params[2],
		(1-params[1])*(1-params[2])*params[4],
		(1-params[1])*(1-params[2])*(1-params[4]))
	VAF = vec_RD[1] / sum(vec_RD)
	# Infer segment status
	stat = get_ALLELE_STAT(params)
	# Infer VC status
	DP = sum(vec_RD)
	LBC = lchoose(DP,vec_RD[1])
	log_vec = rep(NA,4)
	log_vec[1] = log(vec_alpha[1]) - log(DP)
	log_vec[2] = log(vec_alpha[2]) + Rcpp_LL(RD = t(vec_RD),DP = DP,LBC = LBC,pi = 0.5,psi = params[5])
	log_vec[3] = log(vec_alpha[3]) + Rcpp_LL(RD = t(vec_RD),DP = DP,LBC = LBC,pi = params[3],psi = params[5])
	log_vec[4] = log(vec_alpha[4]) + Rcpp_LL(RD = t(vec_RD),DP = DP,LBC = LBC,pi = 1-params[3],psi = params[5])
	
	log_vec = as.numeric(exp(log_vec - Rcpp_logSumExp(log_vec)))
	# which.max(log_vec)
	
	if( is.na(stat) || stat == "noise" ){
		NA
	} else if( stat == "balance" ){
		if( VAF >= 0.5 ){
			0
		} else {
			sum(log_vec[c(1,3,4)])
		}
	} else if( stat == "imbalance" ){
		if( VAF >= max(c(params[3],1-params[3])) ){
			0
		} else {
			sum(log_vec[c(1,2)])
		}
	}
	
}
get_ALLELE_STAT = function(params){
	if( is.na(params[1]) ){
		NA
	} else if(params[1] == 1){
		"noise"
	} else {
		if(params[2] >= 0.5){
			"balance"
		} else {
			"imbalance"
		}
	}
}
get_INFER_geno = function(vec_RD,params){
	if(FALSE){
		ii = which(rownames(vcf) == 7183); ii
		vec_RD = c(vcf$tAD[ii],vcf$tRD[ii]); vec_RD
		params = c(vcf$T_eps[ii],vcf$T_prob_half[ii],
			vcf$T_mBAF[ii],vcf$T_prop_mBAF[ii],vcf$T_PSI[ii]); params
		
		vec_RD = c(603,212)
		params = c(0.3594971,0,0.3511942,1,0.01011822)
		
	}
	
	if( is.na(params[1]) ){
		return(NA)
	}
	
	# sometimes approx 0 values are like -7e-17 or 1.00000 > 1
	params[which(params > 1)] = 1
	params[which(params < 0)] = 0
	
	if( params[1] == 0 ){
		params[1] = 1e-3 
		# Since oxoG variants were pre-filtered, this might have been our "absent noise"
	}
	
	vec_alpha = c(params[1],
		(1-params[1])*params[2],
		(1-params[1])*(1-params[2])*params[4],
		(1-params[1])*(1-params[2])*(1-params[4]))
	VAF = vec_RD[1] / sum(vec_RD)
	# Infer segment status
	stat = get_ALLELE_STAT(params)
	# Infer VC status
	DP = sum(vec_RD)
	LBC = lchoose(DP,vec_RD[1])
	log_vec = rep(NA,4)
	log_vec[1] = log(vec_alpha[1]) - log(DP)
	log_vec[2] = log(vec_alpha[2]) + Rcpp_LL(RD = t(vec_RD),DP = DP,LBC = LBC,pi = 0.5,psi = params[5])
	log_vec[3] = log(vec_alpha[3]) + Rcpp_LL(RD = t(vec_RD),DP = DP,LBC = LBC,pi = params[3],psi = params[5])
	log_vec[4] = log(vec_alpha[4]) + Rcpp_LL(RD = t(vec_RD),DP = DP,LBC = LBC,pi = 1-params[3],psi = params[5])
	log_vec = as.numeric(exp(log_vec - Rcpp_logSumExp(log_vec)))
	# c("noise",
	# which.max(log_vec)
	
	# stop("Update code below")
	
	# New Code
	stat
	prob_geno_cutoff = 0.5
	out_geno = NULL
	if( VAF == 1 ){
		return("1/1")
	} else if( VAF == 0 ){
		return("0/0")
	} else if( is.na(stat) || stat == "noise" ){
		prob_notBAF = NA
		out_geno = "./." # no genotype call
	} else if( stat == "balance" ){
		# prob vc not at BAF = 0.5
		prob_notBAF = sum(log_vec[c(1,3,4)]); prob_notBAF
		if( prob_notBAF > prob_geno_cutoff ){
			if( VAF <= 0.5 ){
				out_geno = "0/0"
			} else {
				out_geno = "1/1"
			}
		} else {
			out_geno = "0/1"
		}
	} else if( stat == "imbalance" ){
		# prob vc not at BAF = mBAF or 1 - mBAF
		prob_notBAF = sum(log_vec[c(1,2)]); prob_notBAF
		if( prob_notBAF > prob_geno_cutoff ){
			small_BAF = min(c(params[3],1-params[3])); small_BAF
			large_BAF = max(c(params[3],1-params[3])); large_BAF
			if( VAF < small_BAF ){
				out_geno = "0/0" # genotype AA or Ref/Ref
			} else if( VAF > large_BAF ){
				out_geno = "1/1" # genotype BB or Alt/Alt
			} else if( VAF > small_BAF && VAF < large_BAF ){
				out_geno = "./." # not sure
			}
		} else {
			out_geno = "0/1" # genotype AB
		}
	}
	
	if( is.null(out_geno) ) stop("Something missing")
	return(out_geno)
	
}
anno_SEG = function(VC,nSEG,tSEG,eps_thres = 0.5,psi_thres = 0.01){
	if(FALSE){
		VC = uniq_vcs; nSEG = nSEG; tSEG = tSEG
		
	}
	
	clust_names = c("eps","prob_half","mBAF","prop_mBAF","PSI","ms")
	
	# Append nSEG results
	message("nANNO ...\n",appendLF = FALSE)
	VC$n_seg = NA
	nSEG$n_seg = NA
	num_nSEG = nrow(nSEG)
	for(ii in seq(num_nSEG)){
		# ii = 1
		if( ii %% 2 == 0 ) message(".",appendLF = FALSE)
		if( ii %% 5e1 == 0 || ii == num_nSEG ) message(sprintf("%s out of %s\n",ii,num_nSEG),appendLF = FALSE)
		
		nSEG$n_seg[ii] = ii
		chr = nSEG$Chr[ii]; chr
		start_pos = nSEG$Start_Position[ii]; start_pos
		end_pos = nSEG$End_Position[ii]; end_pos
		idx = which(VC$Chr == chr
			& VC$Position >= start_pos
			& VC$Position <= end_pos)
		if( length(idx) > 0 ){
			VC$n_seg[idx] = ii
		}
		
		rm(chr,start_pos,end_pos,idx)
	}
	idx_NA = which(is.na(VC$n_seg))
	num_NA_n_seg = length(idx_NA); num_NA_n_seg
	if( num_NA_n_seg > 0 ){
		for(ii in seq(num_NA_n_seg)){
			# ii = 1
			
			chr = VC$Chr[idx_NA[ii]]; chr
			pos = VC$Position[idx_NA[ii]]; pos
			nSEG_chr = nSEG[which(nSEG$Chr == chr),]; nSEG_chr
			
			if( nrow(nSEG_chr) == 0 ) next
			
			# Find closest segment to position
			diff_start_pos = abs(nSEG_chr$Start_Position - pos)
				min_start_diff = diff_start_pos[which.min(diff_start_pos)]
			diff_end_pos = abs(nSEG_chr$End_Position - pos)
				min_end_diff = diff_end_pos[which.min(diff_end_pos)]
			
			if( min_start_diff < min_end_diff ){
				VC$n_seg[idx_NA[ii]] = nSEG_chr$n_seg[which.min(diff_start_pos)]
			} else {
				VC$n_seg[idx_NA[ii]] = nSEG_chr$n_seg[which.min(diff_end_pos)]
			}
			
			rm(chr,pos,nSEG_chr,diff_start_pos,diff_end_pos)
		}
	}
	# smart_table(VC$Chr,!is.na(VC$n_seg))
	VC = smart_rmcols(VC,paste0("N_",clust_names))
	nSEG_2 = nSEG[,c("n_seg",clust_names)]
	names(nSEG_2)[-1] = paste0("N_",names(nSEG_2)[-1])
	# nSEG_2[1:3,]; VC[1:3,]
	VC = smart_merge(VC,nSEG_2,all.x = TRUE)
	VC = smart_rmcols(VC,"n_seg")
	rm(nSEG_2)
	message("\n",appendLF = FALSE)
	
	# Infer H2M
	message("Infer H2M status ...\n",appendLF = FALSE)
	VC$nANNO = ifelse(VC$N_eps >= eps_thres | VC$N_PSI >= psi_thres,"H2M","mappable")

	# Append tSEG results
	message("tANNO ...\n",appendLF = FALSE)
	VC$t_seg = NA
	tSEG$t_seg = NA
	num_tSEG = nrow(tSEG)
	for(ii in seq(num_tSEG)){
		# ii = 1
		if( ii %% 2 == 0 ) message(".",appendLF = FALSE)
		if( ii %% 5e1 == 0 || ii == num_tSEG ) message(sprintf("%s out of %s\n",ii,num_tSEG),appendLF = FALSE)
		
		tSEG$t_seg[ii] = ii
		chr = tSEG$Chr[ii]; chr
		start_pos = tSEG$Start_Position[ii]; start_pos
		end_pos = tSEG$End_Position[ii]; end_pos
		idx = which(VC$Chr == chr
			& VC$Position >= start_pos
			& VC$Position <= end_pos)
		if( length(idx) > 0 ){
			VC$t_seg[idx] = ii
		}
		
		rm(chr,start_pos,end_pos,idx)
	}
	idx_NA = which(is.na(VC$t_seg))
	num_NA_t_seg = length(idx_NA); num_NA_t_seg
	if( num_NA_t_seg > 0 ){
		for(ii in seq(num_NA_t_seg)){
			# ii = 1
			
			chr = VC$Chr[idx_NA[ii]]; chr
			pos = VC$Position[idx_NA[ii]]; pos
			tSEG_chr = tSEG[which(tSEG$Chr == chr),]; tSEG_chr
			
			if( nrow(tSEG_chr) == 0 ) next
			
			# Find closest segment to position
			diff_start_pos = abs(tSEG_chr$Start_Position - pos)
				min_start_diff = diff_start_pos[which.min(diff_start_pos)]
			diff_end_pos = abs(tSEG_chr$End_Position - pos)
				min_end_diff = diff_end_pos[which.min(diff_end_pos)]
			
			if( min_start_diff < min_end_diff ){
				VC$t_seg[idx_NA[ii]] = tSEG_chr$t_seg[which.min(diff_start_pos)]
			} else {
				VC$t_seg[idx_NA[ii]] = tSEG_chr$t_seg[which.min(diff_end_pos)]
			}
			
			rm(chr,pos,tSEG_chr,diff_start_pos,diff_end_pos)
		}
	}
	# smart_table(VC$Chr,!is.na(VC$t_seg))
	VC = smart_rmcols(VC,paste0("T_",clust_names))
	tSEG_2 = tSEG[,c("t_seg",clust_names)]
	names(tSEG_2)[-1] = paste0("T_",names(tSEG_2)[-1])
	# tSEG_2[1:3,]; VC[1:3,]
	VC = smart_merge(VC,tSEG_2,all.x = TRUE)
	VC = smart_rmcols(VC,"t_seg")
	rm(tSEG_2)
	message("\n",appendLF = FALSE)
	
	# Infer BAF status and Posterior probability of BAF
	message("Infer ALLELE_STAT and tANNO ...\n",appendLF = FALSE)
	VC$ALLELE_STAT = apply(VC[,paste0("T_",clust_names)],1,
		function(x) get_ALLELE_STAT(x))
	VC$tANNO = apply(VC[,c("tAD","tRD",paste0("T_",clust_names))],1,
		function(x) get_INFER(vec_RD = x[1:2],params = x[3:8]))
	if(FALSE){
		xx = apply(VC[51:70,c("tAD","tRD",paste0("T_",clust_names))],1,
			function(x) get_INFER(vec_RD = x[1:2],params = x[3:8]))
		
	}
	message("\n",appendLF = FALSE)
	
	# Output
	return(VC)
	
}

# Tumor-only Filtering
proto_simple_TO_pipe = function(npn,rd_thres,ad_thres,qscore_thres,exac_thres,
	prop_call = 0.9,eps_thres = 0.5,psi_thres = 0.01,noise_filter = TRUE,segANNO_filter = FALSE){

	# Want to label VCs as 
	# 1) potential SPMs (R=0)
	# 2) potential Germline VCs (R=1)
	# 3) noise VCs
	
	# Retain TARGET and EXONIC
	npn = npn[which(npn$TARGET == "YES" & npn$EXONIC == "YES"),]
	clust_names = c("eps","prob_half","mBAF","prop_mBAF","PSI","ms")
	
	# Pre-defined thresholds
	# prop_G0 = 0.5
	prop_nonNoise = 0.75
	canon_thres = 10
	my_1000G_thres = 5e-3
	
	min_num_normals = length(unique(npn$STUDYNUMBER))
	npn$nLabel_2 = factor(npn$nLabel,levels = c("G=0","G=1","G=2","noise"))
	mutID_tab = table(npn[,c("mutID","nLabel_2")])
	npn = npn[,names(npn) != "nLabel_2"]
	mutID_df = smart_df(mutID = as.numeric(rownames(mutID_tab)),
		matrix(mutID_tab,ncol=4,byrow=FALSE))
	rm(mutID_tab)
	num_cols = c("G0","G1","G2","noise")
	names(mutID_df)[-1] = num_cols
	
	# Merge in Annotations
	anno_df = smart_uniq_df(npn[,c("mutID","CosmicOverlaps","EXAC_AF","ThsdG_AF",
		"TARGET","EXONIC",paste0("N_",clust_names),"nANNO",
		paste0("T_",clust_names),"ALLELE_STAT","tANNO")],
		vars=c("mutID"))
	anno_df = anno_df[order(anno_df$mutID),]
	mutID_df = cbind(mutID_df,anno_df[,-1])
	
	# Calculate and combine sequencing stats
	mat_RDQ = as.matrix(npn[,c("mutID","tAD","tDP","nDP","Qscore")])
	if( any(is.na(mat_RDQ[,5])) ){
		mat_RDQ[is.na(mat_RDQ[,5]),5] = 0
	}
	DPQ_stat = Rcpp_RDQ(RDQ = mat_RDQ,qq = 0.25)
	DPQ_stat = smart_df(DPQ_stat)
	names(DPQ_stat) = c("mutID","tAD","tDP","nDP_stat","Q_stat")
	DPQ_stat = DPQ_stat[order(DPQ_stat$mutID),]
	mutID_df = cbind(mutID_df,DPQ_stat[,-1])
	mutID_df$decide = ""
	mutID_df = mutID_df[which(mutID_df$EXONIC == "YES" & mutID_df$TARGET == "YES"),]
	mutID_df$G0_rowSums = rowSums(mutID_df[,num_cols])
	
	# Filters
	cosmic_filter 	= mutID_df$CosmicOverlaps >= canon_thres
	qscore_filter 	= mutID_df$Q_stat >= qscore_thres
	nDP_filter 			= mutID_df$nDP_stat >= rd_thres
	tDP_filter 			= mutID_df$tDP >= rd_thres
	tAD_filter 			= mutID_df$tAD >= ad_thres
	exac_filter 		= mutID_df$EXAC_AF < exac_thres
	thsdg_filter 		= mutID_df$ThsdG_AF == 0 | 
		(mutID_df$ThsdG_AF < my_1000G_thres & mutID_df$CosmicOverlaps >= 50)
	propcall_filter = mutID_df$G0_rowSums / min_num_normals >= prop_call
	H2M_filter 			= is.na(mutID_df$N_eps) | ( mutID_df$N_eps < eps_thres & mutID_df$N_PSI < psi_thres)
	germ_filter			= ( is.na(mutID_df$tANNO) | mutID_df$tANNO >= 0.5 )
	# propG0_filter		= mutID_df$G0 / mutID_df$G0_rowSums >= prop_G0
	propNN_filter		= rowSums(mutID_df[,num_cols[-4]]) / mutID_df$G0_rowSums >= prop_nonNoise
	
	# Extract candidate SPMs
	if( noise_filter && segANNO_filter ){
		mutID_df$decide[which( mutID_df$decide == ""
			& ( cosmic_filter 
				| ( 
					# propG0_filter
					propNN_filter
					& H2M_filter & germ_filter
					& qscore_filter & nDP_filter
					& tDP_filter & tAD_filter ) )
			& propcall_filter
			& exac_filter & thsdg_filter
			)] = "cand_SPM"
	} else if( noise_filter && !segANNO_filter ){
		mutID_df$decide[which( mutID_df$decide == ""
			& ( cosmic_filter 
				| ( 
					# propG0_filter
					propNN_filter
					& qscore_filter & nDP_filter
					& tDP_filter & tAD_filter ) )
			& propcall_filter
			& exac_filter & thsdg_filter
			)] = "cand_SPM"
	} else if( !noise_filter && !segANNO_filter ){
		mutID_df$decide[which( mutID_df$decide == ""
			& ( cosmic_filter 
				| ( qscore_filter & nDP_filter
					& tDP_filter & tAD_filter ) )
			& propcall_filter
			& exac_filter & thsdg_filter
			)] = "cand_SPM"
	}
	mutID_df = mutID_df[which(mutID_df$decide == "cand_SPM"),]
	
	out_fNPN = smart_uniq_df(npn[which(npn$mutID %in% mutID_df$mutID),],vars="mutID")
	rm_vars = c("FORMAT","NORMAL","TUMOR","TARGET","EXONIC","Onco_Effect","STUDYNUMBER",
		"nAD","nRD","nDP","nVAF","nLabel","nInterpret","Qscore","ID")
	out_fNPN = smart_merge(out_fNPN,mutID_df[,c("mutID","nDP_stat","Q_stat")])
	out_fNPN = out_fNPN[,which(!(names(out_fNPN) %in% rm_vars))]
	out_fNPN
}
TO_pipe_add_prop_nlabel_ART = function(npn,rd_thres,ad_thres,qscore_thres,exac_thres,
	prop_call = 0.9,eps_thres = 0.5,psi_thres = 0.01,noise_filter = TRUE,segANNO_filter = FALSE){
	if(FALSE){
		npn = vcf_TE;
		prop_call = 0.9; noise_filter = TRUE; segANNO_filter = FALSE
		
	}
	
	fNPN = proto_simple_TO_pipe(npn = npn,
		rd_thres = rd_thres,
		ad_thres = ad_thres,
		qscore_thres = qscore_thres,
		exac_thres = exac_thres,
		prop_call = prop_call,
		eps_thres = eps_thres,
		psi_thres = psi_thres,
		noise_filter = noise_filter,
		segANNO_filter = segANNO_filter)
	
	fNPN$prop_nLabel = NA
	fNPN$OXOG = NA
	fNPN$FFPE = NA
	fNPN$ARTI = NA
	for(ii in seq(nrow(fNPN))){
		# ii = 1
		tmp_index = which(npn$mutID == fNPN$mutID[ii])
		tab = smart_table(npn$nLabel[tmp_index])
		tmp_value = paste0(names(tab),"(",as.numeric(tab),",",round(as.numeric(tab/sum(tab))*100,2),"%)")
		fNPN$prop_nLabel[ii] = paste(tmp_value,collapse=";")
		
		# Properly annotate OXOG column
		tmp_tab = smart_table(npn$OXOG[tmp_index])
		tmp_value = paste0(names(tmp_tab),"(",as.numeric(tmp_tab),",",round(as.numeric(tmp_tab/sum(tmp_tab))*100,2),"%)")
		fNPN$OXOG[ii] = paste(tmp_value,collapse=";")
		
		# Properly annotate FFPE column
		tmp_tab = smart_table(npn$FFPE[tmp_index])
		tmp_value = paste0(names(tmp_tab),"(",as.numeric(tmp_tab),",",round(as.numeric(tmp_tab/sum(tmp_tab))*100,2),"%)")
		fNPN$FFPE[ii] = paste(tmp_value,collapse=";")
		
		# Properly annotate ARTI column
		tmp_tab = smart_table(npn$ARTI[tmp_index])
		tmp_value = paste0(names(tmp_tab),"(",as.numeric(tmp_tab),",",round(as.numeric(tmp_tab/sum(tmp_tab))*100,2),"%)")
		fNPN$ARTI[ii] = paste(tmp_value,collapse=";")
		
	}
	
	# Output
	fNPN
}
plot_TO = function(outdir,tumorID,all_vcs,fNPN,genome,hg,all_tVAF_segs){
	png(filename = file.path(outdir,paste0("TO_",tumorID,".png")),units = 'px',
			height = 1500,width = 3000,res = 250,type = 'cairo')
		par(mar = c(5,4,2,7),xpd = TRUE)
		all_vcs = unique(all_vcs[,c("mutID","x_coor","tVAF","tDP")])
		all_vcs$tDP_cex = 0.5*log10(all_vcs$tDP)
		genome_plot(all_vcs,var = "tVAF",cex = all_vcs$tDP_cex,
			how_color = rgb(0,1,0,0.4),
			pch = 16,hg = hg,main = tumorID,genome = genome)
		apply(all_tVAF_segs[which(all_tVAF_segs$eps < 1),
				c("start_x_coor","end_x_coor","iBAF")],1,
			function(x) segments(x0=x[1],x1=x[2],y0=c(x[3],1-x[3]),
			col="darkorange1",lwd=4))
		fNPN2 = fNPN
		fNPN2$base_indel_status = ifelse(nchar(fNPN2$Ref) == 1 
			& nchar(fNPN2$Alt) == 1,"BS","INDEL")
		points(fNPN2[,c("x_coor","tVAF")],pch = 16,cex = 1)
		my_cex = 1.5
		if( length(which(fNPN2$base_indel_status=="INDEL")) > 0 ){ # Label Indels
			points(fNPN2[which(fNPN2$base_indel_status=="INDEL"),c("x_coor","tVAF")],
				pch = 6,cex = my_cex,lwd = 2,col = "red")
		}
		if( length(which(fNPN2$CosmicOverlaps >= 10)) > 0 ){ # Label cosmic
			points(fNPN2[which(fNPN2$CosmicOverlaps >= 10),c("x_coor","tVAF")],
				pch = 5,cex = my_cex,lwd = 2,col = "magenta")
		}
		if( length(which(fNPN2$nANNO == "H2M")) > 0 ){ # Label H2M
			points(fNPN2[which(fNPN2$nANNO == "H2M"),c("x_coor","tVAF")],
				pch = 4,cex = my_cex,lwd = 2,col = "deeppink")
		}
		leg_df = smart_df(text = c("Candidate SPM","COSMIC >= 10","INDEL","H2M"),
			lwd = 3,col = c("black","magenta","red","deeppink"),
			pch = c(16,5,6,4))
		legend("right",inset = c(-0.12,0),legend = leg_df$text,
			cex = 0.7,col = leg_df$col,bty = "n",pch = leg_df$pch,
			title="Variant Status")
		par(xpd = FALSE)
		abline(h = seq(0.1,0.9,0.1),lty = 2,lwd = 0.5,col = rgb(0,0,0,0.2))
		par(mar = c(5,4,4,2)+0.1)
	dev.off()
}


#' @title UNMASC_definitions
#' @description This function prints the column header definitions of 
#'	the output file \code{tumorOnly_VCs.tsv}.
#' @param headers A character vector of column header names. By default
#'	\code{headers = NULL}, in which all available header definitions 
#'	are printed. If a subset of headers are specified, only the subset 
#'	of header definitions are printed.
#' @return Null value, prints definitions to console.
#' @export
UNMASC_definitions = function(headers = NULL){
	readme_VC = c()
	readme_VC = rbind(readme_VC,smart_df(VAR="mutID",DEF="Sample-specific mutation ID number"))
	readme_VC = rbind(readme_VC,smart_df(VAR="Chr",DEF="Chromosome/Contig"))
	readme_VC = rbind(readme_VC,smart_df(VAR="Position",DEF="Genomic position"))
	readme_VC = rbind(readme_VC,smart_df(VAR="Ref",DEF="Reference base(s)"))
	readme_VC = rbind(readme_VC,smart_df(VAR="Alt",DEF="Alternate base(s)"))
	readme_VC = rbind(readme_VC,smart_df(VAR="GeneName",DEF="Gene name/ Hugo Symbol"))
	readme_VC = rbind(readme_VC,smart_df(VAR="CosmicOverlaps",DEF="Total number of COSMIC mutations at variant site"))
	readme_VC = rbind(readme_VC,smart_df(VAR="ThsdG_AF",DEF="1000 Genomes population allele frequency"))
	readme_VC = rbind(readme_VC,smart_df(VAR="EXAC_AF",DEF="ExAC population allele frequency"))
	readme_VC = rbind(readme_VC,smart_df(VAR="snpEff",DEF="snpEff annotation"))
	readme_VC = rbind(readme_VC,smart_df(VAR="IMPACT",DEF="worst case snpEff Impact Prediction"))
	readme_VC = rbind(readme_VC,smart_df(VAR="gc",DEF="GC Content"))
	readme_VC = rbind(readme_VC,smart_df(VAR="ref_context",DEF="Genomic sequence at variant locus with additional 10 bp of flanking sequence on either side"))
	readme_VC = rbind(readme_VC,smart_df(VAR="tRD",DEF="Tumor Reference read depth"))
	readme_VC = rbind(readme_VC,smart_df(VAR="tAD",DEF="Tumor Alternate read depth"))
	readme_VC = rbind(readme_VC,smart_df(VAR="tDP",DEF="Tumor total read depth (tRD+tAD)"))
	readme_VC = rbind(readme_VC,smart_df(VAR="tVAF",DEF="Tumor Variant Allele Frequency (tAD/tDP)"))
	readme_VC = rbind(readme_VC,smart_df(VAR="x_coor",DEF="Coordinate for plotting variants"))
	readme_VC = rbind(readme_VC,smart_df(VAR=paste0("N_",c("eps","prob_half","ms")),
		DEF=paste0("Local normal ",c("noise","proportion of points near 0.5","model size")," parameter")))
	readme_VC = rbind(readme_VC,smart_df(VAR="nANNO",DEF="Local normal conclusion: H2M = Hard to map"))
	readme_VC = rbind(readme_VC,smart_df(VAR=paste0("T_",c("eps","prob_half","mBAF","prop_mBAF","ms")),
		DEF=paste0("Local tumor ",c("noise","proportion of points near 0.5",
			"tVAF cluster mean","proportion of points near tVAF","model size")," parameter")))
	readme_VC = rbind(readme_VC,smart_df(VAR="ALLELE_STAT",DEF="Local tumor conclusion: balance, imbalance, noise = no clear BAF pattern"))
	readme_VC = rbind(readme_VC,smart_df(VAR="tANNO",DEF="Posterior probability candidate variant is not BAF-like. If NA, tVAF clusters aren't distinguishable."))
	readme_VC = rbind(readme_VC,smart_df(VAR=c("nDP_stat","Q_stat"),DEF=paste0("25th percentile normal ",c("total read depth","Qscore")," across all UMNs")))
	readme_VC = rbind(readme_VC,smart_df(VAR="prop_nLabel",DEF="Proportion of variant's corresponding genotype classification"))
	readme_VC = rbind(readme_VC,smart_df(VAR="INDEL",DEF="TRUE=variant is an indel,FALSE=variant is a base substitution"))
	readme_VC = rbind(readme_VC,smart_df(VAR="RefAlt",DEF="Ref>Alt for base substitution,NA for indels"))
	readme_VC = rbind(readme_VC,smart_df(VAR="TRINUC",DEF="NA=variant is an indel, otherwise reporting the trinucleotide base substition"))
	readme_VC = rbind(readme_VC,smart_df(VAR="OXOG",DEF=paste0("If RefAlt == 'C>A', RefAlt == 'G>T', or has ",
		"a similar tVAF, variant considered 'signal'(aka artifact), 'signal_like' (aka predicted artifact), else 'noise'(aka real).")))
	readme_VC = rbind(readme_VC,smart_df(VAR="FFPE",DEF=paste0("If RefAlt == 'C>T', RefAlt == 'G>A', an indel,",
		" or has a similar tVAF, variant considered 'signal'(aka artifact), else 'noise'(aka real).")))
	readme_VC = rbind(readme_VC,smart_df(VAR="ARTI",DEF=paste0("If RefAlt == 'A>C', RefAlt == 'T>G', or has ",
		"a similar tVAF, variant considered 'signal'(aka artifact), 'signal_like' (aka predicted artifact), else 'noise'(aka real).")))
	readme_VC = rbind(readme_VC,smart_df(VAR="Ref_for",DEF="Forward read counts harboring reference base"))
	readme_VC = rbind(readme_VC,smart_df(VAR="Ref_rev",DEF="Reverse read counts harboring reference base"))
	readme_VC = rbind(readme_VC,smart_df(VAR="Alt_for",DEF="Forward read counts harboring alternate base"))
	readme_VC = rbind(readme_VC,smart_df(VAR="Alt_rev",DEF="Reverse read counts harboring alternate base"))
	readme_VC = rbind(readme_VC,smart_df(VAR="Strand_bias_P",DEF="Fisher Exact Test pvalue of forward/reverse and reference/alternate read counts"))
	readme_VC = rbind(readme_VC,smart_df(VAR="LABEL",DEF=smart_sprintf("Characterizing each variant based on available UNMASC columns. 
		'H2M' = variant located in H2M region. 
		'FFPE'/'OXOG'/'ARTI' = variant's allele frequency classifies to artifact cluster.
		'intron/syn' = variant's highest negative impact is intronic or synonymous.
		'cosmic>=10' = called in COSMIC at least ten times.
		'modify' = snpEff's impact is MODIFIER.
		'strand' = read counts harbor strand bias.
		'BAF-like' = variant's allele frequency classifies to BAF cluster, may be private germline or somatic mutation in high purity sample.
		'no_label' = no COSMIC, artifact (H2M,FFPE,OXOG,ARTI), low impact, or BAF flags.")))
	
	if( is.null(headers) ){
		out = readme_VC
		print(readme_VC,right=FALSE)
	} else if( !is.null(headers) && all(headers %in% readme_VC$VAR) ){
		out = readme_VC[which(readme_VC$VAR %in% headers),]
		print(out,right=FALSE)
	} else {
		miss_headers = headers[which(!(headers %in% readme_VC$VAR))]
		miss_headers = paste0("'",miss_headers,"'")
		miss_headers = paste(miss_headers,collapse=", ")
		stop(sprintf("Headers %s \n don't have a definition or hasn't been defined.\n",miss_headers))
	}
	
	return(NULL)
}

make_uniq_vc_anno = function(all_vcs,strand,all_nVAF_segs,all_tVAF_segs){
	# Combines strand bias, segmentation, FFPE, OXOG, variant call frequency, prop_nLabel
	stop("Function not done yet!!")
	uniq_vcs = smart_uniq_df(input_df = all_vcs,vars = "mutID")
	uniq_vcs = anno_SEG(VC = uniq_vcs,nSEG = all_nVAF_segs,tSEG = all_tVAF_segs)
	tmp_vars = c("mutID","N_eps","N_prob_half","N_mBAF","N_prop_mBAF","N_ms","nANNO",
		"T_eps","T_prob_half","T_mBAF","T_prop_mBAF","T_ms","ALLELE_STAT","tANNO")
	uniq_vcs = smart_merge(uniq_vcs,strand,all.x = TRUE)
}

# New functions for UNMASC

#' @title run_UNMASC
#' @description This function implements the UNMASC workflow from 
#'		processing formatted dataframes of variants toward outputting annotated variants.
#' 
#' @param tumorID A character string unique to a tumor's output 
#'	directory and files.
#' @param outdir A character string of the full path and working directory specifying 
#'	where intermediate and final files are stored.
#' @param vcf If \code{run_UNMASC()} has been run to create an \code{image.rds}
#'	file, this can be set to \code{NULL}. Otherwise, it is a dataframe
#'	containing "Chr" (e.g. 'chr1'), "Position" (e.g. 1000), "Ref" (e.g. "A"), 
#'	"Alt" (e.g. "T"), "Qscore" (e.g. 30), "GeneName" (e.g. 'TP53'),
#'	"TARGET" (e.g. 'YES', 'NO'), "EXONIC" (e.g. 'YES', 'NO'), 
#'	"IMPACT" (e.g. 'MODIFIER','LOW','MODERATE','HIGH'),
#'	"STUDYNUMBER" (e.g. "normal_1"), "CosmicOverlaps" (e.g. 10), "ThsdG_AF" (e.g. 0.5), 
#'	"EXAC_AF" (e.g. 0.1), "nAD" (e.g. 10), "nRD" (e.g. 20), 
#'	"tAD" (e.g. 10), "tRD" (e.g. 20) columns corresponding to 
#'	chromosome, position, reference allele, alternate allele,
#'	variant quality score, Hugo symbol, on/off target status,
#'	exonic/intronic status, snpEff impact, normal control ID, number of COSMIC overlaps,
#'	1000 Genomes population allele frequency, ExAC population allele frequency,
#'	control alternate depth, control reference depth, tumor alternate depth,
#'	tumor reference depth, respectively.
#' @param tBAM_fn A character string specifying the full path 
#'	to the tumor's BAM file. Can be set to \code{NULL} if strand.rds already exists.
#' @param bed_centromere_fn Centromere regions filename. This should
#'	be tab delimited without headers containing columns contig (e.g. 'chr1'),
#'	start position (e.g. 100), and end position (e.g. 100000).
#' @param dict_chrom_fn Chromosome lengths file. This can be
#'	constructed from the output of 'samtools view -H' applied
#'	to a bam file.
#' @param qscore_thres A numeric value specifying a minimum 
#'	allowed Qscore or QUAL value for variant calls.
#' @param exac_thres A numeric value specifying a maximum 
#'	allowed ExAC allele frequency for variant calls.
#' @param ad_thres A numeric value specifying a minimum
#'	allowed number of alternate read counts for variant calls.
#' @param rd_thres A numeric value specifying a minimum
#'	allowed total read depth for variant calls.
#' @param cut_BAF A numeric value specifying the BAF threshold 
#'	to remove variants before running tumor VAF segmentation.
#'	By default, \code{cut_BAF} is set to 0.05.
#' @param minBQ An integer for minimum base quality.
#' @param minMQ An integer for minimum mapping quality.
#' @param eps_thres A numeric value specifying the threshold 
#'	for determining a H2M segment.
#' @param psi_thres A numeric value specifying the threshold 
#'	for determining a H2M segment.
#' @param hg A character string for the human genome. This is 
#'	used for labeling output plots.
#' @param binom Boolean for whether or not to model read counts as 
#'	binomial or beta-binomial distributed.
#' @param gender A single string specifying the subject's gender. 
#'	Valid inputs are \code{"MALE"} and \code{"FEMALE"}.
#'	By default tumor variant read counts on chromosomes 1 thru 
#'	22 are segmented. If \code{gender = "FEMALE"}, chromosome X 
#'	is also segmented.
#' @param flag_samp_depth_thres A vector of two integers thresholds.
#'	The first is a count of unique loci achieving a higher total tumor 
#'	read depth than the second integer specified. For example, if
#'	\code{flag_samp_depth_thres = c(1000,50)}, then if less than 1000
#'	unique loci in the tumor have total depths greater than 50, the sample
#'	is flagged and UNMASC exits.
#' @param ncores A positive integer for the number of threads to 
#'	use for calculating strand-specific read counts.
#' @return Null from function. Outputs UNMASC results to files.
#' @export
run_UNMASC = function(tumorID,outdir,vcf = NULL,tBAM_fn = NULL,
	bed_centromere_fn,dict_chrom_fn,qscore_thres = 30,exac_thres = 5e-3,
	ad_thres = 5,rd_thres = 10,cut_BAF = 5e-2,minBQ = 13,minMQ = 40,
	eps_thres = 0.5,psi_thres = 0.02,hg = "19",binom = TRUE,
	gender = NA,flag_samp_depth_thres = c(1e3,50),ncores = 1){
	
	if(FALSE){
		tumorID = tumorID; outdir = outdir; vcf = vcf
		tBAM_fn = tBAM_fn; bed_centromere_fn = my_dirs$bed_cent_fn
		dict_chrom_fn = my_dirs$dict_chrom_fn; binom = binom; gender = gender
		ncores = ncores
		
		qscore_thres = 30; exac_thres = 5e-3; ad_thres = 3; 
		rd_thres = 10; cut_BAF = 5e-2; minBQ = 13; minMQ = 40; 
		eps_thres = 0.5; psi_thres = 0.02; hg = "19"; 
		flag_samp_depth_thres = c(1e3,50); ncores = 1
		
	}
	
	message("% ------------------------------- %\n",appendLF = FALSE)
	message("% Welcome to the UNMASC workflow! %\n",appendLF = FALSE)
	message("% ------------------------------- %\n",appendLF = FALSE)
	smart_mkdir(outdir)
	
	image_fn = file.path(outdir,"image.rds")
	strand_fn = file.path(outdir,"strand.rds")
	
	if( !file.exists(image_fn) ){
		# Checking required inputs
		req_COLS = c("Chr","Position","Ref","Alt","Qscore",
			"GeneName","TARGET","EXONIC","IMPACT","STUDYNUMBER",
			"CosmicOverlaps","ThsdG_AF","EXAC_AF",
			"nAD","nRD","tAD","tRD")
		smart_header(OBJ = vcf,req_COLS = req_COLS)
		if( !all(unique(vcf$TARGET) %in% c("YES","NO")) ) stop("Check TARGET")
		if( !all(unique(vcf$EXONIC) %in% c("YES","NO")) ) stop("Check EXONIC")
		if( !( is.na(gender) || gender %in% c("FEMALE","MALE") ) ){
			stop("gender should be NA, 'FEMALE', or 'MALE'!")
		}
		
		# Filter based on depth and contig, create mutID
		message(sprintf("%s: Calculate mutID and light filtering ...\n",date()),appendLF = FALSE)
		vcf$nDP = rowSums(vcf[,c("nAD","nRD")])
		vcf$tDP = rowSums(vcf[,c("tAD","tRD")])
		vcf = vcf[which(vcf$tDP >= 5 & vcf$nDP >= 5 & vcf$Qscore >= 5),] # lowest tolerated thresholds
		
		# Flagging low QC samples: Method 1
		if( nrow(vcf) < 1e3 ){
			message(sprintf("%s: LowQCSample b/c low variant count after base filtering ...\n",date()),appendLF = FALSE)
			return(NULL)
		}
		
		vcf = vcf[which(vcf$Chr %in% paste0("chr",c(1:22,"X","Y"))),]
		uvcf = unique(vcf[,c("Chr","Position","Ref","Alt")])
		uvcf$nChr = factor(uvcf$Chr,levels = paste0("chr",c(1:22,"X","Y")))
		uvcf = uvcf[order(uvcf$nChr,uvcf$Position),]
		uvcf = uvcf[,names(uvcf) != "nChr"]
		uvcf$mutID = seq(nrow(uvcf))
		vcf = smart_merge(vcf,uvcf)
		rm(uvcf)
		
		# Calculate VAFs
		message(sprintf("%s: Calculate nVAFs and tVAFs ...\n",date()),appendLF = FALSE)
		if( !("tVAF" %in% names(vcf)) ){
			vcf$tVAF = ifelse(vcf$tDP == 0,0,vcf$tAD/vcf$tDP)
		}
		if( !("nVAF" %in% names(vcf)) ){
			vcf$nVAF = ifelse(vcf$nDP == 0,0,vcf$nAD/vcf$nDP)
		}
		
		# Get chromosome lengths and centromere info
		message(sprintf("%s: Import genome files ...\n",date()),appendLF = FALSE)
		genome = prep_genome(bed_centromere_fn = bed_centromere_fn,
			dict_chrom_fn = dict_chrom_fn)
		
		# Create xCoor
		message(sprintf("%s: Calculate xCoor ...\n",date()),appendLF = FALSE)
		vcf$x_coor = create_xCoor(genome = genome,x = vcf)$x_coor
		
		# Cluster normal counts
		message(sprintf("%s: Normal read count clustering ...\n",date()),appendLF = FALSE)
		out_nCLUST = run_nCLUST(vcf = vcf,ID = "STUDYNUMBER",
			rd_thres = rd_thres,ad_thres = ad_thres,qscore_thres = qscore_thres,
			hg = hg,binom = binom,outdir = outdir,genome = genome,verbose = TRUE)
		vcf = out_nCLUST$out
		SE = out_nCLUST$SE
		
		# Strand
		if( !file.exists(strand_fn) ){
			if( is.null(tBAM_fn) ) stop("Specify tBAM_fn!")
			strand = new_STRAND(BAM_fn = tBAM_fn,vcs = vcf,
				minBQ = minBQ,minMQ = minMQ,NT = ncores)
			saveRDS(strand,strand_fn)
		} else {
			strand = readRDS(strand_fn)
		}
		
		# Save image
		message(sprintf("%s: Save image ...\n",date()),appendLF = FALSE)
		saveRDS(list(vcf = vcf,strand = strand,
			genome = genome,SE = SE,qscore_thres = qscore_thres,
			exac_thres = exac_thres,ad_thres = ad_thres,
			rd_thres = rd_thres,gender = gender,hg = hg),image_fn)
	}
	
	message(sprintf("%s: Import image ...\n",date()),appendLF = FALSE)
	rds 		= readRDS(image_fn)
	vcf 		= rds$vcf
	strand 	= rds$strand
	genome 	= rds$genome
	rm(rds)
	
	# Flagging low QC samples: Method 2
	tab_tumor  = table(CONTROL = vcf$STUDYNUMBER,Tumor_Depth = vcf$tDP)
	tab_normal = table(CONTROL = vcf$STUDYNUMBER,Normal_Depth = vcf$nDP)
	if( ncol(tab_tumor) <= 50 ){
		print(tab_tumor)
		message(sprintf("%s: LowQCSample b/c low variability in tumor depth\n",date()),appendLF = FALSE)
		return(NULL)
	}
	
	# Flagging low QC samples: Method 3
	num_uniq_IDs = length(unique(vcf$mutID[vcf$tDP >= flag_samp_depth_thres[2]]))
	if( num_uniq_IDs < flag_samp_depth_thres[1] ){
		message(sprintf("%s: LowQCSample b/c %s unique tumor loci, which is less than %s, achieve total read depth >= %s\n",date(),num_uniq_IDs,flag_samp_depth_thres[1],flag_samp_depth_thres[2]),appendLF = FALSE)
		return(NULL)
	}
	
	# nSEG
	nSEG_fn = file.path(outdir,"anno_nSEG.rds")
	if( !file.exists(nSEG_fn) ){
		UNMASC_nSEG(outdir = outdir,vcs = vcf,
			tumorID = tumorID,hg = hg,genome = genome,
			binom = binom,eps_thres = eps_thres,
			psi_thres = psi_thres)
	}
	nSEG = readRDS(nSEG_fn)$all_nVAF_segs
	
	# Create RefAlt/INDEL
	vcf = smart_RefAlt(vcs = vcf)
	
	# Artifact hunting
	# art_binom = TRUE
	arti_fn = file.path(outdir,"ARTI.rds")
	art_binom = binom
	arti_oxog = ARTI_find(vcs = vcf,strand = strand,
		nSEG = nSEG,qscore_thres = qscore_thres,rd_thres = rd_thres,
		ad_thres = ad_thres,binom = art_binom,ARTI = "OXOG",
		eps_thres = eps_thres,psi_thres = psi_thres)
		vcf = arti_oxog$vcs
	arti_ffpe = ARTI_find(vcs = vcf,strand = strand,
		nSEG = nSEG,qscore_thres = qscore_thres,rd_thres = rd_thres,
		ad_thres = ad_thres,binom = art_binom,ARTI = "FFPE",
		eps_thres = eps_thres,psi_thres = psi_thres)
		vcf = arti_ffpe$vcs
	arti_arti = ARTI_find(vcs = vcf,strand = strand,
		nSEG = nSEG,qscore_thres = qscore_thres,rd_thres = rd_thres,
		ad_thres = ad_thres,binom = art_binom,ARTI = "ARTI",
		eps_thres = eps_thres,psi_thres = psi_thres)
		vcf = arti_arti$vcs
	saveRDS(list(OXOG = arti_oxog$clust,FFPE = arti_ffpe$clust,
		ARTI = arti_arti$clust),arti_fn)
	rm(arti_oxog,arti_ffpe,arti_arti)
	
	# tSEG
	tSEG_fn = file.path(outdir,"anno_tSEG.rds")
	if( !file.exists(tSEG_fn) ){
		UNMASC_tSEG(outdir = outdir,tumorID = tumorID,
			vcs = vcf,strand = strand,gender = gender,
			cut_BAF = cut_BAF,hg = hg,genome = genome,
			rd_thres = rd_thres,ad_thres = ad_thres,
			qscore_thres = qscore_thres,binom = binom)
	}
	tSEG = readRDS(tSEG_fn)$all_tVAF_segs
	
	# Add segmentation annotations
	vcf_TE = vcf[which(vcf$TARGET == "YES" & vcf$EXONIC == "YES"),]
	uniq_vcs = smart_uniq_df(input_df = vcf_TE,vars = "mutID")
	message(paste0("nrow of uniq_vcs = ",nrow(uniq_vcs),"\n"),appendLF = FALSE)
	uniq_vcs2 = anno_SEG(VC = uniq_vcs,nSEG = nSEG,tSEG = tSEG,
		eps_thres = eps_thres,psi_thres = psi_thres)
	clust_names = c("eps","prob_half","mBAF","prop_mBAF","PSI","ms")
	tmp_vars = c("mutID",paste0("N_",clust_names),"nANNO",
		paste0("T_",clust_names),"ALLELE_STAT","tANNO")
	vcf_TE = smart_merge(vcf_TE,uniq_vcs2[,tmp_vars])
	
	# Add Strand Columns
	vcf_TE = smart_merge(vcf_TE,strand,all.x = TRUE)
	
	# Filter
	fNPN = TO_pipe_add_prop_nlabel_ART(npn = vcf_TE,
		rd_thres = rd_thres,ad_thres = ad_thres,
		qscore_thres = qscore_thres,exac_thres = exac_thres,
		prop_call = 0.9,eps_thres = eps_thres,psi_thres = psi_thres,
		noise_filter = TRUE,segANNO_filter = FALSE)
	
	# Plot: Germline with candidate somatic variants, label indel and canonical
	plot_TO(outdir = outdir,tumorID = tumorID,all_vcs = vcf,
		fNPN = fNPN,genome = genome,hg = hg,all_tVAF_segs = tSEG)
	
	# Label variants
	fNPN$LABEL = ""
	cosmic_label 			= fNPN$CosmicOverlaps >= 10
	intron_syn_label	= rep(FALSE,nrow(fNPN))
	# intron_syn_label[which(grepl("INTRON",fNPN$snpEff_pt1))] = TRUE
	# intron_syn_label[which(grepl("SYNONYMOUS",fNPN$snpEff_pt1))] = TRUE
	# table(intron_syn_label)
	# intron_syn_label 	= fNPN$O_Effect2 %in% c("INTRON","SYNONYMOUS_CODING")
	mod_label					= fNPN$IMPACT %in% c("MODIFIER")
	h2m_label 				= !is.na(fNPN$N_eps) & ( fNPN$N_eps >= eps_thres | fNPN$N_PSI >= psi_thres )
	strand_label 			= !is.na(fNPN$Strand_bias_P) & fNPN$Strand_bias_P < 0.05
	germ_label 				= !is.na(fNPN$tANNO) & fNPN$tANNO < 0.5
	ffpe_label				= !sapply(fNPN$FFPE,artifact_filter,USE.NAMES = !TRUE) 
	oxog_label				= !sapply(fNPN$OXOG,artifact_filter,USE.NAMES = !TRUE)
	arti_label				= !sapply(fNPN$ARTI,artifact_filter,USE.NAMES = !TRUE)
	
	fNPN$LABEL[cosmic_label] 			= paste0(fNPN$LABEL[cosmic_label],"cosmic>=10;")
	fNPN$LABEL[intron_syn_label] 	= paste0(fNPN$LABEL[intron_syn_label],"intron/syn;")
	fNPN$LABEL[mod_label] 				= paste0(fNPN$LABEL[mod_label],"modify;")
	fNPN$LABEL[h2m_label] 				= paste0(fNPN$LABEL[h2m_label],"H2M;")
	fNPN$LABEL[strand_label] 			= paste0(fNPN$LABEL[strand_label],"strand;")
	fNPN$LABEL[germ_label] 				= paste0(fNPN$LABEL[germ_label],"BAF-like;")
	fNPN$LABEL[ffpe_label] 				= paste0(fNPN$LABEL[ffpe_label],"FFPE;")
	fNPN$LABEL[oxog_label] 				= paste0(fNPN$LABEL[oxog_label],"OXOG;")
	fNPN$LABEL[arti_label] 				= paste0(fNPN$LABEL[arti_label],"ARTI;")
	fNPN$LABEL[fNPN$LABEL==""] 		= "no_label"
	smart_table(fNPN$LABEL)
	
	# Output inferred genotypes
	TO_genotype(vcf = vcf,outdir = outdir,hg = hg,genome = genome,
		strand = strand,rd_thres = rd_thres,qscore_thres = qscore_thres,
		eps_thres = eps_thres,psi_thres = psi_thres)
	
	# Output variant calls
	smart_WT(fNPN,file.path(outdir,"tumorOnly_VCs.tsv"),sep = '\t')
	
	return(NULL)
}


###


