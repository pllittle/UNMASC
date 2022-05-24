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

# Mapping Quality functions
get_pileup = function(BAM_fn,CHR,POS,REF,ALT,minBQ,FASTA_fn){
	indel_status = ifelse(nchar(REF) == 1 & nchar(ALT) == 1,"BS",
		ifelse(nchar(REF) > 1,"DEL","INS"))
	out = system(smart_sprintf("samtools view -h %s '%s:%s-%s' \
		| samtools mpileup -Q %s -d 8000 -f %s - \
		| awk '($2==%s)'",
		BAM_fn,CHR,POS,POS+1,minBQ,FASTA_fn,POS),intern = TRUE)
	out = strsplit(out,"\t")[[1]]
	out[5] = toupper(out[5])
	if( indel_status == "INS" ){
		indel_seq = paste(strsplit(ALT,"")[[1]][-1],collapse="")
		indel_length = nchar(indel_seq)
		out[5] = gsub(paste0("\\+",indel_length,indel_seq),"D",out[5])
	} else if( indel_status == "DEL" ){
		indel_seq = paste(strsplit(REF,"")[[1]][-1],collapse="")
		indel_length = nchar(indel_seq)
		out[5] = gsub(paste0("\\-",indel_length,indel_seq),"D",out[5])
	}
	out
}
test_BQ = function(BAM_fn,CHR,POS,REF,ALT,minBQ,FASTA_fn,PERMS=1e3){
	if(FALSE){
		minBQ = 0; PERMS = 1e3
		ii = 31
		CHR = vcs$Chr[ii]; POS = vcs$Position[ii]; REF = vcs$Ref[ii]; ALT = vcs$Alt[ii];
		sprintf("%s:%s:%s:%s",CHR,POS,REF,ALT)
		
	}
	
	chrposrefalt = sprintf("%s:%s:%s:%s",CHR,POS,REF,ALT)
	indel_status = ifelse(nchar(REF) == 1 & nchar(ALT) == 1,"BS",
		ifelse(nchar(REF) > 1,"DEL","INS"))
	
	if( indel_status %in% c("DEL","INS") ){
		return(list(asymp_pvalue = NA,perm_pvalue = NA,exact_pvalue = NA))
	}
	
	aa = get_pileup(BAM_fn = BAM_fn,CHR = CHR,POS = POS,
		REF = REF,ALT = ALT,minBQ = minBQ,FASTA_fn = FASTA_fn)
	if(FALSE){
		aa = rep(NA,6)
		# aa[5] = "...,,,+3ACT,.,."
		# aa[5] = "...,,,-3ACT,.,.,"
	}
	if( length(grep("\\*",aa[5])) > 0 ) aa[5] = gsub("\\^\\*","",aa[5])
	if( length(grep("\\,",aa[5])) > 0 ) aa[5] = gsub("\\^\\,","",aa[5])
	if( length(grep("\\.",aa[5])) > 0 ) aa[5] = gsub("\\^\\.","",aa[5])
	if( length(grep("N",aa[5])) > 0 ) aa[5] = gsub("\\^\\N","",aa[5])
	if( length(grep("+",aa[5])) > 0 ) aa[5] = gsub("\\^\\+","",aa[5])
	if( length(grep("-",aa[5])) > 0 ) aa[5] = gsub("\\^\\-","",aa[5])
	
	if( length(grep("\\^A",aa[5])) > 0 ) aa[5] = gsub("\\^\\A","",aa[5])
	if( length(grep("\\^C",aa[5])) > 0 ) aa[5] = gsub("\\^\\C","",aa[5])
	if( length(grep("\\^G",aa[5])) > 0 ) aa[5] = gsub("\\^\\G","",aa[5])
	if( length(grep("\\^T",aa[5])) > 0 ) aa[5] = gsub("\\^\\T","",aa[5])
	
	bases = c("A","C","G","T")
	
	# If indel exists, replace with single letter "I"
	ss = strsplit(toupper(aa[5]),"")[[1]]
	indel_idx = which(ss %in% c("+","-")); indel_idx
	if( length(indel_idx) > 0 ){
		nums = seq(0,9)
		for(ii in indel_idx){
			# ii = indel_idx[1]
			jj = 1
			while(TRUE){
				if( ss[ii+jj] %in% nums ){
					ss[ii] = paste0(ss[ii],ss[ii+jj])
					ss[ii+jj] = ""
					jj = jj + 1
				} else {
					# Get indel sequence
					indel_status = ifelse(grepl("-",ss[ii]),"-","+"); indel_status
					len_indel = gsub("\\+","",ss[ii]); len_indel = gsub("\\-","",len_indel)
					len_indel = as.numeric(len_indel); len_indel
					indel_seq = paste0(c(indel_status,len_indel,ss[seq(ii+jj,ii+jj+len_indel-1)]),collapse="")
					aa[5] = gsub(indel_seq,"I",aa[5])
					break
				}
			}
		}
		# ss = ss[which(!is.na(ss))]
	}
	
	ss = strsplit(toupper(aa[5]),"")[[1]]
	ss = ss[which(ss != REF)]
	# Notes:
	# 	. base matches reference forward
	# 	, bases matches reference reverse
	# 	$ base at the end of read
	# 	^] base was at the start of read
	
	## how to handle *??
	# ss = ss[which(ss %in% c("*",".",",","D",0,"N",bases))]
	# ss = ss[which(ss %in% c("*",".",",","D","N",bases))] # current code
	ss = ss[which(ss %in% c("*",".",",","N",bases))] # newest code
	ss[which(ss %in% c("*",".",","))] = REF
	# ss = ss[which(ss %in% c(".",",","D",0,"N",bases))]
	# ss[which(ss %in% c(".",","))] = REF
	
	ascii_qual = char2ascii(c = aa[6])
	length(ss); length(ascii_qual)
	
	if( length(ss) > length(ascii_qual) && any("N" %in% ss) ){
		ss = ss[which(ss != "N")]
	}
	if( length(ss) != length(ascii_qual) ){
		if( length(ss) > length(ascii_qual) ){
			tmp_tab = table(ss)
			tmp_diff = length(ss) - length(ascii_qual)
			if( any(tmp_tab == tmp_diff) ){
				extra_char = names(which(tmp_tab == tmp_diff))
				ss = ss[ss != extra_char]
			} else {
				stop(sprintf("Unknown MQ issue: %s",chrposrefalt))
			}
		} else {
			stop(sprintf("Error in number of bases and number of BQ scores: %s",chrposrefalt))
		}
	}
	
	out = smart_df(BASE = ss,BQ = char2ascii(c = aa[6]) - 33)
	out$RefAlt = ifelse(out$BASE == REF,"Ref","Alt")
	
	# Binned mapping quality
	out$QQ = NA
	for(th in rev(seq(0,30,10))){
		out$QQ[which(out$BQ >= th & is.na(out$QQ))] = th + 10
	}
	out = out[which(out$BASE != "0"),]
	# table(out$QQ,useNA="ifany")
	# out[1:5,,drop=FALSE]
	if( length(unique(out$RefAlt)) > 1 ){
		unperm_tab = table(out$RefAlt,out$QQ); # tab
		if( ncol(unperm_tab) == 1 ){
			return(list(samtools_out = aa,tab1 = table(out$BASE,out$QQ),
				tab2 = table(out$RefAlt,out$QQ),asymp_pvalue = NA,
				perm_pvalue = NA,exact_pvalue = NA))
		}
		unperm_res = suppressWarnings(chisq.test(unperm_tab))
		
		# Permutation result: Permute the RefAlt column
		pstat = rep(NA,PERMS+1); nn = nrow(out)
		pstat[1] = unperm_res$statistic
		for(pp in seq(PERMS)){
			pstat[pp+1] = suppressWarnings(chisq.test(
				table(sample(out$RefAlt,nn,replace=FALSE),out$QQ)))$statistic
		}
		perm_pvalue = mean(pstat >= pstat[1])
		exact_pvalue = tryCatch({
			fisher.test(unperm_tab)$p.value
			}, error = function(error_condition){
				NA
			})
		
		list(samtools_out = aa,tab1 = table(out$BASE,out$QQ),
			tab2 = unperm_tab,asymp_pvalue = unperm_res$p.value,
			perm_pvalue = perm_pvalue,exact_pvalue = exact_pvalue)
	} else {
		list(samtools_out = aa,tab1 = table(out$BASE,out$QQ),
			tab2 = table(out$RefAlt,out$QQ),asymp_pvalue = NA,
			perm_pvalue = NA,exact_pvalue = NA)
	}
}
test_parallel_MQ = function(vcs,BAM_fn,FASTA_fn,minBQ=0,PERMS=1e3,NT = 1){
	if(FALSE){
		vcs = vcs[14:15,]
		#BAM_fn = file.path(freeze_dir,samp,"output","tumor.abra.bam")
		#FASTA_fn = FASTA_fn
		minBQ = 0; PERMS = 1e3; NT = 5
	}
	
	if( !file.exists(BAM_fn) ){
		stop(sprintf("%s file does not exist!",BAM_fn))
	}
	
	cat(sprintf("Run MQ feature: Sample has %s variants ...\n",nrow(vcs)))
	cl = makeCluster(NT) # NT = num threads
	registerDoParallel(cl)
	ii = NULL
	xx = foreach(ii=seq(nrow(vcs)),.combine="rbind",
		.export = c("smart_df","smart_sprintf","get_pileup","test_BQ")) %dopar% {
	
		# ii = 2
		out1 = test_BQ(BAM_fn = BAM_fn,CHR = vcs$Chr[ii],POS = vcs$Position[ii],
			REF = vcs$Ref[ii],ALT = vcs$Alt[ii],minBQ = minBQ,FASTA_fn = FASTA_fn,PERMS=PERMS)
		c(out1$asymp_pvalue,out1$perm_pvalue,out1$exact_pvalue)
	}
	stopCluster(cl)
	
	rownames(xx) = NULL
	xx = smart_df(xx)
	names(xx) = c("MQ_asymp","MQ_perm","MQ_exact")
	vcs = smart_df(vcs,xx)
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
		if( ii %% 2 == 0 ) cat(".")
		if( ii %% 20 == 0 || ii == TRIALS ) cat(sprintf("%s out of %s\n",ii,TRIALS))
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
#'
#' @export
run_hsMetric_clustRank = function(DATA,VARS,ID,TRIALS=50,note=""){
	if( class(DATA) != "data.frame" ) stop("DATA should be a data.frame")
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
	cat("Top ranked samples ...\n")
	print(head(out))
	cat("Low ranked samples ...\n")
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
		cat(sprintf("%s: Finding oxoG artifacts ...\n",date()))
	} else if( ARTI == "FFPE" ){
		cat(sprintf("%s: Finding FFPE artifacts ...\n",date()))
	} else if( ARTI == "ARTI" ){
		cat(sprintf("%s: Finding ARTI artifacts ...\n",date()))
	}
	vcs = vcs[,names(vcs) != ARTI]
	vcs$INFER = NA
	
	# Append strand bias column
	cat(sprintf("%s: Merge strand info ...\n",date()))
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
		cat(paste0(id," "))
		if( cnt %% 3 == 0 ) cat("\n")
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
	cat("\n")
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
	
	cat(sprintf("%s: Finding FFPE artifacts ...\n",date()))
	vcs = vcs[,names(vcs) != "FFPE"]
	vcs$INFER = NA
	
	# Append strand bias column
	cat(sprintf("%s: Merge strand info ...\n",date()))
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
		cat(paste0(id," "))
		if(cnt %% 5 == 0) cat("\n")
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
	cat("\n")
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
	if(verbose) cat("Getting target info = ")
	input_df$TARGET = NA
	for(one_contig in sort(unique(input_df$Chr))){
		# one_contig = unique(input_df$Chr)[1]
		if(verbose) cat(paste0(gsub("chr","",one_contig)," "))
		sub_target_bed = target_bed[which(target_bed$chrom == one_contig),]
		one_index = which(input_df$Chr == one_contig)
		input_df$TARGET[one_index] = Rcpp_get_target(intervals = as.matrix(sub_target_bed[,c("start","end")]),
			positions = input_df$Position[one_index],pad = pad)
	}
	if(verbose) cat("\n")
	
	input_df
}
import_bed = function(bed_fn){
	cat(sprintf("Import bed file: %s\n",bed_fn))
	bed_df = fread(file = bed_fn,
		sep = '\t',header = TRUE,
		data.table = FALSE,
		showProgress = FALSE)
	bed_df = bed_df[,c("chrom","start","end")]
	names(bed_df) = c("chrom","start","end")
	bed_df
}

# Segmentation functions
BAF_EM_VAF = function(MAT,binom=TRUE,show_dots=TRUE,clust_tVAF=NULL){
	if(FALSE){
		MAT = as.matrix(sub_data[,c("AD","RD","DP","LBC","log_DP")])
		binom = binom
		show_dots = TRUE
		clust_tVAF = FALSE
	}
	
	# Constants
	init_pi = matrix(c(
		1,0,0,0,	# pure noise
		0,1,0,0,	# pure VAF=0.5
		1,1,0,0,	# mix of noise and VAF=0.5
		1,0,1,1,	# mix of noise and VAF=mm and VAF=1-mm
		0,1,1,1,	# mix of VAF=0.5 and VAF=mm and VAF=1-mm
		1,1,1,1),	# mix of noise and VAF=0.5 and VAF=mm and VAF=1-mm
		ncol=4,byrow=TRUE)
	
	my_max_iter = 3e2
	eps = 1e-3
	
	vec_pp = c(1/20,seq(4)/10)
	if( !is.null(clust_tVAF) ){
		if( clust_tVAF == FALSE ){
			vec_pp = 0.25
			init_pi = init_pi[1:3,,drop=FALSE]
		}
	}
	init_pi = init_pi / rowSums(init_pi)

	# Optimize
	
	# iPSI = ifelse(binom == TRUE,0,0.005)
	iPSI_vec = c(0,0.005)
	if( binom == TRUE ) iPSI_vec = 0
	
	em_out = NULL
	
	for(iPSI in iPSI_vec){
	for(jj in seq(nrow(init_pi))){
	for(mm0 in vec_pp){
		tmp_em = Rcpp_BAF_clust(RD = MAT[,c("AD","RD")],DP = MAT[,"DP"],
			log_DP = MAT[,"log_DP"],LBC = MAT[,"LBC"],props0 = init_pi[jj,],
			mm0 = mm0,psi0 = iPSI,max_iter = my_max_iter,eps = eps,show = FALSE)
		if( is.null(em_out) || ( tmp_em$converge == "yes" && tmp_em$BIC > em_out$BIC ) ){
			# If clust_tVAF, there should be at least 2 distributions mixed together
			if( is.null(clust_tVAF) ){
				em_out = tmp_em
			} else if( clust_tVAF == TRUE ){
				noise_dist = 1 * ( tmp_em$props[1] > 0 )
				half_dist = 1 * ( tmp_em$props[2] > 0 )
				mm_dist = 1 * ( sum(tmp_em$props[3:4]) > 0 )
				num_dist = noise_dist + half_dist + mm_dist
				if( noise_dist == 1 || num_dist >= 2 ){
					em_out = tmp_em
				}
			} else {
				em_out = tmp_em
			}
		}
		rm(tmp_em)
	}}}
	
	em_out
}

#' @title run_AF_SEG
#' @description This function performs normal and tumor variant
#'	allele frequency segmentation. Rows with allele frequencies 
#'	near 0 and 1 should be excluded before segmentation.
#' @param DATA a R data.frame containing columns \code{RD}, 
#'	\code{AD}, and \code{index} corresponding to reference depth,
#'	alternate depth, and the ordering and grouping of \code{DATA}'s 
#'	rows.
#' @param num_divides A positive integer to specify the initial 
#'	number of equally spaced segments prior to segmentation. This is
#'	used if \code{init_seg} is set to \code{NULL}.
#' @param init_seg A vector specifying the user's initial segmentation.
#' @param binom Boolean value set to \code{TRUE} to model read counts
#'	with binomial distribution. Otherwise model read counts with
#'	beta-binomial distribution.
#' @param clust_tVAF If \code{TRUE}, \code{DATA} is assumed to 
#'	originate from tumor read counts. Otherwise, \code{DATA} is 
#'	assumed to originate from normal read counts.
#' @param show_plot If \code{TRUE}, the function will generate 
#'	a sequence of plots to visualize the points and ongoing segmentation.
#' @param show_mess If \code{TRUE}, the function will output detailed 
#'	segmentation messages. By default \code{show_mess = FALSE}.
#' @param min_adj_seg An integer for the minimum number of indices
#'	to constitute a segment.
#' @param min_nUniqIdx_pSeg An integer for the minimum number of unique
#'	indices per segment.
#' @export
run_AF_SEG = function(DATA,num_divides,init_seg = NULL,binom = TRUE,
	clust_tVAF = TRUE,show_plot = FALSE,show_mess = FALSE,min_adj_seg = 20,
	min_nUniqIdx_pSeg = 10){
	
	if(FALSE){
		DATA = input_data
		num_divides = 40
		init_seg = NULL
		binom = binom
		clust_tVAF = !TRUE
		show_plot = TRUE
		show_mess = show_mess
	}
	
	# Functions
	# Plan: Write mini-functions for parts of full segmentation
	get_min_adjacent = function(SEG){
		# SEG = c(1,10,55,100)
		dat = c()
		begin = SEG[1]
		for(ii in SEG[-1]){
			finish = ii
			dat = rbind(dat,smart_df(start = begin,end = finish))
			begin = finish + 1
		}
		dat$length = dat$end - dat$start
		min(dat$length)
	}
	get_START_END = function(SEG,SEG_INT){
		# SEG_INT = takes values between 1 and length(SEG) - 1
		
		# Get START and END
		END = SEG[SEG_INT + 1]
		if(SEG_INT == 1){
			START = SEG[SEG_INT]
		} else {
			START = SEG[SEG_INT] + 1
		}
		# print(paste0("START(",START,") <= INDEX <= END(",END,")"))
		c(START,END)
	}
	get_new_seg_BIC = function(DATA,SEG,input_HIST){
		if(FALSE){
			# DATA = input_data
			# SEG = SEG + kk*jj*vec_shift
			# SEG = c(1,14,24,283,308,334,385,411,487,513)
			SEG = curr_seg
			# SEG = new_SEG
			# SEG = curr_seg[-ii]
			# SEG = curr_seg + pos_neg*increment*e_vec
			input_HIST = hist_LL
		}
		
		start_pos = SEG[1]
		curr_LL = 0
		curr_k = 0
		curr_BIC = 0
		for(ii in SEG[-1]){
			# ii = SEG[-1][1]
			end_pos = ii
			# print(c(start_pos,end_pos))

			# Code for BAF
			tmp_index = which(input_HIST$start == start_pos & input_HIST$end == end_pos)
			tmp_idx = which(DATA$index >= start_pos & DATA$index <= end_pos)
			if( is.null(input_HIST) || length(tmp_index) == 0 ){
				sub_data = DATA[tmp_idx,,drop = FALSE]
				# rownames(sub_data) = NULL
				tmp_int = curr_BAF_EM_VAF(MAT = as.matrix(sub_data[,
					c("AD","RD","DP","LBC","log_DP"),drop = FALSE]))
				
				input_HIST = rbind(input_HIST,smart_df(start = start_pos,
					end = end_pos,LL = tmp_int$LL,ms = tmp_int$ms))
				
			} else if( length(tmp_index) == 1 ){
				tmp_int = list()
				tmp_int$LL = input_HIST$LL[tmp_index]
				tmp_int$ms = input_HIST$ms[tmp_index]
			}
			curr_LL = curr_LL + tmp_int$LL
			curr_k = curr_k + tmp_int$ms
			curr_npts = length(tmp_idx)
			#curr_BIC = curr_BIC + 
			#	2 * tmp_int$LL - log(curr_npts) * tmp_int$ms
			start_pos = end_pos + 1
		}
		
		BIC = 2 * curr_LL - log(nrow(DATA)) * curr_k
		# BIC = curr_BIC
		list(BIC = round(BIC,2),HIST = input_HIST)
	}

	find_new_DIV = function(DATA,SEG,SEG_INT,input_HIST){
		if(FALSE){
			DATA = input_data
			SEG = curr_seg
			# SEG_INT = sample(x = c(2,seq(2,length(SEG))),size = 1) 
			SEG_INT = global_SEG_INT
			input_HIST = hist_LL
		}

		SE = get_START_END(SEG,SEG_INT)

		# Check num unique indices
		thres = 9 #9
		min_num_idx_per_seg = min_adj_seg # minimum end - start index among segmentation
		num_uniq_idx_per_seg = min_nUniqIdx_pSeg
		uniq_indices = unique(DATA$index[which(DATA$index >= SE[1] + 1 & DATA$index <= SE[2] - 1)])
		SEG_BIC = get_new_seg_BIC(DATA,SEG,input_HIST)$BIC
		
		if( length(uniq_indices) >= num_uniq_idx_per_seg ){
			new_INDICES = unique(as.numeric(round(quantile(uniq_indices,seq(thres-1)/thres))))
			vec_BIC = c()
			for(ii in new_INDICES){
				# ii = new_INDICES[2]
				new_SEG = unique(sort(c(SEG,ii)))
				min_adjacent = get_min_adjacent(new_SEG)
				if( min_adjacent >= min_num_idx_per_seg ){
					tmp_list = get_new_seg_BIC(DATA,new_SEG,input_HIST)
					vec_BIC = c(vec_BIC,tmp_list$BIC)
					input_HIST = tmp_list$HIST
				} else {
					vec_BIC = c(vec_BIC,SEG_BIC)
				}
			}
			max_index = which.max(vec_BIC)
			out_SEG = SEG
			if( vec_BIC[max_index] > SEG_BIC){
				out_SEG = sort(c(SEG,new_INDICES[max_index]))
			}
		} else {
			out_SEG = SEG
		}
		
		NEW = ifelse(length(SEG) < length(out_SEG),1,0)
		list(SEG=out_SEG,HIST=input_HIST,NEW=NEW)
	}
	find_new_DEL = function(DATA,SEG,input_HIST){

		curr_BIC = get_new_seg_BIC(DATA,SEG,input_HIST)$BIC
		seg_DEL = 0
		out_SEG_INT = NA
		for(ii in seq(2,length(SEG)-1)){
			tmp_list = get_new_seg_BIC(DATA,SEG[-ii],input_HIST)
			input_HIST = tmp_list$HIST
			if( tmp_list$BIC >= curr_BIC ){
				SEG = SEG[-ii]
				seg_DEL = 1
				out_SEG_INT = ii - 1
				break
			}
		}

		list(SEG=SEG,DEL=seg_DEL,HIST=input_HIST,SEG_INT=out_SEG_INT)
	}
	find_new_ADJ = function(DATA,SEG,input_HIST,vec_INTS,show=TRUE,show_try=!TRUE){
		if(FALSE){
			# DATA = input_data
			SEG = curr_seg
			input_HIST = hist_LL
			vec_INTS = vec_INTS
			show = TRUE
		}

		all_SHIFT = 2^seq(0,4)
		seg_ADJ = 0
		for(ii in vec_INTS){
			# ii = vec_INTS[1]
			vec_shift = rep(0,length(SEG))
			vec_shift[ii] = 1
			seg_ADJ_2 = 0
			while(TRUE){
				# Get max distance for adjustment
				seg_LEFT = SEG[ii]-SEG[ii-1]
				seg_RIGHT = SEG[ii+1]-SEG[ii]
				seg_SHIFT = min(abs(c(seg_LEFT,seg_RIGHT)))
				if( seg_SHIFT < 5 ) break
				all_SHIFT2 = rev(all_SHIFT[which(all_SHIFT < seg_SHIFT)])
				tmp_list = get_new_seg_BIC(DATA,SEG,input_HIST)
				tmp_BIC = tmp_list$BIC
				input_HIST = tmp_list$HIST
				for(jj in all_SHIFT2){
					# jj = all_SHIFT2[3]
					vec_BIC = c()
					for(kk in c(-1,1)){
						if(FALSE){
							kk = -1
							kk = 1 # something wrong here
							SEG + kk*jj*vec_shift
						}
						if(show_try) print(SEG + kk*jj*vec_shift)
						min_adjacent = get_min_adjacent(SEG + kk*jj*vec_shift)
						if( min_adjacent >= 5 ){
							tmp_list = get_new_seg_BIC(DATA,SEG + kk*jj*vec_shift,input_HIST)
							input_HIST = tmp_list$HIST
							vec_BIC = c(vec_BIC,tmp_list$BIC)
						} else {
							vec_BIC = c(vec_BIC,tmp_BIC)
						}
					}
					if( max(vec_BIC) > tmp_BIC ){
						tmp_BIC = max(vec_BIC)
						SEG = SEG + c(-1,1)[which.max(vec_BIC)]*jj*vec_shift
						# get_new_seg_BIC(DATA,SEG,input_HIST)$BIC
						if(show) cat("Adj found!")
						# if(show) print(SEG)
						seg_ADJ = 1
						seg_ADJ_2 = 1
						break
					}
				}
				if( seg_ADJ_2 == 0 || (seg_ADJ_2 == 1 && jj == 1) ) break
			}
		}

		list(SEG=SEG,seg_ADJ=seg_ADJ,HIST=input_HIST)
	}

	plot_SEG = function(DATA,SEG){
		plot(DATA[,c("index","VAF")],pch=16,cex=DATA$pt_cex,col=DATA$pt_col,ylim=c(0,1))
		abline(v=SEG,lwd=2,lty=2,col="red")
		abline(h=0.5,lwd=1,lty=2)
	}
	infer_SEG = function(DATA,SEG,show_plot=TRUE){
		if(FALSE){
			DATA = DATA
			SEG = curr_seg
			show_plot = TRUE
		}
		
		if(show_plot) plot_SEG(DATA,SEG)
		out_infer = c()
		START = SEG[1]
		seg = 1
		# eps_thres = ifelse(clust_tVAF,1,0.5)
		eps_thres = 1
		for(ii in seq(2,length(SEG))){
			# ii = 13
			END = SEG[ii]
			# print(c(START,END))
			my_MAT = as.matrix(DATA[which(DATA$index >= START & DATA$index <= END),c("AD","RD","DP","LBC","log_DP")])
			em_out = curr_BAF_EM_VAF(MAT = my_MAT)
			# em_out[c("converge","BIC","ms","params")]
			maf = NA
			if(em_out$params[1] < eps_thres){
				if( em_out$params[2] >= 0.5 ){
					maf = 0.5
				} else {
					maf = em_out$params[3]
					if( em_out$params[4] < 1 ){
						maf = c(maf,1-maf)
					}
				}
				if(show_plot) segments(x0=START,y0=maf,x1=END,col="magenta",lwd=4,lty=1)
			}
			out_infer = rbind(out_infer,smart_df(seg,START,END,t(em_out$params),ms=em_out$ms))
			START = END + 1
			seg = seg + 1
		}
		
		names(out_infer) = c("seg","start","end","eps","prob_half","mBAF","prop_mBAF","PSI","ms")
		out_infer
	}

	# 0) Check certain data.frame names exist in DATA
	req_cols = c("AD","RD","index")
	if( !all(req_cols %in% names(DATA)) ){
		miss_cols = req_cols[!(req_cols %in% names(DATA))]
		stop(paste0("DATA missing fields: ",paste(miss_cols,collapse=", ")))
	}
	
	# Make auxiliary variables (DP,LBC,VAF,pt_cex)
	DATA$DP = rowSums(DATA[,c("AD","RD")])
	DATA$log_DP = log(DATA$DP)
	DATA$LBC = lchoose(DATA$DP,DATA$AD)
	DATA$VAF = DATA$AD / DATA$DP
	DATA$pt_cex = 0.5*log10(DATA$DP)
	DATA$pt_col = rgb(0,0,0,0.5)
	
	# 1) Determine which VAFs are being clustered
	curr_BAF_EM_VAF = function(...){
		BAF_EM_VAF(...,binom = binom,show_dots = FALSE,clust_tVAF = clust_tVAF)
	}

	# 2) Initialization
	hist_LL = c()
	# DATA = truth_data$all_data
	if( is.null(init_seg) ){
		tmp_range = range(DATA$index)
		curr_seg = as.numeric(round(quantile(seq(tmp_range[1],tmp_range[2]),
			seq(0,num_divides)/num_divides)))
		if( get_min_adjacent(curr_seg) < 5 ){
			curr_seg = tmp_range
		}
	} else {
		curr_seg = init_seg
	}
	curr_list = get_new_seg_BIC(DATA,curr_seg,hist_LL)
	hist_LL = curr_list$HIST
	hist_BIC = curr_list$BIC

	# 3) Run Segmentation
	global_SEG_INT = 1; iter = 0; prop_thres = 0
	while( global_SEG_INT <= length(curr_seg) ){
		# if(iter %% 2 == 0) cat(paste0("Iter = ",iter,"\n"))
		
		if(show_plot){
			plot_SEG(DATA,curr_seg)
			abline(v=curr_seg[c(global_SEG_INT,global_SEG_INT+1)],col="darkgreen",lty=2,lwd=2)
		}
		if(show_mess) print(curr_seg)
		if( !is.na(curr_seg[global_SEG_INT+1]) && curr_seg[global_SEG_INT+1] < max(curr_seg) ){
			tmp_num = curr_seg[global_SEG_INT+1] / max(curr_seg)
			while(TRUE){
				if( tmp_num >= prop_thres ){
					cat(paste0(prop_thres*100,"% "))
					prop_thres = prop_thres + 0.05
				} else {
					break
				}
			}
		}
		
		# ADD SEGMENT
		seg_ADD = 0
		if(TRUE){
			if(show_mess) cat("Try ADD: ")
			list_DIV = find_new_DIV(DATA,curr_seg,global_SEG_INT,hist_LL)
			hist_LL = list_DIV$HIST
			curr_seg = list_DIV$SEG
			seg_ADD = list_DIV$NEW
			if( seg_ADD == 1 ){
				if(show_mess) cat("Added segment")
				global_SEG_INT = max(c(1,global_SEG_INT-1))
			}
			if(show_mess) cat("\n")
		}
		
		# REMOVE SEGMENT
		seg_DEL = 0
		if( seg_ADD == 0 && length(curr_seg) > 3 ){
			if(show_mess) cat("Try DEL: ")
			list_DEL = find_new_DEL(DATA,curr_seg,hist_LL)
			hist_LL = list_DEL$HIST
			curr_seg = list_DEL$SEG
			seg_DEL = list_DEL$DEL
			if(seg_DEL == 1){
				if(show_mess) cat("Removed segment")
				# global_SEG_INT = max(c(1,global_SEG_INT-2))
				global_SEG_INT = max(c(1,min(c(global_SEG_INT-2,list_DEL$SEG_INT-2))))
				# global_SEG_INT = 1
			}
			if(show_mess) cat("\n")
		}
		
		# ADJUST SEGMENTS
		seg_ADJ = 0
		if( seg_ADD == 0 && seg_DEL == 0 && length(curr_seg) > 3 ){
			if(show_mess) cat("Try ADJ: ")
			vec_INTS = seq(global_SEG_INT - 1,global_SEG_INT + 1)
			vec_INTS = vec_INTS[which(vec_INTS > 1 & vec_INTS < length(curr_seg))]
			list_ADJ = find_new_ADJ(DATA,curr_seg,hist_LL,vec_INTS,show=show_mess)
			hist_LL = list_ADJ$HIST
			curr_seg = list_ADJ$SEG
			seg_ADJ = list_ADJ$seg_ADJ
			if(seg_ADJ == 1){
				global_SEG_INT = max(c(1,global_SEG_INT-2))
			}
			if(show_mess) cat("\n")
		}
		
		# Move to next interval if no updates
		if( seg_ADD == 0 && seg_DEL == 0 && seg_ADJ == 0 ){
			global_SEG_INT = global_SEG_INT + 1
		}
		
		iter = iter + 1
	}
	cat("100%\n")
	
	out_INFER = infer_SEG(DATA = DATA,SEG = curr_seg,show_plot = show_plot)
	out_BIC = get_new_seg_BIC(DATA,curr_seg,hist_LL)$BIC

	list(INFER=out_INFER,BIC=out_BIC,SEG=curr_seg)
}
segment_nVAF = function(DATA,chr,binom=TRUE,eps_thres=0.5,psi_thres=0.01,show_plot=FALSE,show_mess=FALSE){
	if(FALSE){
		DATA = vcs; chr = chr; binom = binom
		show_plot = FALSE; show_mess = FALSE
	}
	
	input_data = DATA[which(DATA$Chr == paste0("chr",chr)
		& DATA$nLabel %in% c("noise","G=1")),
		c("Chr","Position","GeneName","TARGET","x_coor","mutID","nAD","nRD","nDP","nVAF")]
	names(input_data) = c("Chr","Position","GeneName","TARGET","x_coor","mutID","AD","RD","DP","VAF")
	input_data = input_data[order(input_data$Position),]
	# input_data = input_data[which(input_data$Position >= 150e6),]
	tmp_df = smart_df(Position = sort(unique(input_data$Position)))
	tmp_df$index = seq(nrow(tmp_df))
	input_data = smart_merge(input_data,tmp_df); rm(tmp_df)
	input_data$pt_cex = 0.5*log10(input_data$DP)

	cat(paste0("\tchr",chr,": "))
	out_AF_SEG = run_AF_SEG(
		DATA = input_data,
		num_divides = 5e1,
		binom = binom,
		clust_tVAF = !TRUE,
		show_plot = show_plot,
		show_mess = show_mess)

	infer_out = out_AF_SEG$INFER
	infer_out$start_x_coor = NA; infer_out$end_x_coor = NA
	infer_out$Chr = paste0("chr",chr)
	infer_out$Start_Position = NA; infer_out$End_Position = NA
	infer_out$Genes = NA; input_data$pt_col = NA
	infer_out$nVAF_mean = NA; infer_out$nVAF_sd = NA; infer_out$N = NA

	for(ii in seq(nrow(infer_out))){
		# ii = 8
		tmp_index1 = which(input_data$index >= infer_out$start[ii]
			& input_data$index <= infer_out$end[ii])
		infer_out[ii,c("start_x_coor","end_x_coor")] = range(input_data$x_coor[tmp_index1])
		infer_out[ii,c("Start_Position","End_Position")] = range(input_data$Position[tmp_index1])
		tmp_index2 = which(input_data$Position >= infer_out$Start_Position[ii]
			& input_data$Position <= infer_out$End_Position[ii])
		tmp_index3 = which(input_data$Position >= infer_out$Start_Position[ii]
			& input_data$Position <= infer_out$End_Position[ii]
			& !is.na(input_data$GeneName)
			& input_data$TARGET == "YES")
		tmp_tab = smart_table(unique(input_data[tmp_index3,c("GeneName","mutID")])[,1])
		if(length(tmp_tab) > 0){
		infer_out$Genes[ii] = ifelse(infer_out$eps[ii] >= eps_thres,
			paste(paste0(names(tmp_tab),"(",tmp_tab,")"),collapse=";"),
			NA)
		}
		input_data$pt_col[tmp_index2] = ifelse(infer_out$eps[ii] >= eps_thres | infer_out$PSI[ii] >= psi_thres,
			rgb(1,0,0,0.5),rgb(0,0,0,0.5))
		infer_out[ii,c("nVAF_mean","nVAF_sd","N")] = c(mean(input_data$VAF[tmp_index1]),
			sd(input_data$VAF[tmp_index1]),length(tmp_index1))
	}
	# infer_out

	# Output annotated input_data
	list(annoDATA = input_data,INFER = infer_out)
}
segment_tVAF = function(DATA,chr,binom=TRUE,show_plot=FALSE,show_mess=FALSE,BAF_thres){
	input_dat_vars = c("Chr","Position","x_coor","mutID","Qscore",
		"tAD","tRD","tDP","tVAF","GeneName","TARGET")
	
	input_data = DATA[which(DATA$Chr == paste0("chr",chr)
		& DATA$tVAF > BAF_thres & DATA$tVAF < 1 - BAF_thres
		& DATA$Qscore >= 10
		),input_dat_vars]

	input_data = name_change(input_data,"tAD","AD")
	input_data = name_change(input_data,"tRD","RD")
	input_data = name_change(input_data,"tDP","DP")
	input_data = name_change(input_data,"tVAF","VAF")
	input_data$LBC = lchoose(input_data$DP,input_data$AD)
	input_data = input_data[order(input_data$Position),]
	
	input_data = smart_uniq_df(input_df = input_data,vars=c("Chr","Position","mutID"))
	input_data = input_data[order(input_data$Position),]
	input_data$pt_cex = 0.5*log10(input_data$DP)
	input_data$index = seq(nrow(input_data))
	
	# Run segmentation!!!
	cat(paste0("\tchr",chr,": "))
	out_AF_SEG = run_AF_SEG(
		DATA = input_data,
		num_divides = 5e1,
		binom = binom,
		clust_tVAF = TRUE,
		show_plot = show_plot,
		show_mess = show_mess)
	
	# After segmentation, infer
	infer_out = out_AF_SEG$INFER
	infer_out$start_x_coor = NA; infer_out$end_x_coor = NA
	infer_out$Chr = paste0("chr",chr)
	infer_out$Start_Position = NA; infer_out$End_Position = NA
	infer_out$iBAF = NA; infer_out$GeneName = ""
	eps_thres = 1
	for(ii in seq(nrow(infer_out))){
		# ii = 1
		tmp_index = which(input_data$index >= infer_out$start[ii]
			& input_data$index <= infer_out$end[ii])
		infer_out[ii,c("start_x_coor","end_x_coor")] = range(input_data$x_coor[tmp_index])
		infer_out[ii,c("Start_Position","End_Position")] = range(input_data$Position[tmp_index])
		tmp_index3 = which(input_data$Position >= infer_out$Start_Position[ii]
			& input_data$Position <= infer_out$End_Position[ii]
			& !is.na(input_data$GeneName)
			& input_data$TARGET == "YES")
		tmp_tab = smart_table(unique(input_data[tmp_index3,c("GeneName","mutID")])[,1])
		
		if( infer_out$eps[ii] == eps_thres ){
			if( length(tmp_tab) > 0 ) infer_out$GeneName[ii] = paste(paste0(names(tmp_tab),"(",tmp_tab,")"),collapse=";")
		} else if( infer_out$eps[ii] < eps_thres ){
			if( infer_out$prob_half[ii] >= 0.5 ){
				infer_out$iBAF[ii] = 0.5
			} else {
				infer_out$iBAF[ii] = min(c(infer_out$mBAF[ii],1-infer_out$mBAF[ii]))
			}
		}
	}

	list(annoDATA=input_data,AF_SEG=out_AF_SEG,INFER=infer_out)
}

UNMASC_nSEG = function(outdir,vcs,tumorID,hg,genome,binom=TRUE,eps_thres=0.5,psi_thres=0.01){
	if(FALSE){
		outdir = outdir; vcs = vcf[which(vcf$Chr == "chr1"),]
		tumorID = tumorID; hg = hg; genome = genome
		binom = binom
	}
	
	cat("normal VAF segmentation ...\n")
	anno_all_vcs_nVAF = c()
	all_nVAF_segs = c()
	nSEG_dir = file.path(outdir,"nSEG")
	smart_mkdir(nSEG_dir)
	for(chr in 1:22){
		# chr = 1
		if( length(which(vcs$Chr == paste0("chr",chr))) == 0 ) next
		tmp_list = segment_nVAF(DATA = vcs,chr = chr,binom = binom,
			eps_thres = eps_thres,psi_thres = psi_thres,show_plot = FALSE,
			show_mess = FALSE)
		
		png(file.path(nSEG_dir,paste0("nSEG_chr",chr,".png")),units='px',
			height=1800,width=3000,res=250,type='cairo')
		par(mfrow=c(2,1),mar=c(5,4,1,6),xpd=TRUE)
		genome_plot(tmp_list$annoDATA,var="VAF",hg=hg,chr=chr,pch=16,
			cex=tmp_list$annoDATA$pt_cex,how_color="pt_col",genome=genome)
		legend("bottomright",inset = c(-0.1,0),legend = c(10^seq(3)),
			pt.cex = 0.5*log10(10^seq(3)),pch=rep(16,3),title="DP")
		plot(tmp_list$annoDATA[,c("index","VAF")],pch = 16,ylab = "nVAF",
			cex = tmp_list$annoDATA$pt_cex,col = tmp_list$annoDATA$pt_col,bty = "n",
			xaxs = "i",ylim = c(0,1))
		legend("bottomright",inset = c(-0.1,0),legend = c("H2M","nonH2M"),
			pt.cex = 1.5,pch = 15,col = c("red","black"))
		par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1,xpd=FALSE)
		dev.off()
		
		anno_all_vcs_nVAF = rbind(anno_all_vcs_nVAF,tmp_list$annoDATA)
		all_nVAF_segs = rbind(all_nVAF_segs,tmp_list$INFER)
		rm(tmp_list)
	}
	
	saveRDS(list(anno_all_vcs_nVAF=anno_all_vcs_nVAF,all_nVAF_segs=all_nVAF_segs),
		file.path(outdir,"anno_nSEG.rds"))
	png(file.path(nSEG_dir,"nSEG.png"),units='px',height=1500,width=3000,res=250,type='cairo')
	par(mar=c(5,4,1,6),xpd=TRUE)
	genome_plot(anno_all_vcs_nVAF,var="VAF",cex=anno_all_vcs_nVAF$pt_cex,pch=16,
		hg=hg,how_color="pt_col",genome=genome,main=tumorID)
	legend("bottomright",inset = c(-0.1,0),legend = c(10^seq(3)),
		pt.cex = 0.5*log10(10^seq(3)),pch=rep(16,3),title="DP")
	par(mar=c(5,4,4,2)+0.1,xpd=FALSE)
	dev.off()
	
	NULL
}
UNMASC_tSEG = function(outdir,tumorID,vcs,strand,gender,cut_BAF,hg,genome,
	rd_thres,ad_thres,qscore_thres,binom=TRUE){
	
	if(FALSE){
		
		outdir = outdir; tumorID = tumorID
		vcs = vcf; strand = strand; gender = gender
		cut_BAF = cut_BAF; hg = hg; genome = genome
		rd_thres = rd_thres; ad_thres = ad_thres
		qscore_thres = qscore_thres; binom = binom
		
	}
	
	cat("Subset strand to variants with depth and qscore thresholds and no strand bias ...\n")
	vcs = vcs[which(vcs$tDP >= rd_thres & vcs$tAD >= ad_thres
		& vcs$Qscore >= qscore_thres),]
	mutIDs = unique(vcs$mutID)
	strand = strand[which(strand$mutID %in% mutIDs),]
	strand = strand[which(is.na(strand$Strand_bias_P) 
		| strand$Strand_bias_P > 0.05),]
	
	cat("Determine OXOG,FFPE,ARTI status ...\n")
	num_uniq_vars = nrow(strand)
	cat(paste0("Number of unique variants = ",num_uniq_vars,"\n"))
	strand$OXOG = NA
	strand$FFPE = NA
	strand$ARTI = NA
	for(ii in seq(num_uniq_vars)){
		# ii = 1
		if(ii %% 5e2 == 0) cat(".")
		if(ii %% 4e4 == 0 || ii == num_uniq_vars ) cat(paste0(ii,"\n"))
		mutid_index = which(vcs$mutID == strand$mutID[ii])
		
		tmp_tab = smart_table(vcs$OXOG[mutid_index])
		tmp_value = paste0(names(tmp_tab),"(",as.numeric(tmp_tab),",",
			round(as.numeric(tmp_tab/sum(tmp_tab))*100,2),"%)")
		strand$OXOG[ii] = paste(tmp_value,collapse=";")
		
		tmp_tab = smart_table(vcs$FFPE[mutid_index])
		tmp_value = paste0(names(tmp_tab),"(",as.numeric(tmp_tab),",",
			round(as.numeric(tmp_tab/sum(tmp_tab))*100,2),"%)")
		strand$FFPE[ii] = paste(tmp_value,collapse=";")
		
		tmp_tab = smart_table(vcs$ARTI[mutid_index])
		tmp_value = paste0(names(tmp_tab),"(",as.numeric(tmp_tab),",",
			round(as.numeric(tmp_tab/sum(tmp_tab))*100,2),"%)")
		strand$ARTI[ii] = paste(tmp_value,collapse=";")
	}
	strand$OXOG_artifact = !sapply(strand$OXOG,artifact_filter,USE.NAMES=FALSE)
	strand$FFPE_artifact = !sapply(strand$FFPE,artifact_filter,USE.NAMES=FALSE)
	strand$ARTI_artifact = !sapply(strand$ARTI,artifact_filter,USE.NAMES=FALSE)
	
	# Segment non-OXOG, non-FFPE, non-strand bias, non-ARTI variants
	strand = strand[which(strand$OXOG_artifact == FALSE
		& strand$FFPE_artifact == FALSE
		& strand$ARTI_artifact == FALSE
		),]
	vcs = vcs[which(vcs$mutID %in% strand$mutID),]
	
	cat("tumor VAF segmentation on variants...\n")
	anno_all_vcs_tVAF = c()
	all_tVAF_segs = c()
	tumor_chrs = 1:22
	if( !is.na(gender) && gender == "FEMALE" ){
		tumor_chrs = c(1:22,"X")
	}
	
	tSEG_dir = file.path(outdir,"tSEG")
	smart_mkdir(tSEG_dir)
	for(chr in tumor_chrs){
		if( length(which(vcs$Chr == paste0("chr",chr))) == 0 ) next
		
		if(FALSE){
			vv = vcs[which(vcs$tVAF > cut_BAF & vcs$tVAF < 1-cut_BAF
				# & vcs$nANNO == "mappable"
				),]
			dim(vv)
			smart_table(vv$RefAlt,vv$tVAF < 0.2)
			vv$RefAlt2 = vv$RefAlt
				vv$RefAlt2[vv$RefAlt %in% c("A>C","T>G")] = "A>C"
				vv$RefAlt2[vv$RefAlt %in% c("A>G","T>C")] = "A>G"
				vv$RefAlt2[vv$RefAlt %in% c("A>T","T>A")] = "A>T"
				vv$RefAlt2[vv$RefAlt %in% c("C>A","G>T")] = "C>A"
				vv$RefAlt2[vv$RefAlt %in% c("C>G","G>C")] = "C>G"
			smart_table(vv$RefAlt2,vv$tVAF < 0.2)
			plot(vv[,c("x_coor","tVAF")],cex = 0.5*log10(vv$tDP),
				col = factor(vv$RefAlt2),pch = 16,ylim = c(0,1))
			vv[which(vv$mutID == 32),]
			
			plot(vv[,c("x_coor","tVAF")],cex = 0.5*log10(vv$tDP),pch = 16,
				col = factor(vv$ThsdG_AF > 0 | vv$EXAC_AF > 0),ylim = c(0,1))
			
		}
		
		# Run segmentation on all non-candidate tumor-only variant calls
		tmp_list = segment_tVAF(DATA = vcs,
			chr = chr,binom = binom,show_plot = FALSE,show_mess = FALSE,
			BAF_thres = cut_BAF)

		png(file.path(tSEG_dir,paste0("tSEG_chr",chr,".png")),units='px',
			height=1800,width=3000,res=250,type='cairo')
		par(mfrow=c(2,1),mar=c(5,4,1,6),xpd=TRUE)
		# 1)
		genome_plot(tmp_list$annoDATA,var="VAF",hg=hg,chr=chr,pch=16,
			cex = tmp_list$annoDATA$pt_cex,how_color = rgb(0,0,0,0.5),
			genome = genome)
		legend("bottomright",inset = c(-0.1,0),legend = c(10^seq(3)),
			pt.cex = 0.5*log10(10^seq(3)),pch=rep(16,3),title="DP")
		apply(tmp_list$INFER[which(!is.na(tmp_list$INFER$iBAF)),c("iBAF","start_x_coor","end_x_coor")],1,
			function(x) segments(x0=x[2],y0=c(x[1],1-x[1]),x1=x[3],lwd=3,col="magenta"))
		# 2)
		plot(tmp_list$annoDATA$VAF,pch = 16,ylab = "tVAF",cex = tmp_list$annoDATA$pt_cex,
			col = rgb(0,0,0,0.5),ylim = c(0,1))
		apply(tmp_list$INFER[which(!is.na(tmp_list$INFER$iBAF)),c("iBAF","start","end")],1,
			function(x) segments(x0=x[2],y0=c(x[1],1-x[1]),x1=x[3],lwd=3,col="magenta"))
		par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1,xpd=FALSE)
		dev.off()
		
		anno_all_vcs_tVAF = rbind(anno_all_vcs_tVAF,tmp_list$annoDATA)
		all_tVAF_segs = rbind(all_tVAF_segs,tmp_list$INFER)
		rm(tmp_list)
	}
	
	saveRDS(list(anno_all_vcs_tVAF=anno_all_vcs_tVAF,all_tVAF_segs=all_tVAF_segs),
		file.path(outdir,"anno_tSEG.rds"))
	
	png(file.path(tSEG_dir,"tSEG.png"),units='px',height=1500,width=3000,res=250,type='cairo')
	par(mar=c(5,4,1,6),xpd=TRUE)
	genome_plot(anno_all_vcs_tVAF,var = "VAF",cex = anno_all_vcs_tVAF$pt_cex,
		pch = 16,hg = hg,how_color = rgb(0,0,0,0.5),genome = genome,
		main = tumorID)
	legend("bottomright",inset = c(-0.1,0),legend = c(10^seq(3)),
		pt.cex = 0.5*log10(10^seq(3)),pch=rep(16,3),title="DP")
	apply(all_tVAF_segs[which(!is.na(all_tVAF_segs$iBAF)),c("iBAF","start_x_coor","end_x_coor")],1,
		function(x) segments(x0=x[2],y0=c(x[1],1-x[1]),x1=x[3],lwd=3,col="magenta"))
	par(mar=c(5,4,4,2)+0.1,xpd=FALSE)
	dev.off()
	
	NULL
}
anno_nSEG = function(DAT,nSEG,eps_thres=0.5,psi_thres=0.01){
	DAT$seg = NA
	for(ii in seq(nrow(nSEG))){
		# ii = 1
		idx = which(DAT$Chr == nSEG$Chr[ii]
			& DAT$Position >= nSEG$Start_Position[ii]
			& DAT$Position <= nSEG$End_Position[ii])
		if( length(idx) > 0 ){
			DAT$seg[idx] = nSEG$seg[ii]
		}
		rm(idx)
	}
	dim(DAT)
	DAT = smart_merge(smart_rmcols(DAT,c("N_eps","N_PSI")),
		nSEG[,c("seg","Chr","eps","PSI")],all.x = TRUE)
	DAT = name_change(DAT,"eps","N_eps")
	DAT = name_change(DAT,"PSI","N_PSI")
	DAT = smart_rmcols(DAT,"seg")
	DAT$nANNO = ifelse(DAT$N_eps >= eps_thres 
		| DAT$N_PSI >= psi_thres,"H2M","mappable")
	
	DAT
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
	cat(sprintf("%s: Total number of loci for strand = %s\n",date(),nrow(vcs)))
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
	cat("nANNO ...\n")
	VC$n_seg = NA
	nSEG$n_seg = NA
	num_nSEG = nrow(nSEG)
	for(ii in seq(num_nSEG)){
		# ii = 1
		if( ii %% 2 == 0 ) cat(".")
		if( ii %% 5e1 == 0 || ii == num_nSEG ) cat(sprintf("%s out of %s\n",ii,num_nSEG))
		
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
	cat("\n")
	
	# Infer H2M
	cat("Infer H2M status ...\n")
	VC$nANNO = ifelse(VC$N_eps >= eps_thres | VC$N_PSI >= psi_thres,"H2M","mappable")

	# Append tSEG results
	cat("tANNO ...\n")
	VC$t_seg = NA
	tSEG$t_seg = NA
	num_tSEG = nrow(tSEG)
	for(ii in seq(num_tSEG)){
		# ii = 1
		if( ii %% 2 == 0 ) cat(".")
		if( ii %% 5e1 == 0 || ii == num_tSEG ) cat(sprintf("%s out of %s\n",ii,num_tSEG))
		
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
	cat("\n")
	
	# Infer BAF status and Posterior probability of BAF
	cat("Infer ALLELE_STAT and tANNO ...\n")
	VC$ALLELE_STAT = apply(VC[,paste0("T_",clust_names)],1,
		function(x) get_ALLELE_STAT(x))
	VC$tANNO = apply(VC[,c("tAD","tRD",paste0("T_",clust_names))],1,
		function(x) get_INFER(vec_RD = x[1:2],params = x[3:8]))
	if(FALSE){
		xx = apply(VC[51:70,c("tAD","tRD",paste0("T_",clust_names))],1,
			function(x) get_INFER(vec_RD = x[1:2],params = x[3:8]))
		
	}
	cat("\n")
	
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

# Simulate functions
gen_genome = function(outdir){

	cat("Create hg19 chromosome lengths and centromere positions.\n")
	genome = smart_df(Chr = paste0("chr",c(1:22,"X","Y")))
	genome$nChr = seq(24)
	genome$length = as.integer(c(249250621,243199373,198022430,
		191154276,180915260,171115067,159138663,146364022,141213431,
		135534747,135006516,133851895,115169878,107349540,102531392,
		90354753,81195210,78077248,59128983,63025520,48129895,
		51304566,155270560,59373566))
	genome$centro_start = as.integer(c(121535434,92326171,90504854,
		49660117,46405641,58830166,58054331,43838887,47367679,39254935,
		51644205,34856694,16000000,16000000,17000000,35335801,22263006,
		15460898,24681782,26369569,11288129,13000000,58632012,10104553))
	genome$centro_end = as.integer(c(124535434,95326171,93504854,
		52660117,49405641,61830166,61054331,46838887,50367679,42254935,
		54644205,37856694,19000000,19000000,20000000,38335801,25263006,
		18460898,27681782,29369569,14288129,16000000,61632012,13104553))
	# genome
	
	## Centromeres
	cent_fn = file.path(outdir,"hg19_centromeres.bed")
	smart_WT(genome[,c("Chr","centro_start","centro_end")],
		cent_fn,sep="\t",col.names=FALSE)
	
	# Lengths
	genome$Chr2 = paste0("SN:",genome$Chr)
	genome$length2 = paste0("LN:",genome$length)
	genome$V1 = "@SQ"
	genome$chr_prop = 1 / (10 + genome$nChr)
	genome$chr_prop[24] = 0
	genome$chr_prop = genome$chr_prop / sum(genome$chr_prop)
	length_fn = file.path(outdir,"hg19.dict")
	smart_WT(genome[,c("V1","Chr2","length2")],length_fn,
		sep="\t",col.names=FALSE)

	# Target regions
	bed = genome[,c("Chr","length")]
	bed$chrom = bed$Chr
	bed$start = 1
	bed$end = bed$length
	bed = bed[,c("chrom","start","end")]
	# dim(bed); bed[1:10,]
	smart_WT(bed,file.path(outdir,"target.bed"),sep="\t")

	genome
}
gen_tumor_state = function(genome,cn_state,purity,max_PSI){
	cat("Generate underlying tumor copy number.\n")
	tumor_segs = c()
	half_bounds = 0.95
	noise_bound = 0.1
	for(chr in paste0("chr",c(1:22,"X"))){
		# chr = "chr1"
		num_segs = sample(seq(2,8),1)
		tmp_length = genome$length[which(genome$Chr == chr)]
		# segment start/stop
		tmp_seg = c(0,sort(sample(seq(tmp_length),num_segs-1)),tmp_length)
		for(seg in seq(num_segs)){
			# underlying CN state
			tmp_row = sample(seq(nrow(cn_state)),1,prob=cn_state$prob)
			tCN = cn_state$tCN[tmp_row]
			CN_A = cn_state$CN_A[tmp_row]
			CN_B = cn_state$CN_B[tmp_row]
			# allelic status
			prop_mBAF = 0.5
			if( CN_A == CN_B ){
				prop_half = runif(1,half_bounds,1)
				mBAF = 0.5
			} else {
				prop_half = runif(1,0,1-half_bounds)
				mBAF = ( min(c(CN_A,CN_B)) * purity + (1 - purity) ) /
					(tCN * purity + 2*(1-purity))
			}
			# add noise
			T_eps = runif(1,0,noise_bound)
			T_PSI = runif(1,0,max_PSI)
			tmp_df = smart_df(Chr = chr,Start = tmp_seg[seg]+1,End = tmp_seg[seg+1],
				tCN = tCN,CN_A = CN_A,CN_B = CN_B,T_eps = T_eps,T_prop_half = prop_half,
				T_mBAF = mBAF,T_prop_mBAF = prop_mBAF,T_PSI = T_PSI)
			tumor_segs = rbind(tumor_segs,tmp_df)
		}
	}
	
	tumor_segs
}
gen_normal_state = function(genome,max_PSI){
	cat("Generate underlying normal genome H2M regions.\n")
	normal_segs = c(); h2m_prob = 0.1; max_h2m_interval = 40e6
	for(chr in paste0("chr",c(1:22,"X"))){
		# chr = "chr1"
		num_segs = sample(seq(20),1)
		tmp_length = genome$length[which(genome$Chr == chr)]
		tmp_seg = c(0,sort(sample(seq(tmp_length),num_segs-1)),tmp_length)
		for(seg in seq(num_segs)){
			H2M = sample(c(TRUE,FALSE),1,prob=c(h2m_prob,1-h2m_prob))
			seg_length = tmp_seg[seg+1] - tmp_seg[seg] + 1
			if(chr == "chrX" || seg_length > max_h2m_interval) H2M = FALSE
			if( H2M ){
				which_H2M = rbinom(1,1,0.5)
				if( which_H2M == 1 ){
					N_PSI = 0
					if( max_PSI > 0.05 ) N_PSI = runif(1,0.05,max_PSI)
					N_eps = runif(1,0,0.5)
				} else {
					N_PSI = runif(1,0,min(c(max_PSI,0.05)))
					N_eps = runif(1,0.5,1)
				}
			} else {
				N_PSI = runif(1,0,min(c(max_PSI,0.05)))
				N_eps = runif(1,0,0.05)
			}
			tmp_df = smart_df(Chr = chr,Start = tmp_seg[seg]+1,
				End = tmp_seg[seg+1],N_eps = N_eps,N_PSI = N_PSI)
			normal_segs = rbind(normal_segs,tmp_df)
		}
	}
	
	normal_segs
}
gen_umn_counts = function(num_normals,uniq_vc,normal_segs,mean_DP,AA_vaf){
	cat("Generate UMN read counts: ")
	nVAF = list()
	probs = c(8,8,2,1,1); probs = probs / sum(probs)
	
	gen_AD = function(dp,maf,psi){
		if(psi == 0){
			rbinom(1,dp,maf)
		} else {
			rbetabinom(n = 1,prob = maf,size = dp,theta = 1/psi)
		}
	}
	
	for(nn in seq(num_normals)){
		cat(".")
		tmp_df = smart_df(uniq_vc[,c("mutID","Chr","Position")],
			STUDYNUMBER = paste0("N_",nn),nAD = NA,nDP = NA)
		for(ii in seq(nrow(normal_segs))){
			# ii = 1
			seg_index = which( uniq_vc$Chr == normal_segs$Chr[ii]
				& uniq_vc$Position >= normal_segs$Start[ii]
				& uniq_vc$Position <= normal_segs$End[ii] )
			tmp_num_loci = length(seg_index)
			if( tmp_num_loci > 0 ){
				tmp_df$nDP[seg_index] = rnbinom(tmp_num_loci,mu=mean_DP,size=1/0.1) + 5
				N_eps = normal_segs$N_eps[ii]
				N_PSI = normal_segs$N_PSI[ii]
				geno_prob = c(N_eps,(1-N_eps)*probs)
				tmp_class = sample(seq(6),tmp_num_loci,replace=TRUE,prob=geno_prob)
				tmp_df$nAD[seg_index][tmp_class == 1] = as.numeric(sapply(tmp_df$nDP[seg_index][tmp_class == 1],
					function(xx) sample(seq(xx),1),USE.NAMES = FALSE))
				tmp_df$nAD[seg_index][tmp_class == 2] = 0
				tmp_df$nAD[seg_index][tmp_class == 3] = as.numeric(sapply(tmp_df$nDP[seg_index][tmp_class == 3],
					function(xx) gen_AD(dp = xx,maf = AA_vaf,psi = N_PSI),USE.NAMES = FALSE))
				tmp_df$nAD[seg_index][tmp_class == 4] = as.numeric(sapply(tmp_df$nDP[seg_index][tmp_class == 4],
					function(xx) gen_AD(dp = xx,maf = 0.5,psi = N_PSI),USE.NAMES = FALSE))
				tmp_df$nAD[seg_index][tmp_class == 5] = as.numeric(sapply(tmp_df$nDP[seg_index][tmp_class == 5],
					function(xx) gen_AD(dp = xx,maf = 1-AA_vaf,psi = N_PSI),USE.NAMES = FALSE))
				tmp_df$nAD[seg_index][tmp_class == 6] = tmp_df$nDP[seg_index][tmp_class == 6]
			}
		}
		tmp_df$nRD = tmp_df$nDP - tmp_df$nAD
		tmp_df$nVAF = round(tmp_df$nAD / tmp_df$nDP,5)
		nVAF[[nn]] = tmp_df; rm(tmp_df)
	}
	cat("\n")
	
	nVAF
}
gen_vc_status = function(uniq_vc,tumor_segs,mean_DP,purity,ffpe_vaf,oxog_vaf){
	
	# Somatic status, artifact status
	cat("Generate variant call statuses\n")
	num_loci = nrow(uniq_vc)
	uniq_vc$EXAC_AF = 0; uniq_vc$ThsdG_AF = 0; uniq_vc$CosmicOverlaps = 0
	uniq_vc$tAD = NA; uniq_vc$status = NA; uniq_vc$Qscore = 50
	uniq_vc$tDP = rnbinom(num_loci,mu=mean_DP,size=1/0.1) + 5
	uniq_vc$status[which(uniq_vc$nGeno == 1)] = "germline_hetero"
	uniq_vc$status[which(uniq_vc$nGeno == 2)] = "germline_homozygous"
	uniq_vc$tAD[which(uniq_vc$nGeno == 2)] = uniq_vc$tDP[which(uniq_vc$nGeno == 2)]
	bases = c("A","C","G","T")
	uniq_vc$Ref = sample(bases,num_loci,replace=TRUE)
	uniq_vc$Alt = sapply(uniq_vc$Ref,function(xx) 
		sample(bases[bases != xx],1),USE.NAMES=FALSE)
	uniq_vc$RefAlt = paste0(uniq_vc$Ref,">",uniq_vc$Alt)
	uniq_vc$gc_content = round(runif(num_loci,0,1),3)
	uniq_vc$ENSG = "ENST00000000000"
	uniq_vc$variant_classification = "Missense_Mutation"
	uniq_vc$Onco_Effect = "missense"
	uniq_vc$EFF = "NON_SYNONYMOUS_CODING(MODERATE|MISSENSE||p.ABCDEFG||protein_coding|CODING|ENST00000000000|1)"
	uniq_vc$left_context = sapply(seq(num_loci),function(xx) 
		paste(sample(bases,10,replace=TRUE),collapse=""),USE.NAMES=FALSE)
	uniq_vc$right_context = sapply(seq(num_loci),function(xx) 
		paste(sample(bases,10,replace=TRUE),collapse=""),USE.NAMES=FALSE)
	uniq_vc$ref_context = paste0(uniq_vc$left_context,uniq_vc$Ref,uniq_vc$right_context)
	uniq_vc = uniq_vc[,!(names(uniq_vc) %in% c("left_context","right_context"))]
	uniq_vc$Strand_bias_P = 1
	uniq_vc$For_Ref = 0; uniq_vc$Rev_Ref = 0
	uniq_vc$For_Alt = 0; uniq_vc$Rev_Alt = 0
	uniq_vc$dbNSFP_MutationTaster = NA
	uniq_vc$CGC_SomMut = NA; uniq_vc$CGC_GermMut = NA
	
	aa_index = which( uniq_vc$nGeno == 0 
		& (is.na(uniq_vc$N_eps) | uniq_vc$N_eps < 0.5) )
	num_soma = 30; num_cosm = 5
	soma_index = sample(aa_index,num_soma)
	cosm_index = sample(soma_index,num_cosm)
	uniq_vc$CosmicOverlaps[cosm_index] = sample(seq(10,500),num_cosm,replace=TRUE)
	uniq_vc$status[soma_index] = "somatic"
	uniq_vc[which(uniq_vc$status == "somatic"),]
	uniq_vc$mult = NA
	for(ii in soma_index){
		# ii = soma_index[1]
		one_index = which(tumor_segs$Chr == uniq_vc$Chr[ii]
			& tumor_segs$Start <= uniq_vc$Position[ii]
			& tumor_segs$End >= uniq_vc$Position[ii])
		tumor_segs[one_index,]
		# select allele
		allele = sample(c("A","B"),1)
		poss_mult = unique(c(1,tumor_segs[one_index,paste0("CN_",allele)]))
		poss_mult = poss_mult[poss_mult > 0]
		uniq_vc$mult[ii] = sample(poss_mult,1)
		tvaf = purity * uniq_vc$mult[ii] / (tumor_segs$tCN[one_index] 
			* purity + 2*(1-purity))
		uniq_vc$tAD[ii] = rbinom(1,uniq_vc$tDP[ii],tvaf)
	}

	prop_ffpe = 0.1; prop_ingerm = 0.4
	miss_status = which(is.na(uniq_vc$status) 
		& uniq_vc$RefAlt %in% c("C>T","G>A"))
	num_ffpe = round(length(miss_status)*prop_ffpe)
	ffpe_index = sample(miss_status,num_ffpe)
	uniq_vc$status[ffpe_index] = "FFPE"
	uniq_vc$tAD[ffpe_index] = rbinom(num_ffpe,uniq_vc$tDP[ffpe_index],ffpe_vaf)
	ffpe_ingerm_index = sample(ffpe_index,round(prop_ingerm*length(ffpe_index)))
	uniq_vc$EXAC_AF[ffpe_ingerm_index] = runif(length(ffpe_ingerm_index),0,1)

	prop_oxog = 0.1; prop_ingerm = 0.4
	miss_status = which(is.na(uniq_vc$status) 
		& uniq_vc$RefAlt %in% c("C>A","G>T"))
	num_oxog = round(length(miss_status)*prop_oxog)
	oxog_index = sample(miss_status,num_oxog)
	uniq_vc$status[oxog_index] = "OXOG"
	uniq_vc$tAD[oxog_index] = rbinom(num_oxog,uniq_vc$tDP[oxog_index],oxog_vaf)
	oxog_ingerm_index = sample(oxog_index,round(prop_ingerm*length(oxog_index)))
	uniq_vc$EXAC_AF[oxog_ingerm_index] = runif(length(oxog_ingerm_index),0,1)

	miss_status = which(is.na(uniq_vc$status))
	uniq_vc$status[miss_status] = "no_mutation"
	uniq_vc$tAD[miss_status] = 0
	
	uniq_vc
}
gen_tumor_counts = function(tumor_segs,uniq_vc){
	# Tumor germline read counts
	cat("Generate tumor germline read counts\n")
	prop_ingerm = 0.8
	
	gen_AD = function(dp,maf,psi){
		if(psi == 0){
			rbinom(1,dp,maf)
		} else {
			rbetabinom(n = 1,prob = maf,size = dp,theta = 1/psi)
		}
	}
	
	for(ii in seq(nrow(tumor_segs))){
		# ii = 1
		seg_index = which( uniq_vc$Chr == tumor_segs$Chr[ii]
			& uniq_vc$Position >= tumor_segs$Start[ii]
			& uniq_vc$Position <= tumor_segs$End[ii] 
			& uniq_vc$nGeno == 1 )
		# tumor_segs[ii,]
		tmp_num_loci = length(seg_index)
		if( tmp_num_loci > 0 ){
			T_eps = tumor_segs$T_eps[ii]
			T_PSI = tumor_segs$T_PSI[ii]
			T_prop_half = tumor_segs$T_prop_half[ii]
			T_mBAF = tumor_segs$T_mBAF[ii]
			T_prop_mBAF = tumor_segs$T_prop_mBAF[ii]
			probs = c(T_eps,(1-T_eps)*c(T_prop_half,
				(1-T_prop_half)*c(T_prop_mBAF,1-T_prop_mBAF)))
			# probs
			tmp_class = sample(seq(4),tmp_num_loci,replace=TRUE,prob=probs)
			uniq_vc$tAD[seg_index][tmp_class == 1] = as.numeric(sapply(uniq_vc$tDP[seg_index][tmp_class == 1],
				function(xx) sample(seq(xx),1),USE.NAMES = FALSE))
			uniq_vc$tAD[seg_index][tmp_class == 2] = as.numeric(sapply(uniq_vc$tDP[seg_index][tmp_class == 2],
				function(xx) gen_AD(dp = xx,maf = 0.5,psi = T_PSI),USE.NAMES = FALSE))
			uniq_vc$tAD[seg_index][tmp_class == 3] = as.numeric(sapply(uniq_vc$tDP[seg_index][tmp_class == 3],
				function(xx) gen_AD(dp = xx,maf = T_mBAF,psi = T_PSI),USE.NAMES = FALSE))
			uniq_vc$tAD[seg_index][tmp_class == 4] = as.numeric(sapply(uniq_vc$tDP[seg_index][tmp_class == 4],
				function(xx) gen_AD(dp = xx,maf = 1-T_mBAF,psi = T_PSI),USE.NAMES = FALSE))
			num_ingerm = round(tmp_num_loci*prop_ingerm)
			ingerm_index = sample(seg_index,num_ingerm)
			uniq_vc$EXAC_AF[ingerm_index] = runif(num_ingerm,0,0.5)
		}
	}
	uniq_vc$tRD = uniq_vc$tDP - uniq_vc$tAD
	uniq_vc$tVAF = round(uniq_vc$tAD / uniq_vc$tDP,5)
	uniq_vc$pt_cex = 0.5*log10(uniq_vc$tDP)
	uniq_vc$pt_col = NA
	uniq_vc$pt_col[which(uniq_vc$status == "germline_hetero")] = rgb(0,1,0,0.25)
	uniq_vc$pt_col[which(uniq_vc$status == "germline_homozygous")] = rgb(0,0,1,0.25)
	uniq_vc$pt_col[which(uniq_vc$status == "somatic")] = rgb(1,0,0,0.5)
	uniq_vc$pt_col[which(uniq_vc$status == "no_mutation")] = rgb(0,0,0,0.25)
	uniq_vc$pt_col[which(uniq_vc$status == "FFPE")] = rgb(1,0,1,0.5)
	uniq_vc$pt_col[which(uniq_vc$status == "OXOG")] = rgb(1,0.5,0,0.5)
	uniq_vc$pt_pch = 1
	uniq_vc$pt_pch[which(uniq_vc$status 
		%in% c("somatic","FFPE","OXOG","germline_hetero"))] = 16
	
	uniq_vc
}
gen_make_vcfs = function(outdir,num_normals,uniq_vc,nVAF){
	# Make VCFs
	cat("Making VCFs: ")
	uniq_vc$EXAC_AF = round(uniq_vc$EXAC_AF,4)
	for(nn in seq(num_normals)){
		# nn = 1
		# dim(uniq_vc); dim(nVAF[[nn]])
		cat(".")
		tmp_vcf = cbind(uniq_vc,nVAF[[nn]][,c("nAD","nRD","nDP","STUDYNUMBER")])
		tmp_vcf$ID = "."; tmp_vcf$FILTER = "PASS"; tmp_vcf$FORMAT = "GT:AD:DP"
		tmp_vcf$CHROM = gsub("chr","",tmp_vcf$Chr)
		tmp_vcf$POS = tmp_vcf$Position; tmp_vcf$REF = tmp_vcf$Ref
		tmp_vcf$ALT = tmp_vcf$Alt; tmp_vcf$QUAL = tmp_vcf$Qscore
		tmp_vcf$INFO = paste0(
			"ExAC_AF=",tmp_vcf$EXAC_AF,";",
			"1000gp3_AF=",tmp_vcf$ThsdG_AF,";",
			"ExAC_ESP_AF_GLOBAL=",tmp_vcf$EXAC_AF,";",
			"gc_content=",tmp_vcf$gc_content,";",
			"gene=Unknown;",
			"COSMIC_n_overlapping_mutations=",tmp_vcf$CosmicOverlaps,";",
			"HGNC_Ensembl_Gene_ID=",tmp_vcf$ENSG,";",
			"variant_classification=",tmp_vcf$variant_classification,";",
			"Ensembl_so_term=",tmp_vcf$Onco_Effect,";",
			"ref_context=",tmp_vcf$ref_context,";",
			"EFF=",tmp_vcf$EFF,";"
			)
		tmp_vcf$NORMAL = paste0("./.:",tmp_vcf$nRD,"|",tmp_vcf$nAD,":",tmp_vcf$nDP)
		tmp_vcf$TUMOR = paste0("./.:",tmp_vcf$tRD,"|",tmp_vcf$tAD,":",tmp_vcf$tDP)
		tmp_vcf$Position = as.integer(tmp_vcf$Position)
		names(tmp_vcf)[names(tmp_vcf) == "CHROM"] = "#CHROM"
		# tmp_vcf[1:3,]
		
		tmp_vcf = tmp_vcf[,c("#CHROM","POS","ID","REF","ALT",
			"QUAL","FILTER","INFO","FORMAT","NORMAL","TUMOR")]
		# tmp_vcf[1:5,]
		smart_WT(tmp_vcf,file.path(outdir,sprintf("UMN_%s.vcf",nn)),sep="\t")
		rm(tmp_vcf)
	}
	cat("\n")
}

#' @title gen_simulated_data
#' @description The function simulates one sample's underlying 
#'	somatic point mutations and copy number aberrations and creates 
#'	vcf files formatted for UNMASC. Additional files created include
#'	bed files defining chromosome lengths, centromere regions, and
#'	target regions.
#' @param outdir Character string specifying the output directory.
#' @param purity Tumor purity. By default, \code{purity = NULL} 
#'	will randomly generate a tumor purity.
#' @param cn_state A R data.frame containing copy number states 
#'	with columns \code{tCN}, \code{CN_A}, \code{CN_B}, and \code{prob} 
#'	corresponding to the total integer copy number, the minor allelic 
#'	integer copy number, the major allelic integer copy number, and 
#'	probability of selecting a copy number state, respectively. 
#'	By default, \code{cn_state = NULL}, will generate a default 
#'	\code{cn_state} object.
#' @param num_loci A positive integer for the number of tumor only 
#'	variant loci simulated.
#' @param num_normals A positive integer for the number of unmatched 
#'	normal controls to simulate.
#' @param AA_vaf The mean non-zero variant allele frequency for 
#'	homozygous reference positions.
#' @param mean_DP A positve numeric value for the mean depth of 
#'	variant calls.
#' @param oxog_vaf A positive numeric value for the expected 
#'	allele frequency of oxoG variants.
#' @param ffpe_vaf A positive numeric value for the expected 
#'	allele frequency of FFPE variants.
#' @param max_PSI A positive numeric value for the maximum overdispersion in
#'	the beta-binomial distribution.
#' @param show_plot By default, set to \code{TRUE}. If \code{TRUE},
#'	the variant calls are plotted. This will include germline, artifact, 
#'	and somatic variants.
#' @param seed A numeric seed value. If \code{NULL}, a default seed is set.
#' @export
gen_simulated_data = function(outdir,purity = NULL,cn_state = NULL,
	num_loci = 1e4,num_normals = 20,AA_vaf = 0.002,mean_DP = 500,
	oxog_vaf = 0.01,ffpe_vaf = 0.1,max_PSI = 0.02,show_plot = TRUE,
	seed = NULL){
	
	if(FALSE){
		outdir = outdir; purity = 0.8; cn_state = NULL;
		num_loci = 1e4; num_normals = 20; AA_vaf = 0.002; mean_DP = 500;
		oxog_vaf = 0.01; ffpe_vaf = 0.1; show_plot = TRUE; seed = 1
	}
	
	smart_mkdir(outdir)
	genome = gen_genome(outdir)
	
	if( is.null(seed) ){
		set.seed(1)
	} else {
		set.seed(seed)
	}
	
	# Check inputs
	if( is.null(purity) ){
		purity = runif(1,0.2,0.9)	
	} else {
		if( !( length(purity) == 1 || is.numeric(purity) ) ){
			stop("Not a valid purity.")
		}
	}
	if( is.null(cn_state) ){
		cn_state = c()
		for(cc in seq(1,5)){
		for(aa in seq(0,floor(cc/2))){
			cn_state = rbind(cn_state,smart_df(tCN=cc,CN_A=aa,CN_B=cc-aa))
		}
		}
		cn_state$prob = 2
		cn_state$prob[cn_state$CN_A != cn_state$CN_B] = 1
		cn_state$prob = cn_state$prob / sum(cn_state$prob)
	} else {
		chk_df = is.data.frame(cn_state)
		chk_nm = all(names(cn_state) %in% c("tCN","CN_A","CN_B","prob"))
		chk_cols = ( is.integer(cn_state$tCN) 
			&& is.integer(cn_state$CN_A)
			&& is.integer(cn_state$CN_B)
			&& sum(cn_state$prob) == 1 )
		if( !( chk_df || chk_nm || chk_cols ) ){
			stop("Issue with cn_state.")
		}
	}
	
	# Simulate H2M regions
	normal_segs = gen_normal_state(genome = genome,max_PSI = max_PSI)

	cat("Allelic copy number states to sample from ...\n")
	print(cn_state,right=FALSE)
	Sys.sleep(2)
	tumor_segs = gen_tumor_state(genome = genome,
		cn_state = cn_state,purity = purity,max_PSI = max_PSI)
	
	# Initialize tumor variant calls
	uniq_vc = smart_df(mutID = seq(num_loci))

	# Chromosome
	uniq_vc$nChr = sort(sample(genome$nChr,num_loci,replace=TRUE,
		prob=genome$chr_prop))
	uniq_vc = smart_merge(uniq_vc,genome[,c("Chr","nChr")])
	uniq_vc = uniq_vc[order(uniq_vc$mutID),]

	# Position
	uniq_vc$Position = NA
	for(chr in paste0("chr",c(1:22,"X"))){
		# chr = "chr1"
		chr_index = which(uniq_vc$Chr == chr)
		tmp_length = genome$length[which(genome$Chr == chr)]
		uniq_vc$Position[chr_index] = sort(sample(seq(tmp_length),
			length(chr_index),replace=FALSE))
	}
	
	# Matched normal genotypes
	geno_prob = c(5,3,1) # Prob AA, AB, BB across loci
	geno_prob = geno_prob / sum(geno_prob)
	uniq_vc$nGeno = sample(c(0,1,2),num_loci,replace=TRUE,prob=geno_prob)
	
	# Generate UMN read counts
	nVAF = gen_umn_counts(num_normals = num_normals,
		uniq_vc = uniq_vc,normal_segs = normal_segs,
		mean_DP = mean_DP,AA_vaf = AA_vaf)

	# Identify H2M regions
	uniq_vc$N_eps = NA
	for(ii in seq(nrow(normal_segs))){
		tmp_index = which(uniq_vc$Chr == normal_segs$Chr[ii]
			& uniq_vc$Position >= normal_segs$Start[ii]
			& uniq_vc$Position <= normal_segs$End[ii])
		N_eps = normal_segs$N_eps[ii]
		if( length(tmp_index) > 0 ){
			uniq_vc$N_eps[tmp_index] = N_eps
		}
	}

	# Somatic status, artifact status
	uniq_vc = gen_vc_status(uniq_vc = uniq_vc,
		tumor_segs = tumor_segs,mean_DP = mean_DP,
		purity = purity,ffpe_vaf = ffpe_vaf,
		oxog_vaf = oxog_vaf)

	# Tumor germline read counts
	uniq_vc = gen_tumor_counts(tumor_segs = tumor_segs,
		uniq_vc = uniq_vc)

	# Plot
	if( show_plot ){
		cat("Plotting ...\n")
		def_par = par()
		par(mfrow=c(2,2),mar=c(4,4,0.2,0.2),bty="n")
		plot(nVAF[[1]]$nVAF,cex=uniq_vc$pt_cex,
			col=rgb(0,0,0,0.5),xlab="Ordered Positions",
			ylab="nVAF",main="")
		plot(uniq_vc$tVAF,cex=uniq_vc$pt_cex,col=uniq_vc$pt_col,
			pch=uniq_vc$pt_pch,xlab="Ordered Positions",ylab="tVAF",
			main="")
		smart_hist(uniq_vc$tVAF[which(uniq_vc$tAD >= 5
			& uniq_vc$tVAF < 0.4 & uniq_vc$EXAC_AF == 0 
			& uniq_vc$RefAlt %in% c("C>A","G>T"))],breaks=40,
			main="",xlim=c(0,0.4),xlab="OXOG distribution")
		smart_hist(uniq_vc$tVAF[which(uniq_vc$tAD >= 5
			& uniq_vc$tVAF < 0.4 & uniq_vc$EXAC_AF == 0 
			& uniq_vc$RefAlt %in% c("C>T","G>A"))],breaks=40,
			main="",xlim=c(0,0.4),xlab="FFPE distribution")
		par(mfrow=def_par$mfrow,mar=def_par$mar,bty=def_par$bty)
	}
	
	# Make strand data.frame
	strand = uniq_vc[,c("mutID","Chr","Position","Ref","Alt",
		"For_Ref","For_Alt","Rev_Ref","Rev_Alt","Strand_bias_P")]
	saveRDS(strand,file.path(outdir,"strand.rds"))

	# Make VCFs
	gen_make_vcfs(outdir = outdir,num_normals = num_normals,
		uniq_vc = uniq_vc,nVAF = nVAF)
	
	# Output arguments for main function
	out = list(vcf_files = file.path(outdir,paste0("UMN_",seq(num_normals),".vcf")),
		bed_regions = file.path(outdir,"target.bed"),
		bed_centromere = file.path(outdir,"hg19_centromeres.bed"),
		dict_chrom = file.path(outdir,"hg19.dict"),
		gender = "FEMALE",simulate = TRUE)
	
	out
}

#' @title simulate_BAF
#' @description This functions simulates B allele frequencies and 
#'	illustrates how UNMASC processes normal and tumor read counts 
#'	for segmentation.
#' @param num_segs An integer specifying the number of segments to simulate.
#' @param mean_depth An integer specifying a mean total depth for loci.
#' @param max_PSI The largest beta-binomial overdispersion term a simulated segment can be.
#' @param tumor To generate tumor read counts, set to \code{TRUE}, 
#'	otherwise set to \code{FALSE}.
#' @param show_plots Set to \code{TRUE} to plot the simulated output.
#' @export
simulate_BAF = function(num_segs = 2,mean_depth = 300,max_PSI = 0,
	tumor = TRUE,show_plots = TRUE){

	all_data = c()
	gen_AD = function(cc,DP,pp,psi){
		if(cc == 1){
			sample(seq(0,DP),1)
		} else if(cc == 2){
			if( psi == 0 ){
				rbinom(1,DP,0.5)
			} else {
				rbetabinom(n = 1,prob = 0.5,size = DP,theta = 1/psi)
			}
		} else if(cc == 3){
			if( psi == 0 ){
				rbinom(1,DP,pp)
			} else {
				rbetabinom(n = 1,prob = pp,size = DP,theta = 1/psi)
			}
		} else if(cc == 4){
			if( psi == 0 ){
				rbinom(1,DP,1-pp)
			} else {
				rbetabinom(n = 1,prob = 1-pp,size = DP,theta = 1/psi)
			}
		}
	}
	
	for(seg in seq(num_segs)){

		# For each segment
		N = sample(seq(5e1,5e2,10),1)
		
		# Parameters governing the VAF mixture distribution
		if( tumor ){ 	# For tumors
			true_eps = ifelse(rbinom(1,1,0.2) == 1,
				runif(1,0.25,0.5),runif(1,0,0.25))
			true_pp = runif(1,0.1,0.9)
			true_pi_half = ifelse(rbinom(1,1,0.5) == 1,1,0)
			true_pi_pp = runif(1,0.3,0.7)
		} else {			# For normals
			true_eps = ifelse(rbinom(1,1,0.2) == 1,
				runif(1,0.8,1),runif(1,0,0.2))
			true_pp = 0.5
			true_pi_half = 1
			true_pi_pp = 0.5
		}
		
		true_params = c(true_eps,true_pi_half,true_pp,true_pi_pp)
		true_probs = c(true_eps,(1-true_eps)*true_pi_half,
			(1-true_eps)*(1-true_pi_half)*true_pi_pp,
			(1-true_eps)*(1-true_pi_half)*(1-true_pi_pp))
		true_PSI = runif(1,0,max_PSI)

		obs_data = smart_df(seg = seg,
			eps = true_eps,prob_half = true_pi_half,
			mBAF = true_pp,prop_mBAF = true_pi_pp,
			DP = rnbinom(N,mu = mean_depth,size = 1.5) + 10,
			class = sample(seq(4),N,replace=TRUE,prob=true_probs),
			PSI = true_PSI)
		obs_data$AD = apply(obs_data[,c("class","DP","mBAF")],1,function(x)
			gen_AD(x[1],x[2],x[3],true_PSI))
		obs_data$VAF = obs_data$AD / obs_data$DP
		obs_data$RD = obs_data$DP - obs_data$AD
		obs_data$pt_cex = 0.5*log10(obs_data$DP)
		all_data = rbind(all_data,obs_data)

	}

	if( show_plots ){
		par(mar=c(5,4,2,1),bty="n")
		plot(all_data$VAF,cex=all_data$pt_cex,
			pch=16,ylim=c(0,1),col=all_data$seg,
			xlab="Ordered Positions",
			ylab=ifelse(tumor,"Tumor VAF","Normal VAF"))
			abline(h=0.5,lty=2)
		par(mar=c(5,4,4,1)+0.1,bty="o")
	}
	
	all_data$index = seq(nrow(all_data))
	true_seg = smart_df(t(sapply(seq(num_segs),function(x) range(all_data$index[which(all_data$seg %in% x)]))))
	names(true_seg) = c("start","end")
	true_seg = round(cbind(unique(all_data[,c("seg","eps","prob_half","mBAF",
		"prop_mBAF","PSI")]),true_seg),2)
	true_seg = true_seg[,c("seg","start","end","eps","prob_half","mBAF","prop_mBAF","PSI")]
	rownames(true_seg) = NULL
	list(truth=true_seg,true_seg=c(1,true_seg$end),all_data=all_data)
}

#' @title UNMASC_definitions
#' @description This function prints the column header definitions of 
#'	the output file \code{tumorOnly_VCs.tsv}.
#' @param headers A character vector of column header names. By default
#'	\code{headers = NULL}, in which all available header definitions 
#'	are printed. If a subset of headers are specified, only the subset 
#'	of header definitions are printed.
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
run_nCLUST_opt = function(mat_RD,binom,verbose = TRUE){
	if(FALSE){
		mat_RD = as.matrix(tmp_vcf[idx_clust,c("nAD","nRD")])
		
		vcf2 = vcf[which(vcf$nDP >= 10 & vcf$Qscore >= 30
			# & vcf$tDP >= 10 
			& vcf$nVAF > 0 & vcf$nVAF < 1),]
		vcf2_fn = file.path(type_dir,"vcf2.rds")
		saveRDS(vcf2,vcf2_fn)
		
		vcf2 = readRDS(vcf2_fn)
		mat_RD = as.matrix(vcf2[,c("nAD","nRD")])
		
		binom = FALSE; verbose = TRUE
		dim(mat_RD)
	}
	
	DP 			= rowSums(mat_RD)
	LBC 		= lchoose(DP,mat_RD[,1])
	log_DP 	= log(DP)
	
	max_iter = 4e3; eps = 1e-4
	iPSI_vec = 0
	if( binom == FALSE ){
		# iPSI_vec = c(0,0.005)
		iPSI_vec = c(0,0.05)
	}
	
	em_out = NULL
	props0_mat = matrix(c(1,0,0,0,
		0,1,1,1,
		1,1,1,1),
		ncol = 4,byrow = TRUE)
	props0_mat = props0_mat / rowSums(props0_mat)
	maf0_vec = seq(5)/100
	
	for(iPSI in iPSI_vec){
	for(jj in seq(nrow(props0_mat))){
	for(mm0 in maf0_vec){
		tmp_em = Rcpp_BAF_clust(RD = mat_RD,DP = DP,log_DP = log_DP,
			LBC = LBC,props0 = props0_mat[jj,],mm0 = mm0,
			psi0 = iPSI,max_iter = max_iter,eps = eps,show = FALSE)
		if( tmp_em$converge == "yes" ){
			if( is.null(em_out) ){
				em_out = tmp_em
			} else if( tmp_em$BIC > em_out$BIC ){
				em_out = tmp_em
			}
			if( verbose ){
				print(em_out[c("iter","converge","props",
					"LL","ms","BIC","maf","psi","params")])
			}
		}
	}}}
	
	if( verbose ){
		print(em_out[c("iter","converge","props",
			"LL","ms","BIC","maf","psi","params")])
	}
	
	em_out
}
run_nCLUST = function(vcf,ID = "STUDYNUMBER",rd_thres,ad_thres,
	qscore_thres,hg,binom = TRUE,outdir,genome,verbose = TRUE){
	
	if(FALSE){
		vcf = vcs
		ID = "STUDYNUMBER"; binom = !TRUE; verbose = TRUE
		
	}
	
	vcf = vcf[,!(names(vcf) %in% c("nInterpret","nLabel"))]
	IDs = unique(vcf[[ID]])
	out = c(); SE = c(); cnt = 1
	if( verbose ){
		if( binom ){
			cat("Clustering normal read counts assuming binomial mixture ...\n")
		} else {
			cat("Clustering normal read counts assuming beta-binomial mixture ...\n")
		}
	}
	
	nCLUST_dir = file.path(outdir,"nCLUST")
	smart_mkdir(nCLUST_dir)
	rd_filter = vcf$nDP >= rd_thres
	qs_filter = vcf$Qscore >= qscore_thres
	
	for(id in IDs){
		# id = IDs[1]; id
		if( verbose ){
			cat(sprintf("%s ",id))
			if( cnt > 0 && cnt %% 5 == 0 ) cat("\n")
		}
		sn_filter = vcf[[ID]] == id
		
		idx = which(sn_filter
			& rd_filter
			& qs_filter)
		tmp_vcf = vcf[idx,]
		tmp_vcf$orig_order = seq(nrow(tmp_vcf))
		
		idx_clust = which(tmp_vcf$nVAF > 0 & tmp_vcf$nVAF < 1
			& tmp_vcf$nAD >= ad_thres)
		mat_RD = as.matrix(tmp_vcf[idx_clust,c("nAD","nRD")]); # head(mat_RD)
		
		opt_list = run_nCLUST_opt(mat_RD = mat_RD,binom = binom,verbose = FALSE)
		# names(opt_list)
		# opt_list[c("iter","converge","props","LL","ms","BIC","maf","psi")]
		
		tmp_vcf$infer_class = NA
		tmp_vcf$infer_class[idx_clust][opt_list$class == 1] = 1 # noise
		tmp_vcf$infer_class[idx_clust][opt_list$class == 2] = 4 # MAF = 0.5
		tmp_vcf$infer_class[idx_clust][opt_list$class == 3] = 3 # MAF near 0
		tmp_vcf$infer_class[idx_clust][opt_list$class == 4] = 5 # MAF near 1
		tmp_vcf$infer_class[tmp_vcf$nVAF == 0] = 2
		tmp_vcf$infer_class[is.na(tmp_vcf$infer_class) & tmp_vcf$nAD < ad_thres] = 3
		tmp_vcf$infer_class[tmp_vcf$nVAF == 1] = 6
		# smart_table(tmp_vcf$infer_class)
		# smart_table(tmp_vcf$infer_class[idx_clust],opt_list$class)
		
		val_df = smart_df(nInterpret = c("noise","G0,VAF=0","G0,VAF>0",
				"G1","G2,VAF<1","G2,VAF=1"),
			nLabel = c("noise","G=0","G=0","G=1","G=2","G=2"),
			infer_class = seq(6))
		val_df$col_test = factor(val_df$nInterpret)
		nSeqError = smart_df(t(c(opt_list$BIC,opt_list$props,opt_list$maf,opt_list$psi)))
		names(nSeqError) = paste0("n_",c("BIC","p_eps","p_G1","p_G0","p_G2","maf","psi"))
		# nSeqError
		
		tmp_vcf = smart_merge(tmp_vcf,val_df)
		tmp_vcf = tmp_vcf[order(tmp_vcf$orig_order),]
		tmp_vcf = tmp_vcf[,which(!(names(tmp_vcf) %in% c("infer_class","orig_order")))]
		
		# Plot
		png(file.path(nCLUST_dir,sprintf("nCLUST_%s.png",id)),
			units = 'px',height = 1500,width = 3000,res = 250,type = "cairo")
		par(mar = c(5,4,1,6),xpd = TRUE)
		genome_plot(DAT = tmp_vcf,var = "nVAF",how_color = "col_test",hg = hg,
			pch = 16,cex = 0.5*log10(tmp_vcf$nDP+1),genome = genome,main = id)
		legend("topright",inset = c(-0.1,0),legend = val_df$nInterpret,
			pch = 15,pt.cex = 1.5,col = val_df$col_test,title="Cluster")
		vec_DP = 10^seq(0,3)
		legend("bottomright",inset = c(-0.1,0),legend = vec_DP,
			pt.cex = 0.5*log10(vec_DP+1),pch=rep(16,3),title="DP")
		par(mar = c(5,4,4,2)+0.1,xpd = FALSE)
		dev.off()
		
		tmp_vcf = tmp_vcf[,which(!(names(tmp_vcf) %in% c("col_test")))]
		out = rbind(out,tmp_vcf); rm(tmp_vcf,idx_clust,mat_RD)
		SE = rbind(SE,smart_df(STUDYNUMBER = id,nSeqError))
		cnt = cnt + 1
	}
	if( verbose ) cat("\n")
	
	if( ID != "STUDYNUMBER" ) SE = name_change(SE,"STUDYNUMBER",ID)
	
	out = smart_merge(vcf,out[,c(ID,"mutID","nInterpret","nLabel")],all.x = TRUE)
	list(out = out,SE = SE)
}
TO_genotype = function(vcf,outdir,hg,genome,strand,
	rd_thres,qscore_thres,eps_thres = 0.5,psi_thres = 0.02){
	
	if(FALSE){
		ID = "UTseq0547"
		TYPE = "full"
		rd_thres = 10; qscore_thres = 30
		eps_thres = 0.5; psi_thres = 0.02
		
		my_dirs = setdirs()
		outdir = file.path(my_dirs$full_dir,ID,TYPE)
		
		rds_fn = file.path(outdir,"image.rds")
		rds = readRDS(rds_fn)
		names(rds)
		vcf = rds$vcf		
		strand = rds$strand
		genome = rds$genome
		qscore_thres = rds$qscore_thres
		rd_thres = rds$rd_thres
		hg = rds$hg
		rm(rds)
		
	}
	
	# Extract HQ unique loci
	dim(vcf)
	vcf = vcf[which(vcf$tDP >= rd_thres & vcf$Qscore >= qscore_thres),]
	dim(vcf)
	vcf = vcf[!duplicated(vcf$mutID),]
	dim(vcf)
	vcf = smart_rmcols(vcf,"Strand_bias_P")
	vcf = smart_merge(vcf,strand[,c("mutID","Strand_bias_P")])
	dim(vcf)
	vcf = vcf[which(vcf$Strand_bias_P > 0.05),]
	dim(vcf); vcf[1:3,]
	
	# Append H2M and BAF clustering
	nSEG_fn = file.path(outdir,"anno_nSEG.rds")
	nSEG 		= readRDS(nSEG_fn)$all_nVAF_segs
	tSEG_fn = file.path(outdir,"anno_tSEG.rds")
	tSEG 		= readRDS(tSEG_fn)$all_tVAF_segs
	
	dim(vcf)
	vcf = anno_SEG(VC = vcf,nSEG = nSEG,tSEG = tSEG,
		eps_thres = eps_thres,psi_thres = psi_thres)
	dim(vcf); vcf[1:3,]
	
	# Infer tumor genotypes
	req_cols = c("tAD","tRD","T_eps","T_prob_half",
		"T_mBAF","T_prop_mBAF","T_PSI")
	vcf$tGENO = apply(vcf[,req_cols],1,function(xx)
		get_INFER_geno(vec_RD = xx[1:2],params = xx[-c(1,2)]))
	smart_table(vcf$tGENO)
	smart_table(vcf$nANNO)
	vcf$tDP_cex = 0.5 * log10(vcf$tDP)
	
	leg_df = smart_df(LAB = c("N/A","H2M","AA","AB","BB"),
		COL = c("gray","black","deeppink2","green3","deepskyblue"))
	vcf$tGENO_col = leg_df$COL[1]
	vcf$tGENO_col[which(vcf$nANNO == "H2M")] = leg_df$COL[2]
	vcf$tGENO_col[which(vcf$tGENO == "0/0" & vcf$nANNO == "mappable")] = leg_df$COL[3]
	vcf$tGENO_col[which(vcf$tGENO == "0/1" & vcf$nANNO == "mappable")] = leg_df$COL[4]
	vcf$tGENO_col[which(vcf$tGENO == "1/1" & vcf$nANNO == "mappable")] = leg_df$COL[5]
	
	cex_df = smart_df(DP = 10^seq(0,3))
	cex_df$pt_cex = 0.5 * log10(cex_df$DP)
	cex_df
	
	png(file.path(outdir,"tumor_genotype.png"),units = "px",
		height = 2500,width = 5000,res = 250,type = "cairo",pointsize = 20)
	par(mar = c(5,4,1,7),xpd = TRUE)
	genome_plot(DAT = vcf,var = "tVAF",cex = vcf$tDP_cex,
		genome = genome,how_color = "tGENO_col",pch = 1,hg = hg,
		main = "Inferred Genotypes and H2M Loci")
	legend("bottomright",inset = c(-0.1,0),legend = leg_df$LAB,
		col = leg_df$COL,pt.cex = 2,pch = 16,title = "Genotype Key",
		cex = 0.8,bty = "n")
	legend("right",inset = c(-0.1,0),legend = cex_df$DP,
		pt.cex = cex_df$pt_cex,pch = 16,title = "Depth Key",
		cex = 0.8,bty = "n")
	par(mar = c(5,4,4,2)+0.1,xpd = FALSE)
	dev.off()
	
	# Output genotypes
	vcf = smart_rmcols(vcf,c("nVAF","nAD","nRD","nDP",
		"mutID","tGENO_col","STUDYNUMBER","ID","TYPE","x_coor",
		"nInterpret","nLabel","tDP_cex","GENO","GENO_col"))
	# vcf[1:5,]
	
	# Make LABEL
	vcf$LABEL = ""
	vcf$LABEL[which(vcf$LABEL == "" & is.na(vcf$tANNO)
		)] = "no_tumor_seg"
	vcf$LABEL[which(vcf$LABEL == "" & !is.na(vcf$tGENO)
		& vcf$nANNO == "H2M"
		)] = "H2M_region"
	vcf$LABEL[which(vcf$LABEL == "" & !is.na(vcf$tGENO)
		& vcf$nANNO == "mappable"
		& vcf$tGENO %in% c("./.")
		)] = "unknown.genotype"
	vcf$LABEL[which(vcf$LABEL == "" & !is.na(vcf$tGENO)
		& vcf$nANNO == "mappable"
		& vcf$tGENO %in% c("0/0")
		)] = "homo.reference"
	vcf$LABEL[which(vcf$LABEL == "" & !is.na(vcf$tGENO)
		& vcf$nANNO == "mappable"
		& vcf$tGENO %in% c("0/1","1/1")
		& ( vcf$ThsdG_AF > 0 | vcf$EXAC_AF > 5e-5)
		)] = "common.SNP"
	vcf$LABEL[which(vcf$LABEL == "" & !is.na(vcf$tGENO)
		& vcf$nANNO == "mappable"
		& vcf$tGENO %in% c("0/1","1/1")
		& !( vcf$ThsdG_AF > 0 | vcf$EXAC_AF > 5e-5)
		)] = "private.SNP"
	smart_table(vcf[,c("tGENO","nANNO","LABEL")])
	smart_table(vcf$LABEL)
	
	out_fn = file.path(outdir,"tumor_genotype.tsv")
	fwrite(x = vcf,file = out_fn,sep = "\t")
	
	return(NULL)
}

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
#'	to the tumor's BAM file.
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
#' @param ncores A positive integer for the number of threads to 
#'	use for calculating strand-specific read counts.
#' @export
run_UNMASC = function(tumorID,outdir,vcf = NULL,tBAM_fn,bed_centromere_fn,dict_chrom_fn,
	qscore_thres = 30,exac_thres = 5e-3,ad_thres = 5,rd_thres = 10,cut_BAF = 5e-2,
	minBQ = 13,minMQ = 40,eps_thres = 0.5,psi_thres = 0.02,hg = "19",binom = TRUE,
	gender = NA,ncores = 1){
	
	if(FALSE){
		tumorID = tumorID; outdir = outdir; vcf = vcf
		tBAM_fn = tBAM_fn; bed_centromere_fn = my_dirs$bed_cent_fn
		dict_chrom_fn = my_dirs$dict_chrom_fn; binom = binom; gender = gender
		ncores = ncores
		
		qscore_thres = 30; exac_thres = 5e-3; ad_thres = 3; 
		rd_thres = 10; cut_BAF = 5e-2; minBQ = 13; minMQ = 40; 
		eps_thres = 0.5; psi_thres = 0.02; hg = "19"; ncores = 1
		
	}
	
	cat("% ------------------------------- %\n")
	cat("% Welcome to the UNMASC workflow! %\n")
	cat("% ------------------------------- %\n")
	Sys.sleep(1)
	smart_mkdir(outdir)
	
	image_fn = file.path(outdir,"image.rds")
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
		cat(sprintf("%s: Calculate mutID and light filtering ...\n",date()))
		vcf$nDP = rowSums(vcf[,c("nAD","nRD")])
		vcf$tDP = rowSums(vcf[,c("tAD","tRD")])
		vcf = vcf[which(vcf$tDP >= 5 & vcf$nDP >= 5 & vcf$Qscore >= 5),] # lowest tolerated thresholds
		if( nrow(vcf) < 1e3 ){
			cat(sprintf("%s: LowQCSample = Low variant count after base filtering ...\n",date()))
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
		cat(sprintf("%s: Calculate nVAFs and tVAFs ...\n",date()))
		if( !("tVAF" %in% names(vcf)) ){
			vcf$tVAF = ifelse(vcf$tDP == 0,0,vcf$tAD/vcf$tDP)
		}
		if( !("nVAF" %in% names(vcf)) ){
			vcf$nVAF = ifelse(vcf$nDP == 0,0,vcf$nAD/vcf$nDP)
		}
		
		# Get chromosome lengths and centromere info
		cat(sprintf("%s: Import genome files ...\n",date()))
		genome = prep_genome(bed_centromere_fn = bed_centromere_fn,
			dict_chrom_fn = dict_chrom_fn)
		
		# Create xCoor
		cat(sprintf("%s: Calculate xCoor ...\n",date()))
		vcf$x_coor = create_xCoor(genome = genome,x = vcf)$x_coor
		
		# Cluster normal counts
		cat(sprintf("%s: Normal read count clustering ...\n",date()))
		out_nCLUST = run_nCLUST(vcf = vcf,ID = "STUDYNUMBER",
			rd_thres = rd_thres,ad_thres = ad_thres,qscore_thres = qscore_thres,
			hg = hg,binom = binom,outdir = outdir,genome = genome,verbose = TRUE)
		vcf = out_nCLUST$out
		SE = out_nCLUST$SE
		
		# Strand
		strand = new_STRAND(BAM_fn = tBAM_fn,vcs = vcf,
			minBQ = minBQ,minMQ = minMQ,NT = ncores)
		
		# Save image
		cat(sprintf("%s: Save image ...\n",date()))
		saveRDS(list(vcf = vcf,strand = strand,
			genome = genome,SE = SE,qscore_thres = qscore_thres,
			exac_thres = exac_thres,ad_thres = ad_thres,
			rd_thres = rd_thres,gender = gender,hg = hg),image_fn)
	}
	
	cat(sprintf("%s: Import image ...\n",date()))
	rds 		= readRDS(image_fn)
	vcf 		= rds$vcf
	strand 	= rds$strand
	genome 	= rds$genome
	rm(rds)
	
	# Flag low QC samples
	tab_tumor 	= table(CONTROL = vcf$STUDYNUMBER,Tumor_Depth = vcf$tDP)
	tab_normal 	= table(CONTROL = vcf$STUDYNUMBER,Normal_Depth = vcf$nDP)
	if( ncol(tab_tumor) <= 50 ){
		print(tab_tumor)
		cat(sprintf("%s: LowQCSample = Low variability in tumor depth\n",date()))
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
	cat(paste0("nrow of uniq_vcs = ",nrow(uniq_vcs),"\n"))
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


