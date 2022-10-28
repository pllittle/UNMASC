
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
		tmp_em = Rcpp_BAF_clust(RD = MAT[,c("AD","RD"),drop = FALSE],DP = MAT[,"DP"],
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
segment_nVAF = function(DATA,chr,binom = TRUE,
	eps_thres = 0.5,psi_thres = 0.01,show_plot = FALSE,
	show_mess = FALSE){
	
	if(FALSE){
		DATA = vcs; chr = chr; binom = binom
		show_plot = FALSE; show_mess = FALSE
	}
	
	input_data = DATA[which( DATA$Chr %in% paste0("chr",chr)
		& DATA$nLabel %in% c("noise","G=1") ),
		c("Chr","Position","GeneName","TARGET","x_coor","mutID",
			"nAD","nRD","nDP","nVAF"),drop = FALSE]
	names(input_data) = c("Chr","Position","GeneName","TARGET",
		"x_coor","mutID","AD","RD","DP","VAF")
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

###

