# UNMASC modeling ITH

gen_ITH_config = function(){
	
	out = list()
	
	# One subclone aka clone
	out[[1]] = list(matrix(1))
	
	# Two subclones
	tmp_list = list()
	tmp_list[[1]] = matrix(c(1,1,0,1),2,2,byrow = TRUE)
	tmp_list[[2]] = matrix(c(1,1,0,1,1,0),3,2,byrow = TRUE) # first subclone disappears
	tmp_list
	out[[2]] = tmp_list
	
	# Three subclones
	tmp_list = list()
	tmp_list[[1]] = matrix(c(1,1,1,0,1,1,0,0,1),3,3,byrow = TRUE)
	tmp_list[[2]] = matrix(c(1,1,1,0,1,0,0,0,1),3,3,byrow = TRUE)
	tmp_list[[3]] = matrix(c(1,1,1,1,0,0,0,1,1,0,0,1),4,3,byrow = TRUE)
	tmp_list
	out[[3]] = tmp_list
	
	# Four subclones
	tmp_list = list()
	tmp_list[[1]] = matrix(c(1,1,1,1,0,1,1,1,0,0,1,1,0,0,0,1),4,4,byrow = TRUE)
	tmp_list[[2]] = matrix(c(1,1,1,1,0,1,0,0,0,0,1,1,0,0,0,1),4,4,byrow = TRUE)
	tmp_list[[3]] = matrix(c(1,1,1,1,0,1,1,1,0,0,1,0,0,0,0,1),4,4,byrow = TRUE)
	tmp_list
	out[[4]] = tmp_list
	
	out
	return(out)
}
get_SC_parent = function(ITH,SC){
	if(FALSE){
		ITH = gen_ITH_config()
		ITH[[3]]
		ITH = ITH[[2]][[2]]
		ITH
		SC = 2
		
	}
	
	sub_SCs = seq(ncol(ITH))
	num_alloc = colSums(ITH); num_alloc
	sub_SCs = sub_SCs[which(num_alloc < num_alloc[SC])]
	sub_SCs
	
	if( length(sub_SCs) == 0 ) return(0)
	
	vec_dist = apply(ITH[,sub_SCs,drop = FALSE],2,
		function(xx) sum(abs(xx - ITH[,SC])))
	vec_dist
	
	if( length(vec_dist) == 1 ){
		out = sub_SCs
	} else {
		out = sub_SCs[which.min(vec_dist)]
	}
	
	out
	
}
get_ITH_relate = function(ITH){
	out = c()
	for(SC in seq(ncol(ITH))){
		PAR = get_SC_parent(ITH = ITH,SC = SC)
		out = rbind(out,c(SC,PAR))
	}
	
	colnames(out) = c("SC","PAR")
	out = smart_df(out)
	out
	
}
gen_ALLE_states = function(ITH,min_aCN,max_aCN){
	if(FALSE){
		ITH = gen_ITH_config()
		ITH = ITH[[1]][[1]]
		ITH
		
		min_aCN = max_aCN = 1
		
	}
	
	nSC = ncol(ITH); nSC
	CN_states = seq(min_aCN,max_aCN)
	CN_states
	
	if( nSC == 1 ){
		nAllele = as.matrix(CN_states)
		return(nAllele)
	}
	
	# Get subclone-parent relations
	RELA = get_ITH_relate(ITH = ITH); RELA
	
	# Enumerate all copy number states for nA and nB
	nAllele = t(combn(x = rep(CN_states,nSC),m = nSC))
	nAllele = unique(nAllele[do.call(order,
		lapply(1:NCOL(nAllele),function(ii) nAllele[,ii])),,drop = FALSE])
	dim(nAllele); # nAllele
	
	# Remove impossible CN sequence
	idx_keep = apply(nAllele,1,function(xx){
		# xx = nAllele[1,]; xx
		chks_out = TRUE
		for(SC in seq(2,nSC)){
			CN_SC = xx[SC]
			PAR = RELA$PAR[RELA$SC == SC]
			if( PAR == 0 ) next
			CN_PAR = xx[PAR]
			if( CN_PAR == 0 && CN_SC != 0 ){
				chks_out = FALSE
				break
			}
		}
		chks_out
	})
	table(idx_keep)
	# nAllele[idx_keep,][1:10,]
	# nAllele[!idx_keep,][1:10,]
	nAllele = nAllele[idx_keep,,drop = FALSE]
	# dim(nAllele); # nAllele
	
	nAllele
}
get_nA_nB_states = function(ITH,min_aCN,max_aCN){
	if(FALSE){
		ITH = gen_ITH_config()
		ITH = ITH[[1]][[1]]
		ITH
		max_aCN = 3
	}
	
	ALLE_states = gen_ALLE_states(ITH = ITH,
		min_aCN = min_aCN,max_aCN = max_aCN)
	# dim(ALLE_states)
	# if( nrow(ALLE_states) <= 20 ) print(ALLE_states)
	
	# Enumerate all pairs of allelic states
	CN_states = t(combn(x = rep(seq(nrow(ALLE_states)),2),m = 2))
	CN_states = unique(CN_states[do.call(order,
		lapply(1:NCOL(CN_states),function(ii) CN_states[,ii])),,drop = FALSE])
	# dim(CN_states)
	
	list(PAIRS = CN_states,STATES = ALLE_states)
}
exp_BAF = function(N_A,N_B,PURITY,sPROP){
	num = PURITY * sum(N_A * sPROP) + (1 - PURITY)
	den = PURITY * sum((N_A + N_B) * sPROP) + 2 * (1 - PURITY)
	num / den
}
exp_LRR = function(N_A,N_B,PURITY,TAU,sPROP){
	num = PURITY * sum((N_A + N_B) * sPROP) + 2 * (1 - PURITY)
	den = PURITY * TAU + 2 * (1 - PURITY)
	log(num/den)
}
gen_ITH_DAT = function(ITH,NN,PUR,min_aCN = 0,max_aCN,
	SIG2_L,SIG2_B,tot_pts = 1e4,sig2_raw_lrr = 1e-2,
	sig2_raw_baf = 1e-2){
	
	if(FALSE){
		nSC = 1; ITH = gen_ITH_config(); ITH = ITH[[nSC]][[1]]; ITH
		NN = 2e1; PUR = 0.8; min_aCN = 1; max_aCN = 3; SIG2_L = 1e-5; SIG2_B = 1e-5
		
	}
	
	# set.seed(1)
	tSC = ncol(ITH)
	if( tSC > 1 ){
		tPROP = runif(tSC - 1,-2,2)
		tPROP = exp(tPROP) / (1 + sum(exp(tPROP)))
		tPROP = c(tPROP,1 - sum(tPROP))
	} else {
		tPROP = 1
	}
	poss_states = get_nA_nB_states(ITH = ITH,
		min_aCN = min_aCN,max_aCN = max_aCN)
	str(poss_states)
	
	dat = smart_df(SEG = seq(NN),LEN = runif(NN),
		PAIR = sample(nrow(poss_states$PAIRS),NN,replace = TRUE))
	dat$LEN = dat$LEN / sum(dat$LEN)
	# dat# [1:5,]
	
	nAB = t(sapply(dat$PAIR,function(xx){
		# xx = dat$PAIR[1]; xx
		tmp_pair = poss_states$PAIR[xx,]
		nA = poss_states$STATES[tmp_pair[1],]
		nB = poss_states$STATES[tmp_pair[2],]
		c(nA,nB)
	},USE.NAMES = FALSE))
	colnames(nAB) = c(sprintf("nA%s",seq(tSC)),sprintf("nB%s",seq(tSC)))
	# nAB
	
	tmp_mat = smart_df(t(sapply(dat$PAIR,function(pp){
		# pp = dat$PAIR[1]; pp
		tmp_pair = poss_states$PAIRS[pp,]
		nA = sum(poss_states$STATES[tmp_pair[1],] * tPROP)
		nB = sum(poss_states$STATES[tmp_pair[2],] * tPROP)
		tmp_vec = c(nA,nB)
		tmp_vec
	},USE.NAMES = FALSE)))
	names(tmp_mat) = c("nA","nB")
	dat = smart_df(dat,tmp_mat); rm(tmp_mat)
	tPLOIDY = sum(rowSums(dat[,c("nA","nB")]) * dat$LEN) / sum(dat$LEN); tPLOIDY
	TRUTH = list(nAB = nAB,PUR = PUR,PLO = tPLOIDY,scPROP = tPROP); # TRUTH
	
	dat$mean_LRR = apply(dat[,c("nA","nB")],1,function(xx){
		num = PUR * sum(xx) + 2 * (1 - PUR)
		den = PUR * tPLOIDY + 2 * (1 - PUR)
		log(num/den)
	})
	dat$mean_BAF = apply(dat[,c("nA","nB")],1,function(xx){
		num = PUR * xx[2] + (1 - PUR)
		den = PUR * sum(xx) + 2 * (1 - PUR)
		out = num/den
		# out = min(c(out,1-out))
		out
	})
	dat$LRR = rnorm(n = NN,mean = dat$mean_LRR,sd = sqrt(SIG2_L / dat$LEN))
	dat$BAF = rnorm(n = NN,mean = dat$mean_BAF,sd = sqrt(SIG2_B / dat$LEN))
	dat$LRR = dat$LRR - mean(dat$LRR) # center LRR
	dat$BAF = sapply(dat$BAF,function(xx) min(c(xx,1-xx)))
	dat = round(dat,4)
	# dat
	
	# Gen raw LRR/BAF
	win = c()
	for(ii in seq(nrow(dat))){
		# ii = 1
		dat[ii,]
		num_pts = round(dat$LEN[ii] * tot_pts,0)
		tmp_df = smart_df(SEG = ii,
			LRR = rnorm(n = num_pts,mean = dat$mean_LRR[ii],sd = sig2_raw_lrr),
			BAF = NA)
		if( dat$mean_BAF[ii] == 0.5 ){
			tmp_df$BAF = rnorm(n = num_pts,mean = dat$mean_BAF[ii],sd = sig2_raw_baf)
		} else {
			tmp_class = rbinom(n = num_pts,1,0.5)
			tmp_df$BAF[tmp_class == 0] = rnorm(n = sum(tmp_class == 0),
				mean = dat$mean_BAF[ii],sd = sig2_raw_baf)
			tmp_df$BAF[tmp_class == 1] = rnorm(n = sum(tmp_class == 1),
				mean = 1 - dat$mean_BAF[ii],sd = sig2_raw_baf)
		}
		win = rbind(win,tmp_df)
	}
	win$idx = seq(nrow(win)) / nrow(win)
	
	return(list(WIN = win,DAT = dat,
		PS = poss_states,TRUTH = TRUTH))
	
}
gridOpt_ITH = function(ITH,DAT,min_aCN = 0,max_aCN,
	uBAF,sig2_up = 1e-5,vPUR = seq(0,1,0.05),vTAU = seq(1,5,0.1),
	ncores = 1){
	
	if(FALSE){
		nSC = 2; ITH = gen_ITH_config(); ITH = ITH[[nSC]][[1]]; ITH
		DAT = dat; min_aCN = 0; max_aCN = 4
		vPUR = seq(0,1,0.1); vTAU = seq(1,5,0.25)
		ncores = 1
		
	}
	
	# Fix and optimize
	nSC = ncol(ITH); nSC
	PS = get_nA_nB_states(ITH = ITH,
		min_aCN = min_aCN,max_aCN = max_aCN)
	print(str(PS))
	
	sPROP = rep(1,nSC); sPROP = sPROP / sum(sPROP); sPROP
	upPARS = rep(1,4 + nSC - 1)
	upPARS[c(1,2)] = 0
	# upPARS
	
	if( nSC == 1 ){
		aPROPS = smart_df(SC1 = 1)
	} else {
		# upPARS[-c(3,4)] = 0
		vec_INTs = c(1,5,9)
		aPROPS = t(combn(x = rep(vec_INTs,nSC),m = nSC))
		aPROPS = smart_df(unique(t(apply(aPROPS,1,function(xx) xx / sum(xx)))))
		aPROPS = aPROPS[do.call(order,aPROPS),]
		rownames(aPROPS) = NULL
		names(aPROPS) = sprintf("SC%s",seq(nSC))
	}
	aPROPS = as.matrix(aPROPS)
	aPROPS
	upPARS
	
	tot = length(vPUR) * length(vTAU)
	cat(sprintf("%s: Grid size = %s ...\n",date(),tot))
	
	if(FALSE){ # R code
		res = c(); cnt = 1
		for(ii in seq(length(vPUR))){
		for(jj in seq(length(vTAU))){
			# ii = 17; jj = 20
			smart_progress(ii = cnt,nn = tot,iter2 = 1e2)
			
			fin_bb = NULL
			for(kk in seq(nrow(aPROPS))){
				# kk = 1
				bb = Rcpp_ITH_LRR_BAF_clust(LRR = DAT$LRR,BAF = DAT$BAF,
					LEN = DAT$LEN,purity0 = vPUR[ii],tau0 = vTAU[jj],
					sigma2_LRR_0 = sig2_up,sigma2_BAF_0 = sig2_up,
					sPROP_0 = as.numeric(aPROPS[kk,]),
					PAIRS = PS$PAIRS - 1,ALLES = PS$STATES,
					upPARS = upPARS,max_iter = 5e2,eps = 1e-4,show = !TRUE)
				# str(bb)
				if( bb$OPT$converge == "yes" ){
					if( is.null(fin_bb) ){
						fin_bb = bb
					} else if( bb$OPT$BIC > fin_bb$OPT$BIC ){
						fin_bb = bb
					}
				}
				rm(bb)
			}
			# str(fin_bb)
			if( is.null(fin_bb) ){
				cnt = cnt + 1
				next
			}
			opt_sPROP = smart_df(t(fin_bb$PARAM$sPROP))
			names(opt_sPROP) = sprintf("SC%s",seq(nSC))
			
			tmp_df = smart_df(PUR = vPUR[ii],TAU = vTAU[jj],
				PLO = fin_bb$PARAM$ploidy,EPS = fin_bb$PARAM$props[1],
				SIG2_L = fin_bb$PARAM$sigma2_LRR,
				SIG2_B = fin_bb$PARAM$sigma2_BAF,
				CONV = fin_bb$OPT$converge,
				nGRAD = Rcpp_norm(fin_bb$OPT$GRAD),
				opt_sPROP,BIC = fin_bb$OPT$BIC)
			# tmp_df
			res = rbind(res,tmp_df); rm(tmp_df,fin_bb)
			cnt = cnt + 1
		}}
		
		res = res[,names(res) != "CONV"]
		# dim(res); res[1:5,]
	}
	
	if( uBAF ){
		BAF = DAT$BAF
	} else {
		BAF = DAT$BAF * 0 - 1
	}
	res = Rcpp_ITH_GRID(LRR = DAT$LRR,BAF = BAF,LEN = DAT$LEN,
		vPUR = vPUR,vTAU = vTAU,sigma2_LRR_0 = sig2_up,sigma2_BAF_0 = sig2_up,
		m_scPROP = aPROPS,PAIRS = PS$PAIRS - 1,ALLES = PS$STATES,
		upPARS = upPARS,max_iter = 5e2,eps = 1e-5,ncores = ncores)
	colnames(res) = c("PUR","TAU","PLO","EPS","S2L","S2B",
		sprintf("SC%s",seq(nSC)),"nGRAD","BIC")
	res = smart_df(res)
	res = res[!duplicated(round(res[,c("PUR","PLO","EPS","BIC")],4)),]
	res = res[order(-res$BIC),]
	rownames(res) = NULL
	
	return(list(RES = res,ITH = ITH,aPROPS = aPROPS,PS = PS))
}
plot_LRR_BAF = function(WIN,DAT,INFER = NULL){
	if(FALSE){
		WIN = win; DAT = dat
		# INFER = NULL
		INFER = bb
		
	}
	
	WIN = WIN[order(WIN$idx),]
	WIN$LRR = WIN$LRR - mean(WIN$LRR)
	if( is.null(INFER) ){
		num_rows = 2
		my_oma = 0
	} else {
		num_rows = 2 + length(INFER$PARAM$sPROP); num_rows
		tmp_class = apply(INFER$OPT$ZZ,1,which.max)
		my_oma = 2
		# stop("code not done")
	}
	
	par(mfrow = c(num_rows,1),mar = c(4,4,2,0.5),oma = c(0,0,my_oma,0))
	
	plot(WIN[,c("idx","LRR")],ylim = 1.2 * c(-1,1) * max(abs(WIN$LRR)),
		col = rgb(0,0,0,0.25),ylab = "Log R Ratio",bty = "n",
		xlab = "Genomic Index",yaxs = "i")
	abline(h = 0,lty = 2,col = "gray")
	cc = apply(DAT[,c("SEG","LRR")],1,function(xx){
		# xx = DAT[1,]
		rr = range(WIN$idx[which(WIN$SEG == xx[1])])
		segments(x0 = rr[1],x1 = rr[2],y0 = xx[2],
			col = "red",lty = 2,lwd = 4)
	})
	if( !is.null(INFER) ){
		cc = sapply(seq(nrow(DAT)),function(xx){
			# xx = 1
			rr = range(WIN$idx[which(WIN$SEG == xx)])
			if( tmp_class[xx] == 1 ) return(NULL)
			mLRR = INFER$INFER$aggINFER[xx,1]
			segments(x0 = rr[1],x1 = rr[2],y0 = mLRR,
				col = "green",lty = 3,lwd = 4)
		})
	}
	
	plot(WIN[,c("idx","BAF")],ylim = c(0,1),ylab = "B Allele Frequency",
		xlab = "Genomic Index",col = rgb(0,0,0,0.25),bty = "n")
	abline(h = c(0,0.5,1),lty = 2,col = "gray")
	cc = apply(DAT[,c("SEG","BAF")],1,function(xx){
		# xx = DAT[1,]
		rr = range(WIN$idx[which(WIN$SEG == xx[1])])
		segments(x0 = rr[1],x1 = rr[2],y0 = xx[2],
			col = "red",lty = 2,lwd = 4)
		segments(x0 = rr[1],x1 = rr[2],y0 = 1 - xx[2],
			col = "red",lty = 2,lwd = 4)
	})
	if( !is.null(INFER) ){
		cc = sapply(seq(nrow(DAT)),function(xx){
			# xx = 1
			rr = range(WIN$idx[which(WIN$SEG == xx)])
			if( tmp_class[xx] == 1 ) return(NULL)
			mBAF = INFER$INFER$aggINFER[xx,2]
			if( mBAF < 0 ) return(NULL)
			segments(x0 = rr[1],x1 = rr[2],y0 = mBAF,
				col = "green",lty = 3,lwd = 4)
			segments(x0 = rr[1],x1 = rr[2],y0 = 1 - mBAF,
				col = "green",lty = 3,lwd = 4)
		})
	}
	
	# Plot inferred CNs
	if( !is.null(INFER) ){
		nSC = length(INFER$PARAM$sPROP)
		max_tCN = max(sapply(seq(nSC),function(xx){
			# xx = 1
			max(INFER$INFER$infer_nA[,xx] + INFER$INFER$infer_nB[,xx])
		},USE.NAMES = FALSE))
		max_tCN
		
		for(sc in seq(nSC)){
			# sc = 1
			plot(0,0,type = "n",bty = "n",ylab = "Copy Number",
				xlab = "Genomic Index",ylim = c(0,max_tCN),xlim = c(0,1),
				main = sprintf("SC%s, Proportion = %s, Ploidy = %s",
					sc,round(INFER$PARAM$sPROP[sc],3),
					round(INFER$PARAM$sPloidy[sc],3)))
			abline(h = seq(0,max_tCN),lty = 2,lwd = 0.8,col = "gray")
			tmp_mat = cbind(INFER$INFER$infer_nA[,sc],
				INFER$INFER$infer_nB[,sc])
			tCN = rowSums(tmp_mat)
			minorCN = apply(tmp_mat,1,min)
			majorCN = apply(tmp_mat,1,max)
			
			cc = sapply(seq(nrow(DAT)),function(xx){
				# xx = 1
				if( tmp_class[xx] == 1 ) return(NULL)
				rr = range(WIN$idx[which(WIN$SEG == xx)])
				# total CN
				segments(x0 = rr[1],x1 = rr[2],y0 = tCN[xx],
					col = "black",lty = 1,lwd = 9)
				mBAF = INFER$INFER$aggINFER[xx,2]
				if( mBAF < 0 ) return(NULL)
				# major CN
				segments(x0 = rr[1],x1 = rr[2],y0 = majorCN[xx],
					col = adjustcolor("blue",alpha = 0.5),lty = 1,lwd = 7)
				# minor CN
				segments(x0 = rr[1],x1 = rr[2],y0 = minorCN[xx],
					col = adjustcolor("red",alpha = 0.5),lty = 1,lwd = 5)
			})
			
		}
		
		mtext(sprintf("BIC = %s, Purity = %s, Ploidy = %s, EPS = %s",
			round(INFER$OPT$BIC,3),round(INFER$PARAM$purity,3),
			round(INFER$PARAM$ploidy,3),round(INFER$PARAM$props[1],3)),
			outer = TRUE,cex = 1.3)
	}
	
	par(mfrow = c(1,1),mar = c(5,4,4,2) + 0.1,oma = rep(0,4))
	
}
plot_ITH_grid = function(RES){
	# RES = res
	
	themes = theme(text = element_text(size = 15),
		panel.background = element_blank(),
		# panel.spacing = unit(0.5,"lines"),
		# panel.border = element_rect(color = "black",fill = NA,size = 1),
		legend.position = c("none","bottom")[2],
		legend.key.width = unit(4,"line"),
		# legend.title = element_text(size = 18),
		legend.text = element_text(size = 15))
	
	gg = ggplot(data = RES,aes(x = PUR,y = TAU)) +
		geom_tile(aes(fill = BIC)) + 
		scale_fill_viridis(discrete = FALSE) + 
		labs(fill = "Bayesian Information Criteria") +
		themes
	
	return(gg)
}
ITH_sim = function(){
	# source(file.path(git_dir,"UNMASC","R","model_ITH.R"))
	
	library(ggplot2); library(viridis)
	
	# Sim data
	set.seed(1)
	tSC = 1; ITH = gen_ITH_config(); ITH = ITH[[tSC]][[1]]; ITH
	sim_out = gen_ITH_DAT(ITH = ITH,NN = 2e1,SIG2_L = 1e-6,
		SIG2_B = 1e-6,PUR = 0.5,min_aCN = 0,max_aCN = 3,
		tot_pts = 2e4,sig2_raw_lrr = 9e-2,sig2_raw_baf = 4e-2)
	dat 	= sim_out$DAT
	TRUTH = sim_out$TRUTH
	win 	= sim_out$WIN
	# str(sim_out)
	plot_LRR_BAF(WIN = win,DAT = dat)
	dat
	
	# Fix and optimize
	sig2_up = 5e-6
	uBAF = TRUE
	nSC = 2; ITH = gen_ITH_config(); ITH = ITH[[nSC]][[1]]; ITH
	# vPUR = seq(0,1,0.1); vTAU = seq(0.5,5,0.25)
	vPUR = seq(0.4,0.6,0.1); vTAU = seq(2,4,0.25)
	iout = gridOpt_ITH(ITH = ITH,DAT = dat,max_aCN = 5,
		uBAF = uBAF,sig2_up = sig2_up,vPUR = vPUR,vTAU = vTAU)
	res = iout$RES; PS = iout$PS
	
	# Plot grid
	res2 = res[seq(min(c(20,nrow(res)))),]
	rownames(res2) = NULL
	plot_ITH_grid(RES = res)
	
	res2$fBIC = res2$fPUR = res2$fPLO = res2$fEPS = res2$fTAU = NA
	# res2$fS2_L = res2$fS2_B = NA
	for(jj in seq(nSC)){
		res2[[sprintf("fSC%s",jj)]] = NA
		# res2[[sprintf("fPLO%s",jj)]] = NA
	}
	for(jj in seq(nSC)){
		# res2[[sprintf("fSC%s",jj)]] = NA
		res2[[sprintf("fPLO%s",jj)]] = NA
	}
	# res2
	upPARS = rep(1,4 + nSC - 1)
	for(ii in seq(nrow(res2))){
		# ii = 1
		smart_progress(ii = ii,nn = nrow(res2),iter = 1,iter2 = 1e1)
		sPROP = as.numeric(res2[ii,sprintf("SC%s",seq(nSC))])
		sPROP = sPROP / sum(sPROP)
		if( uBAF ){
			BAF = dat$BAF
		} else {
			BAF = dat$BAF * 0 - 1
		}
		bb = Rcpp_ITH_LRR_BAF_clust(LRR = dat$LRR,BAF = BAF,
			LEN = dat$LEN,purity0 = res2$PUR[ii],tau0 = res2$TAU[ii],
			sigma2_LRR_0 = sig2_up,sigma2_BAF_0 = sig2_up,sPROP_0 = sPROP,
			PAIRS = PS$PAIRS - 1,ALLES = PS$STATES,
			upPARS = upPARS,max_iter = 5e2,eps = 1e-5,show = !TRUE)
		
		if(FALSE){
			str(bb)
			
			# round(bb$ZZ,5)
			tmp_mat = smart_df(bb$INFER$infer_nA,bb$INFER$infer_nB)
			colnames(tmp_mat) = c(sprintf("EnA%s",seq(nSC)),sprintf("EnB%s",seq(nSC)))
			# tmp_mat
			INFER = smart_df(bb$INFER$aggINFER)
			colnames(INFER) = c("mLRR","mBAF","aCN","tCN")
			INFER$bCN = INFER$tCN - INFER$aCN
			INFER = INFER[,c("mLRR","mBAF","aCN","bCN")]
			INFER = cbind(INFER,tmp_mat)
			round(cbind(dat,INFER),4)
			cbind(tmp_mat,TRUTH$nAB)
			
		}
		
		res2$fBIC[ii] = bb$OPT$BIC
		res2$fPUR[ii] = bb$PARAM$purity
		res2$fPLO[ii] = bb$PARAM$ploidy
		res2$fEPS[ii]	= bb$PARAM$props[1]
		res2$fTAU[ii] = bb$PARAM$tau
		tmp_class = apply(bb$OPT$ZZ,1,which.max); tmp_class
		
		for(jj in seq(nSC)){
			# jj = 2
			res2[[sprintf("fSC%s",jj)]][ii] = bb$PARAM$sPROP[jj]
			tCN = (bb$INFER$infer_nA[,jj] + bb$INFER$infer_nB[,jj])
			tmp_ploidy = sum(tCN[tmp_class != 1] * 
				dat$LEN[tmp_class != 1]) / sum(dat$LEN[tmp_class != 1])
			# tmp_ploidy
			# stop("calculate ploidy per subclone")
			res2[[sprintf("fPLO%s",jj)]][ii] = tmp_ploidy
		}
		# res2$fS2_L[ii] = bb$sigma2_LRR
		# res2$fS2_B[ii] = bb$sigma2_BAF
		
	}
	res2 = res2[!duplicated(round(res2[,c("fPLO","fPUR","fBIC")],4)),]
	rownames(res2) = NULL
	round(res2[order(-res2$fBIC),],3)
	# res2[order(-res2$fBIC),]
	str(TRUTH)
	
	if(FALSE){
		dim(res); res[1:10,]
		
		upPARS = rep(1,4 + nSC - 1)
		upPARS[c(1,2)] = 0; upPARS
		test_res = Rcpp_ITH_GRID(LRR = dat$LRR,BAF = dat$BAF,LEN = dat$LEN,
			vPUR = vPUR,vTAU = vTAU,sigma2_LRR_0 = sig2_up,sigma2_BAF_0 = sig2_up,
			m_scPROP = iout$aPROPS,PAIRS = PS$PAIRS - 1,ALLES = PS$STATES,
			upPARS = upPARS,max_iter = 5e2,eps = 1e-5,ncores = 1)
		colnames(test_res) = c("PUR","TAU","PLO","EPS",sprintf("SC%s",seq(nSC)),"nGRAD","BIC")
		test_res = smart_df(test_res)
		tres2 = test_res[which(test_res$BIC != 0),]
		tres2 = tres2[!duplicated(tres2[,c("PUR","PLO","EPS","BIC")]),]
		tres2 = tres2[order(-tres2$BIC),]
		dim(tres2); tres2[1:10,]
		
		
	}
	
	if( nSC > 1 ){
		res3 = res2[apply(res2[,sprintf("fSC%s",seq(nSC))],1,min) > 0.05,]
		print(round(res3[order(-res3$fBIC),],3))
		# res3[order(-res3$fBIC),]
	}
	
	# Plot inferred result
	ii = 1
	sPROP = as.numeric(res2[ii,sprintf("SC%s",seq(nSC))])
	sPROP = sPROP / sum(sPROP)
	if( uBAF ){
		BAF = dat$BAF
	} else {
		BAF = dat$BAF * 0 - 1
	}
	bb = Rcpp_ITH_LRR_BAF_clust(LRR = dat$LRR,BAF = BAF,
		LEN = dat$LEN,purity0 = res2$PUR[ii],tau0 = res2$TAU[ii],
		sigma2_LRR_0 = sig2_up,sigma2_BAF_0 = sig2_up,sPROP_0 = sPROP,
		PAIRS = PS$PAIRS - 1,ALLES = PS$STATES,upPARS = upPARS,
		max_iter = 5e2,eps = 1e-5,show = !TRUE)
	plot_LRR_BAF(WIN = win,DAT = dat,INFER = bb)
	# str(bb)
	str(TRUTH)
	
	# Check 'entropy' of result
	pp = bb$PARAM$sPROP; pp
	# pp = c(0.1,10); pp = pp / sum(pp); pp
	-sum(pp * log(pp)) # entropy from proportions
	plo = bb$PARAM$sPloidy; plo = plo / sum(plo); plo
	-sum(plo * log(plo))
	
	
}
test_ideas = function(){
	
	ff = function(varp){
		prop = exp(varp)
		prop = prop / (1 + sum(prop))
		c(prop,1-sum(prop))
	}
	
	varp = c(0)
	prop = ff(varp); prop
	QQ = length(prop); QQ
	t(numDeriv::jacobian(ff,varp))
	mm = diag(prop[seq(QQ-1)],1,1) - prop[seq(QQ-1)] %*% t(prop[seq(QQ-1)])
	mm
	mm2 = matrix(c(1,-1),2,1,byrow = TRUE); t(mm2)
	mm %*% t(mm2)
	
}




###

