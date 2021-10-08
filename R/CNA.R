
# ----------
# CNA estimation
# ----------
gen_possCNs2 = function(min_tCN = 0,max_tCN = 5){
	poss_CNs = c()
	for(tCN in seq(min_tCN,max_tCN)){
		# tCN = 1
	for(minCN in seq(0,tCN)){
		# minCN = 0
		BB = minCN
		AA = tCN - minCN
		if( BB > AA ) next
		poss_CNs = rbind(poss_CNs,c(BB,AA))
	}}
	
	colnames(poss_CNs) = c("minCN","majCN")
	return(poss_CNs)
	
}
sim_CNA = function(nseg = 5e1,purity = NULL,
	sig2_LRR = 1e-6,sig2_BAF = 1e-6,max_tCN = 5,
	prop_num = 0.5,cent_tCN = NULL,cex_mult = 3){
	
	if(FALSE){
		nseg = 3e1; purity = 1
		sig2_LRR = 1e-6; sig2_BAF = 1e-6;
		max_tCN = 5; cent_tCN = 2
		prop_num = 0.5; cex_mult = 3
		
	}
	
	# set.seed(1)
	
	# Simulate
	if( is.null(purity) ){
		true_purity = runif(1)
	} else {
		true_purity = purity
	}
	if( is.null(cent_tCN) ){
		cent_tCN = sample(seq(1,5),1)
	}
	
	true_sig2_LRR = sig2_LRR
	true_sig2_BAF = sig2_BAF
	
	poss_CNs = gen_possCNs2(max_tCN = max_tCN)
	LEN = round(runif(nseg),4); LEN = LEN / sum(LEN)
	num_CNs = nrow(poss_CNs)
	prob_CNs = abs(rowSums(poss_CNs) - cent_tCN)
	prob_CNs = exp(-prop_num * prob_CNs)
	prob_CNs = prob_CNs / sum(prob_CNs)
	# cbind(poss_CNs,prob_CNs)
	
	# Underlying truth
	pick_CNs = sample(seq(num_CNs),nseg,
		replace = TRUE,prob = prob_CNs)
	true_CNs = poss_CNs[pick_CNs,]
	if( purity == 0 ){
		true_CNs = true_CNs * 0 + 1
	}
	true_tCN = rowSums(true_CNs)
	true_ploidy = sum( LEN * true_tCN ) / sum( LEN ); # true_ploidy
	true_LRR = log( ( true_tCN * true_purity + 2 * (1 - true_purity) ) /
		( true_ploidy * true_purity + (2 * (1-true_purity) ) ) )
		true_LRR = round(true_LRR,4)
	
	# If purity == 1 and true_tCN = 0, need to subset out these segments
	if( any(is.infinite(true_LRR)) ){
		idx = !is.infinite(true_LRR)
		true_CNs 		= true_CNs[idx,]
		true_LRR 		= true_LRR[idx]
		true_tCN 		= true_tCN[idx]
		LEN					= LEN[idx]
		true_ploidy = sum( LEN * true_tCN ) / sum( LEN )
		nseg = sum(idx)
	}
	
	true_BAF = ( true_CNs[,1] * true_purity + (1 - true_purity) ) /
		( true_tCN * true_purity + 2 * (1 - true_purity) )
		if( true_purity == 1.0 ) true_BAF[true_tCN == 0] = 0
		true_BAF = round(true_BAF,4)
	
	# Simulate Observed
	oLRR = rnorm(n = nseg,mean = true_LRR,sd = sqrt(true_sig2_LRR / LEN))
	oLRR = round(oLRR - median(oLRR),3)
	# summary(oLRR)
	oBAF = rnorm(n = nseg,mean = true_BAF,sd = sqrt(true_sig2_BAF / LEN))
	oBAF = round(oBAF,3)
	oBAF[oBAF < 0] = 0
	oBAF[oBAF > 0.5] = 0.5
	# summary(oBAF)
	
	# Visualize
	par(mfrow = c(3,1),cex.lab = 1.5,cex.axis = 1.2,
		mar = c(4.2,4.2,1,0.5),oma = c(0,0,2,0))
	vec_cex = cex_mult * log10(LEN / quantile(LEN,0.05))
	plot(oLRR,cex = vec_cex,pch = 16,xlab = "Segment",
		ylab = "Log R Ratio")
		abline(h = 0,lty = 2)
	plot(oBAF,cex = vec_cex,pch = 16,ylim = c(-0.1,1.1),
		xlab = "Segment",ylab = "B Allele Freq")
		points(1 - oBAF,cex = vec_cex,pch = 16)
		abline(h = 0.5,lty = 2)
	plot(true_tCN,ylim = range(c(0,true_tCN)),
		cex = vec_cex,pch = 16,
		xlab = "Segment",ylab = "Copy Number")
		points(true_CNs[,2],cex = vec_cex,pch = 1,
			lwd = 2,col = "red")
		points(true_CNs[,1],cex = vec_cex,pch = 1,
			lwd = 2,col = "blue")
	mtext(sprintf("Purity = %s%%; Ploidy = %s",
		round(true_purity*100,1),round(true_ploidy,2)),
		outer = TRUE,cex = 1.3)
	par(mfrow = c(1,1),cex.lab = 1,cex.axis = 1,
		mar = c(5,4,4,2) + 0.1,oma = rep(0,4))
	
	dat = smart_df(SEG = seq(nseg),LEN = LEN,
		LRR = oLRR,BAF = oBAF,minCN = true_CNs[,1],
		tCN = true_tCN,tLRR = true_LRR,tBAF = true_BAF)
	
	poss_CNs = smart_df(poss_CNs,prob = prob_CNs)
	return(list(DAT = dat,PURITY = true_purity,
		PLOIDY = true_ploidy,POSS_CN = poss_CNs))
	
}
opt_CNA = function(LRR,BAF,LEN,PUR,TAU,upPARS,all_CNs,verbose = FALSE){
	if(FALSE){
		PUR = 1; TAU = 1; upPARS = c(0,0,1,1)
	}
	
	upPARS[4] = ifelse(any(BAF >= 0),1,0)
	lout = Rcpp_LRR_BAF_clust(LRR = LRR,BAF = BAF,
		LEN = LEN,purity0 = PUR,tau0 = TAU,
		sigma2_LRR_0 = 1e-4,sigma2_BAF_0 = 1e-4,
		all_CNs = all_CNs,upPARS = upPARS,
		max_iter = 1e2,eps = 1e-5,show = verbose)
	
	# lout[c("converge","BIC","purity","tau","ploidy","props")]
	
	return(lout)
}
find_gOPT = function(MAT,LRR,BAF,LEN,all_CNs,ss = 3){
	if(FALSE){
		MAT = mBIC
		LRR = LRR; BAF = BAF; LEN = LEN
		all_CNs = all_CNs; ss = ss
		
	}
	
	# dim(MAT); MAT[1:5,1:4]
	nc = ncol(MAT)
	nr = nrow(MAT)
	
	neigbor_IDX = function(nr,nc,ii,jj,ss = 1){
		if(FALSE){
			ii = ii2; jj = jj2
			ss = 1
		}
		
		out = c()
		steps = seq(-ss,ss)
		for(ii2 in steps){
		for(jj2 in steps){
			out = rbind(out,c(ii + ii2,jj + jj2))
		}}
		out = out[which(out[,1] >= 1 & out[,2] >= 1
			& out[,1] <= nr & out[,2] <= nc),]
		out
		
	}
	
	res = c(); cnt = 1; tot = nr*nc; # ss = 5
	for(ii in seq(nr)){
	for(jj in seq(nc)){
		# ii = 17; jj = 7
		smart_progress(ii = cnt,nn = tot,iter = 1e1,iter2 = 4e2)
		ii2 = ii
		jj2 = jj
		
		while(TRUE){
			
			# Current BIC
			curr_BIC = MAT[ii2,jj2]; curr_BIC
			
			# Get neigboring indices
			nidx = neigbor_IDX(nr = nr,nc = nc,
				ii = ii2,jj = jj2,ss = ss)
			tmp_BIC = MAT[nidx]
			tmp_df = smart_df(ii = nidx[,1],
				jj = nidx[,2],BIC = tmp_BIC); tmp_df
			max_idx = which.max(tmp_df$BIC)
			new_BIC = tmp_df$BIC[max_idx]; new_BIC
			ii2 = tmp_df$ii[max_idx]; ii2
			jj2 = tmp_df$jj[max_idx]; jj2
			if( curr_BIC == new_BIC ){
				# cat("Finished!\n")
				break
			}
			
		}
		ii2; jj2
		res = rbind(res,smart_df(ii = ii2,
			jj = jj2,BIC = MAT[ii2,jj2]))
		cnt = cnt + 1
	}
		res = unique(res)
	}
	rownames(res) = NULL
	
	# Append purity/tau/ploidy
	upPARS = rep(1,4)
	res$LL = NA; res$PUR = NA; res$TAU = NA; res$PLO = NA
	res$SIG2_LRR = NA; res$SIG2_BAF = NA; res$MS = NA; res$EPS = NA
	res$converge = NA
	rdnum = 4
	for(aa in seq(nrow(res))){
		# aa = 6
		PUR = rownames(MAT)[res$ii[aa]]
		PUR = as.numeric(strsplit(PUR," = ")[[1]][2])
		PUR
		
		TAU = colnames(MAT)[res$jj[aa]]
		TAU = as.numeric(strsplit(TAU," = ")[[1]][2])
		TAU
		
		lout = opt_CNA(LRR = LRR,BAF = BAF,LEN = LEN,
			PUR = PUR,TAU = TAU,upPARS = upPARS,
			all_CNs = all_CNs)
		
		res$LL[aa] = round(lout$LL,rdnum)
		res$BIC[aa] = round(lout$BIC,rdnum)
		res$PUR[aa] = round(lout$purity,rdnum)
		res$TAU[aa] = round(lout$tau,rdnum)
		res$PLO[aa] = round(lout$ploidy,rdnum)
		res$SIG2_LRR[aa] = round(lout$sigma2_LRR,rdnum*2)
		res$SIG2_BAF[aa] = round(lout$sigma2_BAF,rdnum*2)
		res$MS[aa] = lout$ms
		res$EPS[aa] = lout$props[1]
		res$converge[aa] = lout$converge
		
	}
	
	# res = res[which(res$PUR > 0),]
	# res = res[which(res$PLO > 0),]
	res = res[which(res$converge == "yes"),]
	ucols = c("LL","MS","BIC","EPS","PUR","TAU","PLO",
		"SIG2_LRR","SIG2_BAF")
	res = unique(res[,ucols])
	res = res[order(-res$BIC,res$PUR,res$PLO),]
	rownames(res) = NULL
	res$OPT = seq(nrow(res))
	res = res[,c("OPT",ucols)]
	
	return(res)
}
opt_CNA_grid = function(LRR,BAF,LEN,all_CNs,vPUR,vTAU,ss = 3,
	ncores = 1,verbose = TRUE){
	
	if(FALSE){
		LRR = DATA$fLRR; BAF = tmp_BAF
		LEN = DATA$LEN; all_CNs = all_CNs
		vPUR = vPUR; vTAU = vTAU; ss = 3; ncores = ncores
		verbose = TRUE
		
	}
	
	mBIC = Rcpp_opt_CNA_grid(LRR = LRR,BAF = BAF,LEN = LEN,
		all_CNs = all_CNs,vPUR = vPUR,vTAU = vTAU,ncores = ncores,
		show = verbose)
	rownames(mBIC) = paste0("PUR = ",smart_digits(vPUR,digits = 3))
	colnames(mBIC) = paste0("TAU = ",smart_digits(vTAU,digits = 3))
	
	mBIC[is.infinite(mBIC)] = NA
	min_mBIC = min(mBIC,na.rm = TRUE); min_mBIC
	mBIC[is.na(mBIC)] = min_mBIC
	
	if(FALSE){
		bb = smart_heatmap(MAT = mBIC,DIST = FALSE,clustRC = rep(FALSE,2)); rm(bb)
		sort(c(mBIC),decreasing = TRUE)[1:30]
	}
	
	res = find_gOPT(MAT = mBIC,
		LRR = LRR,BAF = BAF,LEN = LEN,
		all_CNs = all_CNs,ss = ss)
	
	list(RES = res,mBIC = mBIC)
}
plot_mBIC = function(MAT,GRID = NULL,main = "",x_cex = 0.8,
	y_cex = 0.8,x_pur = 0.6,y_tau = 0){
	
	if(FALSE){
		MAT = GRID_out$mBIC; GRID = NULL; main = ""
		x_cex = x_cex; y_cex = y_cex; x_pur = x_pur
		y_tau = 0
		
	}
	
	mm = matrix(c(4,4,0,1,2,3),3,2,byrow = TRUE); mm
	height = c(0.5,0.5,4)
	if( main == "" ){
		mm = mm[-1,]
		height = c(1,4)
	}
	layout(mm,width = c(1,4),height = height)
	tmp_quant = quantile(c(MAT),0.3)
	MAT[MAT < tmp_quant] = tmp_quant
	
	par(mar = rep(0,4))
	plot(0,0,xlim = c(0,1),ylim = c(0,1),
		type = "n",axes = FALSE,xaxs = "i")
	nc = ncol(MAT)
	text(seq(nc)/nc - 1/(nc*2),rep(y_tau,nc),
		labels = colnames(MAT),
		srt = 90,cex = y_cex,adj = c(0,0.5)[1])
	
	par(mar = rep(0,4))
	plot(0,0,xlim = c(0,1),ylim = c(0,1),
		type = "n",axes = FALSE,yaxs = "i")
	nr = nrow(MAT)
	text(rep(x_pur,nr),seq(nr,1)/nr - 1/(nr*2),
		labels = rownames(MAT),cex = x_cex,adj = 0)
	
	par(mar = rep(0,4))
	vCOL = gplots::colorpanel(1.5e2,"deepskyblue","white","red")
	smart_image(MAT,col = vCOL,GRID = GRID)
	
	if( main != "" ){
		par(mar = rep(0,4))
		plot(0,0,type = "n",xlim = c(0,1),
			ylim = c(0,1),axes = FALSE,bty = "n")
		text(0.5,0.5,labels = main,cex = 2.5)
	}
	
	par(mfrow = c(1,1),mar = c(5,4,4,2)+0.1)
	
}
UNMASC_CNA_grid = function(DATA,LRR_only = TRUE,max_tCN = 5,
	vPUR = NULL,vTAU = NULL,ncores = 1){
	
	if(FALSE){
		DATA = out; LRR_only = TRUE
		# vPUR = seq(0,1,0.05),vTAU = seq(1,5,0.2),
		vPUR = seq(0.8,1,0.05); vTAU = seq(1,1.5,0.2)
		max_tCN = 6
		ncores = 1
		
	}
	
	# Check colnames
	req_names = c("LEN","fLRR","fBAF")
	if( !all(req_names %in% names(DATA)) ){
		miss_names = req_names[!(req_names %in% names(DATA))]
		stop(sprintf("Missing columns: %s",paste(miss_names,collapse = ", ")))
	}
	
	# Set GRID of parameters
	if( is.null(vPUR) ) vPUR = seq(0,1,0.1)
	if( is.null(vTAU) ) vTAU = seq(0.5,5,0.1)
	
	# Use BAF or not
	# LRR_only = TRUE
	tmp_BAF = DATA$fBAF
	if( LRR_only ){
		tmp_BAF = DATA$fBAF * 0 - 1
		cat("Use LRR only!\n")
	} else {
		cat("Use LRR and BAF!\n")
	}
	
	# Define LEN
	DATA$LEN = DATA$LEN / sum(DATA$LEN)
	
	# Grid and Optimize
	all_CNs = gen_possCNs2(max_tCN = max_tCN)
	out1 = opt_CNA_grid(LRR = DATA$fLRR,BAF = tmp_BAF,
		LEN = DATA$LEN,all_CNs = all_CNs,
		vPUR = vPUR,vTAU = vTAU,ss = 3,ncores = ncores)
	
	if( LRR_only ){
		tmp_RES = out1$RES[,c("OPT","BIC","EPS","PUR","TAU","PLO","SIG2_LRR")]
		dvars = c("BIC","PUR","PLO","SIG2_LRR")
	} else {
		tmp_RES = out1$RES[,c("OPT","BIC","EPS","PUR","TAU","PLO","SIG2_LRR","SIG2_BAF")]
		dvars = c("BIC","PUR","PLO","SIG2_LRR","SIG2_BAF")
	}
	tmp_RES = tmp_RES[!duplicated(tmp_RES[,dvars]),]
	tmp_RES = tmp_RES[which(tmp_RES$TAU <= max(vTAU)),]
	rownames(tmp_RES) = NULL; tmp_RES$OPT = seq(nrow(tmp_RES))
	
	DATA$fBAF = tmp_BAF
	return(list(all_CNs = all_CNs,mBIC = out1$mBIC,
		RES = tmp_RES,DATA = DATA))
	
}
UNMASC_CNA_showGRID = function(GRID_out,x_cex = 0.4,
	y_cex = 0.7,x_pur = 0.6){
	if(FALSE){
		GRID_out = GRID
		x_cex = 0.7; y_cex = 0.7; x_pur = 0.5
		
	}
	
	plot_mBIC(MAT = GRID_out$mBIC,GRID = NULL,
		x_cex = x_cex,y_cex = y_cex,x_pur = x_pur)
	
}
UNMASC_CNA_opt = function(GRID_out,jj = 1){
	all_CNs = GRID_out$all_CNs
	DATA 		= GRID_out$DATA
	RES 		= GRID_out$RES
	
	fout = opt_CNA(LRR = DATA$fLRR,BAF = DATA$fBAF,
		LEN = DATA$LEN,PUR = RES$PUR[jj],
		TAU = RES$TAU[jj],upPARS = rep(1,4),
		all_CNs = all_CNs); # names(fout)
	# fout[c("purity","ploidy","props","sigma2_LRR","sigma2_BAF","ms","BIC")]
	fout$INFER = smart_df(fout$INFER)
	colnames(fout$INFER) = c("mLRR","mBAF","inf_minCN","inf_tCN")
	fout
	
}
UNMASC_SEG_CNA = function(DATA,GRID,SOLU,tumorID,normalID){
	cat(sprintf("%s: Solution = %s ...\n",date(),SOLU))
	one_res = UNMASC_CNA_opt(GRID_out = GRID,jj = SOLU)
	# one_res[c("BIC","purity","ploidy","tau","props")]
	# one_res$INFER[1:10,]
	
	tmp_infer = cbind(GRID$DATA[,c("Chr","Start","End","fLRR","LEN")],one_res$INFER)
	tmp_infer$inf_minCN[tmp_infer$mBAF == 0] = 0
	# tmp_infer = tmp_infer[which(tmp_infer$inf_tCN != -1),]
	max_range = max(abs(tmp_infer$mLRR)) * 2; max_range
	max_range = max(c(1,max_range)); max_range
	# tmp_infer[tmp_infer$Chr %in% paste0("chr",c(5,6,8,9)),]
	
	# Plot LRR, segments, copy number
	plot(DATA$fLRR,# bty = "n",
		xlim = range(DATA$IDX),ylim = c(-1,1) * max_range,
		col = ifelse(DATA$Chr %in% paste0("chr",seq(1,22,2)),"black","slategray"),
		main = sprintf("T = %s; N = %s; BIC = %s; EPS = %s; PUR = %s; PLO = %s",
			tumorID,normalID,round(one_res$BIC,3),
			round(one_res$props[1],3),
			round(one_res$purity,3),round(one_res$ploidy,3)),
		xlab = "Index",ylab = "LRR",axes = !FALSE)
	aa = sapply(unique(DATA$Chr),function(chr){
		tmp_xx = max(DATA$IDX[DATA$Chr == chr])
		segments(x0 = tmp_xx,y0 = -1,y1 = 1,
			lwd = 0.5,lty = 2)})
	aa = sapply(unique(DATA$Chr),function(chr){
		xx = median(DATA$IDX[DATA$Chr == chr])
		text(x = xx,y = max_range,labels = gsub("chr","",chr),
			cex = 0.75,col = "blue")})
	tmp_df = unique(tmp_infer[,c("mLRR","inf_tCN")])
	tmp_df = tmp_df[which(tmp_df$inf_tCN != -1),]; tmp_df
	axis(4,col = "magenta",at = tmp_df$mLRR,cex.axis = 1,
		labels = tmp_df$inf_tCN,pos = max(DATA$IDX) + 1e2,las = 1)
	for(ii in seq(nrow(tmp_infer))){
		# ii = 17
		tmp_infer[ii,]
		aSEG = DATA$IDX[DATA$Chr == tmp_infer$Chr[ii] & DATA$Start == tmp_infer$Start[ii]]
		bSEG = DATA$IDX[DATA$Chr == tmp_infer$Chr[ii] & DATA$End == tmp_infer$End[ii]]
		# c(aSEG = aSEG,bSEG = bSEG)
		segments(x0 = aSEG,y0 = tmp_infer$fLRR[ii],
			x1 = bSEG,lwd = 3,lty = 3,col = "red")
		if( tmp_infer$inf_tCN[ii] != -1 ){
			segments(x0 = aSEG,y0 = tmp_infer$mLRR[ii],
				x1 = bSEG,lwd = 3,lty = 1,col = "green")
		}
	}
	legend("bottomleft",legend = c("seg.mean","pred.mean"),
		col = c("red","green"),pt.cex = 2,pch = 19,bty = "n")
	
	return(one_res)
}
opt_final = function(){
	# Simulate truth
	set.seed(1); sim = sim_CNA(nseg = 1e1,
		sig2_LRR = 1e-6,sig2_BAF = 1e-6,
		purity = 0.25,prop_num = 0.5,cent_tCN = 3)
	# names(sim)
	sim[c("PURITY","PLOIDY")]
	# sim$POSS_CN
	smart_table(sim$DAT$tCN)
	smart_table(sim$DAT[,c("minCN","tCN")])
	# sim$DAT[1:5,]
	
	sim$DAT$fLRR = sim$DAT$LRR
	sim$DAT$fBAF = sim$DAT$BAF
	
	GRID = UNMASC_CNA_grid(DATA = sim$DAT,
		LRR_only = !TRUE,max_tCN = 6,
		vPUR = seq(0,1,by = 0.025),
		vTAU = seq(0.5,5,0.05))
	names(GRID)
	GRID$RES[seq(min(c(10,nrow(GRID$RES)))),]
	
	# See GRID
	UNMASC_CNA_showGRID(GRID_out = GRID,
		x_cex = 0.7,x_pur = 0.5,y_cex = 0.5)
	
	# See top result
	one_res = UNMASC_CNA_opt(GRID_out = GRID,jj = 1)
	one_res[c("BIC","props","purity","ploidy",
		"sigma2_LRR","sigma2_BAF","INFER")]
	
}



###