# LRR/BAF processing code

# ----------
# LRR pre-processing
# ----------
CNT_correction = function(DATA,TARG){
	
	# Import CNTs
	smart_reqNames(DATA = DATA,REQ = c("Chr","Start","End","GC","CNT"))
	DATA$log_CNT = log(1 + DATA$CNT)
	DATA$midwindow = (DATA$Start + DATA$End) / 2
	DATA$Start2 = DATA$Start / 1e6
	
	# Import target regions
	TARG = TARG[,1:3]
	names(TARG) = c("Chr","Start","End")
	TARG$midpos = (TARG$Start + TARG$End) / 2
	
	# Polish
	DATA = DATA[which(DATA$CNT > 0 & DATA$Chr %in% paste0("chr",c(1:22))),]
	int_contigs = intersect(unique(DATA$Chr),unique(TARG$Chr))
	DATA = DATA[DATA$Chr %in% int_contigs,]
	TARG = TARG[TARG$Chr %in% int_contigs,]
	DATA$ChrCol = ifelse(DATA$Chr %in% paste0("chr",
		c(seq(1,22,2),"X")),"black","slategray")
	DATA$TGT = DATA$nTGT = DATA$DT = NA; nn = nrow(DATA)
	for(ii in seq(nn)){
		# ii = 1
		smart_progress(ii = ii,nn = nn,iter = 1e2,iter2 = 2e3)
		idx1 = TARG$Chr == DATA$Chr[ii]
		idx2 = which(idx1 & TARG$Start >= DATA$Start[ii]
			& TARG$End <= DATA$End[ii])
		# table(idx2)
		if( length(idx2) > 0 ){
			DATA$TGT[ii] = "blue"
		} else {
			DATA$TGT[ii] = "black"
		}
		DATA$nTGT[ii] = length(idx2)
		
		tmp_diff = DATA$midwindow[ii] - TARG$midpos[idx1]
		DATA$DT[ii] = tmp_diff[which.min(abs(tmp_diff))]
		rm(idx1,idx2)
	}
	DATA$ncex = 0.25 * log10(DATA$CNT + 5)
	DATA$TGT_CHR_col = factor(sprintf("%s_%s",DATA$ChrCol,DATA$TGT))
	DATA$RANK = seq(nrow(DATA))
	DATA$DT = DATA$DT / 1e6
	DATA$GC2 = (DATA$GC - mean(DATA$GC))^2
	DATA$GC3 = (DATA$GC - mean(DATA$GC))^3
	DATA$DT2 = (DATA$DT - mean(DATA$DT))^2
	DATA$DT3 = (DATA$DT - mean(DATA$DT))^3
	# DATA$pBC2 = (DATA$pBC - mean(DATA$pBC))^2
	# DATA$pBC3 = (DATA$pBC - mean(DATA$pBC))^3
	DATA$TGT_2 = ifelse(DATA$TGT == DATA$TGT[1],0,1)
	DATA$nTGT_col = sapply(DATA$nTGT,function(xx) ifelse(xx <= 5,xx,5) + 1)
	
	# Regress out batch effects
	gout = lm(log_CNT ~ (nTGT
		# + pBC + pBC2 + pBC3 # This is confounded by total copy number
		+ GC + GC2 + GC3
		+ DT + DT2 + DT3
		)^3,data = DATA)
	print(summary(gout))
	DATA$RESI = as.numeric(resid(gout))
	
	return(DATA)
	
}
CNT_combine = function(tCNT,nCNT,SHOW_PLOT = TRUE){
	# Combine
	nCNT$windowID = sprintf("%s:%s-%s",nCNT$Chr,nCNT$Start,nCNT$End)
	tCNT$windowID = sprintf("%s:%s-%s",tCNT$Chr,tCNT$Start,tCNT$End)
	int_ids = intersect(nCNT$windowID,tCNT$windowID); length(int_ids)
	nCNT = nCNT[nCNT$windowID %in% int_ids,]
	tCNT = tCNT[tCNT$windowID %in% int_ids,]
	tCNT = tCNT[match(nCNT$windowID,tCNT$windowID),]
	if( !all(tCNT$windowID == nCNT$windowID) ) stop("Check ordering")
	
	DAT = nCNT[,c("Chr","Start","End","ChrCol",
		"TGT_CHR_col","TGT","nTGT","DT","CNT","RESI")]
	DAT = name_change(DAT,"RESI","nRESI")
	DAT = name_change(DAT,"CNT","nCNT")
	DAT$Start2 = DAT$Start / 1e6
	DAT$tRESI = tCNT$RESI
	DAT$tCNT = tCNT$CNT
	DAT$LRR = DAT$tRESI - DAT$nRESI
	DAT$LRR = DAT$LRR - median(DAT$LRR)
	DAT$RANK = seq(nrow(DAT)) / (nrow(DAT) + 1)
	DAT[1:3,]
	
	if( SHOW_PLOT ){
		lrr_range = quantile(abs(DAT$LRR),0.999); lrr_range
		plot(DAT[,c("RANK","LRR")],col = DAT$TGT_CHR_col,
			xlim = range(DAT$RANK),ylim = c(-1,1) * lrr_range,
			xaxs = "i")
		aa = sapply(unique(DAT$Chr),function(chr){
			xx = median(DAT$RANK[DAT$Chr == chr])
			yy = -lrr_range * 0.9
			yy = yy + ifelse(chr %in% sprintf("chr%s",seq(1,22,2)),
				0.1,-0.1)
			text(x = xx,y = yy,
				labels = gsub("chr","",chr),
				font = 2,cex = 1,col = "blue")})
	}
	
	return(DAT)
	
}
prep_LRR = function(tCNT_fn,nCNT_fn,TARG_fn){
	
	rds_fn = strsplit(tCNT_fn,"/")[[1]]
	rds_fn = paste(rds_fn[-length(rds_fn)],collapse = "/")
	rds_fn = sprintf("%s/LRR.rds",rds_fn)
	if( file.exists(rds_fn) ){
		return(readRDS(rds_fn))
	}
	
	tcnts = data.table::fread(tCNT_fn,
		sep = "\t",header = FALSE,data.table = FALSE)
	names(tcnts) = c("Chr","Start","End","GC","CNT","BC","WIDTH","COV")
	# dim(tcnts); tcnts[1:5,]
	
	ncnts = data.table::fread(nCNT_fn,
		sep = "\t",header = FALSE,data.table = FALSE)
	names(ncnts) = c("Chr","Start","End","GC","CNT","BC","WIDTH","COV")
	# dim(ncnts); ncnts[1:5,]
	
	# Get target
	targ = data.table::fread(TARG_fn,
		data.table = FALSE,sep = "\t",header = TRUE)
	targ = targ[,1:3]
	names(targ) = c("Chr","Start","End")
	# dim(targ); targ[1:5,]
	
	# Get GC corrected log(CNTs)
	tcnts = CNT_correction(DATA = tcnts,TARG = targ)
	ncnts = CNT_correction(DATA = ncnts,TARG = targ)
	cnts 	= CNT_combine(tCNT = tcnts,nCNT = ncnts,
		SHOW_PLOT = FALSE)
	
	saveRDS(cnts,rds_fn)
	return(cnts)
	
}


# ----------
# LRR Segmentation
# ----------
show_SEG = function(DATA,cSEG = NULL,wSEG = NULL,BICs = NULL){
	# Check colnames
	req_names = c("ChrCol","ncex","LRR")
	smart_reqNames(DATA = DATA,REQ = req_names)
	
	par(mfrow = c(1,1),mar = c(4,4,0.4,0) + 0.1)
	ylims = c(-1,1) * max(abs(DATA$LRR))
	plot(DATA$LRR,col = DATA$ChrCol,ylab = "LRR",bty = "n",
		cex = DATA$ncex,ylim = ylims)
	abline(h = 0,lty = 2)
	if( !is.null(cSEG) ) abline(v = cSEG,lwd = 1,lty = 2,col = "blue")
	if( !is.null(wSEG) ){
		if( wSEG == 1 ){
			aSEG = 1; bSEG = cSEG[2]
		} else if( wSEG < length(cSEG) - 1 ){
			aSEG = cSEG[wSEG]; bSEG = cSEG[wSEG+1]
		} else {
			stop("update while loop code")
		}
		abline(v = c(aSEG,bSEG),lwd = 2,lty = 2,col = "darkgreen")
	}
	if( !is.null(BICs) ){
		for(ii in seq(nrow(BICs))){
			segments(x0 = BICs$aSEG[ii],y0 = BICs$MU[ii],
				x1 = BICs$bSEG[ii],col = "magenta",lwd = 5)
		}
	}
	par(mfrow = c(1,1),mar = c(5,4,4,2) + 0.1)
	
}
mod_SEG_LRR = function(sDATA){
	if(FALSE){
		sDATA = DATA[which(DATA$IDX >= aSEG
			& DATA$IDX <= bSEG),]
		sDATA[1:5,]
		
	}
	
	fcout = NULL
	
	for(iPROP in c(0,1e-2,1e-1)){
		# iPROP = c(0,1e-2,1e-1)[1]; iPROP
		
		cout = Rcpp_LRR_clust(LRR = sDATA$LRR,iPROP = iPROP,
			iMU = median(sDATA$LRR),iSIG2 = mad(sDATA$LRR),
			max_iter = 5e2,eps = 1e-5,show = FALSE)
		cout$converge
		
		if( cout$converge == "yes" ){
			if( is.null(fcout) ){
				fcout = cout
			} else if( fcout$BIC < cout$BIC ){
				fcout = cout
			}
		}
		
	}
	
	fcout[c("prop","LL","ms","PARS")]
}
calc_BIC = function(DATA,cSEG,BICs = NULL){
	if(FALSE){
		DATA = DATA
		cSEG = SEG; cSEG
		BICs = NULL
		
		cSEG = cSEG_2
		BICs = BICs
		
	}
	
	out = c(); update_BIC = FALSE
	for(ii in seq(length(cSEG)-1)){
		# ii = 1
		
		if( ii == 1 ){
			aSEG = 1
		} else {
			aSEG = cSEG[ii] + 1
		}
		bSEG = cSEG[ii + 1]
		# c(aSEG = aSEG,bSEG = bSEG)
		
		if( !is.null(BICs) ){
			# Search, if available, get row, else optimize
			idx = which(BICs$aSEG == aSEG
				& BICs$bSEG == bSEG)
			if( length(idx) == 1 ){
				out = rbind(out,BICs[idx,])
				next
			} 
		}
		
		cout = mod_SEG_LRR(sDATA = DATA[which(DATA$IDX >= aSEG
			& DATA$IDX <= bSEG),])
		cout
		out = rbind(out,smart_df(aSEG = aSEG,
			bSEG = bSEG,ms = cout$ms,LL = cout$LL,
			prop = cout$prop,MU = cout$PARS[1],
			SIG2 = cout$PARS[2]))
		update_BIC = TRUE
	}
	out
	
	if( update_BIC ){
		BICs = rbind(BICs,out)
		BICs = unique(BICs)
	}
	
	BIC = 2 * sum(out$LL) - log(nrow(DATA)) * sum(out$ms)
	list(BIC = BIC,BICs = BICs)
}
add_SEG = function(DATA,wSEG,cSEG,BICs = NULL,minSL = 20,verbose = TRUE){
	if(FALSE){
		DATA = DATA
		wSEG = wSEG
		cSEG = SEG
		BICs = BICs
		minSL = minSL
		verbose = verbose2
		
	}
	
	oBIC = calc_BIC(DATA = DATA,cSEG = cSEG,BICs = BICs)
	BICs = oBIC$BICs
	oBIC$BIC
	
	# For a given SEG, try to add a segment
	if( wSEG == 1 ){
		aSEG = 1; bSEG = cSEG[2]
	} else {
		aSEG = cSEG[wSEG] + 1; bSEG = cSEG[wSEG + 1]
	}
	c(aSEG = aSEG,bSEG = bSEG)
	
	# Create candidate grid new break points
	new_breaks = round(seq(aSEG,bSEG,length.out = 20))
	new_breaks = new_breaks[-c(1,length(new_breaks))]
	new_breaks
	
	# Loop over each to find optimal break
	out = c()
	for(new_break in new_breaks){
		# new_break = new_breaks[1]; new_break
		if( verbose ) cat(".")
		cSEG_2 = sort(c(cSEG,new_break)); cSEG_2
		
		# Check newSEG min segment width sufficient
		if( min(diff(cSEG_2)) < minSL ){
			out = rbind(out,smart_df(SEG = new_break,BIC = -Inf))
			next
		}
		
		oBIC_2 = calc_BIC(DATA = DATA,cSEG = cSEG_2,BICs = BICs)
		BICs = oBIC_2$BICs
		out = rbind(out,smart_df(SEG = new_break,BIC = oBIC_2$BIC))
	}
	if( verbose ) cat("\n")
	out
	
	oBIC$BIC; max(out$BIC)
	if( oBIC$BIC < max(out$BIC) ){
		new_break = new_breaks[which.max(out$BIC)]
		out_SEG = sort(c(cSEG,new_break))
		ADD = TRUE
		if( verbose ) cat("Adding segment!\n")
	} else {
		out_SEG = cSEG
		ADD = FALSE
		if( verbose ) cat("Same segment.\n")
	}
	
	return(list(SEG = out_SEG,
		BICs = BICs,ADD = ADD))
}
rem_SEG = function(DATA,cSEG,BICs = NULL,verbose = TRUE){
	if(FALSE){
		DATA = DATA
		cSEG = SEG; cSEG
		verbose = TRUE
		
	}
	
	oBIC = calc_BIC(DATA = DATA,cSEG = cSEG,BICs = BICs)
	BICs = oBIC$BICs
	oBIC$BIC
	
	# Create candidate grid new break points
	new_breaks = cSEG
	new_breaks = new_breaks[-c(1,length(new_breaks))]
	if( length(new_breaks) == 0 ){
		return(list(SEG = cSEG,
			BICs = BICs))
	}
	new_breaks
	
	# Loop over each to find optimal break to remove
	out = c()
	for(new_break in new_breaks){
		# new_break = new_breaks[1]; new_break
		if( verbose ) cat(".")
		cSEG_2 = sort(cSEG[cSEG != new_break]); cSEG_2
		oBIC_2 = calc_BIC(DATA = DATA,cSEG = cSEG_2,BICs = BICs)
		BICs = oBIC_2$BICs
		out = rbind(out,smart_df(SEG = new_break,BIC = oBIC_2$BIC))
		if( oBIC_2$BIC > oBIC$BIC ) break
	}
	if( verbose ) cat("\n")
	out
	
	oBIC$BIC; max(out$BIC)
	if( oBIC$BIC < max(out$BIC) ){
		new_break = new_breaks[which.max(out$BIC)]
		out_SEG = sort(cSEG[cSEG != new_break])
		if( verbose ){
			cat("Removing segment!\n")
		}
		return(rem_SEG(DATA = DATA,cSEG = out_SEG,
			BICs = BICs,verbose = verbose))
	} else {
		out_SEG = cSEG
		if( verbose ) cat("Same segment.\n")
	}
	
	return(list(SEG = out_SEG,
		BICs = BICs))
}
adj_SEG = function(DATA,wSEG,cSEG,BICs = NULL,minSL = 20,verbose = TRUE){
	if(FALSE){
		DATA = DATA
		wSEG = 1; # aka the wSEG-th segment
		cSEG = SEG; cSEG
		minSL = 20
		verbose = TRUE
		
	}
	
	oBIC = calc_BIC(DATA = DATA,cSEG = cSEG,BICs = BICs)
	BICs = oBIC$BICs
	oBIC$BIC
	
	if( wSEG == length(cSEG) - 1 ){
		return(list(SEG = cSEG,BICs = BICs,ADJ = FALSE))
	}
	
	# Make adjustments for one break
	adj_breaks = seq(-500,500,40)
	adj_breaks = c(adj_breaks,seq(-10,10))
	adj_breaks = sort(unique(adj_breaks))
	preSEG = cSEG[wSEG]; centSEG = cSEG[wSEG + 1]; postSEG = cSEG[wSEG + 2]
	c(preSEG = preSEG,centSEG = centSEG,postSEG = postSEG)
	# centSEG - preSEG
	adj_breaks
	
	out = c()
	for(adj_break in adj_breaks){
		# adj_break = adj_breaks[1]
		if( verbose ) cat(".")
		cSEG_2 = cSEG
		cSEG_2[wSEG + 1] = cSEG_2[wSEG + 1] + adj_break
		if( cSEG_2[wSEG + 1] < cSEG[wSEG] || cSEG_2[wSEG + 1] > cSEG[wSEG + 2] ){
			out = rbind(out,smart_df(ADJ = adj_break,BIC = -Inf))
			next
		}
		cSEG_2 = sort(unique(cSEG_2)); cSEG_2
		
		# Check newSEG min segment width sufficient
		if( min(diff(cSEG_2)) < minSL ){
			out = rbind(out,smart_df(ADJ = adj_break,BIC = -Inf))
			next
		}
		
		oBIC_2 = calc_BIC(DATA = DATA,cSEG = cSEG_2,BICs = BICs)
		BICs = oBIC_2$BICs
		out = rbind(out,smart_df(ADJ = adj_break,BIC = oBIC_2$BIC))
		# if( oBIC_2$BIC > oBIC$BIC ) break
	}
	if( verbose ) cat("\n")
	# out
	
	if( oBIC$BIC < max(out$BIC) ){
		out_SEG = cSEG
		adj_break = out$ADJ[which.max(out$BIC)]; adj_break
		out_SEG[wSEG + 1] = out_SEG[wSEG + 1] + adj_break
		ADJ = TRUE
		if( verbose ) cat("Adjusting segment!\n")
	} else {
		out_SEG = cSEG
		ADJ = FALSE
		if( verbose ) cat("Same segment.\n")
	}
	
	return(list(SEG = out_SEG,
		BICs = BICs,ADJ = ADJ))
	
}
seg_LRR = function(DATA,iSEG = NULL,minSL = 20,
	verbose = TRUE,verbose2 = TRUE,show_plot = TRUE){
	
	if(FALSE){
		DATA = sdat; iSEG = NULL; minSL = minSL
			verbose = TRUE; verbose2 = !TRUE; show_plot = TRUE
		
	}
	
	# Check colnames
	req_names = c("Chr","Start","End","LRR")
	smart_reqNames(DATA = DATA,REQ = req_names)
	
	# Order segments
	chr_names = paste0("chr",c(1:22,"X","Y"))
	int_chrs = intersect(unique(DATA$Chr),chr_names)
	if( length(int_chrs) == 0 ) stop("Check Chr coding")
	chr_names = chr_names[chr_names %in% int_chrs]
	DATA$Chr = factor(as.character(DATA$Chr),levels = chr_names)
	DATA = DATA[order(DATA$Chr,DATA$Start),]
	DATA$Chr = as.character(DATA$Chr)
	DATA$IDX = seq(nrow(DATA))
	
	# Initialize segments
	if( is.null(iSEG) ){
		dim(DATA)
		SEG = round(seq(1,nrow(DATA),length.out = 40))
	} else {
		SEG = sort(unique(iSEG))
	}
	# SEG
	
	# Get initial BICs and BIC
	oBIC = calc_BIC(DATA = DATA,cSEG = SEG,BICs = NULL)
	BICs = oBIC$BICs
	dim(BICs)
	oBIC$BIC
	
	wSEG = 1
	while(TRUE){
		
		wSEG; length(SEG)
		if( wSEG >= length(SEG) - 1 ){
			if( verbose ) cat("Segmentation complete!\n")
			break
		}
		
		# Try remove segment
		# cat("will run rm seg\n")
		if( show_plot ) show_SEG(DATA = DATA,cSEG = SEG,wSEG = wSEG)
		rem_out = rem_SEG(DATA = DATA,cSEG = SEG,
			BICs = BICs,verbose = verbose2)
		# names(rem_out)
		REM = FALSE
		if( length(SEG) != length(rem_out$SEG) ){
			REM = TRUE
		} else if( any(SEG != rem_out$SEG) ){
			REM = TRUE
		}
		SEG = rem_out$SEG; BICs = rem_out$BICs
		rm(rem_out)
		if( REM ){
			wSEG = max(c(1,wSEG - 1))
			cat(sprintf("Shift back segment wSEG = %s ...\n",wSEG))
			next
		}
		
		# Try add segment
		# cat("will run add seg\n")
		if( show_plot ) show_SEG(DATA = DATA,cSEG = SEG,wSEG = wSEG)
		add_out = add_SEG(DATA = DATA,wSEG = wSEG,cSEG = SEG,
			BICs = BICs,minSL = minSL,verbose = verbose2)
		# names(add_out)
		SEG = add_out$SEG; BICs = add_out$BICs
		ADD = add_out$ADD; rm(add_out)
		if( ADD ){
			wSEG = max(c(1,wSEG - 1))
			cat(sprintf("Shift back segment wSEG = %s ...\n",wSEG))
			next
		}
		
		# Try adjusting segment
		# cat("will run adj seg\n")
		if( show_plot ) show_SEG(DATA = DATA,cSEG = SEG,wSEG = wSEG)
		adj_out = adj_SEG(DATA = DATA,wSEG = wSEG,cSEG = SEG,
			BICs = BICs,minSL = minSL,verbose = verbose2)
		SEG = adj_out$SEG; BICs = adj_out$BICs
		ADJ = adj_out$ADJ; rm(adj_out)
		if( ADJ ) next
		
		# Shift to next segment
		wSEG = wSEG + 1
		if( verbose ) cat(sprintf("Shift to next segment wSEG = %s ...\n",wSEG))
		
	}
	
	oBIC = calc_BIC(DATA = DATA,cSEG = SEG,BICs = NULL)
	BICs = oBIC$BICs; BICs
	
	if( show_plot ){
		show_SEG(DATA = DATA,cSEG = SEG,BICs = BICs)
	}
	
	BICs$Start = sapply(BICs$aSEG,function(xx){
		DATA$Start[DATA$IDX == xx]
	})
	BICs$End = sapply(BICs$bSEG,function(xx){
		DATA$End[DATA$IDX == xx]
	})
	# BICs
	
	return(BICs)
	
}


# ----------
# BAF pre-processing
# ----------
BAF_prep = function(BAF,mDP = 10,iPSI = 0){
	if(FALSE){
		BAF = tBAF; mDP = 20; iPSI = 0
		# BAF = BAF[which(BAF$QUAL != "."),]
		# BAF$QUAL = as.numeric(BAF$QUAL)
		
	}
	
	BAF$RefAlt = "INDEL"
	idx_bs = which(nchar(BAF$REF) == 1 & nchar(BAF$ALT) == 1)
	BAF$RefAlt[idx_bs] = sprintf("%s>%s",BAF$REF[idx_bs],BAF$ALT[idx_bs])
	
	if(FALSE){ # Check BAF dist
		BAF[1:3,]
		
		smart_table(BAF$RefAlt)
		smart_hist(BAF$VAF[which(BAF$RefAlt %in% c("C>T","G>A","INDEL"))])
		
		idx = which(BAF$QUAL >= 5)
		plot(BAF$VAF[idx],# col = ifelse(BAF$QUAL <= 5,"red","black"),
			cex = BAF$pt_cex[idx],ylim = c(0,1))
		
	}
	
	# Remove VAF near 1 and artifact cluster
	BAF_clust = function(mat_RD,props0,mm0,psi0){
		DP = rowSums(mat_RD)
		log_DP = log(DP)
		LBC = lchoose(DP,mat_RD[,1])
		props0 = props0 / sum(props0)
		
		bout = Rcpp_BAF_clust(RD = mat_RD,DP = DP,log_DP = log_DP,
			LBC = LBC,props0 = props0,mm0 = mm0,psi0 = psi0,
			max_iter = 4e3,eps = 1e-5,show = FALSE)
		bout
	}
	
	# Remove variant clusters to isolate BAF
	BAF$keepPT = "miss"
	BAF$keepPT[BAF$DP <= mDP] = "lowDP"
	BAF$keepPT[BAF$VAF == 1] = "VAF=1"
	# BAF$keepPT[BAF$QUAL == "."] = "missQUAL"
	# smart_hist(BAF$VAF[is.na(BAF$keepPT)],breaks = 30)
	table(BAF$keepPT)
	
	# Removes spike around VAF = 1.0
	props0 = c(1,0,0,1) # c(eps,half,minBAF,maxBAF)
	mm0 = 1e-3 # minBAF
	idx = which(BAF$keepPT == "miss"
		& BAF$VAF > 0.95)
	length(idx)
	mat_RD = as.matrix(BAF[idx,c("AD","RD")])
	bout = BAF_clust(mat_RD = mat_RD,
		props0 = props0,mm0 = mm0,psi0 = iPSI)
	# str(bout); table(bout$class)
	if( any(bout$class == 4) ){
		BAF$keepPT[idx][bout$class == 4] = "VAF~=1"		
	}
	if(FALSE){
		table(BAF$keepPT)
		idx = which(BAF$keepPT == "miss")
		plot(BAF$VAF[idx],cex = BAF$pt_cex[idx],
			col = factor(BAF$keepPT[idx]),ylim = c(0,1))
		
	}
	
	# Remove variants at low AF
	
	# OXOG signature
	idx = which(BAF$keepPT == "miss"
		& BAF$RefAlt %in% c("C>A","G>T")
		# & BAF$DP >= 20
		# & BAF$DP <= 20
		)
	length(idx)
	if( length(idx) > 3e2 ){
		
		mat_RD = as.matrix(BAF[idx,c("AD","RD")])
		props0 = c(1,1,1,0)
		mm0 = 5e-2
		bout = BAF_clust(mat_RD = mat_RD,
			props0 = props0,mm0 = mm0,psi0 = iPSI)
		
		if(FALSE){
			plot(BAF$VAF[idx],cex = BAF$pt_cex[idx],ylim = c(0,1))
			smart_hist(BAF$VAF[idx],breaks = 50)
			str(bout); table(bout$class)
			plot(BAF$VAF[idx],cex = BAF$pt_cex[idx],pch = 16,
				ylim = c(0,1),col = factor(bout$class))
		}
		
		if( any(bout$class == 3) ){
			ss = t(apply(BAF[idx,c("AD","RD")],1,function(xx)
				BINOM_INFER(vec_RD = xx,bout$params)))
			ss = apply(ss,1,function(xx) ifelse(which.max(xx)
				%in% c(1,2),"miss","signal"))
			ss = as.character(ss)
			table(ss)
			BAF$keepPT[idx] = ss; rm(ss)
			
			idx_pred = which(BAF$keepPT == "miss"
				& !(BAF$RefAlt %in% c("C>A","G>T"))
				# & BAF$DP >= 20
				)
			length(idx_pred)
			
			ss = t(apply(BAF[idx_pred,c("AD","RD")],1,function(xx)
				BINOM_INFER(vec_RD = xx,params = bout$params)))
			ss = apply(ss,1,function(xx) ifelse(which.max(xx)
				%in% c(1,2),"miss","signal_like"))
			ss = as.character(ss)
			table(ss)
			BAF$keepPT[idx_pred] = ss; rm(ss)
		}
		
	}
	
	# FFPE signature
	idx = which(BAF$keepPT == "miss"
		& BAF$RefAlt %in% c("C>T","G>A","INDEL")
		# & BAF$DP >= 20
		)
	length(idx)
	if( FALSE && length(idx) > 1e2 ){
		
		mat_RD = as.matrix(BAF[idx,c("AD","RD")])
		props0 = c(1,1,1,0)
		mm0 = 5e-2
		bout = BAF_clust(mat_RD = mat_RD,
			props0 = props0,mm0 = mm0,psi0 = iPSI)
		
		if(FALSE){
			plot(BAF$VAF[idx],cex = BAF$pt_cex[idx],ylim = c(0,1))
			smart_hist(BAF$VAF[idx],breaks = 50)
			str(bout); table(bout$class)
			plot(BAF$VAF[idx],cex = BAF$pt_cex[idx],pch = 16,
				ylim = c(0,1),col = factor(bout$class))
		}
		
	}
	
	if(FALSE){
		idx = which(BAF$keepPT == "miss")
		# idx = which(BAF$Chr == "chr17")
		plot(BAF$VAF[idx],cex = BAF$pt_cex[idx],
			ylim = c(0,1),col = BAF$ChrCol[idx])
		# smart_hist(BAF$VAF[idx])
		
	}
	
	BAF = BAF[which(BAF$keepPT == "miss"),]
	BAF = smart_rmcols(BAF,c("keepPT","RefAlt"))
	return(BAF)
	
}
mod_SEG_BAF = function(MAT){
	if(FALSE){
		MAT = tmp_mat
		
	}
	
	req_names = c("AD","RD","DP","log_DP","LBC")
	smart_reqNames(DATA = MAT,REQ = req_names)
	num_loci = nrow(MAT)
	
	fin_em = NULL
	init_pi = matrix(c(
		# 1,0,0,0,
		1,1,0,0,
		0,1,0,0,
		1,0,1,0,
		1,0,0,1,
		1,0,1,1),ncol = 4,byrow = TRUE)
	vec_pp0 = c(1/20,seq(4) / 10)
	
	for(pp0 in vec_pp0){
	for(ii in seq(nrow(init_pi))){
		# pp0 = vec_pp0[1]; ii = 3
		tmp_em = Rcpp_BAF_clust(RD = MAT[,c("AD","RD"),drop = FALSE],DP = MAT[,"DP"],
			log_DP = MAT[,"log_DP"],LBC = MAT[,"LBC"],props0 = init_pi[ii,],
			mm0 = pp0,psi0 = 0,max_iter = 2e2,eps = 1e-5,show = FALSE)
		if( tmp_em$converge == "yes" ){
			if( tmp_em$props[1] == 1 ) next
			if( is.null(fin_em) ){
				fin_em = tmp_em
			} else if( tmp_em$BIC > fin_em$BIC ){
				fin_em = tmp_em
			}
		}
	}}
	
	return(fin_em)
}

# ----------
# Joint LRR/BAF segmentation
# ----------
LRR_BAF_combine = function(LRR,BAF,VAR = "IDX"){
	if(FALSE){
		LRR = DAT_L; BAF = DAT_B; VAR = "IDX"
		
	}
	
	lev_contigs = sprintf("chr%s",1:22)
	lev_contigs = lev_contigs[lev_contigs %in% unique(LRR$Chr)]
	LRR$Chr2 = factor(LRR$Chr,levels = lev_contigs)
	LRR = LRR[order(LRR$Chr2,LRR$Start),]
	LRR[[VAR]] = seq(nrow(LRR))
	LRR = smart_rmcols(LRR,c("Chr2","RANK"))
	
	BAF[[VAR]] = NA; nn = nrow(BAF)
	for(ii in seq(nn)){
		# ii = 1
		smart_progress(ii = ii,nn = nn,
			iter = 25,iter2 = 5e2)
		idx = which(LRR$Chr == BAF$Chr[ii]
			& LRR$Start <= BAF$Position[ii]
			& LRR$End >= BAF$Position[ii])
		length(idx)
		if( length(idx) > 0 ){
			BAF[[VAR]][ii] = unique(LRR[[VAR]][idx])
		}
		
	}
	
	# dim(LRR); LRR[1:3,]
	# dim(BAF); BAF[1:3,]
	
	return(list(LRR = LRR,BAF = BAF))
	
}
show_SEG_LB = function(DAT_L,DAT_B,cSEG = NULL,wSEG = NULL){
	if(FALSE){
		DAT_L = DAT_L; DAT_B = DAT_B
		cSEG = SEG; wSEG = wSEG
		
	}
	
	par(mfrow = c(2,1),mar = c(4,4,0.5,0.5))
	tmp_rr = range(DAT_L$IDX)
	plot(DAT_L[,c("IDX","LRR")],col = DAT_L$ChrCol,
		ylab = "LRR",ylim = c(-1,1) * 1.2,xlim = tmp_rr)
	if( !is.null(cSEG) ) abline(v = cSEG,lwd = 1,lty = 2,col = "blue")
	if( !is.null(wSEG) ){
		if( wSEG == 1 ){
			aSEG = 1; bSEG = cSEG[2]
		} else if( wSEG < length(cSEG) - 1 ){
			aSEG = cSEG[wSEG]; bSEG = cSEG[wSEG+1]
		} else {
			stop("update while loop code")
		}
		abline(v = c(aSEG,bSEG),lwd = 2,lty = 2,col = "darkgreen")
	}
	
	plot(DAT_B[,c("IDX","VAF")],col = DAT_B$ChrCol,
		ylab = "BAF",ylim = c(0,1),xlim = tmp_rr,cex = DAT_B$pt_cex)
	if( !is.null(cSEG) ) abline(v = cSEG,lwd = 1,lty = 2,col = "blue")
	if( !is.null(wSEG) ){
		if( wSEG == 1 ){
			aSEG = 1; bSEG = cSEG[2]
		} else if( wSEG < length(cSEG) - 1 ){
			aSEG = cSEG[wSEG]; bSEG = cSEG[wSEG+1]
		} else {
			stop("update while loop code")
		}
		abline(v = c(aSEG,bSEG),lwd = 2,lty = 2,col = "darkgreen")
	}
	
	par(mfrow = c(1,1),mar = c(5,4,4,2) + 0.1)
	
}
mod_SEG_LB = function(sDAT_L,sDAT_B = NULL){
	# Model LRR
	fit_LRR = mod_SEG_LRR(sDATA = sDAT_L)
	# str(fit_LRR)
	
	# Model BAF
	if( !is.null(sDAT_B) ){
		tmp_mat = as.matrix(sDAT_B[,
			c("AD","RD","DP","log_DP","LBC")])
		# class(tmp_mat)
		fit_BAF = mod_SEG_BAF(MAT = tmp_mat)
	} else {
		fit_BAF = list(LL = 0,ms = 0,prop = rep(0,4),
			maf = 0,psi = 0)
	}
	
	list(LRR = fit_LRR,BAF = fit_BAF)
}
calc_BIC_LB = function(DAT_L,DAT_B,cSEG,BICs = NULL){
	if(FALSE){
		DAT_L = DAT_L
		DAT_B = DAT_B
		cSEG = SEG; cSEG
		BICs = NULL
		
		# cSEG = cSEG_2
		# BICs = BICs
		
	}
	
	out = c(); update_BIC = FALSE
	for(ii in seq(length(cSEG)-1)){
		# ii = 1
		
		if( ii == 1 ){
			aSEG = 1
		} else {
			aSEG = cSEG[ii] + 1
		}
		bSEG = cSEG[ii + 1]
		# c(aSEG = aSEG,bSEG = bSEG)
		
		if( !is.null(BICs) ){
			# Search, if available, get row, else optimize
			idx = which(BICs$aSEG == aSEG
				& BICs$bSEG == bSEG)
			if( length(idx) == 1 ){
				out = rbind(out,BICs[idx,])
				next
			} 
		}
		
		# cout = mod_SEG_LRR(sDATA = DATA[which(DATA$IDX >= aSEG
			# & DATA$IDX <= bSEG),])
		sDAT_L = DAT_L[which(DAT_L$IDX >= aSEG
			& DAT_L$IDX <= bSEG),]
		sDAT_B = NULL
		idx_B = which(DAT_B$IDX >= aSEG & DAT_B$IDX <= bSEG)
		if( length(idx_B) > 0 ){
			sDAT_B = DAT_B[idx_B,]
		}
		cout = mod_SEG_LB(sDAT_L = sDAT_L,sDAT_B = sDAT_B)
		
		# str(cout)
		out = rbind(out,smart_df(aSEG = aSEG,bSEG = bSEG,
			ms = cout$LRR$ms + cout$BAF$ms,LL = cout$LRR$LL + cout$BAF$LL,
			nL = nrow(sDAT_L),nB = length(idx_B),
			L_EPS = cout$LRR$prop,L_MU = cout$LRR$PARS[1],L_S2 = cout$LRR$PARS[2],
			B_EPS = cout$BAF$prop[1],B_MAF = cout$BAF$maf,B_PSI = cout$BAF$psi))
		update_BIC = TRUE
	}
	out = round(out,4)
	out
	
	if( update_BIC ){
		BICs = rbind(BICs,out)
		BICs = unique(BICs)
	}
	
	BIC = 2 * sum(out$LL) - log(sum(out$nL) + sum(out$nB)) * sum(out$ms)
	list(BIC = BIC,BICs = BICs)
}
add_SEG_LB = function(DAT_L,DAT_B,wSEG,cSEG,BICs = NULL,minSL = 20,verbose = TRUE){
	if(FALSE){
		wSEG = wSEG
		cSEG = SEG
		BICs = BICs
		minSL = minSL
		verbose = verbose2
		
	}
	
	oBIC = calc_BIC_LB(DAT_L = DAT_L,DAT_B = DAT_B,
		cSEG = cSEG,BICs = BICs)
	BICs = oBIC$BICs
	oBIC$BIC
	
	# For a given SEG, try to add a segment
	if( wSEG == 1 ){
		aSEG = 1; bSEG = cSEG[2]
	} else {
		aSEG = cSEG[wSEG] + 1; bSEG = cSEG[wSEG + 1]
	}
	c(aSEG = aSEG,bSEG = bSEG)
	
	# Create candidate grid new break points
	new_breaks = round(seq(aSEG,bSEG,length.out = 20))
	new_breaks = new_breaks[-c(1,length(new_breaks))]
	new_breaks
	
	# Loop over each to find optimal break
	out = c()
	for(new_break in new_breaks){
		# new_break = new_breaks[1]; new_break
		if( verbose ) cat(".")
		cSEG_2 = sort(c(cSEG,new_break)); cSEG_2
		
		# Check newSEG min segment width sufficient
		if( min(diff(cSEG_2)) < minSL ){
			out = rbind(out,smart_df(SEG = new_break,BIC = -Inf))
			next
		}
		
		oBIC_2 = calc_BIC_LB(DAT_L = DAT_L,DAT_B = DAT_B,
			cSEG = cSEG_2,BICs = BICs)
		BICs = oBIC_2$BICs
		out = rbind(out,smart_df(SEG = new_break,BIC = oBIC_2$BIC))
	}
	if( verbose ) cat("\n")
	out
	
	oBIC$BIC; max(out$BIC)
	if( oBIC$BIC < max(out$BIC) ){
		new_break = new_breaks[which.max(out$BIC)]
		out_SEG = sort(c(cSEG,new_break))
		ADD = TRUE
		if( verbose ) cat("Adding segment!\n")
	} else {
		out_SEG = cSEG
		ADD = FALSE
		if( verbose ) cat("Same segment.\n")
	}
	
	return(list(SEG = out_SEG,BICs = BICs,ADD = ADD))
	
}
rem_SEG_LB = function(DAT_L,DAT_B,cSEG,BICs = NULL,verbose = TRUE){
	if(FALSE){
		# DATA = DATA
		cSEG = SEG; cSEG
		verbose = TRUE
		
	}
	
	oBIC = calc_BIC_LB(DAT_L = DAT_L,DAT_B = DAT_B,
		cSEG = cSEG,BICs = BICs)
	BICs = oBIC$BICs
	oBIC$BIC
	
	# Create candidate grid new break points
	new_breaks = cSEG
	new_breaks = new_breaks[-c(1,length(new_breaks))]
	if( length(new_breaks) == 0 ){
		return(list(SEG = cSEG,
			BICs = BICs))
	}
	new_breaks
	
	# Loop over each to find optimal break to remove
	out = c()
	for(new_break in new_breaks){
		# new_break = new_breaks[1]; new_break
		if( verbose ) cat(".")
		cSEG_2 = sort(cSEG[cSEG != new_break]); cSEG_2
		oBIC_2 = calc_BIC_LB(DAT_L = DAT_L,DAT_B = DAT_B,
			cSEG = cSEG_2,BICs = BICs)
		BICs = oBIC_2$BICs
		out = rbind(out,smart_df(SEG = new_break,BIC = oBIC_2$BIC))
		if( oBIC_2$BIC > oBIC$BIC ) break
	}
	if( verbose ) cat("\n")
	out
	
	oBIC$BIC; max(out$BIC)
	if( oBIC$BIC < max(out$BIC) ){
		new_break = new_breaks[which.max(out$BIC)]
		out_SEG = sort(cSEG[cSEG != new_break])
		if( verbose ){
			cat("Removing segment!\n")
		}
		return(rem_SEG_LB(DAT_L = DAT_L,DAT_B = DAT_B,
			cSEG = out_SEG,BICs = BICs,verbose = verbose))
	} else {
		out_SEG = cSEG
		if( verbose ) cat("Same segment.\n")
	}
	
	return(list(SEG = out_SEG,BICs = BICs))
}
adj_SEG_LB = function(DAT_L,DAT_B,wSEG,cSEG,BICs = NULL,minSL = 20,verbose = TRUE){
	if(FALSE){
		# DATA = DATA
		wSEG = 1; # aka the wSEG-th segment
		cSEG = SEG; cSEG
		minSL = 20
		verbose = TRUE
		
	}
	
	oBIC = calc_BIC_LB(DAT_L = DAT_L,DAT_B = DAT_B,
		cSEG = cSEG,BICs = BICs)
	BICs = oBIC$BICs
	oBIC$BIC
	
	if( wSEG == length(cSEG) - 1 ){
		return(list(SEG = cSEG,BICs = BICs,ADJ = FALSE))
	}
	
	# Make adjustments for one break
	adj_breaks = seq(-500,500,40)
	adj_breaks = c(adj_breaks,seq(-10,10))
	adj_breaks = sort(unique(adj_breaks))
	preSEG = cSEG[wSEG]; centSEG = cSEG[wSEG + 1]; postSEG = cSEG[wSEG + 2]
	c(preSEG = preSEG,centSEG = centSEG,postSEG = postSEG)
	# centSEG - preSEG
	adj_breaks
	
	out = c()
	for(adj_break in adj_breaks){
		# adj_break = adj_breaks[1]
		if( verbose ) cat(".")
		cSEG_2 = cSEG
		cSEG_2[wSEG + 1] = cSEG_2[wSEG + 1] + adj_break
		if( cSEG_2[wSEG + 1] < cSEG[wSEG] || cSEG_2[wSEG + 1] > cSEG[wSEG + 2] ){
			out = rbind(out,smart_df(ADJ = adj_break,BIC = -Inf))
			next
		}
		cSEG_2 = sort(unique(cSEG_2)); cSEG_2
		
		# Check newSEG min segment width sufficient
		if( min(diff(cSEG_2)) < minSL ){
			out = rbind(out,smart_df(ADJ = adj_break,BIC = -Inf))
			next
		}
		
		oBIC_2 = calc_BIC_LB(DAT_L = DAT_L,DAT_B = DAT_B,
			cSEG = cSEG_2,BICs = BICs)
		BICs = oBIC_2$BICs
		out = rbind(out,smart_df(ADJ = adj_break,BIC = oBIC_2$BIC))
		# if( oBIC_2$BIC > oBIC$BIC ) break
	}
	if( verbose ) cat("\n")
	# out
	
	if( oBIC$BIC < max(out$BIC) ){
		out_SEG = cSEG
		adj_break = out$ADJ[which.max(out$BIC)]; adj_break
		out_SEG[wSEG + 1] = out_SEG[wSEG + 1] + adj_break
		ADJ = TRUE
		if( verbose ) cat("Adjusting segment!\n")
	} else {
		out_SEG = cSEG
		ADJ = FALSE
		if( verbose ) cat("Same segment.\n")
	}
	
	return(list(SEG = out_SEG,BICs = BICs,ADJ = ADJ))
	
}
seg_LRR_BAF = function(DAT_L,DAT_B,iSEG = NULL,minSL = 20,
	verbose = TRUE,verbose2 = TRUE,show_plot = TRUE){
	
	if(FALSE){
		chr = "chr7"
		DAT_L = lrr[which(lrr$Chr == chr),]
		DAT_B = baf[which(baf$Chr == chr),]
		iSEG = NULL; minSL = 20; verbose = TRUE;
		verbose2 = TRUE; show_plot = TRUE
		
		
	}
	
	# Check colnames
	req_names = c("Chr","Start","End","LRR")
	smart_reqNames(DATA = DAT_L,REQ = req_names)
	req_names = c("Chr","Position","AD","RD")
	smart_reqNames(DATA = DAT_B,REQ = req_names)
	
	# Order segments
	chr_names = paste0("chr",c(1:22,"X","Y"))
	int_chrs = intersect(unique(DAT_L$Chr),chr_names)
	if( length(int_chrs) == 0 ) stop("Check Chr coding")
	chr_names = chr_names[chr_names %in% int_chrs]
	DAT_L$Chr = factor(as.character(DAT_L$Chr),levels = chr_names)
	DAT_L = DAT_L[order(DAT_L$Chr,DAT_L$Start),]
	DAT_L$Chr = as.character(DAT_L$Chr)
	# DAT_L$IDX = seq(nrow(DAT_L))
	
	# Make BAF data match up
	pDAT = LRR_BAF_combine(LRR = DAT_L,
		BAF = DAT_B,VAR = "IDX")
	DAT_L = pDAT$LRR
	DAT_B = pDAT$BAF
	
	# Prep BAF data
	DAT_B$DP = rowSums(DAT_B[,c("AD","RD")])
	DAT_B$log_DP = log(DAT_B$DP)
	DAT_B$LBC = lchoose(DAT_B$DP,DAT_B$AD)
	
	# Initialize segments
	if( is.null(iSEG) ){
		dim(DAT_L)
		SEG = round(seq(1,nrow(DAT_L),length.out = 40))
	} else {
		SEG = sort(unique(iSEG))
	}
	# SEG
	
	# Get initial BICs and BIC
	oBIC = calc_BIC_LB(DAT_L = DAT_L,DAT_B = DAT_B,
		cSEG = SEG,BICs = NULL)
	BICs = oBIC$BICs
	# dim(BICs); oBIC$BIC
	
	wSEG = 1
	while(TRUE){
		
		wSEG; length(SEG)
		if( wSEG >= length(SEG) - 1 ){
			if( verbose ) cat("Segmentation complete!\n")
			break
		}
		
		# Try remove segment
		# cat("will run rm seg\n")
		if( show_plot ) show_SEG_LB(DAT_L = DAT_L,DAT_B = DAT_B,
			cSEG = SEG,wSEG = wSEG)
		rem_out = rem_SEG_LB(DAT_L = DAT_L,DAT_B = DAT_B,
			cSEG = SEG,BICs = BICs,verbose = verbose2)
		# names(rem_out)
		REM = FALSE
		if( length(SEG) != length(rem_out$SEG) ){
			REM = TRUE
		} else if( any(SEG != rem_out$SEG) ){
			REM = TRUE
		}
		SEG = rem_out$SEG; BICs = rem_out$BICs
		rm(rem_out)
		if( REM ){
			wSEG = max(c(1,wSEG - 1))
			cat(sprintf("Shift back segment wSEG = %s ...\n",wSEG))
			next
		}
		
		# Try add segment
		# cat("will run add seg\n")
		if( show_plot ) show_SEG_LB(DAT_L = DAT_L,DAT_B = DAT_B,
			cSEG = SEG,wSEG = wSEG)
		add_out = add_SEG_LB(DAT_L = DAT_L,DAT_B = DAT_B,
			wSEG = wSEG,cSEG = SEG,BICs = BICs,minSL = minSL,
			verbose = verbose2)
		# names(add_out)
		SEG = add_out$SEG; BICs = add_out$BICs
		ADD = add_out$ADD; rm(add_out)
		if( ADD ){
			wSEG = max(c(1,wSEG - 1))
			cat(sprintf("Shift back segment wSEG = %s ...\n",wSEG))
			next
		}
		
		# Try adjusting segment
		# cat("will run adj seg\n")
		if( show_plot ) show_SEG_LB(DAT_L = DAT_L,DAT_B = DAT_B,
			cSEG = SEG,wSEG = wSEG)
		adj_out = adj_SEG_LB(DAT_L = DAT_L,DAT_B = DAT_B,
			wSEG = wSEG,cSEG = SEG,BICs = BICs,minSL = minSL,
			verbose = verbose2)
		SEG = adj_out$SEG; BICs = adj_out$BICs
		ADJ = adj_out$ADJ; rm(adj_out)
		if( ADJ ) next
		
		# Shift to next segment
		wSEG = wSEG + 1
		if( verbose ) cat(sprintf("Shift to next segment wSEG = %s ...\n",wSEG))
		
	}
	
	oBIC = calc_BIC_LB(DAT_L = DAT_L,DAT_B = DAT_B,
		cSEG = SEG,BICs = NULL)
	BICs = oBIC$BICs; BICs
	
	tmp_start = sapply(BICs$aSEG,function(xx){
		DAT_L$Start[DAT_L$IDX == xx]
	})
	tmp_end = sapply(BICs$bSEG,function(xx){
		DAT_L$End[DAT_L$IDX == xx]
	})
	tmp_chr = unique(DAT_L$Chr)
	BICs = smart_df(Chr = tmp_chr,Start = tmp_start,
		End = tmp_end,BICs)
	# BICs
	
	if( show_plot ) show_SEG_LB(DAT_L = DAT_L,DAT_B = DAT_B,
		cSEG = SEG,wSEG = NULL)
	
	return(BICs)
	
}




###

