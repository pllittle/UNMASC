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


###

