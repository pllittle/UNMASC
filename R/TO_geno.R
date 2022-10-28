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

###