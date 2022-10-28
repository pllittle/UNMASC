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

###
