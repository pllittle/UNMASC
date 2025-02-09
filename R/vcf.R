# Import Files

import_strelkaVCF = function(vcf_fn,FILTER = NULL,ncores = 1){
	if(FALSE){
		vcf_fn = fn; FILTER = NULL; ncores = ncores
		
	}
	
	# Check format
	aa = readLines(vcf_fn,n = 100)
	
	bb = aa[grepl("^##fileformat=",aa)]
	bb = strsplit(bb,"=")[[1]][2]
	if( bb != "VCFv4.2" ) warning("This function may not prepare 
		columns correctly and may need updating")
	
	bb = aa[grepl("^##cmdline=",aa)]
	bb = bb[grepl("configureStrelkaSomaticWorkflow.py",bb)]
	if( length(bb) != 1 ) warning("This function may not prepare 
		columns correctly and may need updating")
	
	# Import
	tmp_df = fread(vcf_fn,sep = "\t",header = TRUE,
		data.table = FALSE,skip = "#CHROM",nThread = ncores)
	if( !is.null(FILTER) ){
		tmp_df = tmp_df[which(tmp_df$FILTER %in% FILTER),]
	}
	
	# Extract columns Qscore,tAD,tRD,nAD,nRD
	tmp_df$Qscore = NA; tmp_df$tAD = NA; tmp_df$tRD = NA
	tmp_df$nAD = NA; tmp_df$nRD = NA
	
	idx_bs = nchar(tmp_df$REF) == 1 & nchar(tmp_df$ALT) == 1
	first_bs = which(idx_bs)[1]
	chk_elem_QSS_NT = grep("^QSS_NT",strsplit(tmp_df$INFO[first_bs],";")[[1]])
	tmp_df$Qscore[idx_bs] = as.integer(gsub("QSS_NT=","",
		sapply(tmp_df$INFO[idx_bs],function(xx) strsplit(xx,";")[[1]][chk_elem_QSS_NT],
		USE.NAMES = FALSE)))
	bases = c("A","C","G","T")
	for(base in bases){
		# base = bases[1]; base
		if( base == "A" ){
			sub_idx = grep("^AU",strsplit(tmp_df$FORMAT[first_bs],":")[[1]])
		} else if( base == "C" ){
			sub_idx = grep("^CU",strsplit(tmp_df$FORMAT[first_bs],":")[[1]])
		} else if( base == "G" ){
			sub_idx = grep("^GU",strsplit(tmp_df$FORMAT[first_bs],":")[[1]])
		} else if( base == "T" ){
			sub_idx = grep("^TU",strsplit(tmp_df$FORMAT[first_bs],":")[[1]])
		}
		idx_alt_base = idx_bs & tmp_df$ALT == base
		idx_ref_base = idx_bs & tmp_df$REF == base
		tmp_df$tAD[idx_alt_base] = as.integer(sapply(tmp_df$TUMOR[idx_alt_base],
			function(xx) strsplit(strsplit(xx,":")[[1]][sub_idx],",")[[1]][1],USE.NAMES = FALSE))
		tmp_df$tRD[idx_ref_base] = as.integer(sapply(tmp_df$TUMOR[idx_ref_base],
			function(xx) strsplit(strsplit(xx,":")[[1]][sub_idx],",")[[1]][1],USE.NAMES = FALSE))
		
		tmp_df$nAD[idx_alt_base] = as.integer(sapply(tmp_df$NORMAL[idx_alt_base],
			function(xx) strsplit(strsplit(xx,":")[[1]][sub_idx],",")[[1]][1],USE.NAMES = FALSE))
		tmp_df$nRD[idx_ref_base] = as.integer(sapply(tmp_df$NORMAL[idx_ref_base],
			function(xx) strsplit(strsplit(xx,":")[[1]][sub_idx],",")[[1]][1],USE.NAMES = FALSE))
	}
	
	idx_indel = !idx_bs
	if( sum(idx_indel) > 0 ){
		first_indel = which(idx_indel)[1]
		get_indel_qscore = function(xx){
			# xx = tmp_df$INFO[idx_indel][1]
			xx2 = strsplit(xx,";")[[1]]
			xx2 = xx2[grepl("^QSI_NT=",xx2)]
			xx2 = as.numeric(gsub("QSI_NT=","",xx2))
			xx2
		}
		tmp_df$Qscore[idx_indel] = sapply(tmp_df$INFO[idx_indel],function(xx) 
			get_indel_qscore(xx),USE.NAMES = FALSE)
		tmp_df[idx_indel,]
		chk_elem_TIR = grep("^TIR",strsplit(tmp_df$FORMAT[first_indel],":")[[1]])
		chk_elem_TAR = grep("^TAR",strsplit(tmp_df$FORMAT[first_indel],":")[[1]])
		
		tmp_df$tAD[idx_indel] = as.integer(sapply(tmp_df$TUMOR[idx_indel],
			function(xx) strsplit(strsplit(xx,":")[[1]][chk_elem_TIR],",")[[1]][1],USE.NAMES = FALSE))
		tmp_df$tRD[idx_indel] = as.integer(sapply(tmp_df$TUMOR[idx_indel],
			function(xx) strsplit(strsplit(xx,":")[[1]][chk_elem_TAR],",")[[1]][1],USE.NAMES = FALSE))
		tmp_df$nAD[idx_indel] = as.integer(sapply(tmp_df$NORMAL[idx_indel],
			function(xx) strsplit(strsplit(xx,":")[[1]][chk_elem_TIR],",")[[1]][1],USE.NAMES = FALSE))
		tmp_df$nRD[idx_indel] = as.integer(sapply(tmp_df$NORMAL[idx_indel],
			function(xx) strsplit(strsplit(xx,":")[[1]][chk_elem_TAR],",")[[1]][1],USE.NAMES = FALSE))
	}
	
	return(tmp_df)
	
}
import_ANNO = function(anno_fn,nlines = 100,ncores = 1){
	if(FALSE){
		anno_fn = anno_fn; nlines = nlines; ncores = ncores
		
	}
	
	# Check format
	aa = readLines(anno_fn,n = nlines)
	aa = aa[grepl("^#",aa)]
	aa = gsub("\"","",aa)
	bb = aa[grepl("^##INFO=<ID=CSQ",aa)]
	if( length(bb) != 1 ) stop("Can't find ID=CSQ line in annotated VCF")
	bb = strsplit(bb," ")[[1]]
	bb = bb[grepl("IMPACT",bb)]
	bb = gsub(">","",bb)
	if( length(bb) != 1 || 
		bb[1] != "IMPACT|Consequence|SYMBOL|HGVSc|HGVSp|AF|gnomAD_AF|COSMIC|COSMIC_CNT|COSMIC_LEGACY_ID" )
		stop("Annotation file may not be correctly processed")
	
	message(sprintf("%s: Import annotation file ...\n",date()),appendLF = FALSE)
	ANNO = fread(anno_fn,sep = "\t",header = TRUE,
		skip = "#CHROM",nThread = ncores,data.table = FALSE)
	
	message(sprintf("%s: Parse annotation ...\n",date()),appendLF = FALSE)
	
	# Extract out annotations
	ANNO = name_change(ANNO,"#CHROM","Chr")
	ANNO = name_change(ANNO,"POS","Position")
	ANNO = name_change(ANNO,"REF","Ref")
	ANNO = name_change(ANNO,"ALT","Alt")
	ANNO = smart_rmcols(ANNO,c("ID","QUAL","FILTER"))
	ANNO = ANNO[which(ANNO$Chr %in% paste0("chr",c(1:22,"X","Y"))),]
	
	# Cosmic/legacy ids
	message(sprintf("%s: Process COSMIC ...\n",date()),appendLF = FALSE)
	ANNO$CosmicOverlaps = ANNO$CosmicIDs = NA
	idx = which(grepl("\\|COSV",ANNO$INFO) 
		| grepl("\\|COSM",ANNO$INFO))
	length(idx)
	
	if( length(idx) == 0 ){
		message(sprintf("%s: No COSMIC hits found ...\n",date()),appendLF = FALSE)
	} else {
		message(sprintf("%s: %s COSMIC hits found ...\n",date(),length(idx)),appendLF = FALSE)
		tmp_mat = t(sapply(ANNO$INFO[idx],function(xx){
			# xx = ANNO$INFO[idx][92430]; xx
			xx2 = strsplit(xx,";")[[1]]
			xx2 = xx2[grepl("^CSQ=",xx2)]
			xx2 = strsplit(xx2,"=")[[1]][2]
			xx2 = strsplit(xx2,",")[[1]]
			# xx2
			
			xx2 = t(sapply(xx2,function(zz){
				strsplit(zz,"[|]")[[1]][8:10] # COSV_ID, CNT, COSM_ID
			},USE.NAMES = FALSE))
			xx2 = unique(xx2)
			if( nrow(xx2) != 1 ) stop("update COSMIC parsing")
			xx2 = c(xx2)
			xx2
			
			# Process COSV_ID
			if( !is.na(xx2[1]) )
				xx2[1] = paste(unique(strsplit(xx2[1],"&")[[1]]),collapse = "&")
			
			# Process CNT
			xx2[2] = max(as.integer(strsplit(xx2[2],"&")[[1]]))
			
			# Process COSM_ID
			if( !is.na(xx2[3]) )
				xx2[3] = paste(unique(strsplit(xx2[3],"&")[[1]]),collapse = "&")
			
			c(xx2[2],paste(xx2[-2][!is.na(xx2[-2])],collapse = ";"))
		},USE.NAMES = FALSE))
		tmp_mat = smart_df(tmp_mat)
		names(tmp_mat) = c("CosmicOverlaps","CosmicIDs")
		tmp_mat$CosmicOverlaps = as.integer(tmp_mat$CosmicOverlaps)
		ANNO[idx,c("CosmicOverlaps","CosmicIDs")] = tmp_mat
		rm(idx,tmp_mat)
	}
	ANNO$CosmicOverlaps[is.na(ANNO$CosmicOverlaps)] = 0
	ANNO$CosmicIDs[is.na(ANNO$CosmicIDs)] = ""
	
	message(sprintf("%s: Process snpEff/Impact ...\n",date()),appendLF = FALSE)
	tmp_mat = t(sapply(ANNO$INFO,function(xx){
		# xx = ANNO$INFO[1]
		xx2 = strsplit(xx,";")[[1]]
		xx2 = xx2[grepl("^CSQ=",xx2)]
		xx2 = strsplit(xx2,"=")[[1]][2]
		xx2 = strsplit(xx2,",")[[1]]
		# xx2
		
		for(EFF in c("HIGH","MODERATE","LOW","MODIFIER")){
			if( any(grepl(EFF,xx2)) ){
				idx = grep(EFF,xx2)
				break
			}
		}
		xx2 = xx2[idx]
		# xx2
		
		# extract impact|consequence|symbol|hgvsc|hgvsp, drop cosmic/exac/1000g
		xx2 = sapply(xx2,function(zz){
			paste(strsplit(zz,"[|]")[[1]][1:5],collapse = "|")
		},USE.NAMES = FALSE)
		# xx2
		
		c(paste(xx2,collapse = ","),EFF)
	},USE.NAMES = FALSE))
	tmp_mat = smart_df(tmp_mat)
	names(tmp_mat) = c("snpEff","IMPACT")
	ANNO = smart_rmcols(ANNO,c("snpEff","IMPACT"))
	ANNO = cbind(ANNO,tmp_mat); rm(tmp_mat)
	smart_table(ANNO$IMPACT)
	
	message(sprintf("%s: Process snpEff_pt1 ...\n",date()),appendLF = FALSE)
	ANNO$snpEff_pt1 = sapply(ANNO$snpEff,function(xx){
		# xx = ANNO$snpEff[1]
		xx2 = strsplit(xx,",")[[1]]
		xx2 = sapply(xx2,function(zz){
			strsplit(zz,"[|]")[[1]][2]
		},USE.NAMES = FALSE)
		xx2 = unique(xx2)
		if( any(grepl("&",xx2)) ){
			xx2 = unlist(sapply(xx2,function(zz)
				strsplit(zz,"&")[[1]],
				USE.NAMES = FALSE))
		}
		xx2 = sort(unique(xx2))
		paste(xx2,collapse = ",")
	},USE.NAMES = FALSE)
	ANNO$snpEff_pt1 = toupper(ANNO$snpEff_pt1)
	tab = names(smart_table(ANNO$snpEff_pt1))
	
	tmp_df = smart_df(EFF = sort(unique(unlist(sapply(tab,
		function(xx) strsplit(xx,",")[[1]],USE.NAMES = FALSE)))))
	tmp_df$EXONIC = NA
	
	# Prioritizing variants: 
	#		https://www.targetvalidation.org/variants
	# 	https://pcingola.github.io/SnpEff/se_inputoutput/#effect-prediction-details
	message(sprintf("%s: Process EXONIC ...\n",date()),appendLF = FALSE)
	poss_exonic = c("CODING_SEQUENCE_VARIANT","INFRAME_DELETION",
		"MISSENSE_VARIANT","SPLICE_ACCEPTOR_VARIANT","SPLICE_DONOR_VARIANT",
		"SPLICE_REGION_VARIANT","START_LOST","STOP_GAINED","STOP_LOST",
		"FRAMESHIFT_VARIANT","INFRAME_INSERTION","PROTEIN_ALTERING_VARIANT",
		"5_PRIME_UTR_PREMATURE_START_CODON_GAIN_VARIANT","STOP_RETAINED_VARIANT",
		"INITIATOR_CODON_VARIANT","RARE_AMINO_ACID_VARIANT","START_RETAINED_VARIANT",
		"GENE_FUSION","BIDIRECTIONAL_GENE_FUSION","EXON_LOSS_VARIANT",
		"5_PRIME_UTR_TRUNCATION","INCOMPLETE_TERMINAL_CODON_VARIANT",
		"SPLICE_DONOR_5TH_BASE_VARIANT","SPLICE_DONOR_REGION_VARIANT",
		"SPLICE_POLYPYRIMIDINE_TRACT_VARIANT")
	poss_nonexonic = c("3_PRIME_UTR_VARIANT","5_PRIME_UTR_VARIANT",
		"DOWNSTREAM_GENE_VARIANT","INTRON_VARIANT","NMD_TRANSCRIPT_VARIANT",
		"NON_CODING_TRANSCRIPT_EXON_VARIANT","NON_CODING_TRANSCRIPT_VARIANT",
		"SYNONYMOUS_VARIANT","UPSTREAM_GENE_VARIANT","INTERGENIC_VARIANT",
		"MATURE_MIRNA_VARIANT","INTERGENIC_REGION","INTRAGENIC_VARIANT")
	ANNO$EXONIC = NA
	for(pp in poss_exonic){
		idx = which(is.na(ANNO$EXONIC) & grepl(pp,ANNO$snpEff_pt1))
		ANNO$EXONIC[idx] = "YES"
		idx = which(is.na(tmp_df$EXONIC) & grepl(pp,tmp_df$EFF))
		tmp_df$EXONIC[idx] = "YES"
		rm(idx)
	}
	for(pp in poss_nonexonic){
		idx = which(is.na(ANNO$EXONIC) & grepl(pp,ANNO$snpEff_pt1))
		ANNO$EXONIC[idx] = "NO"
		idx = which(is.na(tmp_df$EXONIC) & grepl(pp,tmp_df$EFF))
		tmp_df$EXONIC[idx] = "NO"
		rm(idx)
	}
	# ANNO$EXONIC[which(is.na(ANNO$EXONIC) & ANNO$snpEff_pt1 %in% poss_exonic)] = "YES"
	# ANNO$EXONIC[which(is.na(ANNO$EXONIC) & ANNO$snpEff_pt1 %in% poss_nonexonic)] = "NO"
	smart_table(ANNO$EXONIC)
	print(tmp_df)
	if( any(is.na(ANNO$EXONIC)) || any(is.na(tmp_df$EXONIC)) ){
		print("Refer to https://platform.opentargets.org/variants")
		print("Need to include more EXONIC/nonEXONIC classifications!")
		stop("Need to include more EXONIC/nonEXONIC classifications!")
	}
	ANNO = smart_rmcols(ANNO,c("snpEff_pt1"))
	
	message(sprintf("%s: Process GeneName ...\n",date()),appendLF = FALSE)
	ANNO$GeneName = sapply(ANNO$snpEff,function(xx){
		# xx = ANNO$snpEff[1]
		xx2 = strsplit(xx,",")[[1]]
		xx2 = sapply(xx2,function(zz){
			strsplit(zz,"[|]")[[1]][3]
		},USE.NAMES = FALSE)
		xx2 = sort(unique(xx2))
		paste(xx2,collapse = ",")
	},USE.NAMES = FALSE)
	
	message(sprintf("%s: Process ThsdG_AF ...\n",date()),appendLF = FALSE)
	ANNO$ThsdG_AF = sapply(ANNO$INFO,function(xx){
		# xx = ANNO$INFO[1]
		xx2 = strsplit(xx,";")[[1]]
		xx2 = xx2[grepl("^CSQ=",xx2)]
		xx2 = strsplit(xx2,"=")[[1]][2]
		xx2 = strsplit(xx2,",")[[1]]
		xx2 = sapply(xx2,function(zz){
			strsplit(zz,"[|]")[[1]][6]
		},USE.NAMES = FALSE)
		xx2 = unique(xx2)
		xx2
	},USE.NAMES = FALSE)
	ANNO$ThsdG_AF[which(ANNO$ThsdG_AF == "")] = 0
	ANNO$ThsdG_AF = as.numeric(ANNO$ThsdG_AF)
	
	message(sprintf("%s: Process EXAC_AF ...\n",date()),appendLF = FALSE)
	ANNO$EXAC_AF = sapply(ANNO$INFO,function(xx){
		# xx = ANNO$INFO[1]
		xx2 = strsplit(xx,";")[[1]]
		xx2 = xx2[grepl("^CSQ=",xx2)]
		xx2 = strsplit(xx2,"=")[[1]][2]
		xx2 = strsplit(xx2,",")[[1]]
		xx2 = sapply(xx2,function(zz){
			strsplit(zz,"[|]")[[1]][7]
		},USE.NAMES = FALSE)
		xx2 = unique(xx2)
		if( xx2 == "" ) return(xx2)
		xx2 = as.numeric(strsplit(xx2,"&")[[1]])
		xx2 = max(xx2)
		if( length(xx2) != 1 ) stop("update ExAC parsing")
		xx2
	},USE.NAMES = FALSE)
	ANNO$EXAC_AF[which(ANNO$EXAC_AF == "")] = 0
	ANNO$EXAC_AF = as.numeric(ANNO$EXAC_AF)
	
	ANNO = smart_rmcols(ANNO,"INFO")
	
	message(sprintf("%s: Finishing import_ANNO() ...\n",date()),appendLF = FALSE)
	return(ANNO)
	
}
prep_TARGET = function(TARGET,ANNO){
	
	message(sprintf("%s: Append TARGET ...\n",date()),appendLF = FALSE)
	req_names = c("Chr","Start","End")
	smart_reqNames(DATA = TARGET,REQ = req_names)
	
	req_names = c("Chr","Position")
	smart_reqNames(DATA = ANNO,REQ = req_names)
	
	nsegs = nrow(TARGET)
	ANNO$TARGET = "NO"
	for(ii in seq(nsegs)){
		smart_progress(ii = ii,nn = nsegs,iter = 20,iter2 = 1e3)
		idx = which(ANNO$Chr == TARGET$Chr[ii]
			& ANNO$Position >= TARGET$Start[ii]
			& ANNO$Position <= TARGET$End[ii])
		if( length(idx) > 0 ){
			ANNO$TARGET[idx] = "YES"
		}
		rm(idx)
		if( ii %% 1000 == 0 || ii == nsegs ) print(table(ANNO$TARGET))
	}
	
	return(ANNO)
	
}

#' @title import_VCFs
#' @description Imports Strelka VCFs containing SNVs and INDELs
#' @param DAT A data.frame containing column names 'FILENAME' for
#'	the full vcf filename and 'STUDYNUMBER' for a unique string
#'	mapped to the control sample used when calling variants 
#'	against the same tumor.
#' @param FILTER A list containing named elements 'nDP', 'tDP', 'Qscore',
#'	and 'contigs' for minimum normal depth, minimum tumor depth, minimum
#'	Qscore, and a string of contigs to retain, respectively. Default values
#'	are nDP = tDP = 2, Qscore = 3, and contigs = chr1-chr22, chrX, chrY.
#' @inheritParams run_UNMASC
#' @return A R data.frame containing vcf fields
#' @export
import_VCFs = function(DAT,FILTER = NULL,ncores = 1){
	
	if( is.null(FILTER) ){
		FILTER = list(nDP = 2,tDP = 2,Qscore = 3,
			contigs = sprintf("chr%s",c(1:22,"X","Y")))
	}
	if( is.null(FILTER$nDP) || !is.numeric(FILTER$nDP) )
		stop("Issue with FILTER$nDP")
	if( is.null(FILTER$tDP) || !is.numeric(FILTER$tDP) )
		stop("Issue with FILTER$tDP")
	if( is.null(FILTER$Qscore) || !is.numeric(FILTER$Qscore) )
		stop("Issue with FILTER$Qscore")
	if( is.null(FILTER$contigs) || !is.character(FILTER$contigs) )
		stop("Issue with FILTER$contigs")
	
	smart_header(OBJ = DAT,req_COLS = c("FILENAME","STUDYNUMBER"))
	DAT$CNT = seq(nrow(DAT))
	
	vcf = c()
	
	for(CNT in DAT$CNT){
		message(sprintf("%s: Import vcf normal_%s ...\n",date(),CNT),appendLF = FALSE)
		fn = DAT$FILENAME[DAT$CNT == CNT]
		tmp_vcf = import_strelkaVCF(vcf_fn = fn,
			FILTER = NULL,ncores = ncores)
		tmp_vcf = name_change(tmp_vcf,"#CHROM","Chr")
		tmp_vcf = name_change(tmp_vcf,"POS","Position")
		tmp_vcf = name_change(tmp_vcf,"REF","Ref")
		tmp_vcf = name_change(tmp_vcf,"ALT","Alt")
		tmp_vcf = smart_rmcols(tmp_vcf,c("ID","QUAL","FILTER",
			"INFO","FORMAT","NORMAL","TUMOR"))
		
		# Filtering
		tmp_vcf = tmp_vcf[which(tmp_vcf$Qscore >= FILTER$Qscore
			& rowSums(tmp_vcf[,c("nAD","nRD")]) >= FILTER$nDP
			& rowSums(tmp_vcf[,c("tAD","tRD")]) >= FILTER$tDP
			& tmp_vcf$Chr %in% FILTER$contigs
			),]
		
		if( nrow(tmp_vcf) == 0 ){
			message(smart_sprintf("%s: Skipping normal_%s b/c no variants after basic/low filtering\n",date(),CNT),appendLF = FALSE)
			CNT = CNT + 1
			next
		}
		
		tmp_vcf$STUDYNUMBER = DAT$STUDYNUMBER[DAT$CNT == CNT]
		vcf = rbind(vcf,tmp_vcf); rm(tmp_vcf)
		CNT = CNT + 1
		
	}
	
	return(vcf)
}

#' @title prep_UNMASC_VCF
#' @description A comprehensive function to import VCFs, annotations,
#'	and target bed file to prepare the necessary VCF for UNMASC's main 
#'	workflow.
#' @param target_fn A filename for the tab delimited target 
#'	capture bed file. The first three column headers should be 
#'	'Chr', 'Start', and 'End'.
#' @param anno_fn A filename for the unique set of variant loci called
#'	across multiple normals and annotated by VEP.
#' @param nlines An integer specifying the number of initial lines to
#'	read in from the \code{anno_fn} to assess the expected formatting.
#' @inheritParams run_UNMASC
#' @inheritParams import_VCFs
#' @return A R data.frame containing VCF fields and annotation.
#' @export
prep_UNMASC_VCF = function(outdir,DAT,FILTER,target_fn,
	anno_fn,nlines = 100,ncores = 1){
	if(FALSE){
		outdir = type_dir; DAT = ndat
		FILTER = NULL; target_fn = targ_bed_fn
		anno_fn = anno_fn; nlines = 100; ncores = 1
		
	}
	
	# Get VCFs
	vcfs_rds_fn = file.path(outdir,"vcfs.rds")
	if( !file.exists(vcfs_rds_fn) ){
		vcfs = import_VCFs(DAT = DAT,FILTER = FILTER,
			ncores = ncores)
		saveRDS(vcfs,vcfs_rds_fn)
	}
	vcfs = readRDS(vcfs_rds_fn)
	
	# Get ANNO
	anno_rds_fn = file.path(outdir,"anno.rds")
	if( !file.exists(anno_rds_fn) ){
		anno = import_ANNO(anno_fn = anno_fn,
			nlines = nlines,ncores = ncores)
		saveRDS(anno,anno_rds_fn)
	}
	anno = readRDS(anno_rds_fn)
	
	# Get TARGET regions, append to anno
	aa = readLines(target_fn,n = 1)
	aa = strsplit(aa,"\t")[[1]]
	header = FALSE
	if( aa[1] == "Chr" && aa[2] == "Start" && aa[3] == "End" ){
		header = TRUE
	}
	if( !header ) stop("Fix target_fn headers")
	target = fread(target_fn,sep = "\t",header = header,
		data.table = FALSE)
	target = target[,c("Chr","Start","End")]
	dim(target); target[1:10,]
	anno = prep_TARGET(TARGET = target,ANNO = anno)
	
	message(sprintf("%s: Merge vcf and anno ...\n",date()),appendLF = FALSE)
	final_vcf = smart_merge(vcfs,anno)
	
	return(final_vcf)
}

###

