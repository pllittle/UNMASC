# Test R code

## samtools + mpileup + seqTools to get each read's map quality
rm(list=ls())
library(seqTools) # https://bioconductor.org/packages/3.4/bioc/html/seqTools.html
library(foreach)
source("~/github/smart/RR/SOURCE.R") # get smart_RT(), smart_WT(), smart_ncores()

smart_df = function(...){
	data.frame(...,stringsAsFactors=FALSE)
}
smart_sprintf = function(...){
	orig = sprintf(...)
	orig = gsub("\n","",orig)
	orig = gsub("\t","",orig)
	orig
}
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
	
	ascii_qual = seqTools::char2ascii(c = aa[6])
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
	
	out = smart_df(BASE = ss,BQ = seqTools::char2ascii(c = aa[6]) - 33)
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
check_MQ = function(BAM_fn,CHR,POS){
	# sam = system(smart_sprintf("samtools view -h %s '%s:%s-%s' \
		# | grep -v 'SQ' | grep 'chr' | cut -f1-11",BAM_fn,CHR,POS,POS),intern=TRUE)
	# sam = t(sapply(sam,function(xx) strsplit(xx,"\t")[[1]],USE.NAMES=FALSE))
	xx = as.numeric(system(smart_sprintf("samtools view -h %s '%s:%s-%s' \
		| grep -v 'SQ' | grep 'chr' | cut -f5",BAM_fn,CHR,POS,POS),intern=TRUE))
	print(table(xx))
	out = smart_df(MQ = xx,MQ2 = NA)
	for(th in rev(seq(0,30,10))){
		out$MQ2[which(out$MQ >= th & is.na(out$MQ2))] = th + 10
	}
	out
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
	cl = parallel::makeCluster(NT) # NT = num threads
	doParallel::registerDoParallel(cl)
	# utils::globalVariables("ii")
	ii = NULL
	xx = foreach::foreach(ii=seq(nrow(vcs)),.combine="rbind",
		.export = c("smart_df","smart_sprintf","get_pileup","test_BQ")) %dopar% {
	
		# ii = 2
		out1 = test_BQ(BAM_fn = BAM_fn,CHR = vcs$Chr[ii],POS = vcs$Position[ii],
			REF = vcs$Ref[ii],ALT = vcs$Alt[ii],minBQ = minBQ,FASTA_fn = FASTA_fn,PERMS=PERMS)
		c(out1$asymp_pvalue,out1$perm_pvalue,out1$exact_pvalue)
	}
	parallel::stopCluster(cl)
	
	rownames(xx) = NULL
	xx = smart_df(xx)
	names(xx) = c("MQ_asymp","MQ_perm","MQ_exact")
	vcs = smart_df(vcs,xx)
	vcs
}

if(FALSE){ # Testing above functions

FASTA_fn = "/datastore/nextgenout2/share/labs/UNCseq/tools-data/reference/hg19_hs37d5_viral.fa"
freeze_dir = "/datastore/nextgenout2/share/labs/UNCseq/datafreeze/2018-05-01/data-use"

# RUNX1 issue
# samp = paste0("UNCseq",c(1603,1987,2140,2333,2476))[1]
# CHR = "chr21"; POS = 36164610; REF = "T"

# Samples/Variants that are fine
# samp = paste0("UNCseq",1371)[1]; CHR = "chr12"; POS = 115112278; REF = "C"; ALT = "T"
samp = paste0("UNCseq",1438)[1]; CHR = "chr1"; POS = 144918930; REF = "C"; ALT = "A"

# Tumor
BAM_fn = file.path(freeze_dir,samp,"output","tumor.abra.bam")
tt = test_BQ(BAM_fn = BAM_fn,CHR = CHR,POS = POS,REF = REF,ALT = ALT,
	minBQ = 0,FASTA_fn = FASTA_fn,PERMS = 1e3); tt[-1]
mm = check_MQ(BAM_fn = BAM_fn,CHR = CHR,POS = POS)

# Matched normal
BAM_fn = file.path(freeze_dir,samp,"output","normal.abra.bam")
tt = test_BQ(BAM_fn = BAM_fn,CHR = CHR,POS = POS,REF = REF,ALT = ALT,
	minBQ = 0,FASTA_fn = FASTA_fn,PERMS = 1e3); tt[-1]
mm = check_MQ(BAM_fn = BAM_fn,CHR = CHR,POS = POS)

}
if(FALSE){ # Loop through some benchmark tumor-only variants
pl_dir = "/datastore/hayeslab/Current_members/Paul_Little"
bm_dir = file.path(pl_dir,"UNCseq_proj3/step6_real_TO/s1_process/UNCseq_BM_3")
aa = readRDS(file.path(bm_dir,"output","npn.rds"))
bb = aa[which(nchar(aa$Ref)==1 & nchar(aa$Alt) == 1 
	& aa$inPN == "NO" & aa$nVAF < 0.4 & aa$tVAF < 0.4),
	c("STUDYNUMBER","Chr","Position","Ref","Alt","GeneName",
	"tAD","tRD","tDP","tVAF","nVAF","N_eps","tANNO")][11:30,]
bb$BQ = NA; bb$MQ = NA
for(ii in seq(nrow(bb))){
	# ii = 1
	# cat(".")
	samp = bb$STUDYNUMBER[ii]; CHR = bb$Chr[ii]; POS = bb$Position[ii]
	REF = bb$Ref[ii]; ALT = bb$Alt[ii]
	BAM_fn = file.path(freeze_dir,samp,"output","tumor.abra.bam")
	tmp_list = test_BQ(BAM_fn = BAM_fn,CHR = CHR,POS = POS,REF = REF,ALT = ALT,
		minBQ = 0,FASTA_fn = FASTA_fn,PERMS = 1e3)
	bb$BQ[ii] = tmp_list$perm_pvalue
}
}
if(FALSE){ # Run on one sample's variants
###########################
###########################
FASTA_fn = "/datastore/nextgenout2/share/labs/UNCseq/tools-data/reference/hg19_hs37d5_viral.fa"
freeze_dir = "/datastore/nextgenout2/share/labs/UNCseq/datafreeze/2018-05-01/data-use"
pl_dir = "/datastore/hayeslab/Current_members/Paul_Little"
bm_dir = file.path(pl_dir,"UNCseq_proj3/step6_real_TO/s1_process/UNCseq_BM_3")
mand_dir = "/datastore/nextgenout5/share/labs/bioinformatics/mand/mand2115"
s7_dir = file.path(mand_dir,"somatic_TO_v7_2"); s8_dir = file.path(mand_dir,"somatic_TO_v8_2")
s9_dir = file.path(mand_dir,"somatic_TO_v9_2")

# benchmark samples
samp = paste0("UNCseq",c(1149,1161,1202))[2]; samp 
tmp_fn = file.path(bm_dir,samp,"tumorOnly_VCs.tsv")

# real tumorOnly v8 samples
# samp = paste0("UNCseq",c(1603))[1]; samp 
# samp2 = list.files(s8_dir); samp2 = samp2[grep(samp,samp2)]; samp2
# tmp_fn = file.path(s8_dir,samp2,"tumorOnly_VCs.tsv")

vcs = read.table(tmp_fn,header=TRUE,sep="\t",stringsAsFactors=FALSE,quote="")
vcs = vcs[which(vcs$tANNO >= 0.5 & vcs$N_eps < 0.5 
	# & nchar(vcs$Ref) == 1 & nchar(vcs$Alt) == 1 
	& vcs$Strand_bias_P > 0.05),
	c("Chr","Position","Ref","Alt","GeneName","CosmicOverlaps","tAD","tDP","tVAF",
	"OXOG","FFPE","N_eps","tANNO","Strand_bias_P","nDP_stat","Q_stat")]
dim(vcs)
BAM_fn = file.path(freeze_dir,samp,"output","tumor.abra.bam")
out = test_parallel_MQ(vcs = vcs,BAM_fn = BAM_fn,FASTA_fn = FASTA_fn,
	minBQ = 0,PERMS = 1e3,NT = smart_ncores())

}
if(TRUE){  # Run code on benchmark and real samples

FASTA_fn 		= "/datastore/nextgenout2/share/labs/UNCseq/tools-data/reference/hg19_hs37d5_viral.fa"
freeze_dir 	= "/datastore/nextgenout2/share/labs/UNCseq/datafreeze/2018-05-01/data-use"
pl_dir 			= "/datastore/hayeslab/Current_members/Paul_Little"
bm_dir 			= file.path(pl_dir,"UNCseq_proj3/step6_real_TO/s1_process/UNCseq_BM_3")
mand_dir 		= "/datastore/nextgenout5/share/labs/bioinformatics/mand/mand2115"
s7_dir 			= file.path(mand_dir,"somatic_TO_v7_2")
s8_dir 			= file.path(mand_dir,"somatic_TO_v8_2")
s9_dir 			= file.path(mand_dir,"somatic_TO_v9_2")

cohort = "blah_cohort"; samp = "blah_samp"
# cohort = "bm"; samp = "UNCseq1277"
# cohort = "v7"; samp = "UNCseq0517"

if( cohort %in% "bm" ){
	samp_dir = file.path(bm_dir,samp)
	BAM_dir = freeze_dir
} else if( cohort %in% c("v7","v8","v9") ){
	if( cohort == "v7" ){
		tmp_dir = s7_dir
		BAM_dir = file.path(mand_dir,"somatic_runs/somatic_runs_v7")
	} else if( cohort == "v8" ){
		tmp_dir = s8_dir
		BAM_dir = file.path(mand_dir,"somatic_runs/somatic_runs_sge_v8")
	} else {
		tmp_dir = s9_dir
		BAM_dir = file.path(mand_dir,"somatic_runs/somatic_runs_v9")
	}
	
	samp_dir = file.path(tmp_dir,samp)
}

vcs = smart_RT(file.path(samp_dir,"tumorOnly_VCs.tsv"),
	header = TRUE,sep = "\t",quote = "")
dim(vcs)
BAM_fn = file.path(BAM_dir,samp,"output","tumor.abra.bam")
if( !file.exists(BAM_fn) ){
	BAM_fn = list.files(file.path(BAM_dir,samp),recursive = TRUE)
	BAM_fn = BAM_fn[grep("output/tumor.abra.bam$",BAM_fn)]
	BAM_fn = file.path(BAM_dir,samp,BAM_fn)
	BAM_fn
}

mq_out = test_parallel_MQ(
	vcs = vcs,
	# vcs = vcs[30:31,],
	BAM_fn = BAM_fn,FASTA_fn = FASTA_fn,
	minBQ = 0,PERMS = 1e3,NT = smart_ncores())
# vcs[which(vcs$Chr=="chr14" & vcs$Position==105243640),]

smart_WT(mq_out,file.path(samp_dir,"tumorOnly_VCs_MQ.tsv"),sep="\t")

if(FALSE){

mm = mq_out
mm = mm[which(mm$tANNO >= 0.5 & mm$N_eps < 0.5 & mm$Strand_bias_P > 0.05),
	c("Chr","Position","Ref","Alt","tVAF","tDP","GeneName","tANNO","N_eps",
		"CosmicOverlaps","OXOG","FFPE","MQ_asymp","MQ_perm","MQ_exact")]
dim(mm); mm

#


}

q("no")

}
if(FALSE){ # Aggregate variants bm and CHANCE 531 sample cohort

pl_dir 			= "/datastore/hayeslab/Current_members/Paul_Little"
bm_dir 			= file.path(pl_dir,"UNCseq_proj3/step6_real_TO/s1_process/UNCseq_BM_3")
mand_dir 		= "/datastore/nextgenout5/share/labs/bioinformatics/mand/mand2115"
s7_dir 			= file.path(mand_dir,"somatic_TO_v7_2")
s8_dir 			= file.path(mand_dir,"somatic_TO_v8_2")
s9_dir 			= file.path(mand_dir,"somatic_TO_v9_2")
hnsc531_dir = file.path(mand_dir,"output/hnsc531")

hh = smart_RT(file.path(hnsc531_dir,"hnsc531_HJ_DNH_merge.tsv"),
	sep = "\t",header = TRUE,quote = ""); dim(hh)
hh$MQ_asymp = NA; hh$MQ_perm = NA; hh$MQ_exact = NA
tmp_names = paste0("MQ_",c("asymp","perm","exact"))
us = unique(hh[,c("cohort","STUDYNUMBER_1")]); dim(us)
for(ii in seq(nrow(us))){
	# ii = 1
	if(ii %% 5 == 0) cat(".")
	if(ii %% 2e2 == 0) cat(sprintf("%s\n",ii))
	cc = us$cohort[ii]; uu = us$STUDYNUMBER_1[ii]
	if( cc == "v7" ){
		tmp_fn = file.path(s7_dir,uu,"tumorOnly_VCs_MQ.tsv")
	} else if( cc == "v8" ){
		tmp_fn = file.path(s8_dir,uu,"tumorOnly_VCs_MQ.tsv")
	} else if( cc == "v9" ){
		tmp_fn = file.path(s9_dir,uu,"tumorOnly_VCs_MQ.tsv")
	}
	tmp_idx = which(hh$STUDYNUMBER_1 == uu)
	tmp_df = smart_RT(tmp_fn,header = TRUE,sep = "\t",quote = "")
	if( nrow(tmp_df) == length(tmp_idx) ){
		tmp_df = tmp_df[match(tmp_df$mutID,hh$mutID[tmp_idx]),]
	} else {
		dup_idx = which(duplicated(hh$mutID[tmp_idx]))
		hh = hh[-tmp_idx[dup_idx],]
		tmp_idx = which(hh$STUDYNUMBER_1 == uu)
		if( nrow(tmp_df) != length(tmp_idx) ) stop("Issue with sample's variants.")
	}
	
	# dim(tmp_df); dim(hh[which(hh$STUDYNUMBER_1 == uu),])
	hh[tmp_idx,tmp_names] = tmp_df[,tmp_names]
	
	hh[which(hh$STUDYNUMBER_1 == uu & hh$mutID == 3128),]
	
	rm(tmp_df)
}
out_fn = file.path(hnsc531_dir,"hnsc531_HJ_DNH_PL")
write.table(hh,sprintf("%s.tsv",out_fn),sep = "\t",row.names = FALSE)
saveRDS(hh,sprintf("%s.rds",out_fn))

# Remarks
# RUNX1: 					Positions 36164605 and 36164610
# SLC39A6: 				Positions 33706642 and 33706644
# RP11-201K10.3: 	Position 155161773
# AKT1: 					Position 105240430
# SLC16A3: 				Position 80195478
# PGR: 						Position 100998618





}
if(FALSE){ # Check on BM variants

pl_dir 			= "/datastore/hayeslab/Current_members/Paul_Little"
bm_dir 			= file.path(pl_dir,"UNCseq_proj3/step6_real_TO/s1_process/UNCseq_BM_3")
all_samps 	= list.files(bm_dir); all_samps = all_samps[grep("UNCseq",all_samps)]

samp = all_samps[13]; samp
aa = smart_RT(file.path(bm_dir,samp,"fNPN_segANNO_F.tsv"),sep="\t",header=TRUE,quote="")
aa = aa[,names(aa) != "snpEff"]
mm = smart_RT(file.path(bm_dir,samp,"tumorOnly_VCs_MQ.tsv"),sep="\t",header=TRUE,quote="")
mm = mm[,names(mm) != "snpEff"]
notpn_idx = which(aa$inPN == "NO")

########
########
########

aa[notpn_idx,]
mm[which(mm$mutID %in% aa$mutID[notpn_idx]),]



}
#






###

