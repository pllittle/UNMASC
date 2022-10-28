
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

###
