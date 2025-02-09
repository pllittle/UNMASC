# Genome functions
prep_genome = function(bed_centromere_fn,dict_chrom_fn){
	message("Import/Prep genome files, used for plotting genome position on x-axis.\n",appendLF = FALSE)
	
	# Prep centro
	centro = fread(bed_centromere_fn,sep = '\t',
		header = FALSE,data.table = FALSE)
	names(centro) = c("Chr","centro_start","centro_end")
	centro$nChr = gsub("chr","",centro$Chr)
	centro$nChr[centro$nChr == "X"] = "23"
	centro$nChr[centro$nChr == "Y"] = "24"
	centro$nChr = as.numeric(centro$nChr)
	
	tab = table(centro$Chr)
	if( length(tab[tab > 1]) > 0 )
		stop("There are chromosomes with multiple start/end centromere intervals.")
	
	# Prep chr_length
	chr_length = fread(dict_chrom_fn,sep = '\t',
		header = FALSE,data.table = FALSE)
	names(chr_length)[c(2,3)] = c("Chr","length")
	chr_length = chr_length[,c("Chr","length")]
	chr_length = chr_length[grep("chr",chr_length$Chr),]
	chr_length$Chr = gsub("SN:","",chr_length$Chr)
	chr_length = chr_length[which(chr_length$Chr %in% paste0("chr",c(1:22,"X","Y"))),]
	chr_length$length = as.numeric(gsub("LN:","",chr_length$length))
	chr_length$cumul_length = cumsum(chr_length$length)
	
	# Merge and output
	genome = smart_merge(chr_length,centro)
	genome = genome[order(genome$nChr),]
	genome
}
import_bed = function(bed_fn){
	bed_df = fread(file = bed_fn,
		sep = '\t',header = TRUE,
		data.table = FALSE,
		showProgress = FALSE)
	bed_df = bed_df[,c("chrom","start","end")]
	names(bed_df) = c("chrom","start","end")
	bed_df
}
get_genome_coor = function(genome,chr,pos){
	one_index = which(genome$Chr == chr)
	x_coor = pos / genome$length[one_index] + genome$nChr[one_index] - 1
	x_coor
}
create_xCoor = function(genome,x){
	if( !all(c("Chr","Position") %in% names(x)) ){
		stop("Missing Chr and Position columns")
	}
	
	uniq_contigs = unique(x$Chr)
	uniq_contigs = uniq_contigs[which(!(uniq_contigs %in% paste0("chr",c(1:22,"X","Y"))))]
	if( length(uniq_contigs) > 0 ){
		stop("Contigs other than 1-22,X,Y are contained in dataframe")
	}

	x$x_coor = NA
	for(one_chr in unique(x$Chr)){
		one_index = which(x$Chr == one_chr)
		x$x_coor[one_index] = get_genome_coor(genome,one_chr,x$Position[one_index])
	}
	x
}
genome_plot = function(DAT,var,hg="?",chr=NA,horiz=FALSE,
	how_color=NA,my_ylim=NULL,draw_alt_color=TRUE,draw_centromere=TRUE,
	genome,...){
	
	if(FALSE){
		DAT=tmp_list$annoDATA; var="VAF"; hg=hg; chr=chr; pch=16
		cex=tmp_list$annoDATA$pt_cex; how_color="pt_col"; genome=genome
		my_ylim=NULL;draw_alt_color=TRUE;draw_centromere=TRUE
	}
	
	# Range of Chromosomes
	if( is.na(chr) ){
		my_xlim = c(0,24)
		my_xlab = paste0("Genomic Position (hg",hg,")")
	} else {
		if(chr %in% as.character(1:22) ){
			chr = as.numeric(chr)
			my_xlim = c(chr-1,chr)
		} else if( chr %in% c(23,"X") ){
			chr = "X"
			my_xlim = c(22,23)
		} else if(chr %in% c(24,"Y") ){
			chr = "Y"
			my_xlim = c(23,24)
		}
		DAT = DAT[which(DAT$Chr == paste0("chr",chr)),]
		my_xlab = paste0("Position in Chromosome ",chr," in MB (hg",hg,")")
	}
	
	# Point Color
	if( length(how_color)==1 && is.na(how_color) ){
		pt_col = "black"
		DAT = unique(DAT[,c("x_coor",var)])
	} else if( length(how_color) != 1 || (length(how_color)==1 && !(how_color %in% names(DAT))) ){
		pt_col = how_color
		DAT = unique(DAT[,c("x_coor",var)])
	} else if( length(how_color) == 1 && how_color %in% names(DAT) ){
		DAT = unique(DAT[,c("x_coor",var,how_color)])
		pt_col = DAT[,how_color]
	}
	
	# Plot
	if( is.null(my_ylim) ){
		if(var %in% c("tVAF","nVAF","VAF","BAF","gc","UMAP")){
			my_ylim = c(0,1)
		} else if(var %in% c("log2r","LRR")){
			my_ylim = (max(abs(DAT[,var]))+0.2) * c(-1,1)
		} else if(var == "min_tVAF"){
			my_ylim = c(0,0.5)
		} else if(var == "infer_tCN"){
			my_ylim = c(0,max(DAT[,var],na.rm=TRUE)+1)
		} else {
			my_ylim = round(as.numeric(summary(DAT[,var])[c(1,6)]),1)
		}
	}
	
	plot(DAT[,c("x_coor",var)],
		xlab = my_xlab,ylab = var,xlim = my_xlim,ylim = my_ylim,
		xaxs = "i",xaxt = "n",col = pt_col,bty = "n",...)
	
	# Detailing x-axis and vertical labels
	alt_color = "gray"
	if( is.na(chr) ){
		axis(side=1,at=seq(0.5,23.5,1),labels=c(1:22,"X","Y"))
		
		for(cc in seq(0,24)){
			segments(x0=cc,y0=my_ylim[1],x1=cc,y1=my_ylim[2],lwd=0.5,lty=2)
		}
		
		for(cc in seq(1,24,2)){
			#cc = 1
			polygon(x=c(cc-1,cc-1,cc,cc),
				y=c(my_ylim,rev(my_ylim)),
				col=adjustcolor(alt_color,alpha.f=0.5),
				border=NA
				)
		}
	} else {
		if(chr %in% 1:22){
			chr_num = chr
			chr_index = which(genome$nChr == chr_num)
		} else if(chr == "X"){
			chr_num = 23
			chr_index = which(genome$nChr == chr_num)
		} else if(chr == "Y"){
			chr_num = 24
			chr_index = which(genome$nChr == chr_num)
		}
		chr_length = genome$length[chr_index]
		chr_max = chr_length / 1e6
		my_labels = seq(0,chr_max,25)
		my_labels2 = unique(c(my_labels,chr_max))
		my_locate = sapply(my_labels * 1e6,function(x) get_genome_coor(genome,genome$Chr[chr_index],x))
		my_locate2 = sapply(my_labels2 * 1e6,function(x) get_genome_coor(genome,genome$Chr[chr_index],x))
		axis(side=1,at=my_locate,labels=my_labels)
		count = 0
		if( draw_alt_color ){
		for(cc in seq(length(my_locate2))){
			if(count %% 2 != 0){
			polygon(x=c(my_locate2[cc-1],my_locate2[cc-1],my_locate2[cc],my_locate2[cc]),
				y=c(my_ylim,rev(my_ylim)),
				col=adjustcolor(alt_color,alpha.f=0.5),
				border=NA)
			}
			count = count + 1
		}
		}
	}
	
	if(horiz == TRUE) abline(h=seq(0,1,0.1),lwd=0.3,lty=3)
	if( !is.null(var) && var %in% c("tVAF","nVAF","VAF","gc") ){
		segments(x0=my_xlim[1],y0=0.5,x1=my_xlim[2],lwd=0.5,lty=3)
	} else if( !is.null(var) && var %in% c("log2r","LRR")){
		abline(h=0,lwd=0.5,lty=3)
	}
	
	# Draw Centromere regions
	if( draw_centromere ){
		for(cc in genome$Chr){
			# cc = "chr22"
			one_index = which(genome$Chr==cc)
			xx1 = get_genome_coor(genome,cc,genome$centro_start[one_index])
			xx2 = get_genome_coor(genome,cc,genome$centro_end[one_index])
			polygon(x=c(xx1,xx1,xx2,xx2),
				y=c(my_ylim,rev(my_ylim)),
				col=adjustcolor("red",alpha.f=0.7),
				border=NA)
			segments(
				x0=mean(c(xx1,xx2)),y0=my_ylim[1],
				x1=mean(c(xx1,xx2)),y1=my_ylim[2],col=adjustcolor("red",alpha.f=0.7))
		}
	}
	
}


###
