rm(list=ls())
library(data.table)
library(smarter)
library(UNMASC)

outdir = "/scratch/primary/UTHSC/Current_members/Paul_Little/prac/UNMASC"
vcf_dir = file.path(outdir,"prep_UNMASC_VCF/VCF")
DAT = list.files(vcf_dir,pattern = ".vcf$")
DAT = smart_df(FILENAME = sprintf("%s/%s",vcf_dir,DAT))
DAT

mat = t(sapply(DAT$FILENAME,function(xx){
	strsplit(xx,"_")[[1]]
},USE.NAMES = FALSE))
colnames(mat) = sprintf("V%s",seq(ncol(mat)))

for(nm in colnames(mat)){
	if( length(unique(mat[,nm])) == 1 )
		mat = mat[,-which(colnames(mat) == nm),drop = FALSE]
}
DAT$STUDYNUMBER = as.character(mat)
DAT

anno_fn = file.path(outdir,"prep_UNMASC_VCF","anno_fn","cjs.51A.high_allvar.vcf")
target_fn = file.path(outdir,"prep_UNMASC_VCF","target_fn","S07604514_Regions.bed")

highvcf = UNMASC::prep_UNMASC_VCF(
	outdir = outdir,
	DAT = DAT,
	anno_fn = anno_fn,
	target_fn = target_fn,
	FILTER = NULL,
	nlines = 100,
	ncores = 1)

###

