# Make pUNMASC package

# ----------
# Shortcuts, Libraries
# ----------
rm(list=ls())
git_dir = "C:/Users/Admin/Desktop/github"
curr_dir = "C:/Users/Admin/Desktop/packages"

if( FALSE ){

makeRpack_fn = file.path(git_dir,"makeRpack","other","makeRpack_1.0.0.tar.gz")
chk_makeRpack = ("makeRpack" %in% installed.packages()[,"Package"]); chk_makeRpack
# remove.packages("makeRpack")
if( !chk_makeRpack ) install.packages(makeRpack_fn,type = "source",repos = NULL)

}
library(makeRpack)

if( !dir.exists(curr_dir) ) dir.create(curr_dir)
setwd(curr_dir)

pack_name = "pUNMASC"
pack_version = "1.0.0"
pack_dir = file.path(curr_dir,pack_name)
if( dir.exists(pack_dir) ) unlink(pack_dir,recursive = TRUE)
pack_gz = file.path(curr_dir,sprintf("%s_%s.tar.gz",pack_name,pack_version))
if( file.exists(pack_gz) ) unlink(pack_gz)

pack_installed = (pack_name %in% installed.packages()[,"Package"]); pack_installed
if( pack_installed ){
	cat(sprintf("Removing previously installed %s package ...\n",pack_name))
	remove.packages(pack_name)
}

if( !pack_installed ){

R_dir 	= file.path(git_dir,pack_name,"R")
src_dir = file.path(git_dir,pack_name,"src")

makeRpack::makeRpack(pack_name = pack_name,work_dir = curr_dir,
	r_scripts = file.path(R_dir,c("UNMASC.R","details.R")),
	cpp_scripts = file.path(src_dir,"UNMASC.cpp"),
	# vig_scripts = file.path(git_dir,pack_name,"vignettes","UNMASC.Rmd"),
	pack_imports = c("data.table","doParallel","foreach","stringr","parallel",
		"GenomicRanges","IRanges","Rsamtools","emdbook","scales","seqTools"),
	pack_title = "Tumor-Only Variant Calling with Unmatched Normal Control Samples",
	pack_author = "Paul Little",pack_email = "plittle@fredhutch.org",
	pack_desc = "Tumor-only variant calling is conducted through tumor 
		and unmatched normal annotation, filtering, clustering, and segmentation.",
	pack_version = pack_version,
	CHECK = FALSE,pandoc_dir = "C:/Program Files/RStudio/bin/pandoc")

}

# vignette(pack_name)