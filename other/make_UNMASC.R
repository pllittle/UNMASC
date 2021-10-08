# Make UNMASC package

# ----------
# Shortcuts, Libraries
# ----------
rm(list=ls())
git_dir = "C:/Users/Admin/Desktop/github"
drive = strsplit(getwd(),"/")[[1]][1]; drive
# source(file.path(git_dir,"smart","RR","SOURCE.R"))
if(FALSE){
makeRpack_fn = file.path(git_dir,"makeRpack","other","makeRpack_1.0.0.tar.gz")
chk_makeRpack = ("makeRpack" %in% installed.packages()[,"Package"]); chk_makeRpack
# remove.packages("makeRpack")
if( !chk_makeRpack ) install.packages(makeRpack_fn,type = "source",repos = NULL)

}
library(makeRpack)

curr_dir = file.path(drive,"Users/Admin/Desktop/packages")
if( !dir.exists(curr_dir) ) dir.create(curr_dir)
setwd(curr_dir)

pack_name = "UNMASC"; pack_version = "1.0.0"
pack_dir = file.path(curr_dir,pack_name)
if( dir.exists(pack_dir) ) unlink(pack_dir,recursive = TRUE)
pack_gz = file.path(curr_dir,sprintf("%s_%s.tar.gz",pack_name,pack_version))
if( file.exists(pack_gz) ) unlink(pack_gz)

# ----------
# Use my makeRpack package function
# ----------
pack_installed = (pack_name %in% installed.packages()[,"Package"]); pack_installed
if( pack_installed ){
	cat(sprintf("Removing previously installed %s package ...\n",pack_name))
	remove.packages(pack_name)
}

if( !pack_installed ){
makeRpack::makeRpack(pack_name = pack_name,work_dir = curr_dir,
	r_scripts = file.path(file.path(git_dir,pack_name,"R"),c("UNMASC.R","details.R")),
	cpp_scripts = file.path(git_dir,pack_name,"src","UNMASC.cpp"),
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

##

