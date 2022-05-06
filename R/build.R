
#' @importFrom smarter smart_df smart_table smart_merge
#'	smart_mkdir smart_hist smart_rmcols name_change
#'	smart_WT smart_progress smart_reqNames
#' @importFrom Rcpp sourceCpp
#' @importFrom emdbook rbetabinom
#' @importFrom stats rbinom runif density chisq.test 
#'	fisher.test var quantile sd rnbinom
#' @importFrom graphics par plot hist lines axis segments
#'	polygon abline legend points
#' @importFrom utils head tail read.table write.table
#' @importFrom scales hue_pal
#' @importFrom seqTools char2ascii
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom data.table fread fwrite
#' @importFrom grDevices adjustcolor rgb png dev.off
#' @importFrom Rsamtools BamFile PileupParam ScanBamParam pileup
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @import foreach
#' @useDynLib UNMASC
NULL

# Create package
# rm(list=ls()); library(smarter)
# smart_prepPack(pack_dir = "C:/Users/Admin/Desktop/github/UNMASC",
#		pandoc = "C:/Program Files/RStudio/bin/pandoc",
#		make_vign = TRUE,cran = TRUE,
#		build_dir = NULL)
#		build_dir = "C:/Users/Admin/Desktop")

###
