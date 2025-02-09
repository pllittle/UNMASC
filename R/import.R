
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
#' @importFrom methods is
#' @import foreach
#' @useDynLib UNMASC
NULL

###
