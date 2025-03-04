% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{gen_simulated_data}
\alias{gen_simulated_data}
\title{gen_simulated_data}
\usage{
gen_simulated_data(
  outdir,
  purity = NULL,
  cn_state = NULL,
  num_loci = 10000,
  num_normals = 20,
  AA_vaf = 0.002,
  mean_DP = 500,
  oxog_vaf = 0.01,
  ffpe_vaf = 0.1,
  max_PSI = 0.02,
  show_plot = TRUE,
  seed = NULL
)
}
\arguments{
\item{outdir}{Character string specifying the output directory.}

\item{purity}{Tumor purity. By default, \code{purity = NULL} 
will randomly generate a tumor purity.}

\item{cn_state}{A R data.frame containing copy number states 
with columns \code{tCN}, \code{CN_A}, \code{CN_B}, and \code{prob} 
corresponding to the total integer copy number, the minor allelic 
integer copy number, the major allelic integer copy number, and 
probability of selecting a copy number state, respectively. 
By default, \code{cn_state = NULL}, will generate a default 
\code{cn_state} object.}

\item{num_loci}{A positive integer for the number of tumor only 
variant loci simulated.}

\item{num_normals}{A positive integer for the number of unmatched 
normal controls to simulate.}

\item{AA_vaf}{The mean non-zero variant allele frequency for 
homozygous reference positions.}

\item{mean_DP}{A positve numeric value for the mean depth of 
variant calls.}

\item{oxog_vaf}{A positive numeric value for the expected 
allele frequency of oxoG variants.}

\item{ffpe_vaf}{A positive numeric value for the expected 
allele frequency of FFPE variants.}

\item{max_PSI}{A positive numeric value for the maximum overdispersion in
the beta-binomial distribution.}

\item{show_plot}{By default, set to \code{TRUE}. If \code{TRUE},
the variant calls are plotted. This will include germline, artifact, 
and somatic variants.}

\item{seed}{A numeric seed value. If \code{NULL}, a default seed is set.}
}
\value{
A list of simulated data including vcfs, filepaths,
target BED file, dictionary chrom BED file, gender indicator
to reflect the data type and format UNMASC is expecting.
}
\description{
The function simulates one sample's underlying 
somatic point mutations and copy number aberrations and creates 
vcf files formatted for UNMASC. Additional files created include
bed files defining chromosome lengths, centromere regions, and
target regions.
}
