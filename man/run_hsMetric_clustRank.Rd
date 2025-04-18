% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UNMASC.R
\name{run_hsMetric_clustRank}
\alias{run_hsMetric_clustRank}
\title{run_hsMetric_clustRank}
\usage{
run_hsMetric_clustRank(DATA, VARS, ID, TRIALS = 50, note = "")
}
\arguments{
\item{DATA}{A data.frame of sample ID and metrics to cluster on}

\item{VARS}{A character vector of at least two metrics to cluster on
within \code{DATA}}

\item{ID}{A character string specifying the column name in \code{DATA}
that corresponds to sample ID.}

\item{TRIALS}{An integer number of times to re-run the EM algorithm,
more trials help to avoid using a sub-optimal solution.}

\item{note}{A character string to label the plot title}
}
\value{
R data.frame of cluster results and sample labeling.
}
\description{
Clusters sample metrics to identify outlier and 
reliable samples to establish a normal control cohort
}
