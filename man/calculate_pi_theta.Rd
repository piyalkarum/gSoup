% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.R
\name{calculate_pi_theta}
\alias{calculate_pi_theta}
\title{Calculate nucleotide diversity \eqn{\theta[\pi]}}
\usage{
calculate_pi_theta(
  fast_path,
  window = 100,
  pairwise_deletion = TRUE,
  plot = TRUE,
  progress = TRUE,
  ...
)
}
\arguments{
\item{fast_path}{description
must have at least the three colums \code{type, start, end}}

\item{window}{size of the sliding window. \code{default=100}}

\item{pairwise_deletion}{logical. Whether to ignore gaps in the pairwise alignment. \code{default=TRUE}}

\item{plot}{logical. Whether to plot the output}

\item{progress}{logical. Whether to show progress}

\item{...}{any other parameters to pass to \link{plot.default}}
}
\value{
Returns a list with nucleotide diversity \eqn{\pi}, watterson's theta \eqn{\theta[W]},
total number of sites, number of segregating sites and Tajima's D in a sliding window and average
}
\description{
This function calculates nucleotide diversity and Watterson's theta
from a multiple sequnce alignment
}
\examples{
\dontrun{fasta<-paste0(path.package("gSoup"), "/alignment.fasta")
PT<-calculate_pi_theta(fasta)}



}
\author{
Piyal Karunarathne
}
