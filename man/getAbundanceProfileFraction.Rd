% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manageInAndOut.R
\name{getAbundanceProfileFraction}
\alias{getAbundanceProfileFraction}
\title{getAbundanceProfileFraction}
\usage{
getAbundanceProfileFraction(metagenomeDir, part, minMapq)
}
\arguments{
\item{metagenomeDir}{path to the datasystem}

\item{part}{share of the reads to take from every genome}

\item{minMapq}{minMapq of the reads to be cosidered}
}
\description{
This function gets the coverage needed to get the fraction of the reads specified in part from a MGCS datasystem.
}
