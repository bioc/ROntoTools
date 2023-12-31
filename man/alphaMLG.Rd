% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{alphaMLG}
\alias{alphaMLG}
\title{Compute alpha weights}
\usage{
alphaMLG(pv, threshold = max(pv))
}
\arguments{
\item{pv}{vector of p-values}

\item{threshold}{the threshold value that was used to select DE genes}
}
\description{
Transform a vector of p-values into weights.
}
\details{
Computes a set of weights from p-values using the formula \code{-log10(pv/threshold)}.
}
\examples{

load(system.file("extdata/E-GEOD-21942.topTable.RData", package = "ROntoTools"))

head(alphaMLG(top$adj.P.Val))

}
\author{
Calin Voichita and Sorin Draghici
}
\seealso{
\code{\link{pe}}
}

