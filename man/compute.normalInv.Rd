\name{compute.normalInv}
\alias{compute.normalInv}
\title{Combine independent p-values using the normal inversion method}
\usage{
  compute.normalInv(p, eps = 1e-06)
}
\arguments{
  \item{p}{a vector of independent p-values}

  \item{eps}{the minimal p-value considered (all p-values
  smaller will be set to this value)}
}
\description{
  Combine independent p-values using the normal inversion
  method
}
\examples{
p <- c(.1, .01)
compute.normalInv(p)
}
\author{
  Calin Voichita and Sorin Draghici
}
\seealso{
  \code{\link{pe}},\code{\link{compute.fischer}}
}

