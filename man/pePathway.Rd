\name{pePathway}
\alias{pePathway}
\alias{pePathway-class}
\title{Class that encodes the result of Pathway-Express for a single pathway}
\description{
  Class that encodes the result of Pathway-Express for a
  single pathway
}
\section{Slots}{
  \describe{ \item{map}{an object of type graph (e.g.,
  \code{\link{graphNEL}}).} \item{input}{named vector of
  fold changes for genes on this pathway. The names of the
  genes are the orignal IDS used in the analysis}
  \item{ref}{vector of reference IDs on this pathway}
  \item{boot}{an object of class \code{boot} encoding the
  bootstrap information.} \item{PF}{the gene perturbation
  factors for all genes on the pathway, as computed by
  Pathway-Express.} \item{Acc}{the gene accumulations for
  all genes on the pathway, as computed by
  Pathway-Express.} }
}
\author{
  Calin Voichita and Sorin Draghici
}
\seealso{
  \code{\link{pe}} \code{\link{peRes}}
}

