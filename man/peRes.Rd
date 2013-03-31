\name{peRes}
\alias{peRes}
\alias{peRes-class}
\title{Pathway-Express result class}
\description{
  This class is used to encode the results of the pathway
  analysis performed by the function \code{\link{pe}}.
}
\details{
  The slots \code{input} and \code{ref} record global
  information related to the whole analysis, while the
  \code{pathways} slot records the specific results as
  \code{\link{pePathway}} for each one of the pathways used
  in the analysis.
}
\section{Slots}{
  \describe{ \item{\code{pathways}:}{A list of
  \code{\link{pePathway}} objects.}
  \item{\code{input}:}{named vector of fold changes used
  for the analysis. The names of the vector are the IDs
  originaly used.} \item{\code{ref}:}{character vector
  containing the IDs used as reference in the analysis.}
  \item{\code{cutOffFree}:}{boolean value indicating if a
  cut-of-free analysis has been performed.} }
}
\author{
  Calin Voichita and Sorin Draghici
}
\seealso{
  \code{\link{pe}}, \code{\link{pePathway}}
}

