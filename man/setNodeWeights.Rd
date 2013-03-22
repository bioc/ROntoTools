\name{setNodeWeights}
\alias{setNodeWeights}
\title{Set node weights}
\usage{
  setNodeWeights(graphList, weights = NULL,
    defaultWeight = 1)
}
\arguments{
  \item{graphList}{a list of \code{graph} (e.g.,
  \code{\link{graphNEL}}) objects}

  \item{weights}{named vector or matrix; if vector, the
  node is going to have the same weight in all graphs it
  appears; if matrix, the rows represent nodes and columns
  represent graphs and the node will have different weights
  in each pathway}

  \item{defaultWeight}{the default weight for all nodes not
  set by the parameter \code{weights}}
}
\value{
  The \code{graphList} with the node weights set.
}
\description{
  Set node weights
}
\examples{
# load the set of pathways
kpg <- keggPathwayGraphs("hsa")

kpg <- setNodeWeights(kpg)

nodeWeights(kpg[["path:hsa04110"]])
}
\author{
  Calin Voichita and Sorin Draghici
}

