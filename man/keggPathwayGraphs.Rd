\name{keggPathwayGraphs}
\alias{keggPathwayGraphs}
\title{Download and parse KEGG pathway data}
\usage{
  keggPathwayGraphs(organism = "hsa",
    targRelTypes = c("GErel", "PCrel", "PPrel"),
    relPercThresh = 0.9, nodeOnlyGraphs = FALSE,
    updateCache = FALSE, verbose = TRUE)
}
\arguments{
  \item{organism}{organism code as defined by KEGG}

  \item{targRelTypes}{target relation types}

  \item{relPercThresh}{percentage of the number of relation
  types over all possible realtions in the pathway}

  \item{nodeOnlyGraphs}{allow graphs with no edges}

  \item{updateCache}{re-download KEGG data}

  \item{verbose}{show progress of downloading and parsing}
}
\value{
  A list of \code{\link{graphNEL}} objects encoding the
  pathway information.
}
\description{
  Download and parse KEGG pathway data
}
\examples{
# The pathway cache provided as part of the pathway contains only the
# pathways that passed the default filtering. We recommend, re-downloading
# the pathways using the updateCache parameter
kpg <- keggPathwayGraphs("hsa")

# to update the pathway cache for human run:
# kpg <- keggPathwayGraphs("hsa", updateCache = TRUE)
# this is time consuming and depends on the available bandwith.

head(names(kpg))

kpg[["path:hsa04110"]]
head(nodes(kpg[["path:hsa04110"]]))
head(edges(kpg[["path:hsa04110"]]))
}
\author{
  Calin Voichita and Sorin Draghici
}
\seealso{
  \code{\link{keggPathwayNames}}
}

