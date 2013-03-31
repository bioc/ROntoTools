\name{peEdgeRenderInfo}
\alias{peEdgeRenderInfo}
\title{Extract edge render information from a \code{pePathway} object}
\usage{
  peEdgeRenderInfo(x, pos.col = "black", pos.lty = "solid",
    pos.ah = "vee", neg.col = "black", neg.lty = "dashed",
    neg.ah = "tee", zero.col = "lightgray",
    zero.lty = "dotted", zero.ah = "none")
}
\arguments{
  \item{x}{an object of class \code{\link{pePathway}}}

  \item{pos.col}{color of the edges with possitive weight}

  \item{pos.lty}{line type of the edges with possitive
  weight}

  \item{pos.ah}{arrow head of the edges with possitive
  weight}

  \item{neg.col}{color of the edges with negative weight}

  \item{neg.lty}{line type of the edges with negative
  weight}

  \item{neg.ah}{arrow head of the edges with negative
  weight}

  \item{zero.col}{color of the edges with zero weight}

  \item{zero.lty}{color of the edges with zero weight}

  \item{zero.ah}{color of the edges with zero weight}
}
\description{
  Extract edge render information from a \code{pePathway}
  object
}
\examples{
# load experiment
load(system.file("extdata/E-GEOD-21942.topTable.RData", package = "ROntoTools"))
fc <- top$logFC[top$adj.P.Val <= .01]
names(fc) <- top$entrez[top$adj.P.Val <= .01]
ref <- top$entrez

# load the set of pathways
kpg <- keggPathwayGraphs("hsa")
kpg <- setEdgeWeights(kpg)
kpg <- setNodeWeights(kpg, defaultWeight = 1)

# perform the pathway analysis
peRes <- pe(fc, graphs = kpg, ref = ref, nboot = 100, verbose = TRUE)

p <- peRes@pathways[[50]]
g <- layoutGraph(p@map, layoutType = "dot")
graphRenderInfo(g) <- list(fixedsize = FALSE)
edgeRenderInfo(g) <- peEdgeRenderInfo(p)
nodeRenderInfo(g) <- peNodeRenderInfo(p)
# notice the different type of edges in the graph (solid/dashed/dotted)
renderGraph(g)
}
\author{
  Calin Voichita and Sorin Draghici
}
\seealso{
  \code{\link{edgeRenderInfo}},\code{\link{par}}
}

