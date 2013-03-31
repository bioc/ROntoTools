\name{peNodeRenderInfo}
\alias{peNodeRenderInfo}
\title{Extract node render information from a \code{pePathway} object}
\usage{
  peNodeRenderInfo(x, y = "Pert", input.shape = "box",
    default.shape = "ellipse", pos.col = "red",
    neg.col = "blue", zero.col = "white")
}
\arguments{
  \item{x}{an object of class \code{\link{pePathway}}}

  \item{y}{a string representing the factor to be
  represented (\code{Pert, Acc} or \code{input}; see
  \code{\link{pePathway}})}

  \item{input.shape}{shape of nodes that have measured
  expression change}

  \item{default.shape}{shape of all other nodes}

  \item{pos.col}{color of nodes with a positive \code{y}
  factor}

  \item{neg.col}{color of nodes with a negative \code{y}
  factor}

  \item{zero.col}{color of nodes with the \code{y} factor
  equal to zero}
}
\description{
  Extract node render information from a \code{pePathway}
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
# notice the different type of nodes in the graph (box/circle)
# the color of each node represents the perturbation (red = positive, blue = negative)
# the shade represents the stregth of the perturbation
renderGraph(g)

nodeRenderInfo(g) <- peNodeRenderInfo(p, "Acc")
# now, the color of each node represents the accumulation (red = positive, blue = negative)
# notice that square nodes with no parents have no accumulation
renderGraph(g)
}
\author{
  Calin Voichita and Sorin Draghici
}
\seealso{
  \code{\link{nodeRenderInfo}},\code{\link{par}}
}

