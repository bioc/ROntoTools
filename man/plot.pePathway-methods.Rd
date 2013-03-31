\name{plot.pePathway}
\alias{plot,pePathway,character-method}
\alias{plot,pePathway,missing-method}
\alias{plot.pePathway}
\title{Plot pathway level statistics}
\arguments{
  \item{x}{an object of type \code{\link{pePathway}}}

  \item{y}{if provided, the factor to be ploted (either
  \code{Acc} (default) or \code{Pert}; see
  \code{\link{pePathway}})}

  \item{...}{Arguments to be passed to methods, such as
  \code{\link{par}}}

  \item{type}{type of plot (either \code{two.way} (default)
  or \code{boot})}

  \item{eps}{any value smaller than this will be ploted as
  0}
}
\description{
  Display graphical representation of pathway level
  statistic like: i) two way comparison between the
  measured expression change and one of the factors
  computed by Pathway-Express (\code{\link{pe}}) or ii) the
  boostrap statistics of the same factors.
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

# perform the pathway analysis (for more accurate results use nboot = 2000)
peRes <- pe(fc, graphs = kpg, ref = ref, nboot = 100, verbose = TRUE)

plot(peRes@pathways[[50]])

plot(peRes@pathways[[50]], "Pert", main = "Perturbation factor")

plot(peRes@pathways[[50]], type = "boot")

plot(peRes@pathways[[50]], "Pert", type = "boot", main = "Perturbation factor")
}
\author{
  Calin Voichita and Sorin Draghici
}
\seealso{
  \code{\link{pe}}, \code{\link{plot.peRes}},
  \code{\link{peNodeRenderInfo}},
  \code{\link{peEdgeRenderInfo}}
}

