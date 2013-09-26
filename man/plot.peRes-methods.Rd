\name{plot.peRes}
\alias{plot,peRes,character-method}
\alias{plot,peRes,missing-method}
\alias{plot.peRes}
\title{Plot Pathway-Express result}
\arguments{
  \item{x}{an object of type \code{\link{peRes}}}

  \item{y}{vector of two p-values names to be combined
  using \code{comb.pv.func} (default: \code{c("pAcc",
  "pORA")}).}

  \item{...}{Arguments to be passed to methods, such as
  \code{\link{par}}.}

  \item{comb.pv.func}{the function to combine the p-values
  - takes as input a vector of p-values and returns the
  combined p-value (default:
  \code{\link{compute.fisher}}).}

  \item{adjust.method}{the name of the method to adjust the
  p-value (see \code{\link{p.adjust}})}

  \item{threshold}{corrected p-value threshold}

  \item{eps}{any value smaller than this will be considered
  as \code{eps} (default: \code{1e-6}).}
}
\description{
  Display a two-way plot using two of the p-values from the
  Pathway-Express analysis.
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

plot(peRes)

plot(peRes, c("pPert","pORA"), comb.pv.func = compute.normalInv, threshold = .01)
}
\author{
  Calin Voichita and Sorin Draghici
}
\seealso{
  \code{\link{pe}}, \code{\link{Summary.peRes}},
  \code{\link{plot.pePathway}}
}

