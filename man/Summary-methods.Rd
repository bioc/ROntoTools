\name{Summary}
\alias{Summary,peRes-method}
\alias{Summary.peRes}
\title{Summarize the results of a Pathway-Express analysis}
\arguments{
  \item{x}{Pathways-Express result object obtained using
  \code{\link{pe}}}

  \item{pathNames}{named vector of pathway names; the names
  of the vector are the IDs of the pathways}

  \item{totalAcc}{boolean value indicating if the total
  accumulation should be computed}

  \item{totalPert}{boolean value indicating if the total
  perturbation should be computed}

  \item{normalize}{boolean value indicating if
  normalization with regards to the boostrap simulations
  should be performed on totalAcc and totalPert}

  \item{pPert}{boolean value indicating if the significance
  of the total perturbation in regards to the bootstrap
  permutations should be computed}

  \item{pAcc}{boolean value indicating if the significance
  of the total accumulation in regards to the bootstrap
  permutations should be computed}

  \item{pORA}{boolean value indicating if the
  over-represtation p-value should be computed}

  \item{comb.pv}{vector of the p-value names to be combine
  (any of the above p-values)}

  \item{comb.pv.func}{the function to combine the p-values;
  takes as input a vector of p-values and returns the
  combined p-value}

  \item{order.by}{the name of the p-value that is used to
  order the results}

  \item{adjust.method}{the name of the method to adjust the
  p-value (see \link{p.adjust})}
}
\description{
  Summarize the results of a Pathway-Express analysis
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

# obtain summary of results
head(Summary(peRes))

kpn <- keggPathwayNames("hsa")

head(Summary(peRes))

head(Summary(peRes, pathNames = kpn, totalAcc = FALSE, totalPert = FALSE,
             pAcc = FALSE, pORA = FALSE, comb.pv = NULL, order.by = "pPert"))
}
\seealso{
  \code{\link{pe}}
}

