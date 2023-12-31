% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.pDisRes.R
\name{summary.pDisRes}
\alias{summary.pDisRes}
\title{Summarize the results of a primary dis-regulation (pDis) analysis}
\usage{
summary.pDisRes(object, ..., pathNames = NULL, totalpDis = TRUE, normalize = TRUE, 
 ppDis = TRUE, pORA = TRUE, 
 comb.pv = c("ppDis", "pORA"), comb.pv.func = compute.fisher,
 order.by = "pComb", adjust.method = "fdr")
}
\arguments{
\item{object}{pDis analysis result object obtained using \code{\link{pDis}}}

\item{...}{ignored}

\item{pathNames}{named vector of pathway names; the names of the vector are the IDs of the pathways}

\item{totalpDis}{boolean value indicating if the total primary dis-regulation should be computed}

\item{normalize}{boolean value indicating if normalization with regards to the boostrap simulations should be performed on totalpDis}

\item{ppDis}{boolean value indicating if the significance of the total primary dis-regulation in regards to the bootstrap permutations should be computed}

\item{pORA}{boolean value indicating if the over-represtation p-value should be computed}

\item{comb.pv}{vector of the p-value names to be combine (any of the above p-values)}

\item{comb.pv.func}{the function to combine the p-values; takes as input a vector of p-values and returns the combined p-value}

\item{order.by}{the name of the p-value that is used to order the results}

\item{adjust.method}{the name of the method to adjust the p-value (see \link{p.adjust})}
}
\description{
Summarize the results of a primary dis-regulation (pDis) analysis
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
pDisRes <- pDis(fc, graphs = kpg, ref = ref, nboot = 100, verbose = TRUE)

# obtain summary of results
head(summary(pDisRes))

kpn <- keggPathwayNames("hsa")

head(summary(pDisRes))

head(summary(pDisRes, pathNames = kpn, totalpDis = FALSE, 
            pORA = FALSE, comb.pv = NULL, order.by = "pDis"))

}
\seealso{
\code{\link{pDis}}
}

