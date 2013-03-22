\name{pe}
\alias{pe}
\title{Pathway-Express: Pathway analysis of signaling pathways}
\usage{
  pe(x, graphs, ref = NULL, nboot = 2000, verbose = TRUE,
    cluster = NULL, seed = NULL)
}
\arguments{
  \item{x}{named vector of log fold changes for the
  differentially expressed genes; \code{names(x)} must use
  the same id's as \code{ref} and the nodes of the
  \code{graphs}}

  \item{graphs}{list of pathway graphs as objects of type
  \code{graph} (e.g., \code{\link{graphNEL}}); the graphs
  must be weighted graphs (i.e., have an atrribute
  \code{weight} for both nodes and edges)}

  \item{ref}{the reference vector for all genes in the
  analysis; if the reference is not provided or it is
  identical to \code{names(x)} a cut-off free analysis is
  performed}

  \item{nboot}{number of bootstrap iterations}

  \item{verbose}{print progress output}

  \item{cluster}{a cluster object created by makeCluster
  for parallel computations}

  \item{seed}{an integer value passed to set.seed() during
  the boostrap permutations}
}
\value{
  An object of class \code{\link{peRes}}.
}
\description{
  Pathway-Express: Pathway analysis of signaling pathways
}
\details{
  See details in the cited articles.
}
\examples{
# load a multiple sclerosis study (public data available in Array Express
# ID: E-GEOD-21942)
# This file contains the top table, produced by the limma package with
# added gene information. All the probe sets with no gene associate to them,
# have been removed. Only the most significant probe set for each gene has been
# kept (the table is already ordered by p-value)
# The table contains the expression fold change and signficance of each
# probe set in peripheral blood mononuclear cells (PBMC) from 12 MS patients
# and 15 controls.
load(system.file("extdata/E-GEOD-21942.topTable.RData", package = "ROntoTools"))
head(top)

# select differentially expressed genes at 1\% and save their fold change in a
# vector fc and their p-values in a vector pv
fc <- top$logFC[top$adj.P.Val <= .01]
names(fc) <- top$entrez[top$adj.P.Val <= .01]

pv <- top$P.Value[top$adj.P.Val <= .01]
names(pv) <- top$entrez[top$adj.P.Val <= .01]

# alternativly use all the genes for the analysis
# NOT RUN:
# fc <- top$logFC
# names(fc) <- top$entrez

# pv <- top$P.Value
# names(pv) <- top$entrez

# get the reference
ref <- top$entrez

# load the set of pathways
kpg <- keggPathwayGraphs("hsa")

# set the beta information (see the citated documents for meaning of beta)
kpg <- setEdgeWeights(kpg)

# inlcude the significance information in the analysis (see Voichita:2012
# for more information)
# set the alpha information based on the pv with one of the predefined methods
kpg <- setNodeWeights(kpg, weights = alphaMLG(pv), defaultWeight = 1)

# perform the pathway analysis
# in order to obtain accurate results the number of boostraps, nboot, should
# be increase to a number like 2000
peRes <- pe(fc, graphs = kpg, ref = ref, nboot = 100, verbose = TRUE)

# obtain summary of results
head(Summary(peRes))
}
\author{
  Calin Voichita and Sorin Draghici
}
\references{
  Voichita C., Donato M., Draghici S.: "Incorporating gene
  significance in the impact analysis of signaling
  pathways", IEEE Machine Learning and Applications
  (ICMLA), 2012 11th International Conference on, Vol. 1,
  p.126-131, 2012

  Tarca AL., Draghici S., Khatri P., Hassan SS., Kim J.,
  Kim CJ., Kusanovic JP., Romero R.: "A Signaling Pathway
  Impact Analysis for Microarray Experiments", 2008,
  Bioinformatics, 2009, 25(1):75-82.

  Khatri P., Draghici S., Tarca AL., Hassan SS., Romero R.:
  "A system biology approach for the steady-state analysis
  of gene signaling networks". Progress in Pattern
  Recognition, Image Analysis and Applications, Lecture
  Notes in Computer Science. 4756:32-41, November 2007.

  Draghici S., Khatri P., Tarca A.L., Amin K., Done A.,
  Voichita C., Georgescu C., Romero R.: "A systems biology
  approach for pathway level analysis". Genome Research,
  17, 2007.
}
\seealso{
  \code{\link{Summary}} \code{\link{keggPathwayGraphs}}
  \code{\link{setNodeWeights}} \code{\link{setEdgeWeights}}
}

