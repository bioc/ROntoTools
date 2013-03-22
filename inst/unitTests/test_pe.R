test_pe <- function()
{
  load(system.file("extdata/E-GEOD-21942.topTable.RData", package = "ROntoTools"))
  fc <- top$logFC[top$adj.P.Val <= .01]
  names(fc) <- top$entrez[top$adj.P.Val <= .01]
  pv <- top$P.Value[top$adj.P.Val <= .01]
  names(pv) <- top$entrez[top$adj.P.Val <= .01]
  ref <- top$entrez
  
  # load the set of pathways
  kpg <- keggPathwayGraphs("hsa")
  kpg <- setEdgeWeights(kpg)
  kpg <- setNodeWeights(kpg, weights = alphaMLG(pv), defaultWeight = 1)
  
  # perform the pathway analysis
  # in order to obtain accurate results the number of boostraps, nboot, should
  # be increase to a number like 2000
  peRes <- pe(fc, graphs = kpg, ref = ref, nboot = 100, verbose = TRUE)
  
  checkTrue(all(names(peRes@pathways) %in% names(kpg)))
  checkTrue(all(peRes@input %in% fc))
  checkTrue(all(peRes@ref %in% ref))
  checkEquals(peRes@cutOffFree, FALSE)
}