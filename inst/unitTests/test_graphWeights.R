test_setNodeWeights <- function()
{
  # load the set of pathways
  kpg <- keggPathwayGraphs("hsa")
  
  wt <- runif(length(nodes(kpg[["path:hsa04110"]])), 0, 1)
  names(wt) <- nodes(kpg[["path:hsa04110"]])
  
  kpg <- setNodeWeights(kpg, wt)
  
  checkEquals(nodeWeights(kpg[["path:hsa04110"]]), wt)
}

test_setEdgeWeights <- function()
{
  # load the set of pathways
  kpg <- keggPathwayGraphs("hsa")
  
  edgeWeightByType <- list(activation = 1, inhibition = -1, expression = 1, repression = -1)
  kpg <- setEdgeWeights(kpg)
  
  checkEquals(
    unlist(edgeWeights(kpg[["path:hsa04110"]], "hsa:7029")),
    unlist(edgeWeightByType[unlist(edgeData(kpg[["path:hsa04110"]], "hsa:7029", attr = "subtype"))]),
    checkNames = FALSE
  )
  
}