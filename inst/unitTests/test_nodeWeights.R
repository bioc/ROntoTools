test_nodeWeights <- function()
{
  V <- LETTERS[1:4]
  g <- graphNEL(nodes = V, edgemode = "directed")
  
  checkEquals(length(nodeWeights(g)), length(nodes(g)))
  checkEquals(nodeWeights(g, "B"), 1, checkNames=FALSE)
  checkEquals(nodeWeights(g, "A", attr = "WT", default = 3), 3, checkNames=FALSE)
}