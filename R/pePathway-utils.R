# setMethod("plot", signature(x="pePathway", y="ANY"),
#           function(x, y, ...){
#             
#             nShape <- rep("circle", length(nodes(x@map)))
#             names(nShape) <- nodes(x@map)
#             nShape[names(x@input)] <- "box"
#             
#             m <- as(x@map, "matrix")
#             nodeToPlot <- nodes(x@map)[(colSums(abs(m)) != 0) & (rowSums(abs(m)) != 0)]
#             
#             plot(x@map,  nodeAttrs = list(
#               fillcolor = pf2col(x@PF), 
#               fontsize = addNames(340.0, nodes(x@map)),
#               shape = nShape
#               ))
#           })

compute.pAcc <- function(pePath)
  ifelse( !all(pePath@boot$t[,"tAcc"] == 0),
    compute.bootPV(pePath@boot$t0["tAcc"], pePath@boot$t[,"tAcc"]),
          NA)

compute.pPert <- function(pePath)
  ifelse( !all(pePath@boot$t[,"tPert"] == 0),
          compute.bootPV(pePath@boot$t0["tPert"], pePath@boot$t[,"tPert"]),
          NA)


compute.pORA <- function(pePath, inputSize, refSize)
  phyper(q = length(pePath@input)-1,
         m = length(pePath@ref),
         n = refSize-length(pePath@ref),
         k = inputSize,
         lower.tail = FALSE) 

get.totalAcc <- function(pePath)
  as.numeric(pePath@boot$t0["tAcc"])

get.totalAccNorm <- function(pePath)
  ifelse( !all(pePath@boot$t[,"tAcc"] == 0),
          as.numeric((pePath@boot$t0["tAcc"] - mean(pePath@boot$t[,"tAcc"])) / sd(pePath@boot$t[,"tAcc"])),
          NA)

get.totalPert <- function(pePath)
  as.numeric(pePath@boot$t0["tPert"])

get.totalPertNorm <- function(pePath)
  ifelse( !all(pePath@boot$t[,"tPert"] == 0),
          as.numeric((pePath@boot$t0["tPert"] - mean(pePath@boot$t[,"tPert"])) / sd(pePath@boot$t[,"tPert"])),
          NA)

