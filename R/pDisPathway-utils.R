
compute.ppDis <- function(pDisPath)
  ifelse( !all(pDisPath@boot$t[,"tpDis"] == 0),
          compute.bootPV(pDisPath@boot$t0["tpDis"], pDisPath@boot$t[,"tpDis"]),
          NA)


compute.pORA <- function(pDisPath, inputSize, refSize)
  phyper(q = length(pDisPath@input)-1,
         m = length(pDisPath@ref),
         n = refSize-length(pDisPath@ref),
         k = inputSize,
         lower.tail = FALSE) 


get.totalpDis <- function(pDisPath)
  as.numeric(pDisPath@boot$t0["tpDis"])

get.totalpDisNorm <- function(pDisPath)
  ifelse( !all(pDisPath@boot$t[,"tpDis"] == 0),
          as.numeric((pDisPath@boot$t0["tpDis"] - mean(pDisPath@boot$t[,"tpDis"])) / sd(pDisPath@boot$t[,"tpDis"])),
          NA)

