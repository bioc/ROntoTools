
#' @keywords internal
pf2col <- function(pf)
{
  retcol <- rep("white", length(pf))
  
  ### possitive
  i <- pf > 0
  #retcol[i] <- rainbow(256, start = 1/6, end = 2/6)[ceiling(pf[i] / max(pf[i]) * 256)]  
  retcol[i] <- colorRampPalette(c("white", "red"))(256)[ceiling(pf[i] / max(pf[i]) * 256)]
  
  i <- pf < 0
  #retcol[i] <- rainbow(256, start = 1, end = 1/6)[256:1][ceiling(abs(pf[i]) / max(abs(pf[i])) * 256)]
  retcol[i] <- colorRampPalette(c("white", "blue"))(256)[ceiling(abs(pf[i]) / max(abs(pf[i])) * 256)]
  
  names(retcol) <- names(pf)
  return(retcol)
}

#' @keywords internal
addNames <- function(x, nms)
{
  if (length(x) == length(nms))
    names(x) <- nms
  
  x <- rep(x, length(nms))
  names(x) <- nms
  
  return(x)
}

#' @keywords internal
compute.bootPV <- function(real, dist)
  ( sum(abs(dist - mean(dist)) > abs(real - mean(dist))) + 1 ) / (1 + length(dist))

#' @keywords internal
compute.fischer <- function(p)
  {k <- p[1] * p[2]; return(k-k*log(k))}

#' @keywords internal
compute.normalInv <- function(p)
  pnorm( (qnorm(p[1]) + qnorm(p[2])) / sqrt(2) )

#` @keywords internal
graph2ftM <- function(g)
{
  ind <- which(as.vector(as(g, "matrix"))!=0)
  
  from <- (ind-1) %% length(nodes(g)) + 1
  to <- (ind-1) %/% length(nodes(g)) + 1
  
  return(cbind(nodes(g)[from], nodes(g)[to]))
}

#' Compute alpha weights
#' 
#' Transform a vector of p-values into weights.
#' 
#' @details
#' 
#' Computes a set of weights from p-values using the formula \code{1-pv/threshold}.
#' 
#' @param pv vector of p-values
#' @param threshold the threshold value that was used to select DE genes
#' 
#' @author
#' 
#' Calin Voichita and Sorin Draghici
#'
#' @seealso \code{\link{pe}}
#' 
#' @examples
#' 
#' load(system.file("extdata/E-GEOD-21942.topTable.RData", package = "ROntoTools"))
#' 
#' head(alpha1MR(top$adj.P.Val))
#' 
#' @export
alpha1MR <- function(pv, threshold = max(pv))
{
  if (!is.numeric(pv))
    stop("pv is not numeric")
  if (threshold == 0)
    stop("threshold cannot be 0")
  
  return(1-pv/threshold)
}


#' Compute alpha weights
#' 
#' Transform a vector of p-values into weights.
#' 
#' @details
#' 
#' Computes a set of weights from p-values using the formula \code{-log10(pv/threshold)}.
#' 
#' @param pv vector of p-values
#' @param threshold the threshold value that was used to select DE genes
#' 
#' @author
#' 
#' Calin Voichita and Sorin Draghici
#'
#' @seealso \code{\link{pe}}
#' 
#' @examples
#' 
#' load(system.file("extdata/E-GEOD-21942.topTable.RData", package = "ROntoTools"))
#' 
#' head(alphaMLG(top$adj.P.Val))
#' 
#' @export
alphaMLG <- function(pv, threshold = max(pv))
{
  if (!is.numeric(pv))
    stop("pv is not numeric")
  if (threshold == 0)
    stop("threshold cannot be 0")
  
  return(-log10(pv/threshold))
}


