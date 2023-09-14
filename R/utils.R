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

#' Combine independent p-values using the Fisher method
#' 
#' @param p a vector of independent p-values
#' @param eps the minimal p-value considered (all p-values smaller will be set to this value), default value 1e-24
#' 
#' @return the combined p-value
#' 
#' @author Calin Voichita and Sorin Draghici
#' 
#' @references
#' 
#' Tarca AL., Draghici S., Khatri P., Hassan SS., Kim J., Kim CJ., Kusanovic JP., Romero R.: "A Signaling Pathway Impact Analysis for Microarray Experiments", 2008, Bioinformatics, 2009, 25(1):75-82.
#' 
#' @seealso \code{\link{pe}},\code{\link{compute.normalInv}}
#' 
#' @examples
#' 
#' p <- c(.1, .01)
#' compute.fisher(p)
#' 
#' @export
compute.fisher <- function(p, eps = 1e-24) #default eps changed 2023-09-14 from 1e-6 to 1e-24
{
  stopifnot(any(p >= 0 & p<=1))  
  p[p < eps] <- eps
  
  k <- prod(p); 
  return(k-k*log(k))
}

#' Combine independent p-values using the normal inversion method
#' 
#' @param p a vector of independent p-values
#' @param eps the minimal p-value considered (all p-values smaller will be set to this value)
#' 
#' @return the combined p-value
#' 
#' @author Calin Voichita and Sorin Draghici
#' 
#' @references
#' 
#' Tarca AL., Draghici S., Romero R.: "A Mmore Specific Method To Combine Perturbation and Over-representation Evidence in Pathway Analysis", PSB 2010 poster.
#' 
#' @seealso \code{\link{pe}},\code{\link{compute.fisher}}
#' 
#' @examples
#' 
#' p <- c(.1, .01)
#' compute.normalInv(p)
#' 
#' @export
compute.normalInv <- function(p, eps = 1e-6)
{
  stopifnot(any(p >= 0 & p<=1))  
  p[p < eps] <- eps
  return(pnorm(sum(sapply(p, qnorm)) / sqrt(length(p))))
}

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


