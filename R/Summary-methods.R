
#' Summarize the results of a Pathway-Express analysis
#' 
#' @param x Pathway-Express analysis result object obtained using \code{\link{pe}}
#' @param ... see \code{\link{summary.peRes}}
#' @param na.rm ignored
#'
#' @aliases Summary,peRes-method
#'
#' @export
setMethod("Summary", c("x" = "peRes"),
  function(x, ..., na.rm = FALSE) {
    return(summary.peRes(x, ...));
  }
)

#' Summarize the results of a Pathway-Express analysis
#' 
#' @param x Primary dis-regulation analysis result object obtained using \code{\link{pDis}}
#' @param ... see \code{\link{summary.pDisRes}}
#' @param na.rm ignored
#'
#' @aliases Summary,pDisRes-method
#'
#' @export
setMethod("Summary", c("x" = "pDisRes"),
  function(x, ..., na.rm = FALSE) {
    return(summary.pDisRes(x, ...));
  }
)
