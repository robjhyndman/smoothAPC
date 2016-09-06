#' Prints the result of smoothing (object of class \code{sm2D}).
#'
#' @param x Result of smoothing (object of class \code{sm2D}).
#' @param ... Other parameters. They are currently ignored.
#' @examples
#' \dontrun{
#'
#' library(demography)
#' m = log(fr.mort$rate$female[1:30, 150:160])
#' sm = autoDemogSmooth(m)
#' sm
#'
#' }
#' @author Alexander Dokumentov
#' @export

print.sm2D = function(x, ...)
{
  cat("\nObject of class sm2D")
  return(invisible(x))
}
