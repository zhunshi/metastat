#' INT transform
#' a function to permform INT-transformation
#' @param x numeric vector
#'
#' @return
#' @export
#'
#' @examples
INT <- function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}
