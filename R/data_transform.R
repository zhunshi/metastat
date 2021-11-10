#' Title
#'
#' @param x numeric vector. Input
#'
#' @return
#' @export
#'
#' @examples
INT <- function(x){qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
