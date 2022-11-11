#' customize pheatmap
#'
#' @param scc output of scc_self, scc_between, pcc_self, pcc_between
#' @param ReAdjPvalLowerTri if re-adjust pvalue (lower tri)
#' @param maxRange Whether to keep the maximum and minimum values the same, default TRUE
#' @param trans Whether to transpose the matrix
#' @param Psignificance If display sigificance for uncorrected P values
#' @param color fill colors (eg. colorRampPalette(colors = c("navy", "white", "firebrick3"))[50])
#' @param ... other parameters passed to pheatmap
#'
#' @return
#' @export
#'
#' @examples
pheatmap.scc <- function (scc,maxRange=TRUE,ReAdjPvalLowerTri=FALSE,trans=FALSE, color=NULL, Psignificance=FALSE, ...) {
  n <- ncol(scc)/3
  id <- (1:n) * 3
  scc.c <- scc[, id - 2, drop = F]
  scc.p <- scc[, id - 1, drop = F]
  scc.f <- scc[, id, drop = F]
  colnames(scc.c) <- colnames(scc.f) <- colnames(scc.p) <- gsub("[\\. ]rho", "", colnames(scc.c))
  # if re-adjusted pvalues in lower-triangles
  if(ReAdjPvalLowerTri){
    p <- scc.p[lower.tri(scc.p)]
    f <- p.adjust(p,method = "BH")
    scc.f[lower.tri(scc.f)] <- f
    scc.f <- t(scc.f)
    scc.f[lower.tri(scc.f)] <- f
    scc.f <- as.data.frame(scc.f)
  }
  if(Psignificance){
    scc.num <- apply(scc.p, 2, function(x) ifelse(x < 0.05, "#", ""))
    scc.num[scc.f<0.05] <- "*"
  }else{
    scc.num <- apply(scc.f, 2, function(x) ifelse(x < 0.05, "*", ""))
  }

  # colors
  paletteLength <- 50
  if (is.null(color)){
    myColor <- colorRampPalette(colors = c("navy", "white", "firebrick3"))(paletteLength)
  }else{
    myColor <- color
  }
  # if use max vaulues in color bars
  if(maxRange){
    v <- max(abs(scc.c))
    mybreak <- c(
      seq(-v, 0, length.out = ceiling(paletteLength/2) + 1),
      seq(v/paletteLength, v, length.out = floor(paletteLength/2))
    )
  }else{
    mybreak <- c(
      seq(min(scc.c), 0, length.out = ceiling(paletteLength/2) + 1),
      seq(max(scc.c)/paletteLength, max(scc.c), length.out = floor(paletteLength/2))
    )
  }
  # if transpose
  if (trans){
    scc.c <- as.data.frame(t(scc.c))
    scc.num <- as.data.frame(t(scc.num))
  }
  pheatmap::pheatmap(scc.c, display_numbers = scc.num, color = myColor, breaks = mybreak, ...)
}

