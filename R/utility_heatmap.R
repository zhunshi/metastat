#' customize pheatmap
#'
#' @param scc output of scc_self, scc_between, pcc_self, pcc_between
#' @param ReAdjPvalLowerTri if re-adjust pvalue (lower tri)
#' @param FilterRow FilterRow
#' @param FilterCol FilterCol
#' @param cellheight cellheight
#' @param cellwidth cellwidth
#' @param fontsize fontsize
#' @param fontsize_number fontsize_number
#' @param cluster_rows cluster_rows
#' @param cluster_cols cluster_cols
#' @param ... other parameters passed to pheatmap
#'
#' @return
#' @export
#'
#' @examples
mypheatmap <- function(scc,ReAdjPvalLowerTri=FALSE,FilterRow=F,FilterCol=F,cellheight = 15,cellwidth = 20,
                    fontsize = 10,fontsize_number=10,cluster_rows = T,cluster_cols = T,...){
  n <- ncol(scc)/3
  id <- (1:n)*3
  scc.c <- scc[,id-2,drop=F]
  scc.p <- scc[,id-1,drop=F]
  scc.f <- scc[,id,drop=F]
  colnames(scc.c) <- colnames(scc.f) <- colnames(scc.p) <- gsub("[- ]rho","",colnames(scc.c))

  if(ReAdjPvalLowerTri){
    p <- scc.p[lower.tri(scc.p)]
    f <- p.adjust(p,method = "BH")
    scc.f[lower.tri(scc.f)] <- f
    scc.f <- t(scc.f)
    scc.f[lower.tri(scc.f)] <- f
    scc.f <- as.data.frame(scc.f)
  }

  scc.num <- apply(scc.f, 2, function(x) ifelse(x<0.05,"*",""))
  id1 <- apply(scc.num, 1, function(x) all(x==""))
  id2 <- apply(scc.num, 2, function(x) all(x==""))
  if(FilterRow){
    scc.num <- scc.num[!id1,,drop=F]
  }
  if(FilterCol){
    scc.num <- scc.num[,!id2,drop=F]
  }
  scc.c <- scc.c[rownames(scc.num),colnames(scc.num),drop=F]

  paletteLength <- 50
  #myColor <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(paletteLength)
  myColor <- colorRampPalette(colors = c("navy","white","firebrick3"))(paletteLength)
  mybreak <- c(seq(min(scc.c), 0, length.out=ceiling(paletteLength/2) + 1),
               seq(max(scc.c)/paletteLength, max(scc.c), length.out=floor(paletteLength/2)))
  pheatmap::pheatmap(scc.c,display_numbers = scc.num,
           #cellheight = cellheight,cellwidth = cellwidth,
           #fontsize = fontsize,fontsize_number = fontsize_number,
           #cluster_rows = cluster_rows,cluster_cols = cluster_cols,
           color = myColor,
           breaks = mybreak,...)
}
