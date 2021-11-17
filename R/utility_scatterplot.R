#' scatterplot
#' customized scatterplot
#' @param dat dataframe, input
#' @param x string, lab for x
#' @param y string, lab for y
#' @param group string, lab for grouping, default NULL
#'
#' @return
#' @export
#'
#' @examples
scatterplot <- function (dat, x, y, group = NULL)
{
  dat <- dat[!is.na(dat[, x]), , drop = F]
  s0 <- cor.test(dat[, x], dat[, y], method = "s")
  if (is.null(group)) {
    lab <- paste0("scc rho = ", round(s0$estimate, 3), "; p = ",
                  formatC(s0$p.value, digits = 2))
    p <- ggplot(dat, aes_string(x, y)) +
      geom_point(size = 0.8, alpha = 0.5) +
      geom_smooth(method = "lm", se = F) +
      annotate("text", x = -Inf, y = Inf, vjust = 1.2,
               hjust = 0, label = lab, size = 3) + theme_bw()
  }
  else {
    lst_levels <- levels(as.factor(dat[, group]))
    lab <- c()
    for(l in lst_levels){
      a <- cor.test(dat[dat[,group]==l,x],dat[dat[,group]==l,y],method = "s")
      lab <- c(lab,paste0(l," scc rho=",round(a$estimate,3),"; p=",formatC(a$p.value, digits = 2)))
    }
    lab <- paste(lab,collapse = "\n")
    p <- ggplot(dat, aes_string(x, y, color = group)) +
      geom_point(size = 0.8, alpha = 0.5) +
      geom_smooth(method = "lm",se = F) +
      geom_smooth(data = dat, method = "lm",se = F, aes_string(x, y), color = "black") +
      annotate("text", x = -Inf, y = Inf, vjust = 1.2, hjust = 0, label = lab,size = 3) +
      theme_bw()
  }
  p
}
