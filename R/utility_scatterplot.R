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
scatterplot <- function(dat,x,y,group=NULL){
  dat <- dat[!is.na(dat[,x]),,drop=F]
  s0 <- cor.test(dat[,x],dat[,y],method = "s")
  if(is.null(group)){
    lab <- paste0("scc rho = ",round(s0$estimate,3),"; p = ",formatC(s0$p.value,digits = 2))
    p <- ggplot(dat,aes_string(x,y))+
      geom_point(size=0.8,alpha=0.5)+
      geom_smooth(method = "lm",se=F)+
      annotate("text",x=-Inf,y=Inf,vjust=1.2,hjust=0,label=lab,size=3)+
      theme_bw()
  }else{
    lst_levels <- levels(as.factor(dat[,group]))
    s1 <- cor.test(dat[dat[,group]==lst_levels[1],x],dat[dat[,group]==lst_levels[1],y],method = "s")
    s2 <- cor.test(dat[dat[,group]==lst_levels[2],x],dat[dat[,group]==lst_levels[2],y],method = "s")
    lab <- paste0(
      "Total scc rho = ",round(s0$estimate,3),"; p = ",formatC(s0$p.value,digits = 2),"\n",
      paste0(lst_levels[1]," scc rho = "),round(s1$estimate,3),"; p = ",formatC(s1$p.value,digits = 2),"\n",
      paste0(lst_levels[2]," scc rho = "),round(s2$estimate,3),"; p = ",formatC(s2$p.value,digits = 2)
    )
    p <-
      ggplot(dat,aes_string(x,y,color=group))+
      geom_point(size=0.8,alpha=0.5)+
      geom_smooth(method = "lm",se=F)+
      geom_smooth(data=dat,method = "lm",se=F,aes_string(x,y),color="black")+
      annotate("text",x=-Inf,y=Inf,vjust=1.2,hjust=0,label=lab,size=3)+
      theme_bw()
  }
  p
}
