#' Title
#'
#' @param PC output of prcomp function
#' @param x PC1 to display
#' @param y PC2 to display
#' @param Group dataframe with one column to group
#' @param main main title
#' @param point_size size of point
#' @param point_alpha alpha of point
#' @param n_label number of labels to display in every PCs
#'
#' @return ggplot object
#' @export
#'
#' @examples
PCbiplot <- function(PC, x="PC1", y="PC2",Group=NULL,main=NULL,point_size=3,point_alpha=1,n_label=5) {
  # points
  data <- data.frame(obsnames=row.names(PC$x), PC$x)
  # labels
  loading <- as.data.frame(PC$rotation[,1:2])
  loading <- abs(loading)
  dis_name <- rownames(loading)[order(loading$PC1,decreasing = T)][1:n_label]
  dis_name2 <- rownames(loading)[order(loading$PC2,decreasing = T)][1:n_label]
  dis_name <- union(dis_name,dis_name2)
  loading <- PC$rotation[dis_name,1:2]
  datapc <- data.frame(varnames=rownames(loading), loading)
  # arrows
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(
    datapc,
    v1 = .8 * mult * (get(x)),
    v2 = .8 * mult * (get(y))
  )
  # eig values
  eig <- PC$sdev ^ 2
  eig1 <- eig[1]/sum(eig)
  eig2 <- eig[2]/sum(eig)
  eig1 <- paste("PC1(",round(eig1*100,2),"%)",sep = "")
  eig2 <- paste("PC2(",round(eig2*100,2),"%)",sep = "")
  # plot
  if(!is.null(Group)){
    colnames(Group) <- "Group"
    Group <- Group[data$obsnames,,drop=F]
    data <- cbind(data,Group)
    p <- ggplot(data, aes_string(x=x, y=y))
    p <- p + geom_point(aes(color=Group),size=point_size,alpha=point_alpha)
  }else{
    p <- ggplot(data, aes_string(x=x, y=y))
    p <- p + geom_point(size=point_size,alpha=point_alpha)
  }
  p <- p +
    ggrepel::geom_text_repel(
      data=datapc,
      aes(x=v1, y=v2, label=varnames),
      size = 5, color="black"
    )+
    geom_segment(
      data=datapc,
      aes(x=0, y=0, xend=v1, yend=v2),
      arrow=arrow(length=unit(0.2,"cm")),
      alpha=0.7, color="black",size=1
    ) +
    scale_color_brewer(palette = "Set1")+
    guides(colour = guide_legend(override.aes = list(size=3)))+
    theme_bw()+
    theme(
      axis.title = element_text(size = 13,color = "black"),
      axis.text = element_text(size = 12,color = "black"),
      panel.grid = element_blank(),
      legend.title = element_text(size = 13,colour = "black"),
      legend.text = element_text(size = 12,color = "black"),
      plot.title = element_text(size = 13,color = "black",hjust=0.5))+
    xlab(eig1)+ylab(eig2)
  if(!is.null(main)){
    p <- p + labs(title=main)
  }
  p
}

