#' Split metaphlan2 profile by given rank
#'
#' Split metaphlan2 profile by given rank
#'
#' @param pro     input combined metaphlan2 profile
#' @param rank    input a given rank from [a,p,c,o,f,g,s],default a for all ranks
#'
#' @return dataframe or list
#'
#' @export
SplitMetaphlan2Profile <- function(pro,rank="a"){
  rn <- rownames(pro)
  lst <- list(
    Kingdom = grep("\\|",rn,invert = T),
    Phylum = intersect(grep("p__",rn), grep("c__",rn,invert = T)),
    Class = intersect(grep("c__",rn), grep("o__",rn,invert = T)),
    Order = intersect(grep("o__",rn), grep("f__",rn,invert = T)),
    Family = intersect(grep("f__",rn), grep("g__",rn,invert = T)),
    Genus = intersect(grep("g__",rn), grep("s__",rn,invert = T)),
    Species = intersect(grep("s__",rn), grep("t__",rn,invert = T))
  )

  lst <- lapply(lst,function(x,rn) rn[x],rn=rn)
  datlst <- lapply(lst,function(x,pro) pro[x,],pro=pro)

  if(rank=="a"){
    out <- datlst
  }else{
    out <- datlst[[rank]]
  }
  return(out)
}

#' Filter metaphlan2 profile
#'
#' Filter metaphlan2 profile
#'
#' @param pro                  input splited metaphlan2 profile in a rank, e.g. genus
#' @param include              any kingdom included, default k__Bacteria and K__Archaea. The prefix "K__" is needed.
#' @param remove_unclassified  if remove unclasssified features in profiles.
#' @param remove_noname        if remove features with noname.
#' @param normalize            if normalization after above process.
#'
#' @return dataframe
#'
#' @export
FilterMetaphlan2Profile <- function(pro,include=c("k__Bacteria","k__Archaea"),remove_unclassified=TRUE,remove_noname=TRUE,normalize=TRUE){
  name_split <- strsplit(rownames(pro),"\\|")
  name_k <- sapply(name_split,function(x) x[1])
  pro <- pro[name_k%in%include,]
  rownames(pro) <- gsub(".*__","",rownames(pro))

  # normalization
  if(normalize){
    pro <- sweep(pro,2,apply(pro,2,sum),"/")
  }

  if(remove_unclassified){
    pro <- pro[!grepl("unclassified",rownames(pro)),,drop=F]
  }
  if(remove_noname){
    pro <- pro[!grepl("noname",rownames(pro)),,drop=F]
  }
  pro
}

#' Standarize species name in  metaphlan2 profile
#'
#' A function to standarize species names in metaphlan2 profile
#'
#' @param x           vector of species name
#' @param italic      if use italic face
#'
#' @return value      vector of standarized species name
#'
#' @export
standarizeSpeciesName <- function(x){
  x <- gsub("s__","",x)
  x <- gsub("_"," ",x)
  a <- x
  b <- gsub("(\\d) (\\d)","\\1_\\2",a)
  while(!all(a==b)){
    a <- b
    b <- gsub("(\\d) (\\d)","\\1_\\2",a)
  }
  return(b)
}


#' Barplot for humann2 profile
#'
#' A function to implement barplot in humann2
#'
#' @param dat        full combined profile
#' @param group      group information with one column
#' @param pathwayID  pathway ID
#'
#' @return value      vector of standarized species name
#'
#' @export
Barplot_Pathway2Species <- function(dat,group,pathwayID){
  # Extract data for given pathwayID
  dat <- dat[grep(pathwayID,rownames(dat)),]
  FullName = rownames(dat)[1]
  dat <- dat[grep("\\|",rownames(dat)),]

  # Rename contributors (species level)
  rownames(dat) <- gsub(".*s__","",rownames(dat))
  rownames(dat) <- gsub("_"," ",rownames(dat))

  # order features
  b <- data.frame(
    Max = apply(dat,1,mean)
  ) %>%
    rownames_to_column(var="Species") %>%
    arrange(desc(Max))
  order_feature <- b$Species[1:7]
  if(nrow(b)>7){
    dat <- dat %>%
      rownames_to_column(var="Species") %>%
      mutate(
        Species = replace(Species, !Species %in% order_feature, "Other")
      ) %>%
      group_by(Species)  %>%
      summarise_all(sum) %>%
      column_to_rownames(var="Species")
    order_feature <- c(order_feature, "Other")
  }

  # order samples
  a <- data.frame(Sum = colSums(dat))
  group <- group[rownames(a),,drop=F]
  a <- cbind(a,group) %>% rownames_to_column(var="SampleID")
  a <- a[order(a$Sum,decreasing = T),]
  a <- a[order(a$Group),]
  order_sampels <- as.character(a$SampleID)
  a$y = ""
  a$SampleID <- factor(a$SampleID,levels = order_sampels)

  # Dataframe transform
  dat2 <- dat %>%
    rownames_to_column(var="Species") %>%
    gather(Samples,Abundance,-Species) %>%
    mutate(
      Species = factor(Species,levels = order_feature),
      Samples = factor(Samples,levels = order_sampels)
    )

  # barplot
  Colors <- c("#000080","#0029FF","#00D5FF","#7DFF7A","#FFE600","#FF4700","#800000","#808080")
  p1 <- ggplot(dat2,aes(Samples,Abundance,fill=Species))+
    geom_bar(stat = "identity",color="black")+
    scale_y_continuous(expand = expand_scale(mult = c(0,0.1)))+
    scale_fill_manual(values = Colors)+
    xlab("") + ylab("Relative abundance")+
    labs(title = FullName)+
    theme_bw()+
    theme(
      axis.title = element_text(size = 14,color = "black"),
      axis.text = element_text(size = 13,color = "black"),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text(size = 14,color = "black"),
      legend.text = element_text(size = 13,color = "black",face = "italic"),
      panel.grid = element_blank(),
      plot.title = element_text(size = 14,color = "black",hjust = 0.5)
    )

  p2 <- ggplot(a,aes(SampleID,y,fill=Group))+
    geom_tile()+
    xlab("Samples")+ylab("")+
    scale_y_discrete(expand = c(0,0))+
    theme_bw()+
    theme(
      axis.title = element_text(size = 14,color = "black"),
      axis.text = element_text(size = 13,color = "black"),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text(size = 14,color = "black"),
      legend.text = element_text(size = 13,color = "black"),
      panel.grid = element_blank()
    )

  ggarrange(p1,p2,nrow = 2,ncol = 1,align = "v",heights = c(5,1),legend = "right")
}
