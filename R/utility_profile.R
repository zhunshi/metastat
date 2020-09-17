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


