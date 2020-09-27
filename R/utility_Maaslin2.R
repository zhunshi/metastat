#' read maaslin2 results and re-adjusted(BH)
#'
#' read maaslin2 results and re-adjusted(BH) in a specified directory
#'
#' @param path        specify a directory of the maaslin results
#' @return dataframe
#'
#' @export
ReAdjustBHMaaslin2 <- function(Path){
  dirs <- dir(Path)
  files <- sapply(dirs,function(x) paste0(Path,"/",x,"/all_results.tsv"))
  datlst <- lapply(files,function(x) read.table(x,header = T,check.names = F,sep = "\t"))
  datlst <- lapply(
    datlst,function(dat){dat %>% group_by(metadata) %>% mutate(BH=p.adjust(pval,method = "BH"))})
  names(datlst) <- dirs
  return(datlst)
}
