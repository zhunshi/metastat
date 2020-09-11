#' Combine Dataframes
#'
#' Inplementation of combine tables function in plyr
#'
#' @param datalst        input list of dataframes
#' @param by             which column used to combine
#' @param method         combine method in dplyr packages used to combine dataframes. Default is full_join. Supported methods including full_join,left_join,right_join,inner_jion
#' @param use_list_name  if use list name replace original names.
#' @return dataframe
#'
#' @export
CombineTidyTables <- function(datlst,by,method= dplyr::full_join,use_list_name=TRUE){
  dat1 <- datlst[[1]]
  for(i in 2:length(datlst)){
    dat1 <- dat1 %>%
      method(datlst[[i]],by=by)
  }
  if(use_list_name){
    colnames(dat1) <- c(by,names(datlst))
  }
  return(dat1)
}
