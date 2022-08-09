
#' KW test
#'
#' This function perform Kruskal-Walis test among multiple gorups
#'
#' @param dat , input data
#' @param grp , input groups with one column
#'
#' @return dataframe
#'
#' @export

KW_phenotypes <- function(dat,grp){
  inte <- intersect(rownames(dat),rownames(grp))
  dat <- dat[inte,]
  grp <- grp[inte,,drop=F]

  #variables
  num <- ncol(dat)
  fr <- grp[, 1]
  group <- levels(fr)
  #output
  len <- length(group)
  num2 <- len*len/2+3.5*len + 2
  out <- matrix(NA, num, num2)
  out <- as.data.frame(out)
  #kurskal test
  ktp <- apply(dat, 2, function(x,fr){kruskal.test(x ~ fr)$p.value},fr=fr)
  #post hoc dunn test
  #library(PMCMRplus)
  for (i in 1:num) {
    #    print(i)
    index1 <- is.na(fr)
    index2 <- is.na(dat[,i])
    index <- index1 | index2
    fr1 <- fr[!index]
    rk  <- rank(dat[,i][!index])
    dat1 <- dat[,i][!index]
    dat1 <- as.numeric(dat1)
    res <- c(ktp[i],
             tapply(dat1, fr1, median),
             tapply(dat1,fr1,mean),
             tapply(!is.na(dat[,i][!index1]), fr[!index1],sum))
    dtp <- PMCMRplus::kwAllPairsDunnTest(dat1, fr1, p.adjust.method = "BH")$p.value
    dtp <- cbind(dtp, rep(NA, len - 1))
    dtp <- rbind(rep(NA, len), dtp)
    dtp[upper.tri(dtp)] <- t(dtp)[upper.tri(dtp)]
    rownames(dtp)[1] <- colnames(dtp)[1]
    colnames(dtp)[len] <- rownames(dtp)[len]

    mean_rank <- tapply(rank(dat1),fr1,mean)
    res <- c(res,dtp[lower.tri(dtp)], mean_rank)
    #
    conclude <- rep(0,2*len-1)
    or <- order(mean_rank)
    conclude[2*(1:len)-1] <- group[or]
    op <- rep(1,len-1)
    for(j in 1:(len-1)){op[j] <- dtp[or[j],or[j+1]]}
    symbol <- rep("=",len-1)
    symbol[!is.na(op) & op <= 0.05] <- "<"
    symbol[is.na(op) | op == 1] <- "<=>"
    for(x in 1:(len-1)){
      if(symbol[x]=="<"){
        p_tmp <- c()
        for(y in 1:x){
          for(z in (x+1):len){
            p_tmp <- c(p_tmp,dtp[or[y],or[z]])
          }
        }
        if(any(p_tmp>0.05)){symbol[x] <- "="}
      }
    }

    conclude[(1:(len - 1)) * 2] <- symbol
    tmp <- paste(conclude, collapse = " ")
    #res <- c(res, paste(conclude, collapse = " "))
    if(length(res)==(ncol(out)-1)){
      out[i, 1:length(res)] <- res
      out[i, ncol(out)] <- tmp
    }else{print (res)}
  }
  rownames(out) <- colnames(dat)
  cn <- c("kw.p", paste("median", group[1:len], sep = "_"))
  cn <- c(cn,paste("mean", group[1:len], sep = "_"))
  cn <- c(cn, paste("or", group[1:len], sep = "_"))
  cn <- c(cn, paste0("p", "_", group[row(dtp)[lower.tri(dtp)]], "_", group[col(dtp)[lower.tri(dtp)]]))
  cn <- c(cn, paste("mean_rank",group[1:len], sep = "_"))
  cn <- c(cn, "nearby")
  colnames(out) <- cn
  out$FDR <- p.adjust(out$kw.p,method = "BH")
  out <- out[,c(1,ncol(out),2:(ncol(out)-1))]
  return(out)
}


#' customized wilcoxon test
#'
#' This function perform wilcox test between gorups
#'
#' @param dat , input data
#' @param grp , input groups with one column
#'
#' @return dataframe
#'
#' @export

wilcox.customized <- function(dat,grp,type="phenotype"){
  grp <- na.omit(grp)
  inte <- intersect(rownames(dat),rownames(grp))
  dat <- dat[inte,]
  grp <- grp[inte,,drop=F]
  if(min(dat,na.rm = T)>=0 & max(dat,na.rm = T)<1 & all(apply(dat,2,function(x) class(x)=="numeric"))){
    type <- "profile"
  }
  cat("The input data type is: ",type,"\n")
  grp[,1] <- as.factor(grp[,1])
  grp.level <- levels(grp[,1])
  grp.level.n <- length(grp.level)
  if(grp.level.n!=2){
    stop("Number of Group levels must be 2")
  }
  out.cn <- paste(rep(c("median","mean","SD","mean_rank","occ_rate","n"),each=2),
                  rep(grp.level,6),sep = ".")
  out.cn <- c("pvalue","FDR",out.cn,"Enrichment.FDR","Effectsize")
  out <- matrix(NA,ncol(dat),length(out.cn))
  out <- as.data.frame(out)
  colnames(out) <- out.cn
  rownames(out) <- colnames(dat)

  out$pvalue <- apply(dat,2,function(x) wilcox.test(x~grp[,1])$p.value)
  out$FDR <- p.adjust(out$pvalue,method = "BH")

  for(i in 1:nrow(out)){
    out[i,3:4] <- tapply(dat[,i],grp[,1],function(y) median(y,na.rm = T))
    out[i,5:6] <- tapply(dat[,i],grp[,1],function(y) mean(y,na.rm = T))
    out[i,7:8] <- tapply(dat[,i],grp[,1],function(y) sd(y,na.rm = T))
    tmp <- dat[,i]
    tmp2 <- tmp[!is.na(tmp)]
    out[i,9:10] <- tapply(rank(tmp2),grp[,1][!is.na(tmp)],function(y) mean(y))
    out[i,11:12] <- tapply(dat[,i],grp[,1],
                           function(y) ifelse(type=="phenotype",sum(!is.na(y))/length(y),sum(y>0)/length(y)))
    out[i,13:14] <- tapply(dat[,i],grp[,1],
                           function(y) ifelse(type=="phenotype",sum(!is.na(y)),sum(y>0)))
    out[i,16] <- rstatix::wilcox_effsize(dat[,i]~grp[,1])$effsize
  }
  out$Enrichment.FDR <- ifelse(out[,9]>out[,10],grp.level[1],grp.level[2])
  out$Enrichment.pval <- out$Enrichment.FDR
  out$Enrichment.FDR[out$FDR>0.05 | is.na(out$FDR)] <- "NONE"
  out$Enrichment.pval[out$pvalue>0.05 | is.na(out$pvalue)] <- "NONE"
  return(out)
}


#' customized paired wilcoxon-rank sum test
#'
#' The function to perform paired wilcoxon-rank sum test between gorups
#'
#' @param dat     input data
#' @param grp     groups with two column, "paired" for patient ID and "group" for treatment groups
#' @param paired  which column used as patiend ID
#' @param group   which column used as treatment groups
#'
#' @return dataframe
#'
#' @export
wilcox.Paired <- function(dat,grp,paired,group){
  grp <- na.omit(grp)
  inte <- intersect(rownames(dat),rownames(grp))
  grp <- grp[inte,,drop=F]
  dat <- dat[inte,,drop=F]

  out.cn <- paste(rep(c("median","mean","SD","mean_rank","N"),rep(2,5)),
                  rep(levels(factor(grp[,group])),5),sep = ".")
  out.cn <- c("pvalue","FDR",out.cn,"Enrichment.Direction")
  out <- matrix(NA,ncol(dat),length(out.cn))
  out <- as.data.frame(out)
  colnames(out) <- out.cn
  rownames(out) <- colnames(dat)

  for(i in 1:nrow(out)){
    tmp <- cbind(dat[,i,drop=F],grp)
    tmp2 <- tmp %>% spread(get(group),get(colnames(dat)[i]))
    tmp2 <- na.omit(tmp2)
    out[i,1] <- wilcox.test(tmp2[,2],tmp2[,3],paired = T)$p.value
    out[i,3:4] <- c(median(tmp2[,2],na.rm = T),median(tmp2[,3],na.rm=T))
    out[i,5:6] <- c(mean(tmp2[,2],na.rm = T),mean(tmp2[,3],na.rm=T))
    out[i,7:8] <- c(sd(tmp2[,2],na.rm = T),sd(tmp2[,3],na.rm=T))
    rk <- rank(c(tmp2[,2],tmp2[,3]))
    out[i,9:10] <- c(mean(rk[1:nrow(tmp2)]),mean(rk[(nrow(tmp2)+1):(nrow(tmp2)*2)]))
    out[i,11:12] <- c(sum(tmp2[,2]>0),sum(tmp2[,3]>0))
  }
  out$FDR <- p.adjust(out$pvalue,method = "BH")

  l <- levels(factor(grp[,group]))
  out$Enrichment.Direction <- ifelse(out[,9]>out[,10],l[1],l[2])
  out$Enrichment.Direction[out$FDR>=0.05 | is.na(out$FDR)] <- "NONE"
  return(out)
}
