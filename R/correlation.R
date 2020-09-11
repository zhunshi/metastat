#' Spearman's correlation
#'
#' This function perform spearman's correlation between two dataframe
#'
#' @param dat1  input data1
#' @param dat2  input data2
#'
#' @return dataframe
#'
#' @export
#' @example

scc_between <- function(dat1,dat2){
  inte <- intersect(rownames(dat1),rownames(dat2))
  dat1 <- dat1[inte,,drop=F]
  dat2 <- dat2[inte,,drop=F]
  out <- as.data.frame(matrix(NA,ncol(dat1),ncol(dat2)*3))
  rownames(out) <- colnames(dat1)
  colnames(out) <- paste(rep(colnames(dat2),rep(3,ncol(dat2))),rep(c("rho","pvalue","FDR"),ncol(dat2)))
  for(i in 1:ncol(dat2)){
    #print(i)
    x <- as.numeric(dat2[,i])
    for(j in 1:ncol(dat1)){
      y <- as.numeric(dat1[,j])
      res <- cor.test(x,y,method = "spearman")
      out[j,i*3-2] <- res$estimate
      out[j,i*3-1] <- res$p.value
    }
    out[,i*3] <- p.adjust(out[,i*3-1],method = "BH")
  }
  return(out)
}

#' Partial Spearman's correlation
#'
#' This function perform partial spearman's correlation between two dataframe
#'
#' @param dat1     input data1
#' @param dat2     input data2
#' @param dat_adj  input adjustment data
#'
#' @return dataframe
#'
#' @export
#' @example
pcc_between <- function(dat1,dat2,dat_adj){
  inte <- intersect(rownames(dat1),rownames(dat2))
  dat1 <- dat1[inte,,drop=F]
  dat2 <- dat2[inte,,drop=F]
  dat_adj <- dat_adj[inte,,drop=F]
  out <- as.data.frame(matrix(NA,ncol(dat1),ncol(dat2)*3))
  rownames(out) <- colnames(dat1)
  colnames(out) <- paste(rep(colnames(dat2),rep(3,ncol(dat2))),rep(c("rho","pvalue","FDR"),ncol(dat2)))
  for(i in 1:ncol(dat2)){
    #print(i)
    x <- dat2[,i,drop=F]
    for(j in 1:ncol(dat1)){
      y <- dat1[,j,drop=F]
      a <- cbind(x,y,dat_adj)
      a <- na.omit(a)
      res <- pcor.test(a[,1,drop=F],a[,2,drop=F],dat_adj[rownames(a),,drop=F],method = "spearman")
      out[j,i*3-2] <- res$estimate
      out[j,i*3-1] <- res$p.value
    }
    out[,i*3] <- p.adjust(out[,i*3-1],method = "BH")
  }
  return(out)
}

#' Spearman's correlation
#'
#' This function perform spearman's correlation in one dataframe
#'
#' @param dat   input data
#'
#' @return dataframe
#'
#' @export
#' @example
scc_self <- function(dat){
  out <- as.data.frame(matrix(NA,ncol(dat),ncol(dat)*3))
  rownames(out) <- colnames(dat)
  colnames(out) <- paste(rep(colnames(dat),rep(3,ncol(dat))),rep(c("rho","pvalue","FDR"),ncol(dat)))
  for(i in 1:ncol(dat)){
    #print(i)
    x <- as.numeric(dat[,i])
    for(j in 1:ncol(dat)){
      y <- as.numeric(dat[,j])
      res <- cor.test(x,y,method = "spearman")
      out[j,i*3-2] <- res$estimate
      out[j,i*3-1] <- res$p.value
    }
    out[,i*3] <- p.adjust(out[,i*3-1],method = "BH")
  }
  return(out)
}

#' Partial Spearman's correlation
#'
#' This function perform partial spearman's correlation in one dataframe
#'
#' @param dat       input data
#' @param data_adj  input adjustment
#'
#' @return dataframe
#'
#' @export
#' @example
pcc_self <- function(dat,dat_adj){
  dat_adj <- dat_adj[rownames(dat),,drop=F]
  out <- as.data.frame(matrix(NA,ncol(dat),ncol(dat)*3))
  rownames(out) <- colnames(dat)
  colnames(out) <- paste(rep(colnames(dat),rep(3,ncol(dat))),rep(c("rho","pvalue","FDR"),ncol(dat)))
  for(i in 1:ncol(dat)){
    x <- dat[,i,drop=F]
    for(j in 1:ncol(dat)){
      if(i==j){
        out[j,i*3-2] <- 1
        out[j,i*3-1] <- 0
        next
      }
      y <- dat[,j,drop=F]
      a <- cbind(x,y,dat_adj)
      a <- na.omit(a)
      res <- pcor.test(a[,1,drop=F],a[,2,drop=F],dat_adj[rownames(a),,drop=F],method = "spearman")
      out[j,i*3-2] <- res$estimate
      out[j,i*3-1] <- res$p.value
    }
    out[,i*3] <- p.adjust(out[,i*3-1],method = "BH")
  }
  return(out)
}


#' Customized mediation analysis
#'
#' Customized mediation analysis.
#' fit1: m ~ x + adj
#' fit2: y ~ x + m + adj
#'
#' @param data   input data
#' @param x      x
#' @param y      y
#' @param m      mediation variable
#' @param adj    adjust variable
#'
#' @return dataframe
#'
#' @export
mediation_analysis <- function(data,x,y,m,adj){
  #print(c(x,y,m))
  data <- na.omit(data)
  f1 <- as.formula(paste0(m,"~",x,"+",paste(adj,collapse = "+")))
  f2 <- as.formula(paste0(y,"~",x,"+",m,"+",paste(adj,collapse = "+")))
  fit1 <- glm(f1,data=data)
  fit2 <- glm(f2,data=data)
  fit3 <- mediation::mediate(fit1,fit2,treat = x,mediator = m)
  med.out <- fit3
  lines <- summary(med.out) %>%
    capture.output() %>%
    discard(`==`,"")
  i <- which(grepl("^ ", lines))
  res <- lines[i:(i+4)] %>%
    sub("^       ", "med.out", .) %>%
    gsub("95% CI Lower", "CI_95%_Lower", .) %>%
    gsub("95% CI Upper", "CI_95%_Upper", .) %>%
    sub(" \\(", "_(", .) %>%
    sub("p-", "p_", .) %>%
    sub("\\*+","",.) %>%
    sub("<","",.) %>%
    sub("\\. *$","",.) %>%
    sub(" ", "_", ., fixed=TRUE)  %>%
    textConnection() %>%
    read.table(header=TRUE) %>%
    setNames(sub("_$", "", colnames(.))) %>%
    dplyr::mutate(med.out=sub("\\.|_$", "", med.out),
                  med.out=gsub("_", " ", med.out))
  out <- data.frame(
    X = x,Y = y, mediator = m,
    ACME.estimate = res$Estimate[1],
    ACME.p = res$p_value[1],
    ADE.estimate = res$Estimate[2],
    ADE.p = res$p_value[2],
    Total.estimate = res$Estimate[3],
    Total.p = res$p_value[3],
    Prop.mediate = res$Estimate[4],
    Prop.p = res$p_value[4]
  )
  return(out)
}
