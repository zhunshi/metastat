% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_heatmap.R
\name{pheatmap.scc}
\alias{pheatmap.scc}
\title{customize pheatmap}
\usage{
pheatmap.scc(
  scc,
  maxRange = TRUE,
  ReAdjPvalLowerTri = FALSE,
  trans = FALSE,
  color = NULL,
  Psignificance = FALSE,
  ...
)
}
\arguments{
\item{scc}{output of scc_self, scc_between, pcc_self, pcc_between}

\item{maxRange}{Whether to keep the maximum and minimum values the same, default TRUE}

\item{ReAdjPvalLowerTri}{if re-adjust pvalue (lower tri)}

\item{trans}{Whether to transpose the matrix}

\item{color}{fill colors (eg. colorRampPalette(colors = c("navy", "white", "firebrick3"))[50])}

\item{Psignificance}{If display sigificance for uncorrected P values}

\item{...}{other parameters passed to pheatmap}
}
\value{

}
\description{
customize pheatmap
}
