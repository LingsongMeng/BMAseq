% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BMAseq.multi.DEG.R
\name{BMAseq.multi.DEG}
\alias{BMAseq.multi.DEG}
\title{BMAseq.multi.DEG}
\usage{
BMAseq.multi.DEG(postprob.output, ind.incl.add = NULL, cut.FDR = 0.05)
}
\arguments{
\item{postprob.output}{The output from BMAseq.multi.postprob function, a list.}

\item{ind.incl.add}{Indices of additional included model.}

\item{cut.FDR}{False discovery rate criterion for identifying DE genes, default value is 0.05.}
}
\value{
A list consisting of
\item{post.incl.modelprob}{Posterior inclusive model probability for each gene associated with main effect of each variables if interaction is NULL; with main, interaction, main or interaction effects if interaction not NULL..}
\item{post.incl.modelprob.JointMain}{Posterior inclusive model probability for each gene associated with all possible joint effect of variables, will be output if interaction is NULL.}
\item{post.incl.modelprob.Add}{Posterior inclusive model probability for each gene associated with additional inclusive models, will be output if ind.incl.add is not Null.}
\item{index.incl}{Indices of the inclusive models.}
\item{eFDR}{Estimated FDR for each gene associated with main effect of each variables if interaction is NULL; with main, interaction, main or interaction effects if interaction is not NULL.}
\item{eFDR.JointMain}{Estimated FDR for each gene associated with all possible joint effect of variables, will be output if interaction is NULL.}
\item{eFDR.Add}{Estimated FDR for each gene associated with additional inclusive models, will be output if ind.incl.add is not Null.}
\item{summary.nDEG}{A summary table of the number of identified DE gene associated with main effect of each variable if interaction is NULL; with main, interaction, main or interaction effects if interaction is not NULL.}
\item{summary.nDEG.JointMain}{A summary table of the number of identified DE gene associated with all possible joint effcts of variables, will be output if interaction is NULL.}
\item{summary.nDEG.Add}{A summary table of the number of identified DE gene associated with additional inclusive models, will be output if ind.incl.add is not Null.}
\item{DEG.bestmodel}{DE genes associated with main effect of each variable if interaction is NULL; with main, interaction, main or interaction effects of each variable if interaction is not NULL; and the best model used to identify each DE gene.}
}
\description{
DE genes identification with obtained posterior model probability for each gene in multivariate analysis using BMAseq approach
}
\details{
DE genes identification with obtained posterior model probability for each gene in multivariate analysis on RNA-seq count data using BMAseq approach
}
\author{
Lingsong Meng
}
