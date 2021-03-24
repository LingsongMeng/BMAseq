##' Calculation of posterior model probability for each gene in multivariate analysis using BMAseq approach
##'
##' Calculation of posterior model probability for each gene in multivariate analysis on RNA-seq count data using BMAseq approach
##' @title BMAseq.multi.postprob
##' @param dat.expr.counts RNA-seq count data matrix (rows for genes and columns for subjects).
##' @param dat.pheno Phenotypic data matrix (rows for subjects and columns for variables).
##' @param var.pool Variables of interest, a vector.
##' @param max.nvar The maximum number of variables in a model.
##' @param interaction Variables considered in interaction terms, a vector, default is NULL.
##' @param cut.BF Bayes factor criterion, default value is 1.
##' @return A list consisting of
##' \item{dat.expr.logcpm}{Normalized RNA-seq data matrix (rows for genes and columns for subjects).}
##' \item{weights}{Estimated voom weights.}
##' \item{dat.pheno.new}{Phenotypic data matrix including new interaction variables, will be same as input dat.pheno if interaction=NULL.}
##' \item{model.space}{Model space including all possible models.}
##' \item{post.modelprob}{Posterior model probability for each gene.}
##' @export
##' @author Lingsong Meng


BMAseq.multi.postprob.Laplace <- function(dat.expr.counts, dat.pheno, var.pool, max.nvar, interaction = NULL, cut.BF = 1) {
    
    
    if (sum(colnames(dat.expr.counts) != rownames(dat.pheno)) > 0) 
        stop("Two datasets dat.expr.counts and dat.pheno must be matched by subjects.")
    if (is.null(var.pool) || length(var.pool) < 2) 
        stop("Must provide at least two variables of interest.")
    if (sum(var.pool %in% colnames(dat.pheno)) != length(var.pool)) 
        stop("Variables of interest must be in datasets dat.pheno.")
    if (is.null(max.nvar)) 
        stop("Must provide max.nvar.")
    
    
    n.var <- length(var.pool)
    n.gene <- nrow(dat.expr.counts)
    
    
    # estimate voom weights
    y0 <- voom(dat.expr.counts)
    input.dat.expr <- y0$E
    input.weights <- y0$weights
    
    
    # build model space
    result.modelspace <- Modelspace(dat.pheno, var.pool, max.nvar, interaction = interaction)
    input.model.space <- result.modelspace$model.space
    input.dat.pheno <- result.modelspace$dat.pheno.new
    
    
    # calculate Bayes factor
    out.logbf <- vapply(1:length(input.model.space), function(i) {
        print(paste0(i, ". ", input.model.space[i]))
        design <- model.matrix(formula(input.model.space[i]), data = input.dat.pheno)
        colnames(design)[1] <- "Int"
        vapply(1:n.gene, function(k) Bayesfactor.Laplace(design, input.dat.expr[k, ], input.weights[k, ]), FUN.VALUE = as.double(1))
    }, FUN.VALUE = as.double(1:n.gene))
    
    
    # obtain prior model probability
    logbf.max <- t(apply(out.logbf, 1, function(x) ifelse(x == max(x) & x > log(cut.BF), 1, 0)))
    prior.modelprob <- apply(logbf.max, 2, sum)/n.gene
    pm.est <- prior.modelprob
    n.model <- length(prior.modelprob)
    pm.prior0 <- pm.prior1 <- rep(0.5/(n.model - 1), n.model - 1)
    alpha <- 0.2
    for (j in 1:30) {
        denom <- t(apply(out.logbf[, -1], 1, function(y) exp(y - max(y) - log((1 - sum(pm.prior0))/max(exp(y)) + sum(exp(y + log(pm.prior0) - max(y)))))))
        pm.prior1 <- pm.est[-1]/apply(denom, 2, mean)
        pm.prior0 <- pm.prior0 + (pm.prior1 - pm.prior0) * alpha
    }
    prior.modelprob <- c(1 - sum(pm.prior0), pm.prior0)
    
    
    # calculate posterior model probability
    post.modelprob <- t(apply(out.logbf, 1, function(y) exp(y + log(prior.modelprob) - max(y) - log(sum(exp(y + log(prior.modelprob) - max(y)))))))
    rownames(post.modelprob) <- rownames(dat.expr.counts)
    colnames(post.modelprob) <- input.model.space
    
    
    
    return(list(dat.expr.logcpm = input.dat.expr, weights = input.weights, model.space = input.model.space, dat.pheno.new = input.dat.pheno, post.modelprob = post.modelprob))
}



