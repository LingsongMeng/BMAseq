##' Analysis using single model approach
##'
##' Analysis on RNA-seq count data using single model approach
##' @title BMAseq.multi
##' @param dat.expr.counts RNA-seq count data matrix (rows for genes and columns for subjects).
##' @param dat.pheno Phenotypic data matrix (rows for subjects and columns for variables).
##' @param var.pool Variables of interest, a vector.
##' @param max.nvar The maximum number of variables in a model.
##' @param interaction Interactions.
##' @return A list consisting of
##' \item{dat.expr.logcpm}{Normalized RNA-seq data matrix (rows for genes and columns for subjects).}
##' \item{weights}{Estimated voom weights.}
##' \item{summary.lmFit}{Output results of single linear models}
##' @export
##' @author Lingsong Meng


lmFit.single <- function(dat.expr.counts, dat.pheno, var.pool, max.nvar, interaction = NULL) {


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
    result.modelspace <- Modelspace(dat.pheno, var.pool, max.nvar, interaction)
    if (!is.null(interaction)) {
        input.model.space <- result.modelspace$model.space
        input.dat.pheno <- result.modelspace$dat.pheno.new
    } else {
        input.model.space <- result.modelspace$model.space
        input.dat.pheno <- dat.pheno
    }

    out.lmFit <- list()
    for (i in 1:length(input.model.space)) {
        print(paste0(i, ". ", input.model.space[i]))
        design <- model.matrix(formula(input.model.space[i]), data = input.dat.pheno)
        out.lmFit[[i]] <- my.lmFit(input.dat.expr, design, input.weights)
    }
    names(out.lmFit) <- input.model.space


    return(list(dat.expr.logcpm = input.dat.expr, weights = input.weights, summary.lmFit = out.lmFit))
}

