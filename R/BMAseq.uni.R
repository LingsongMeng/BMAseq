##' Univariate analysis using BMAseq approach
##'
##' Univariate analysis on RNA-seq count data using BMAseq approach
##' @title BMAseq.uni
##' @param dat.expr.counts RNA-seq count data matrix (rows for genes and columns for subjects).
##' @param dat.pheno Phenotypic data matrix (rows for subjects and columns for variables).
##' @param var.pool Variables of interest, a vector.
##' @param cut.BF Bayes factor criterion, default value is 1.
##' @param cut.FDR False discovery rate criterion for identifying DE genes, default value is 0.25.
##' @param laplace Use Laplace approximation to calculate Bayes factor or not, default is FALSE.}
##' @return A list consisting of
##' \item{dat.expr.logcpm}{Normalized RNA-seq data matrix (rows for genes and columns for subjects).}
##' \item{weights}{Estimated voom weights.}
##' \item{eFDR}{Estimated false discovery rate matrix (rows for genes and columns for variables of interest).}
##' \item{nDEG}{The number of DE genes associated with each variable of interest.}
##' \item{DEG}{DE genes associated with each variable of interest.}
##' @export
##' @author Lingsong Meng


BMAseq.uni <- function(dat.expr.counts, dat.pheno, var.pool, cut.BF = 1, cut.FDR = 0.25, laplace = FALSE) {


    if (sum(colnames(dat.expr.counts) != rownames(dat.pheno)) > 0)
        stop("Two datasets dat.expr.counts and dat.pheno must be matched by subjects.")
    if (is.null(var.pool))
        stop("Must provide at least one variable of interest.")
    if (sum(var.pool %in% colnames(dat.pheno)) != length(var.pool))
        stop("Variables of interest must be in datasets dat.pheno.")

    n.var <- length(var.pool)
    n.gene <- nrow(dat.expr.counts)

    # estimate voom weights
    y0 <- voom(dat.expr.counts)
    input.dat.expr <- y0$E
    input.weights <- y0$weights

    eFDR <- eFDR.notOrder <- ind.order.gene <- rank.gene <- matrix(NA, nrow = n.gene, ncol = n.var)
    colnames(eFDR) <- var.pool
    colnames(eFDR.notOrder) <- var.pool
    rownames(eFDR.notOrder) <- rownames(dat.expr.counts)

    for (i.var in 1:n.var) {

        print(i.var)
        # build model space
        result.modelspace <- Modelspace(dat.pheno, var.pool[i.var], max.nvar = 1, interaction = NULL)
        input.model.space <- result.modelspace$model.space
        input.dat.pheno <- dat.pheno

        if (!laplace) {

            # calculate Bayes factor
            out.bf <- NULL
            for (i in 1:length(input.model.space)) {
                print(input.model.space[i])
                design <- model.matrix(formula(input.model.space[i]), data = input.dat.pheno)
                colnames(design)[1] <- "Int"
                x <- NULL
                for (k in 1:nrow(input.dat.expr)) {
                  x <- c(x, Bayesfactor(design, input.dat.expr[k, ], input.weights[k, ]))
                }
                out.bf <- cbind(out.bf, x)
            }

            # obtain prior model probability
            bf.max <- t(apply(out.bf, 1, function(x) ifelse(x == max(x) & x > cut.BF, 1, 0)))
            prior.modelprob <- apply(bf.max, 2, sum)/nrow(input.dat.expr)
            pm.est <- prior.modelprob
            n.model <- length(prior.modelprob)
            pm.prior0 <- pm.prior1 <- rep(0.5/(n.model - 1), n.model - 1)
            alpha <- 0.2
            for (j in 1:30) {
                if (is.vector(out.bf[, -1])) {
                  pm.prior1 <- pm.est[-1]/mean(out.bf[, -1]/(1 - pm.prior0 + out.bf[, -1] * pm.prior0))
                } else pm.prior1 <- pm.est[-1]/apply(t(apply(out.bf[, -1], 1, function(x) x/(1 - sum(pm.prior0) + sum(x * pm.prior0)))),
                  2, mean)
                pm.prior0 <- pm.prior0 + (pm.prior1 - pm.prior0) * alpha
            }
            prior.modelprob <- c(1 - sum(pm.prior0), pm.prior0)

            # calculate posterior model probability
            post.modelprob <- t(apply(out.bf, 1, function(x) x * prior.modelprob/sum(x * prior.modelprob)))

        } else {

            # calculate Bayes factor
            out.logbf <- NULL
            for (i in 1:length(input.model.space)) {
                print(input.model.space[i])
                design <- model.matrix(formula(input.model.space[i]), data = input.dat.pheno)
                colnames(design)[1] <- "Int"
                x <- NULL
                for (k in 1:nrow(input.dat.expr)) {
                  x <- c(x, Bayesfactor.Laplace(design, input.dat.expr[k, ], input.weights[k, ]))
                }
                out.logbf <- cbind(out.logbf, x)
            }

            # obtain prior model probability
            logbf.max <- t(apply(out.logbf, 1, function(x) ifelse(x == max(x) & x > log(cut.BF), 1, 0)))
            prior.modelprob <- apply(logbf.max, 2, sum)/nrow(input.dat.expr)
            pm.est <- prior.modelprob
            n.model <- length(prior.modelprob)
            pm.prior0 <- pm.prior1 <- rep(0.5/(n.model - 1), n.model - 1)
            alpha <- 0.2
            for (j in 1:30) {
                if (is.vector(out.logbf[, -1])) {
                  pm.prior1 <- pm.est[-1]/mean(exp(-log((1 - pm.prior0)/exp(out.logbf[, -1]) + pm.prior0)))
                } else {
                  denom <- t(apply(out.logbf[, -1], 1, function(y) exp(y - max(y) - log((1 - sum(pm.prior0))/max(exp(y)) +
                    sum(exp(y + log(pm.prior0) - max(y)))))))
                  pm.prior1 <- pm.est[-1]/apply(denom, 2, mean)
                }
                pm.prior0 <- pm.prior0 + (pm.prior1 - pm.prior0) * alpha
            }
            prior.modelprob <- c(1 - sum(pm.prior0), pm.prior0)

            # calculate posterior model probability
            post.modelprob <- t(apply(out.logbf, 1, function(y) exp(y + log(prior.modelprob) - max(y) - log(sum(exp(y + log(prior.modelprob) -
                max(y)))))))
        }

        rownames(post.modelprob) <- rownames(dat.expr.counts)

        # calculate estimated FDR
        p.excl <- 1 - post.modelprob[, 2]
        ntop <- seq(1, nrow(input.dat.expr), by = 1)
        ind.order.gene[, i.var] <- order(p.excl)
        rank.gene[, i.var] <- rank(p.excl)
        eFDR[, i.var] <- apply(matrix(1 - p.excl), 2, function(x) cumsum(1 - x[order(1 - x)])[ntop]/ntop)
        eFDR.notOrder[, i.var] <- eFDR[rank.gene[, i.var], i.var]
    }

    nDEG <- apply(eFDR, 2, function(x) sum(x < cut.FDR))
    nDEG.o <- order(nDEG, decreasing = T)
    nDEG.dec <- nDEG[nDEG.o]
    DEG <- list()
    for (i in 1:n.var) {
        index.DEG <- ind.order.gene[, i][eFDR[, i] < cut.FDR]
        name.DEG <- rownames(dat.expr.counts)[index.DEG]
        DEG[[i]] <- cbind(index.DEG, name.DEG)
    }
    names(DEG) <- var.pool


    return(list(dat.expr.logcpm = input.dat.expr, weights = input.weights, eFDR = eFDR.notOrder, nDEG = nDEG.dec, DEG = DEG))
}

