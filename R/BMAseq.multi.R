##' Multivariate analysis using BMAseq approach
##'
##' Multivariate analysis on RNA-seq count data using BMAseq approach
##' @title BMAseq.multi
##' @param dat.expr.counts RNA-seq count data matrix (rows for genes and columns for subjects).
##' @param dat.pheno Phenotypic data matrix (rows for subjects and columns for variables).
##' @param var.pool Variables of interest, a vector.
##' @param max.nvar The maximum number of variables in a model.
##' @param interaction Specific interaction terms with input format 'A&B', default is Null.
##' @param ind.incl.add Indices of additional included model.
##' @param cut.BF Bayes factor criterion, default value is 1.
##' @param cut.FDR False discovery rate criterion for identifying DE genes, default value is 0.05.
##' @param laplace Use Laplace approximation to calculate Bayes factor or not, default is FALSE.}
##' @return A list consisting of
##' \item{dat.expr.logcpm}{Normalized RNA-seq data matrix (rows for genes and columns for subjects).}
##' \item{weights}{Estimated voom weights.}
##' \item{dat.pheno.new}{Phenotypic data matrix including new interaction variables, will be same as input dat.pheno if interaction is NULL.}
##' \item{model.space}{Model space including all possible models.}
##' \item{post.modelprob}{Posterior model probability for each gene.}
##' \item{post.incl.modelprob}{Posterior inclusive model probability for each gene associated with main effect of each variables if interaction is NULL; with main, interaction, main or interaction effects if interaction not NULL.}
##' \item{post.incl.modelprob.JointMain}{Posterior inclusive model probability for each gene associated with joint main effects of variables, will be output if interaction is NULL.}
##' \item{post.incl.modelprob.Add}{Posterior inclusive model probability for each gene associated with additional inclusive models, will be output if ind.incl.add is not Null.}
##' \item{index.incl}{Indices of the inclusive models.}
##' \item{eFDR}{Estimated FDR for each gene associated with main effect of each variables if interaction is NULL; with main, interaction, main or interaction effects if interaction is not NULL.}
##' \item{eFDR.JointMain}{Estimated FDR for each gene associated with all possible joint effect of variables, will be output if interaction is NULL.}
##' \item{eFDR.Add}{Estimated FDR for each gene associated with additional inclusive models, will be output if ind.incl.add is not Null.}
##' \item{summary.nDEG}{A summary table of the number of identified DE gene associated with main effect of each variable if interaction is NULL; with main, interaction, main or interaction effects if interaction is not NULL.}
##' \item{summary.nDEG.JointMain}{A summary table of the number of identified DE gene associated with joint main effcts of variables, will be output if interaction is NULL.}
##' \item{summary.nDEG.Add}{A summary table of the number of identified DE gene associated with additional inclusive models, will be output if ind.incl.add is not Null.}
##' \item{DEG.bestmodel}{DE genes associated with main effect of each variable if interaction is NULL; with main, interaction, main or interaction effects of each variable if interaction is not NULL; and the best model used to identify each DE gene.}
##' @export
##' @author Lingsong Meng


BMAseq.multi <- function(dat.expr.counts, dat.pheno, var.pool, max.nvar, interaction = NULL, ind.incl.add = NULL, cut.BF = 1,
    cut.FDR = 0.05, laplace = FALSE) {

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

    if (!laplace) {

        # calculate Bayes factor
        out.bf <- NULL
        for (i in 1:length(input.model.space)) {
            print(paste0(i, ". ", input.model.space[i]))
            design <- model.matrix(formula(input.model.space[i]), data = input.dat.pheno)
            colnames(design)[1] <- "Int"
            x <- NULL
            for (k in 1:n.gene) {
                x <- c(x, Bayesfactor(design, input.dat.expr[k, ], input.weights[k, ]))
            }
            out.bf <- cbind(out.bf, x)
        }

        # obtain prior model probability
        bf.max <- t(apply(out.bf, 1, function(x) ifelse(x == max(x) & x > cut.BF, 1, 0)))
        prior.modelprob <- apply(bf.max, 2, sum)/n.gene
        pm.est <- prior.modelprob
        n.model <- length(prior.modelprob)
        pm.prior0 <- pm.prior1 <- rep(0.5/(n.model - 1), n.model - 1)
        alpha <- 0.2
        for (j in 1:30) {
            pm.prior1 <- pm.est[-1]/apply(t(apply(out.bf[, -1], 1, function(x) x/(1 - sum(pm.prior0) + sum(x * pm.prior0)))),
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
            print(paste0(i, ". ", input.model.space[i]))
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
            denom <- t(apply(out.logbf[, -1], 1, function(y) exp(y - max(y) - log((1 - sum(pm.prior0))/max(exp(y)) + sum(exp(y +
                log(pm.prior0) - max(y)))))))
            pm.prior1 <- pm.est[-1]/apply(denom, 2, mean)
            pm.prior0 <- pm.prior0 + (pm.prior1 - pm.prior0) * alpha
        }
        prior.modelprob <- c(1 - sum(pm.prior0), pm.prior0)

        # calculate posterior model probability
        post.modelprob <- t(apply(out.logbf, 1, function(y) exp(y + log(prior.modelprob) - max(y) - log(sum(exp(y + log(prior.modelprob) -
            max(y)))))))
    }

    rownames(post.modelprob) <- rownames(dat.expr.counts)
    colnames(post.modelprob) <- input.model.space

    # calculate posterior model inclusion probability main effect, interaction effect, main or interaction effect of each
    # variable
    output1 <- Post.incl.modelprob(input.dat.pheno, input.model.space, post.modelprob, var.pool, interaction, ind.incl.add)
    post.incl.prob <- output1$post.incl.modelprob.Main
    if (!is.null(interaction)) {
        post.incl.prob <- output1$post.incl.modelprob.LargeTable
    } else {
        post.incl.prob.joint <- output1$post.incl.modelprob.Joint
    }
    if (!is.null(ind.incl.add))
        post.incl.prob.add <- output1$post.incl.modelprob.add
    index.incl <- output1$ind.incl

    # calculate estimated FDR and number of DE genes associated with main effect, joint main effect, interaction effect,
    # main or interaction effect of each variable
    output2 <- eFDR.DEG(post.incl.prob, cut.FDR)
    eFDR.notOrder <- output2$eFDR
    indicator.eFDR <- output2$indicator.eFDR
    id.DEgene <- output2$ind.DEG
    summary.nDEG <- rbind(output2$nDEG)
    rownames(summary.nDEG) <- "nDEG"

    if (is.null(interaction)) {
        output2.joint <- eFDR.DEG(post.incl.prob.joint, cut.FDR)
        eFDR.notOrder.joint <- output2.joint$eFDR
        indicator.eFDR.joint <- output2.joint$indicator.eFDR
        summary.nDEG.joint <- output2.joint$nDEG
    }

    if (!is.null(ind.incl.add)) {
        output2.add <- eFDR.DEG(post.incl.prob.add, cut.FDR)
        eFDR.notOrder.add <- output2.add$eFDR
        indicator.eFDR.add <- output2.add$indicator.eFDR
        summary.nDEG.add <- output2.add$nDEG
    }

    # identified DE genes associated with main effect, joint effects, interaction effect, main or interaction effect of
    # each variable
    output3 <- DEG.bestmodel(post.modelprob, input.model.space, id.DEgene)
    DEG.bestmodel <- output3$DEG.bestmodel

    if (is.null(ind.incl.add)) {
        if (!is.null(interaction)) {
            return(list(dat.expr.logcpm = input.dat.expr, weights = input.weights, dat.pheno.new = input.dat.pheno, model.space = input.model.space,
                post.modelprob = post.modelprob, post.incl.modelprob = post.incl.prob, index.incl = index.incl, eFDR = eFDR.notOrder,
                summary.nDEG = summary.nDEG, DEG.bestmodel = DEG.bestmodel))
        } else {
            return(list(dat.expr.logcpm = input.dat.expr, weights = input.weights, dat.pheno.new = input.dat.pheno, model.space = input.model.space,
                post.modelprob = post.modelprob, post.incl.modelprob = post.incl.prob, post.incl.modelprob.JointMain = post.incl.prob.joint,
                index.incl = index.incl, eFDR = eFDR.notOrder, eFDR.JointMain = eFDR.notOrder.joint, summary.nDEG = summary.nDEG,
                summary.nDEG.JointMain = summary.nDEG.joint, DEG.bestmodel = DEG.bestmodel))
        }
    } else {
        if (!is.null(interaction)) {
            return(list(dat.expr.logcpm = input.dat.expr, weights = input.weights, dat.pheno.new = input.dat.pheno, model.space = input.model.space,
                post.modelprob = post.modelprob, post.incl.modelprob = post.incl.prob, post.incl.modelprob.Add = post.incl.prob.add,
                index.incl = index.incl, eFDR = eFDR.notOrder, eFDR.Add = eFDR.notOrder.add, summary.nDEG = summary.nDEG, summary.nDEG.Add = summary.nDEG.add,
                DEG.bestmodel = DEG.bestmodel))
        } else {
            return(list(post.incl.modelprob = post.incl.prob, post.incl.modelprob.JointMain = post.incl.prob.joint, post.incl.modelprob.Add = post.incl.prob.add,
                index.incl = index.incl, eFDR = eFDR.notOrder, eFDR.JointMain = eFDR.notOrder.joint, eFDR.Add = eFDR.notOrder.add,
                summary.nDEG = summary.nDEG, summary.nDEG.JointMain = summary.nDEG.joint, summary.nDEG.Add = summary.nDEG.add,
                DEG.bestmodel = DEG.bestmodel))
        }
    }
}

