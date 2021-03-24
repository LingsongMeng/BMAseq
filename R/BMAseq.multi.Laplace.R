##' Multivariate analysis using BMAseq approach
##'
##' Multivariate analysis on RNA-seq count data using BMAseq approach
##' @title BMAseq.multi
##' @param dat.expr.counts RNA-seq count data matrix (rows for genes and columns for subjects).
##' @param dat.pheno Phenotypic data matrix (rows for subjects and columns for variables).
##' @param var.pool Variables of interest, a vector.
##' @param max.nvar The maximum number of variables in a model.
##' @param interaction Specific interaction terms with input format "A&B", default is Null.
##' @param cut.BF Bayes factor criterion, default value is 1.
##' @param cut.FDR False discovery rate criterion for identifying DE genes, default value is 0.05.
##' @return A list consisting of
##' \item{dat.expr.logcpm}{Normalized RNA-seq data matrix (rows for genes and columns for subjects).}
##' \item{weights}{Estimated voom weights.}
##' \item{dat.pheno.new}{Phenotypic data matrix including new interaction variables, will be same as input dat.pheno if interaction is NULL.}
##' \item{model.space}{Model space including all possible models.}
##' \item{post.modelprob}{Posterior model probability for each gene.}
##' \item{post.incl.modelprob.Main}{Posterior inclusive model probability for each gene associated with main effect of each varibles.}
##' \item{post.incl.modelprob.Interaction}{Posterior inclusive model probability for each gene associated with interaction effect of each varibles, will be output if interaction is not NULL.}
##' \item{post.incl.modelprob.MainInteraction}{Posterior inclusive model probability for each gene associated with main or interaction effect of each varibles, will be output if interaction is not NULL.}
##' \item{summary.nDEG}{A summary table of the number of identified DE gene associated with each variable.}
##' \item{DEG.bestmodel.Main}{DE genes associated with main effect of each variables and the best model used to identify each DE gene.}
##' \item{DEG.bestmodel.Interaction}{DE genes associated with interaction effect of each variables and the best model used to identify each DE gene, will be output if interaction is not NULL.}
##' \item{DEG.bestmodel.MainInteraction}{DE genes associated with main or interaction effect of each variables and the best model used to identify each DE gene, will be output if interaction is not NULL.}
##' \item{bestmodel.DEG.Main}{The best models used to identify DE genes associated with main effect of each variable and the DE gene identified by each model.}
##' \item{bestmodel.DEG.Interaction}{The best models used to identify DE genes associated with interaction effect of each variable and the DE gene identified by each model, will be output if interaction is not NULL.}
##' \item{bestmodel.DEG.MainInteraction}{The best models used to identify DE genes associated with main and interaction effect of each variable and the DE gene identified by each model, will be output if interaction is not NULL.}
##' @export
##' @author Lingsong Meng


BMAseq.multi.Laplace <- function(dat.expr.counts, dat.pheno, var.pool, max.nvar, interaction = NULL, cut.BF = 1, cut.FDR = 0.05) {
    
    
    if (sum(colnames(dat.expr.counts) != rownames(dat.pheno)) > 0) 
        stop("Two datasets dat.expr.counts and dat.pheno must be matched by subjects.")
    if (is.null(var.pool) || length(var.pool) < 2) 
        stop("Must provide at least two variables of interest.")
    if (sum(var.pool %in% colnames(dat.pheno)) != length(var.pool)) 
        stop("Variables of interest must be in datasets dat.pheno.")
    if (is.null(max.nvar)) 
        stop("Must provide max.nvar.")
    
    if (is.null(interaction)) 
        interact = FALSE else interact = TRUE
    
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
        denom <- t(apply(out.logbf[, -1], 1, function(y) exp(y - max(y) - log((1 - sum(pm.prior0))/max(exp(y)) + sum(exp(y + log(pm.prior0) - max(y)))))))
        pm.prior1 <- pm.est[-1]/apply(denom, 2, mean)
        pm.prior0 <- pm.prior0 + (pm.prior1 - pm.prior0) * alpha
    }
    prior.modelprob <- c(1 - sum(pm.prior0), pm.prior0)
    
    
    # calculate posterior model probability
    post.modelprob <- t(apply(out.logbf, 1, function(y) exp(y + log(prior.modelprob) - max(y) - log(sum(exp(y + log(prior.modelprob) - max(y)))))))
    rownames(post.modelprob) <- rownames(dat.expr.counts)
    colnames(post.modelprob) <- input.model.space
    
    
    # calculate posterior model inclusion probability main effect, interaction effect, main or interaction effect of each variable
    output1 <- Post.incl.modelprob(input.dat.pheno, input.model.space, post.modelprob, var.pool, interact)
    post.incl.prob.main <- output1$post.incl.modelprob.Main
    if (interact) {
        post.incl.prob.int <- output1$post.incl.modelprob.Interaction
        post.incl.prob.mainInt <- output1$post.incl.modelprob.MainInteraction
    }
    
    
    # calculate estimated FDR and number of DE genes associated with main effect, interaction effect, main or interaction effect of each variable
    output2.main <- eFDR.DEG(post.incl.prob.main, cut.FDR)
    id.DEgene1 <- output2.main$ind.DEG
    summary.nDEG <- rbind(output2.main$nDEG)
    rownames(summary.nDEG) <- c("nDEG.MainEffect")
    if (interact) {
        output2.int <- eFDR.DEG(post.incl.prob.int, cut.FDR)
        id.DEgene2 <- output2.int$ind.DEG
        output2.mainInt <- eFDR.DEG(post.incl.prob.mainInt, cut.FDR)
        id.DEgene3 <- output2.mainInt$ind.DEG
        summary.nDEG <- rbind(output2.main$nDEG, output2.int$nDEG, output2.mainInt$nDEG)
        rownames(summary.nDEG) <- c("nDEG.MainEffect", "nDEG.InteractionEffect", "nDEG.MainInteractionEffect")
    }
    
    
    # identified DE genes associated with main effect, interaction effect, main or interaction effect of each variable
    output3.main <- DEG.bestmodel(post.modelprob, input.model.space, id.DEgene1)
    DEG.bestmodel.main <- output3.main$DEG.bestmodel
    bestmodel.nDEG.main <- output3.main$bestmodel.nDEG
    bestmodel.DEG.main <- output3.main$bestmodel.DEG
    if (interact) {
        output3.int <- DEG.bestmodel(post.modelprob, input.model.space, id.DEgene2)
        DEG.bestmodel.int <- output3.int$DEG.bestmodel
        bestmodel.nDEG.int <- output3.int$bestmodel.nDEG
        bestmodel.DEG.int <- output3.int$bestmodel.DEG
        output3.mainInt <- DEG.bestmodel(post.modelprob, input.model.space, id.DEgene3)
        DEG.bestmodel.mainInt <- output3.mainInt$DEG.bestmodel
        bestmodel.nDEG.mainInt <- output3.mainInt$bestmodel.nDEG
        bestmodel.DEG.mainInt <- output3.mainInt$bestmodel.DEG
    }
    
    
    
    if (interact) {
        return(list(dat.expr.logcpm = input.dat.expr, weights = input.weights, dat.pheno.new = input.dat.pheno, model.space = input.model.space, post.modelprob = post.modelprob, 
            post.incl.modelprob.Main = post.incl.prob.main, post.incl.modelprob.Interaction = post.incl.prob.int, post.incl.modelprob.MainInteraction = post.incl.prob.mainInt, 
            summary.nDEG = summary.nDEG, DEG.bestmodel.Main = DEG.bestmodel.main, DEG.bestmodel.Interaction = DEG.bestmodel.int, DEG.bestmodel.MainInteraction = DEG.bestmodel.mainInt, 
            bestmodel.DEG.Main = bestmodel.DEG.main, bestmodel.DEG.Interaction = bestmodel.DEG.int, bestmodel.DEG.MainInteraction = bestmodel.DEG.mainInt))
    } else {
        return(list(dat.expr.logcpm = input.dat.expr, weights = input.weights, dat.pheno.new = input.dat.pheno, model.space = input.model.space, post.modelprob = post.modelprob, 
            post.incl.modelprob.Main = post.incl.prob.main, summary.nDEG = summary.nDEG, DEG.bestmodel.Main = DEG.bestmodel.main, bestmodel.DEG.Main = bestmodel.DEG.main))
    }
}
