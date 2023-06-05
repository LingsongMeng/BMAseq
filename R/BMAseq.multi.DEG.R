##' DE genes identification with obtained posterior model probability for each gene in multivariate analysis using BMAseq approach
##'
##' DE genes identification with obtained posterior model probability for each gene in multivariate analysis on RNA-seq count data using BMAseq approach
##' @title BMAseq.multi.DEG
##' @param postprob.output The output from BMAseq.multi.postprob function, a list.
##' @param ind.incl.add Indices of additional included model.
##' @param cut.FDR False discovery rate criterion for identifying DE genes, default value is 0.05.
##' @return A list consisting of
##' \item{post.incl.modelprob}{Posterior inclusive model probability for each gene associated with main effect of each variables if interaction is NULL; with main, interaction, main or interaction effects if interaction not NULL..}
##' \item{post.incl.modelprob.JointMain}{Posterior inclusive model probability for each gene associated with all possible joint effect of variables, will be output if interaction is NULL.}
##' \item{post.incl.modelprob.Add}{Posterior inclusive model probability for each gene associated with additional inclusive models, will be output if ind.incl.add is not Null.}
##' \item{index.incl}{Indices of the inclusive models.}
##' \item{eFDR}{Estimated FDR for each gene associated with main effect of each variables if interaction is NULL; with main, interaction, main or interaction effects if interaction is not NULL.}
##' \item{eFDR.JointMain}{Estimated FDR for each gene associated with all possible joint effect of variables, will be output if interaction is NULL.}
##' \item{eFDR.Add}{Estimated FDR for each gene associated with additional inclusive models, will be output if ind.incl.add is not Null.}
##' \item{summary.nDEG}{A summary table of the number of identified DE gene associated with main effect of each variable if interaction is NULL; with main, interaction, main or interaction effects if interaction is not NULL.}
##' \item{summary.nDEG.JointMain}{A summary table of the number of identified DE gene associated with all possible joint effcts of variables, will be output if interaction is NULL.}
##' \item{summary.nDEG.Add}{A summary table of the number of identified DE gene associated with additional inclusive models, will be output if ind.incl.add is not Null.}
##' \item{DEG.bestmodel}{DE genes associated with main effect of each variable if interaction is NULL; with main, interaction, main or interaction effects of each variable if interaction is not NULL; and the best model used to identify each DE gene.}
##' @export
##' @author Lingsong Meng


BMAseq.multi.DEG <- function(postprob.output, ind.incl.add = NULL, cut.FDR = 0.05) {

    if (is.null(postprob.output))
        stop("Must give input for argument postprob.output.")

    dat.pheno <- postprob.output$dat.pheno.new
    model.space <- postprob.output$model.space
    post.modelprob <- postprob.output$post.modelprob
    var.pool <- postprob.output$var.pool
    interaction <- postprob.output$interaction

    # calculate posterior model inclusion probability main effect, interaction effect, main or interaction effect of each
    # variable
    output1 <- Post.incl.modelprob(dat.pheno, model.space, post.modelprob, var.pool, interaction, ind.incl.add)
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
    output3 <- DEG.bestmodel(post.modelprob, model.space, id.DEgene)
    DEG.bestmodel <- output3$DEG.bestmodel

    if (is.null(ind.incl.add)) {
        if (!is.null(interaction)) {
            return(list(post.incl.modelprob = post.incl.prob, index.incl = index.incl, eFDR = eFDR.notOrder, summary.nDEG = summary.nDEG,
                DEG.bestmodel = DEG.bestmodel))
        } else {
            return(list(post.incl.modelprob = post.incl.prob, post.incl.modelprob.JointMain = post.incl.prob.joint, index.incl = index.incl,
                eFDR = eFDR.notOrder, eFDR.JointMain = eFDR.notOrder.joint, summary.nDEG = summary.nDEG, summary.nDEG.JointMain = summary.nDEG.joint,
                DEG.bestmodel = DEG.bestmodel))
        }

    } else {
        if (!is.null(interaction)) {
            return(list(post.incl.modelprob = post.incl.prob, post.incl.modelprob.Add = post.incl.prob.add, index.incl = index.incl,
                eFDR = eFDR.notOrder, eFDR.Add = eFDR.notOrder.add, summary.nDEG = summary.nDEG, summary.nDEG.Add = summary.nDEG.add,
                DEG.bestmodel = DEG.bestmodel))
        } else {
            return(list(post.incl.modelprob = post.incl.prob, post.incl.modelprob.JointMain = post.incl.prob.joint, post.incl.modelprob.Add = post.incl.prob.add,
                index.incl = index.incl, eFDR = eFDR.notOrder, eFDR.JointMain = eFDR.notOrder.joint, eFDR.Add = eFDR.notOrder.add,
                summary.nDEG = summary.nDEG, summary.nDEG.JointMain = summary.nDEG.joint, summary.nDEG.Add = summary.nDEG.add,
                DEG.bestmodel = DEG.bestmodel))
        }
    }
}
