##' DE genes identification with obtained posterior model probability for each gene in multivariate analysis using BMAseq approach
##'
##' DE genes identification with obtained posterior model probability for each gene in multivariate analysis on RNA-seq count data using BMAseq approach
##' @title BMAseq.multi.DEG
##' @param dat.pheno Phenotypic data matrix with new interaction variables if interact=T.
##' @param model.space The model space including all possible models.
##' @param post.modelprob Posterior model probability for each gene.
##' @param var.pool Variables of interest, a vector.
##' @param interact Whether interaction terms are considered.
##' @param ind.incl.add Indices of additional included model.
##' @param cut.FDR False discovery rate criterion for identifying DE genes, default value is 0.05.
##' @return A list consisting of
##' \item{post.incl.modelprob.Main}{Posterior inclusive model probability for each gene associated with main effect of each variables.}
##' \item{post.incl.modelprob.Interaction}{Posterior inclusive model probability for each gene associated with interaction effect of each variables, will be output if interact=TRUE.}
##' \item{post.incl.modelprob.MainInteraction}{Posterior inclusive model probability for each gene associated with main or interaction effect of each variables, will be output if interact=TRUE.}
##' \item{post.incl.modelprob.add}{Posterior inclusive model probability for each gene associated with additional inclusive models, will be output if ind.incl.add != Null.}
##' \item{index.incl}{Indices of the inclusive models.}
##' \item{eFDR.Main}{Estimated FDR for each gene associated with main effect of each variables.}
##' \item{eFDR.Interaction}{Estimated FDR for each gene associated with interaction effect of each variables, will be output if interact=TRUE.}
##' \item{eFDR.MainInteraction}{Estimated FDR for each gene associated with main or interaction effect of each variables, will be output if interact=TRUE.}
##' \item{eFDR.add}{Estimated FDR for each gene associated with additional inclusive models, will be output if ind.incl.add != Null.}
##' \item{indicator.eFDR.Main}{Indicators for DE genes associated with main effect of each variables.}
##' \item{indicator.eFDR.Interaction}{Indicators for DE genes associated with interaction effect of each variables, will be output if interact=TRUE.}
##' \item{indicator.eFDR.MainInteraction}{Indicators for DE genes associated with mean or interaction effect of each variables, will be output if interact=TRUE.}
##' \item{indicator.eFDR.add}{Indicators for DE genes associated with additional inclusive models, will be output if ind.incl.add != Null.}
##' \item{summary.nDEG}{A summary table of the number of identified DE gene associated with each variable.}
##' \item{summary.nDEG.add}{A summary table of the number of identified DE gene associated with additional inclusive models, will be output if ind.incl.add != Null.}
##' \item{DEG.bestmodel.Main}{DE genes associated with main effect of each variables and the best model used to identify each DE gene.}
##' \item{DEG.bestmodel.Interaction}{DE genes associated with interaction effect of each variables and the best model used to identify each DE gene, will be output if interact=TRUE.}
##' \item{DEG.bestmodel.MainInteraction}{DE genes associated with main or interaction effect of each variables and the best model used to identify each DE gene, will be output if interact=TRUE.}
##' \item{bestmodel.DEG.Main}{The best models used to identify DE genes associated with main effect of each variable and the DE gene identified by each model.}
##' \item{bestmodel.DEG.Interaction}{The best models used to identify DE genes associated with interaction effect of each variable and the DE gene identified by each model, will be output if interact=TRUE.}
##' \item{bestmodel.DEG.MainInteraction}{The best models used to identify DE genes associated with main and interaction effect of each variable and the DE gene identified by each model, will be output if interact=TRUE.}
##' @export
##' @author Lingsong Meng


BMAseq.multi.DEG <- function(dat.pheno, model.space, post.modelprob, var.pool, interact = FALSE, ind.incl.add = NULL, cut.FDR = 0.05) {
    
    # calculate posterior model inclusion probability main effect, interaction effect, main or interaction effect of each variable
    output1 <- Post.incl.modelprob(dat.pheno, model.space, post.modelprob, var.pool, interact, ind.incl.add)
    post.incl.prob.main <- output1$post.incl.modelprob.Main
    if (interact) {
        post.incl.prob.int <- output1$post.incl.modelprob.Interaction
        post.incl.prob.mainInt <- output1$post.incl.modelprob.MainInteraction
    }
    if (is.null(ind.incl.add) == FALSE) 
        post.incl.prob.add <- output1$post.incl.modelprob.add
    index.incl <- output1$ind.incl
    
    
    # calculate estimated FDR and number of DE genes associated with main effect, interaction effect, main or interaction effect of each variable
    output2.main <- eFDR.DEG(post.incl.prob.main, cut.FDR)
    eFDR.notOrder1 <- output2.main$eFDR
    indicator.eFDR1 <- output2.main$indicator.eFDR
    id.DEgene1 <- output2.main$ind.DEG
    summary.nDEG <- rbind(output2.main$nDEG)
    rownames(summary.nDEG) <- c("nDEG.MainEffect")
    if (interact) {
        output2.int <- eFDR.DEG(post.incl.prob.int, cut.FDR)
        eFDR.notOrder2 <- output2.int$eFDR
        indicator.eFDR2 <- output2.int$indicator.eFDR
        id.DEgene2 <- output2.int$ind.DEG
        output2.mainInt <- eFDR.DEG(post.incl.prob.mainInt, cut.FDR)
        eFDR.notOrder3 <- output2.mainInt$eFDR
        indicator.eFDR3 <- output2.mainInt$indicator.eFDR
        id.DEgene3 <- output2.mainInt$ind.DEG
        summary.nDEG <- rbind(output2.main$nDEG, output2.int$nDEG, output2.mainInt$nDEG)
        rownames(summary.nDEG) <- c("nDEG.MainEffect", "nDEG.InteractionEffect", "nDEG.MainInteractionEffect")
        if (is.null(ind.incl.add) == FALSE) {
            output2.add <- eFDR.DEG(post.incl.prob.add, cut.FDR)
            eFDR.notOrder4 <- output2.add$eFDR
            indicator.eFDR4 <- output2.add$indicator.eFDR
            summary.nDEG.add <- output2.add$nDEG
        }
    }
    
    
    # identified DE genes associated with main effect, interaction effect, main or interaction effect of each variable
    output3.main <- DEG.bestmodel(post.modelprob, model.space, id.DEgene1)
    DEG.bestmodel.main <- output3.main$DEG.bestmodel
    bestmodel.nDEG.main <- output3.main$bestmodel.nDEG
    bestmodel.DEG.main <- output3.main$bestmodel.DEG
    if (interact) {
        output3.int <- DEG.bestmodel(post.modelprob, model.space, id.DEgene2)
        DEG.bestmodel.int <- output3.int$DEG.bestmodel
        bestmodel.nDEG.int <- output3.int$bestmodel.nDEG
        bestmodel.DEG.int <- output3.int$bestmodel.DEG
        output3.mainInt <- DEG.bestmodel(post.modelprob, model.space, id.DEgene3)
        DEG.bestmodel.mainInt <- output3.mainInt$DEG.bestmodel
        bestmodel.nDEG.mainInt <- output3.mainInt$bestmodel.nDEG
        bestmodel.DEG.mainInt <- output3.mainInt$bestmodel.DEG
    }
    
    
    if (is.null(ind.incl.add)) {
        if (interact) {
            return(list(post.incl.modelprob.Main = post.incl.prob.main, post.incl.modelprob.Interaction = post.incl.prob.int, post.incl.modelprob.MainInteraction = post.incl.prob.mainInt, 
                eFDR.Main = eFDR.notOrder1, eFDR.Interaction = eFDR.notOrder2, eFDR.MainInteraction = eFDR.notOrder3, indicator.eFDR.Main = indicator.eFDR1, indicator.eFDR.Interaction = indicator.eFDR2, 
                indicator.eFDR.MainInteraction = indicator.eFDR3, summary.nDEG = summary.nDEG, DEG.bestmodel.Main = DEG.bestmodel.main, DEG.bestmodel.Interaction = DEG.bestmodel.int, 
                DEG.bestmodel.MainInteraction = DEG.bestmodel.mainInt, bestmodel.DEG.Main = bestmodel.DEG.main, bestmodel.DEG.Interaction = bestmodel.DEG.int, bestmodel.DEG.MainInteraction = bestmodel.DEG.mainInt))
        } else return(list(post.incl.modelprob.Main = post.incl.prob.main, eFDR.Main = eFDR.notOrder1, indicator.eFDR.Main = indicator.eFDR1, summary.nDEG = summary.nDEG, 
            DEG.bestmodel.Main = DEG.bestmodel.main, bestmodel.DEG.Main = bestmodel.DEG.main))
    } else {
        if (interact) {
            return(list(post.incl.modelprob.Main = post.incl.prob.main, post.incl.modelprob.Interaction = post.incl.prob.int, post.incl.modelprob.MainInteraction = post.incl.prob.mainInt, 
                post.incl.modelprob.add = post.incl.prob.add, index.incl = index.incl, eFDR.Main = eFDR.notOrder1, eFDR.Interaction = eFDR.notOrder2, eFDR.MainInteraction = eFDR.notOrder3, 
                eFDR.add = eFDR.notOrder4, indicator.eFDR.Main = indicator.eFDR1, indicator.eFDR.Interaction = indicator.eFDR2, indicator.eFDR.MainInteraction = indicator.eFDR3, 
                indicator.eFDR.add = indicator.eFDR4, summary.nDEG = summary.nDEG, summary.nDEG.add = summary.nDEG.add, DEG.bestmodel.Main = DEG.bestmodel.main, DEG.bestmodel.Interaction = DEG.bestmodel.int, 
                DEG.bestmodel.MainInteraction = DEG.bestmodel.mainInt, bestmodel.DEG.Main = bestmodel.DEG.main, bestmodel.DEG.Interaction = bestmodel.DEG.int, bestmodel.DEG.MainInteraction = bestmodel.DEG.mainInt))
        } else return(list(post.incl.modelprob.Main = post.incl.prob.main, post.incl.modelprob.add = post.incl.prob.add, index.incl = index.incl, eFDR.Main = eFDR.notOrder1, 
            eFDR.add = eFDR.notOrder4, indicator.eFDR.Main = indicator.eFDR1, indicator.eFDR.add = indicator.eFDR4, summary.nDEG = summary.nDEG, summary.nDEG.add = summary.nDEG.add, 
            DEG.bestmodel.Main = DEG.bestmodel.main, bestmodel.DEG.Main = bestmodel.DEG.main))
    }
}

