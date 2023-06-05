##' Obtain specific joint main effects output
##'
##' Obtain specific joint main effects output
##' @title get.jointMain
##' @param postprob.output The output from BMAseq.multi.postprob function, a list.
##' @param joint specific joint main effects with format 'A.B', a vector.
##' @param cut.FDR False discovery rate criterion for identifying DE genes, default value is 0.05.
##' @return A list consisting of
##' \item{post.incl.modelprob.JointMain}{Posterior inclusive model probability for specific joint effect of variables, only use if interaction is NULL.}
##' \item{eFDR.JointMain}{Estimated FDR for specific joint effect, only use if interaction is NULL.}
##' \item{summary.nDEG.JointMain}{A summary table of the number of identified DE gene for specific joint effcts of variables, only use if interaction is NULL.}
##' \item{ind.incl.JointMain}{Indices of the inclusive models for specific joint effcts of variables, only use if interaction is NULL.}
##' @export
##' @author Lingsong Meng

get.jointMain <- function(postprob.output, joint = NULL, cut.FDR = 0.05) {

    if (is.null(postprob.output))
        stop("Must give input for argument postprob.output.")

    model.space <- postprob.output$model.space
    post.modelprob <- postprob.output$post.modelprob
    var.pool <- postprob.output$var.pool
    interaction <- postprob.output$interaction
    if (!is.null(interaction))
        stop("Please only use this function when interaction is NULL.")

    if (is.null(joint)) {
        stop("Please give an input for argument joint.")
    } else {
        if (sum(grep(".", joint)) != sum(1:length(joint)))
            stop("Each two variables in joint input must be seperated by .")
        joint.spec.split.list <- strsplit(joint, split = ".", fixed = T)
        if (sum(lapply(joint.spec.split.list, function(x) length(x)) < 2) > 0)
            stop("Each joint effect in joint input must include two or more than two variables.")
        ind.joint.spec <- rep(0, length(joint))
    }

    # calculate posterior model inclusive probability by joint main effect if no interaction effect is considered
    n.var <- length(var.pool)
    n.gene <- nrow(post.modelprob)
    Genes <- rownames(post.modelprob)
    model.space.split1 <- strsplit(model.space, split = "+", fixed = T)
    var.joint.list <- list()
    for (j in 1:(n.var - 1)) {
        var.sub <- subsets(n.var, j + 1)
        if (!is.matrix(var.sub))
            var.sub <- matrix(var.sub, nrow = 1)
        for (i in 1:nrow(var.sub)) {
            var.joint.list[[length(var.joint.list) + 1]] <- var.pool[var.sub[i, ]]
        }
    }

    name.joint <- NULL
    ind.incl.joint <- list()
    post.incl.prob.joint <- matrix(NA, nrow = n.gene, ncol = length(var.joint.list))
    for (i in 1:length(var.joint.list)) {
        var.joint <- var.joint.list[[i]]
        name.tmp <- paste0(var.joint, collapse = ".")
        ind.incl.tmp <- which(unlist(lapply(model.space.split1, function(x) length(intersect(var.joint, x)) == length(var.joint))) ==
            TRUE)
        if (length(ind.incl.tmp) == 1)
            post.incl.prob.tmp <- post.modelprob[, ind.incl.tmp] else post.incl.prob.tmp <- apply(post.modelprob[, ind.incl.tmp], 1, sum)
        name.joint <- c(name.joint, name.tmp)
        ind.incl.joint[[i]] <- ind.incl.tmp
        post.incl.prob.joint[, i] <- post.incl.prob.tmp

        if (is.null(joint) == FALSE) {
            check.joint.spec <- unlist(lapply(joint.spec.split.list, function(x) length(intersect(var.joint, x)) == max(length(var.joint),
                length(x))))
            ind.joint.spec[which(check.joint.spec == TRUE)] <- i
        }
    }
    if (sum(ind.joint.spec == 0) != 0)
        stop("Variables in joint must be included in var.pool.")

    colnames(post.incl.prob.joint) <- name.joint
    rownames(post.incl.prob.joint) <- Genes
    names(ind.incl.joint) <- name.joint

    post.incl.prob.joint.spec <- NULL
    ind.incl.joint.spec <- NULL
    if (is.null(joint) == FALSE) {
        post.incl.prob.joint.spec <- post.incl.prob.joint[, ind.joint.spec]
        ind.incl.joint.spec <- ind.incl.joint[ind.joint.spec]
    }

    # Calculate eFDR and summary of nDEG
    tmp.joint <- eFDR.DEG(post.incl.prob.joint.spec, cut.FDR)
    eFDR.joint.spec <- tmp.joint$eFDR
    summary.nDEG.joint.spec <- tmp.joint$nDEG

    if (length(ind.joint.spec) == 1) {
        post.incl.prob.joint.spec <- as.matrix(post.incl.prob.joint.spec, ncol = 1)
        eFDR.joint.spec <- as.matrix(eFDR.joint.spec, ncol = 1)
        colnames(post.incl.prob.joint.spec) <- name[ind.joint.spec]
        colnames(eFDR.joint.spec) <- name[ind.joint.spec]
    }

    return(list(post.incl.modelprob.JointMain = post.incl.prob.joint.spec, eFDR.JointMain = eFDR.joint.spec, summary.nDEG.JointMain = summary.nDEG.joint.spec,
        ind.incl.JointMain = ind.incl.joint.spec))
}

