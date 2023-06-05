eFDR.DEG <- function(post.incl.modelprob, cut.FDR) {

    post.incl.modelprob <- as.matrix(post.incl.modelprob)
    Genes <- rownames(post.incl.modelprob)
    name.incl <- colnames(post.incl.modelprob)
    n.incl <- ncol(post.incl.modelprob)
    ngene <- nrow(post.incl.modelprob)
    ntop <- c(1:ngene)
    p <- 1 - post.incl.modelprob
    eFDR <- apply(1 - p, 2, function(x) cumsum(1 - x[order(1 - x)])[ntop]/ntop)
    ind.order.gene <- apply(1 - p, 2, function(x) order(1 - x))
    eFDR.notOrder <- matrix(nrow = ngene, ncol = n.incl)
    for (i in 1:n.incl) {
        eFDR.notOrder[, i] <- eFDR[, i][match(Genes, Genes[ind.order.gene[, i]])]
    }
    rownames(eFDR.notOrder) <- Genes
    colnames(eFDR.notOrder) <- name.incl
    indicator.eFDR <- as.matrix(apply(eFDR.notOrder, 2, function(x) ifelse(x < cut.FDR, 1, 0)))
    rownames(indicator.eFDR) <- Genes
    colnames(indicator.eFDR) <- name.incl
    nDEG <- apply(eFDR, 2, function(x) sum(x < cut.FDR))
    ind.DEG <- list()
    for (i in 1:n.incl) {
        ind.DEG[[i]] <- ind.order.gene[, i][eFDR[, i] < cut.FDR]
    }
    names(ind.DEG) <- name.incl

    return(list(eFDR = eFDR.notOrder, indicator.eFDR = indicator.eFDR, nDEG = nDEG, ind.DEG = ind.DEG))
}
