DEG.bestmodel <- function(post.modelprob, model.space, ind.DEG) {

    var <- names(ind.DEG)
    n.var <- length(ind.DEG)
    id.model.max1 <- model.max1 <- list()
    for (i in 1:n.var) {
        if (length(ind.DEG[[i]]) == 1) {
            id.model.max1[[i]] <- which.max(post.modelprob[ind.DEG[[i]], ])
        } else id.model.max1[[i]] <- as.vector(apply(post.modelprob[ind.DEG[[i]], ], 1, function(x) which.max(x)))
        model.max1[[i]] <- model.space[id.model.max1[[i]]]
    }

    # DE genes associated with each variable and the best model for each DE gene
    DEgene.modelmax1 <- list()
    for (i in 1:n.var) {
        DEgene.modelmax1[[i]] <- cbind(index.DEG = ind.DEG[[i]], name.DEG = rownames(post.modelprob)[ind.DEG[[i]]], best.model = model.max1[[i]])
    }
    names(DEgene.modelmax1) <- var

    model.max.list1 <- list()
    for (i in 1:n.var) {
        model.max.list1[[i]] <- sort(table(model.max1[[i]]), decreasing = T)
    }

    # DE genes identified by each best model and associated with each variable
    DEgene.max.list1 <- list()
    for (i in 1:n.var) {
        DEgene.max.list0 <- list()
        if (length(model.max.list1[[i]]) == 0) {
            DEgene.max.list0 <- NA
        } else {
            for (j in 1:length(model.max.list1[[i]])) {
                name.DEG <- rownames(post.modelprob)[ind.DEG[[i]][which(model.max1[[i]] == names(model.max.list1[[i]])[j])]]
                index.DEG <- ind.DEG[[i]][which(model.max1[[i]] == names(model.max.list1[[i]])[j])]
                DEgene.max.list0[[j]] <- cbind(index.DEG, name.DEG)
            }
        }
        names(DEgene.max.list0) <- names(model.max.list1[[i]])
        DEgene.max.list1[[i]] <- DEgene.max.list0
    }
    names(DEgene.max.list1) <- var


    return(list(DEG.bestmodel = DEgene.modelmax1, bestmodel.DEG = DEgene.max.list1))
}

