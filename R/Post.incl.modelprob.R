Post.incl.modelprob <- function(dat.pheno, model.space, post.modelprob, var.pool, interact, ind.incl.add = NULL) {
    
    n.var <- length(var.pool)
    n.gene <- nrow(post.modelprob)
    Genes <- rownames(post.modelprob)
    
    # calculate posterior model inclusion probability by main effect
    model.space.split1 <- strsplit(model.space, split = "+", fixed = T)
    ind.incl.main <- list()
    for (i in 1:n.var) {
        ind.incl.main[[i]] <- which(unlist(lapply(model.space.split1, function(x) var.pool[i] %in% x)) == TRUE)
    }
    names(ind.incl.main) <- var.pool
    post.incl.prob.main <- matrix(NA, nrow = n.gene, ncol = n.var)
    for (i in 1:n.var) {
        if (length(ind.incl.main[[i]]) == 1) 
            post.incl.prob.main[, i] <- post.modelprob[, ind.incl.main[[i]]] else post.incl.prob.main[, i] <- apply(post.modelprob[, ind.incl.main[[i]]], 1, sum)
    }
    colnames(post.incl.prob.main) <- var.pool
    rownames(post.incl.prob.main) <- Genes
    ind.incl.model <- list(ind.incl.main = ind.incl.main)
    
    if (interact) {
        # calculate posterior model inclusion probability by interaction effect
        model.space.split2 <- lapply(model.space.split1, function(x) x[-c(1, which(x %in% var.pool))])
        model.space.split3 <- lapply(model.space.split2, function(x) unlist(strsplit(x, split = ".", fixed = T)))
        ind.incl.int <- list()
        for (i in 1:n.var) {
            var.level <- paste0(var.pool[i], levels(dat.pheno[, colnames(dat.pheno) == var.pool[i]]))
            ind.incl.int[[i]] <- which(unlist(lapply(model.space.split3, function(x) sum(var.level %in% x) > 0)) == TRUE)
        }
        names(ind.incl.int) <- var.pool
        post.incl.prob.int <- matrix(NA, nrow = n.gene, ncol = n.var)
        for (i in 1:n.var) {
            if (length(ind.incl.int[[i]]) == 1) 
                post.incl.prob.int[, i] <- post.modelprob[, ind.incl.int[[i]]] else post.incl.prob.int[, i] <- apply(post.modelprob[, ind.incl.int[[i]]], 1, sum)
        }
        colnames(post.incl.prob.int) <- var.pool
        rownames(post.incl.prob.int) <- Genes
        
        # calculate posterior model inclusion probability by main effect or interaction effect
        ind.incl.mainInt <- list()
        for (i in 1:n.var) {
            ind.incl.mainInt[[i]] <- union(ind.incl.main[[i]], ind.incl.int[[i]])
        }
        names(ind.incl.mainInt) <- var.pool
        post.incl.prob.mainInt <- matrix(NA, nrow = n.gene, ncol = n.var)
        for (i in 1:n.var) {
            if (length(ind.incl.mainInt[[i]]) == 1) 
                post.incl.prob.mainInt[, i] <- post.modelprob[, ind.incl.mainInt[[i]]] else post.incl.prob.mainInt[, i] <- apply(post.modelprob[, ind.incl.mainInt[[i]]], 1, sum)
        }
        colnames(post.incl.prob.mainInt) <- var.pool
        rownames(post.incl.prob.mainInt) <- Genes
        ind.incl.model <- list(ind.incl.Main = ind.incl.main, ind.incl.Interaction = ind.incl.int, ind.incl.MainInteraction = ind.incl.mainInt)
    }
    
    if (is.null(ind.incl.add) == FALSE) {
        if (is.list(ind.incl.add) == FALSE) 
            ind.incl.add <- list(ind.incl.add)
        post.incl.modelprob <- lapply(ind.incl.add, function(x) apply(as.matrix(post.modelprob[, x]), 1, sum))
        post.incl.modelprob <- matrix(unlist(post.incl.modelprob), ncol = length(ind.incl.add))
        rownames(post.incl.modelprob) <- Genes
        colnames(post.incl.modelprob) <- names(ind.incl.add)
        ind.incl.model[[length(ind.incl.model) + 1]] <- ind.incl.add
        names(ind.incl.model)[length(ind.incl.model)] <- "ind.incl.add"
    }
    
    
    if (is.null(ind.incl.add)) {
        if (interact) 
            return(list(post.incl.modelprob.Main = post.incl.prob.main, post.incl.modelprob.Interaction = post.incl.prob.int, post.incl.modelprob.MainInteraction = post.incl.prob.mainInt, 
                ind.incl = ind.incl.model)) else return(list(post.incl.modelprob.Main = post.incl.prob.main, ind.incl = ind.incl.model))
    } else {
        if (interact) 
            return(list(post.incl.modelprob.Main = post.incl.prob.main, post.incl.modelprob.Interaction = post.incl.prob.int, post.incl.modelprob.MainInteraction = post.incl.prob.mainInt, 
                post.incl.modelprob.add = post.incl.modelprob, ind.incl = ind.incl.model)) else return(list(post.incl.modelprob.Main = post.incl.prob.main, post.incl.modelprob.add = post.incl.modelprob, ind.incl = ind.incl.model))
    }
}

