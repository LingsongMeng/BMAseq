Post.incl.modelprob <- function(dat.pheno, model.space, post.modelprob, var.pool, interaction, ind.incl.add = NULL) {

    n.var <- length(var.pool)
    n.gene <- nrow(post.modelprob)
    Genes <- rownames(post.modelprob)

    # calculate posterior model inclusion probability by main effect for each variable
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

    if (!is.null(interaction)) {
        # calculate posterior model inclusion probability by interaction effect for each variable
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

        # calculate posterior model inclusion probability by interaction effect for specific interaction terms
        ind.incl.int.spec <- list()
        n.int <- length(interaction)
        for (i in 1:n.int) {
            vars.int <- strsplit(interaction[i], split = "&", fixed = T)[[1]]
            var1.level <- paste0(vars.int[1], levels(dat.pheno[, colnames(dat.pheno) == vars.int[1]]))
            var2.level <- paste0(vars.int[2], levels(dat.pheno[, colnames(dat.pheno) == vars.int[2]]))
            ind.incl.int.spec[[i]] <- which(unlist(lapply(model.space.split3, function(x) (sum(var1.level %in% x) > 0) & (sum(var2.level %in%
                x) > 0))) == TRUE)
        }
        names(ind.incl.int.spec) <- interaction
        post.incl.prob.int.spec <- matrix(NA, nrow = n.gene, ncol = n.int)
        for (i in 1:n.int) {
            if (length(ind.incl.int.spec[[i]]) == 1)
                post.incl.prob.int.spec[, i] <- post.modelprob[, ind.incl.int.spec[[i]]] else post.incl.prob.int.spec[, i] <- apply(post.modelprob[, ind.incl.int.spec[[i]]], 1, sum)
        }
        colnames(post.incl.prob.int.spec) <- interaction
        rownames(post.incl.prob.int.spec) <- Genes

        # calculate posterior model inclusion probability by main effect or interaction effect for each gene
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
        ind.incl.model <- list(ind.incl.Main = ind.incl.main, ind.incl.Interaction = ind.incl.int, ind.incl.Interaction.specific = ind.incl.int.spec,
            ind.incl.MainInteraction = ind.incl.mainInt)

        # a large table with posterior model inclusion probability by main, interaction, main or interaction effect.
        for (i in 1:length(interaction)) {
            vars.int <- strsplit(interaction[i], split = "&", fixed = T)[[1]]
            if (i == 1) {
                tmp.mainInt <- post.incl.prob.main[, vars.int] + post.incl.prob.int.spec[, i]
                colnames(tmp.mainInt) <- c(paste0(vars.int[1], ".", interaction[i]), paste0(vars.int[2], ".", interaction[i]))
            } else {
                tmp <- post.incl.prob.main[, vars.int] + post.incl.prob.int.spec[, i]
                colnames(tmp) <- c(paste0(vars.int[1], ".", interaction[i]), paste0(vars.int[2], ".", interaction[i]))
                tmp.mainInt <- cbind(tmp.mainInt, tmp)
            }
        }
        post.incl.prob.LargeTable <- cbind(post.incl.prob.main, post.incl.prob.int.spec, tmp.mainInt)


    } else {

        # calculate posterior model inclusion probability by joint main effect if no interaction effect is considered
        # only consdier joint effect for two variabels and three variables.
        model.space.split1 <- strsplit(model.space, split = "+", fixed = T)
        var.joint.list <- list()
        var.sub2 <- subsets(n.var, 2)
        if (!is.matrix(var.sub2))
            var.sub2 <- matrix(var.sub2, nrow = 1)
        for (i in 1:nrow(var.sub2)) {
            var.joint.list[[i]] <- var.pool[var.sub2[i, ]]
        }
        if (n.var >= 3) {
            var.sub3 <- subsets(n.var, 3)
            if (!is.matrix(var.sub3))
                var.sub3 <- matrix(var.sub3, nrow = 1)
            for (i in 1:nrow(var.sub3)) {
                var.joint.list[[nrow(var.sub2) + i]] <- var.pool[var.sub3[i, ]]
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
        }
        colnames(post.incl.prob.joint) <- name.joint
        rownames(post.incl.prob.joint) <- Genes
        names(ind.incl.joint) <- name.joint
        ind.incl.model <- list(ind.incl.Main = ind.incl.main, ind.incl.Joint = ind.incl.joint)
    }

    # calculate posterior model inclusion probability with additional input ind.incl.add
    if (!is.null(ind.incl.add)) {
        if (!is.list(ind.incl.add))
            ind.incl.add <- list(ind.incl.add)
        post.incl.modelprob.add <- lapply(ind.incl.add, function(x) apply(as.matrix(post.modelprob[, x]), 1, sum))
        post.incl.modelprob.add <- matrix(unlist(post.incl.modelprob.add), ncol = length(ind.incl.add))
        rownames(post.incl.modelprob.add) <- Genes
        colnames(post.incl.modelprob.add) <- names(ind.incl.add)
        ind.incl.model[[length(ind.incl.model) + 1]] <- ind.incl.add
        names(ind.incl.model)[length(ind.incl.model)] <- "ind.incl.add"
    }


    if (is.null(ind.incl.add)) {
        if (!is.null(interaction)) {
            return(list(post.incl.modelprob.LargeTable = post.incl.prob.LargeTable, post.incl.modelprob.Main = post.incl.prob.main,
                post.incl.modelprob.Interaction = post.incl.prob.int, post.incl.modelprob.Interaction.specific = post.incl.prob.int.spec,
                post.incl.modelprob.MainInteraction = post.incl.prob.mainInt, ind.incl = ind.incl.model))
        } else {
            return(list(post.incl.modelprob.Main = post.incl.prob.main, post.incl.modelprob.Joint = post.incl.prob.joint, ind.incl = ind.incl.model))
        }

    } else {
        if (!is.null(interaction)) {
            return(list(post.incl.modelprob.LargeTable = post.incl.prob.LargeTable, post.incl.modelprob.Main = post.incl.prob.main,
                post.incl.modelprob.Interaction = post.incl.prob.int, post.incl.modelprob.Interaction.specific = post.incl.prob.int,
                post.incl.modelprob.MainInteraction = post.incl.prob.mainInt, post.incl.modelprob.add = post.incl.modelprob.add,
                ind.incl = ind.incl.model))
        } else {
            return(list(post.incl.modelprob.Main = post.incl.prob.main, post.incl.modelprob.Joint = post.incl.prob.joint, post.incl.modelprob.add = post.incl.modelprob.add,
                ind.incl = ind.incl.model))
        }
    }
}

