Modelspace <- function(dat.pheno, var.pool, max.nvar, interaction = NULL) {
    
    id.var <- match(var.pool, colnames(dat.pheno))
    dat.pheno <- dat.pheno[id.var]
    n.var <- length(var.pool)
    # model space with no interaction
    model.space <- "~1"
    for (i in 1:min(n.var, max.nvar)) {
        var.sub <- subsets(n.var, i)
        if (!is.matrix(var.sub)) 
            var.sub <- matrix(var.sub, nrow = 1)
        for (j in 1:nrow(var.sub)) {
            if (ncol(var.sub) == 1) 
                model.space <- c(model.space, paste("~1", var.pool[var.sub[j, ]], sep = "+")) else model.space <- c(model.space, paste("~1", paste(var.pool[var.sub[j, ]], collapse = "+"), sep = "+"))
        }
    }
    
    
    # Interactions
    if (!is.null(interaction)) {
        
        if (sum(grep("&", interaction)) != sum(1:length(interaction))) 
            stop("Two variables in interactions must be seperated by &.")
        int.split.list <- strsplit(interaction, split = "&")
        if (sum(lapply(int.split.list, function(x) length(x)) != 2) > 0) 
            stop("Interactions must include two variables.")
        int.split <- matrix(unlist(int.split.list), ncol = 2, byrow = T)
        n.var.sub <- nrow(int.split)
        var.int <- cbind(rep(int.split[, 1], each = 5), rep(int.split[, 2], each = 5))
        var.int.match <- match(var.int, var.pool)
        if (sum(is.na(var.int.match)) > 0) 
            stop("Variables in interactions must be included in var.pool.")
        ind.var.int <- matrix(var.int.match, ncol = 2)
        n.sub <- nrow(dat.pheno)
        
        # Create new interaction variables
        pheno.int <- matrix(nrow = n.sub, ncol = 5 * n.var.sub)
        name.int <- NULL
        for (i in 1:n.var.sub) {
            pheno.int[, 5 * i - 4] <- ifelse(dat.pheno[, ind.var.int[5 * i - 4, 1]] == levels(dat.pheno[, ind.var.int[5 * i - 4, 1]])[1] & dat.pheno[, ind.var.int[5 * 
                i - 4, 2]] == levels(dat.pheno[, ind.var.int[5 * i - 4, 2]])[1], 1, 0)
            pheno.int[, 5 * i - 3] <- ifelse(dat.pheno[, ind.var.int[5 * i - 3, 1]] == levels(dat.pheno[, ind.var.int[5 * i - 3, 1]])[1] & dat.pheno[, ind.var.int[5 * 
                i - 3, 2]] == levels(dat.pheno[, ind.var.int[5 * i - 3, 2]])[2], 1, 0)
            pheno.int[, 5 * i - 2] <- ifelse(dat.pheno[, ind.var.int[5 * i - 2, 1]] == levels(dat.pheno[, ind.var.int[5 * i - 2, 1]])[2] & dat.pheno[, ind.var.int[5 * 
                i - 2, 2]] == levels(dat.pheno[, ind.var.int[5 * i - 2, 2]])[1], 1, 0)
            pheno.int[, 5 * i - 1] <- ifelse(dat.pheno[, ind.var.int[5 * i - 1, 1]] == levels(dat.pheno[, ind.var.int[5 * i - 1, 1]])[2] & dat.pheno[, ind.var.int[5 * 
                i - 1, 2]] == levels(dat.pheno[, ind.var.int[5 * i - 1, 2]])[2], 1, 0)
            pheno.int[, 5 * i] <- ifelse((dat.pheno[, ind.var.int[5 * i, 1]] == levels(dat.pheno[, ind.var.int[5 * i, 1]])[1] & dat.pheno[, ind.var.int[5 * i, 2]] == 
                levels(dat.pheno[, ind.var.int[5 * i, 2]])[1]) | (dat.pheno[, ind.var.int[5 * i, 1]] == levels(dat.pheno[, ind.var.int[5 * i, 1]])[2] & dat.pheno[, ind.var.int[5 * 
                i, 2]] == levels(dat.pheno[, ind.var.int[5 * i, 2]])[2]), 1, 0)
            name.int[5 * i - 4] <- paste0(var.int[5 * i - 4, 1], levels(dat.pheno[, ind.var.int[5 * i - 4, 1]])[1], ".", var.int[5 * i - 4, 2], levels(dat.pheno[, ind.var.int[5 * 
                i - 4, 2]])[1])
            name.int[5 * i - 3] <- paste0(var.int[5 * i - 3, 1], levels(dat.pheno[, ind.var.int[5 * i - 3, 1]])[1], ".", var.int[5 * i - 3, 2], levels(dat.pheno[, ind.var.int[5 * 
                i - 3, 2]])[2])
            name.int[5 * i - 2] <- paste0(var.int[5 * i - 2, 1], levels(dat.pheno[, ind.var.int[5 * i - 2, 1]])[2], ".", var.int[5 * i - 2, 2], levels(dat.pheno[, ind.var.int[5 * 
                i - 2, 2]])[1])
            name.int[5 * i - 1] <- paste0(var.int[5 * i - 1, 1], levels(dat.pheno[, ind.var.int[5 * i - 1, 1]])[2], ".", var.int[5 * i - 1, 2], levels(dat.pheno[, ind.var.int[5 * 
                i - 1, 2]])[2])
            name.int[5 * i] <- paste0(var.int[5 * i, 1], levels(dat.pheno[, ind.var.int[5 * i, 1]])[1], ".", var.int[5 * i, 2], levels(dat.pheno[, ind.var.int[5 * i, 
                2]])[1], ".", var.int[5 * i, 1], levels(dat.pheno[, ind.var.int[5 * i, 1]])[2], ".", var.int[5 * i, 2], levels(dat.pheno[, ind.var.int[5 * i, 2]])[2])
        }
        colnames(pheno.int) <- name.int
        dat.pheno.new <- cbind(dat.pheno, pheno.int)
        for (i in (n.var + 1):ncol(dat.pheno.new)) {
            dat.pheno.new[, i] <- factor(dat.pheno.new[, i])
        }
        
        # Main Effect and Interactions (main effect variables, interactions)
        var.mainInt <- matrix(nrow = 4 * n.var.sub, ncol = 2)
        for (i in 1:n.var.sub) {
            var.mainInt[4 * i - 3, ] <- c(var.int[5 * i - 3, 1], name.int[5 * i - 3])
            var.mainInt[4 * i - 2, ] <- c(var.int[5 * i - 2, 2], name.int[5 * i - 2])
            var.mainInt[4 * i - 1, ] <- c(var.int[5 * i - 1, 1], name.int[5 * i - 1])
            var.mainInt[4 * i, ] <- c(var.int[5 * i - 1, 2], name.int[5 * i - 1])
        }
        # 2 Interactions (interaction1, interaction2)
        var.int2 <- matrix(nrow = 2 * n.var.sub, ncol = 2)
        for (i in 1:n.var.sub) {
            var.int2[2 * i - 1, ] <- c(name.int[5 * i], name.int[5 * i - 2])
            var.int2[2 * i, ] <- c(name.int[5 * i], name.int[5 * i - 1])
        }
        # 2 Main Effects and Interactions (main effect1, main effect2, interaction)
        var.main2Int <- matrix(nrow = n.var.sub, ncol = 3)
        for (i in 1:n.var.sub) {
            var.main2Int[i, ] <- c(var.int[5 * i - 1, 1], var.int[5 * i - 1, 2], name.int[5 * i - 1])
        }
        # model space including models with interactions
        for (i in 1:nrow(var.int)) {
            model.space <- c(model.space, paste("~1", name.int[i], sep = "+"))
        }
        for (i in 1:nrow(var.mainInt)) {
            model.space <- c(model.space, paste("~1", paste(var.mainInt[i, ], collapse = "+"), sep = "+"))
        }
        for (i in 1:nrow(var.int2)) {
            model.space <- c(model.space, paste("~1", paste(var.int2[i, ], collapse = "+"), sep = "+"))
        }
        for (i in 1:nrow(var.main2Int)) {
            model.space <- c(model.space, paste("~1", paste(var.main2Int[i, ], collapse = "+"), sep = "+"))
        }
        
        # remove model including sparse interactions
        min.n.level <- apply(dat.pheno.new[, (n.var + 1):ncol(dat.pheno.new)], 2, function(x) min(table(x)))
        ind.sparse <- which(min.n.level < 10 | min.n.level == n.sub)
        inter.sparse <- names(ind.sparse)
        model.space.split <- strsplit(model.space, split = "+", fixed = T)
        ind.model.rm <- NULL
        for (i in 1:length(ind.sparse)) {
            sparse.true <- unlist(lapply(model.space.split, function(x) inter.sparse[i] %in% x))
            ind.model.rm <- c(ind.model.rm, which(sparse.true == TRUE))
        }
        dup.true <- duplicated(ind.model.rm)
        if (sum(dup.true) > 0) {
            ind.dup <- which(dup.true == TRUE)
            ind.model.rm <- ind.model.rm[-ind.dup]
        }
        if (length(ind.model.rm) > 0) 
            model.space <- model.space[-ind.model.rm]
    }
    
    
    if (!is.null(interaction)) 
        return(list(model.space = model.space, dat.pheno.new = dat.pheno.new)) else return(list(model.space = model.space, dat.pheno.new = dat.pheno))
}


