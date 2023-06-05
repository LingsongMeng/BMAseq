##' Create default model space
##'
##' Create default model space consisting of all possible models
##' @title Modelspace
##' @param dat.pheno Phenotypic data matrix (rows for subjects and columns for variables).
##' @param var.pool Variables of interest, a vector.
##' @param max.nvar The maximum number of variables in a model.
##' @param interaction Specific interaction terms with input format 'A&B', default is Null.
##' @return A list consisting of
##' \item{model.space}{Model space including all possible models.}
##' \item{dat.pheno.new}{Phenotypic data matrix including new interaction variables, will be same as input dat.pheno if interaction=NULL.}
##' \item{var.pool}{Variables of interest, a vector.}
##' \item{interaction}{Specific interaction terms with format 'A&B', default is Null.}
##' @export
##' @author Lingsong Meng


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
            pheno.int[, 5 * i - 4] <- ifelse(dat.pheno[, ind.var.int[5 * i - 4, 1]] == levels(dat.pheno[, ind.var.int[5 * i -
                4, 1]])[1] & dat.pheno[, ind.var.int[5 * i - 4, 2]] == levels(dat.pheno[, ind.var.int[5 * i - 4, 2]])[1], 1,
                0)
            pheno.int[, 5 * i - 3] <- ifelse(dat.pheno[, ind.var.int[5 * i - 3, 1]] == levels(dat.pheno[, ind.var.int[5 * i -
                3, 1]])[1] & dat.pheno[, ind.var.int[5 * i - 3, 2]] == levels(dat.pheno[, ind.var.int[5 * i - 3, 2]])[2], 1,
                0)
            pheno.int[, 5 * i - 2] <- ifelse(dat.pheno[, ind.var.int[5 * i - 2, 1]] == levels(dat.pheno[, ind.var.int[5 * i -
                2, 1]])[2] & dat.pheno[, ind.var.int[5 * i - 2, 2]] == levels(dat.pheno[, ind.var.int[5 * i - 2, 2]])[1], 1,
                0)
            pheno.int[, 5 * i - 1] <- ifelse(dat.pheno[, ind.var.int[5 * i - 1, 1]] == levels(dat.pheno[, ind.var.int[5 * i -
                1, 1]])[2] & dat.pheno[, ind.var.int[5 * i - 1, 2]] == levels(dat.pheno[, ind.var.int[5 * i - 1, 2]])[2], 1,
                0)
            pheno.int[, 5 * i] <- ifelse((dat.pheno[, ind.var.int[5 * i, 1]] == levels(dat.pheno[, ind.var.int[5 * i, 1]])[1] &
                dat.pheno[, ind.var.int[5 * i, 2]] == levels(dat.pheno[, ind.var.int[5 * i, 2]])[1]) | (dat.pheno[, ind.var.int[5 *
                i, 1]] == levels(dat.pheno[, ind.var.int[5 * i, 1]])[2] & dat.pheno[, ind.var.int[5 * i, 2]] == levels(dat.pheno[,
                ind.var.int[5 * i, 2]])[2]), 1, 0)
            name.int[5 * i - 4] <- paste0(var.int[5 * i - 4, 1], levels(dat.pheno[, ind.var.int[5 * i - 4, 1]])[1], ".", var.int[5 *
                i - 4, 2], levels(dat.pheno[, ind.var.int[5 * i - 4, 2]])[1])
            name.int[5 * i - 3] <- paste0(var.int[5 * i - 3, 1], levels(dat.pheno[, ind.var.int[5 * i - 3, 1]])[1], ".", var.int[5 *
                i - 3, 2], levels(dat.pheno[, ind.var.int[5 * i - 3, 2]])[2])
            name.int[5 * i - 2] <- paste0(var.int[5 * i - 2, 1], levels(dat.pheno[, ind.var.int[5 * i - 2, 1]])[2], ".", var.int[5 *
                i - 2, 2], levels(dat.pheno[, ind.var.int[5 * i - 2, 2]])[1])
            name.int[5 * i - 1] <- paste0(var.int[5 * i - 1, 1], levels(dat.pheno[, ind.var.int[5 * i - 1, 1]])[2], ".", var.int[5 *
                i - 1, 2], levels(dat.pheno[, ind.var.int[5 * i - 1, 2]])[2])
            name.int[5 * i] <- paste0(var.int[5 * i, 1], levels(dat.pheno[, ind.var.int[5 * i, 1]])[1], ".", var.int[5 * i,
                2], levels(dat.pheno[, ind.var.int[5 * i, 2]])[1], ".", var.int[5 * i, 1], levels(dat.pheno[, ind.var.int[5 *
                i, 1]])[2], ".", var.int[5 * i, 2], levels(dat.pheno[, ind.var.int[5 * i, 2]])[2])
        }
        colnames(pheno.int) <- name.int
        dat.pheno.new <- cbind(dat.pheno, pheno.int)
        for (i in (n.var + 1):ncol(dat.pheno.new)) {
            dat.pheno.new[, i] <- factor(dat.pheno.new[, i])
        }

        # Reove the models including variables in interaction term (only consider 1 interaction)
        MS.splitList <- sapply(strsplit(model.space, split = "+", fixed = T), function(x) x[-1])
        keep.bool <- sapply(MS.splitList, function(x) length(intersect(x, unlist(int.split.list))) == 0)
        MS.adding <- sapply(MS.splitList, function(x) paste0("+", paste(x, collapse = "+")))
        MS.adding <- MS.adding[keep.bool][-1]
        len.MS.adding <- length(MS.adding)  # number of models adding to the basic 16 patterns


        # Pattern1: Main Effect and Interactions (main effect variables, interactions)
        var.mainInt <- matrix(nrow = 4 * n.var.sub, ncol = 2)
        for (i in 1:n.var.sub) {
            var.mainInt[4 * i - 3, ] <- c(var.int[5 * i - 3, 1], name.int[5 * i - 3])
            var.mainInt[4 * i - 2, ] <- c(var.int[5 * i - 2, 2], name.int[5 * i - 2])
            var.mainInt[4 * i - 1, ] <- c(var.int[5 * i - 1, 1], name.int[5 * i - 1])
            var.mainInt[4 * i, ] <- c(var.int[5 * i - 1, 2], name.int[5 * i - 1])
        }
        # Pattern2: 2 Interactions (interaction1, interaction2)
        var.int2 <- matrix(nrow = 2 * n.var.sub, ncol = 2)
        for (i in 1:n.var.sub) {
            var.int2[2 * i - 1, ] <- c(name.int[5 * i], name.int[5 * i - 2])
            var.int2[2 * i, ] <- c(name.int[5 * i], name.int[5 * i - 1])
        }
        # Pattern3: Add 2 Main Effects and Interactions (main effect1, main effect2, interaction)
        var.main2Int <- matrix(nrow = n.var.sub, ncol = 3)
        for (i in 1:n.var.sub) {
            var.main2Int[i, ] <- c(var.int[5 * i - 1, 1], var.int[5 * i - 1, 2], name.int[5 * i - 1])
        }

        # model space including models with interactions
        MS.basic <- "~1"
        MS.addMain <- NULL
        for (i in 1:nrow(int.split)) {
            MS.addMain <- c(MS.addMain, c(paste0("~1+", int.split[i, 1]), paste0("~1+", int.split[i, 2]), paste0("~1+", int.split[i,
                1], "+", int.split[i, 2])))
        }
        MS.basic <- c(MS.basic, MS.addMain[which(duplicated(MS.addMain) == F)])

        for (i in 1:nrow(var.int)) {
            MS.basic <- c(MS.basic, paste("~1", name.int[i], sep = "+"))
        }
        for (i in 1:nrow(var.mainInt)) {
            MS.basic <- c(MS.basic, paste("~1", paste(var.mainInt[i, ], collapse = "+"), sep = "+"))
        }
        for (i in 1:nrow(var.int2)) {
            MS.basic <- c(MS.basic, paste("~1", paste(var.int2[i, ], collapse = "+"), sep = "+"))
        }
        for (i in 1:nrow(var.main2Int)) {
            MS.basic <- c(MS.basic, paste("~1", paste(var.main2Int[i, ], collapse = "+"), sep = "+"))
        }
        len.MS.before.add <- length(MS.basic)  # 16 patterns related to interactions

        for (i in 1:len.MS.adding) {
            MS.AddNew <- paste0(MS.basic[1:len.MS.before.add], MS.adding[i])
            MS.basic <- c(MS.basic, MS.AddNew)
        }
        model.space <- MS.basic
    }


    if (!is.null(interaction))
        return(list(model.space = model.space, dat.pheno.new = dat.pheno.new, var.pool = var.pool, interaction = interaction)) else return(list(model.space = model.space, dat.pheno.new = dat.pheno, var.pool = var.pool, interaction = interaction))
}
