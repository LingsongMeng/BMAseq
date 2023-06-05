my.lmFit <- function(data, design, weights, interact = F) {
    ## this function augment the lmFit output with additional statistics if interaction term is true, the function
    ## provides an additional estimate i.e beta2+beta4

    out <- lmFit(data, design, weights = weights)
    tstat <- (out$coef/out$stdev.unscaled)/out$sigma
    pval <- 2 * pt(-abs(tstat), df = out$df.resid)
    coef <- out$coef
    if (interact) {
        contrast.matrix <- cbind(var1.2 = c(0, 1, 0, 1), var1.avg = c(0, 1, 0, 0.5))
        out.contr <- contrasts.fit(out, contrasts = contrast.matrix)
        tstat.contr <- (out.contr$coef/out.contr$stdev.unscaled)/out.contr$sigma
        pval.contr <- 2 * pt(-abs(tstat.contr), df = out.contr$df.resid)
        coef <- cbind(out$coef[, c(1, 3, 2)], out.contr$coef, interact = out$coef[, 4])
        tstat <- cbind(tstat[, c(1, 3, 2)], tstat.contr, interact = tstat[, 4])
        pval <- cbind(pval[, c(1, 3, 2)], pval.contr, interact = pval[, 4])
    }
    colnames(tstat) <- paste0("t.", colnames(tstat))
    colnames(pval) <- paste0("pval.", colnames(pval))
    n <- dim(out$design)[1]
    K <- out$rank + 1
    rss <- out$df.resid * (out$sigma)^2
    loglik <- -(n/2) * log(2 * pi * rss/n) - n/2
    summary.lmFit <- data.frame(coef, tstat, pval, RSS = rss, loglik = loglik)

    return(list(out.lmFit = out, summary.lmFit = summary.lmFit))
}
