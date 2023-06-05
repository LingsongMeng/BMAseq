Bayesfactor <- function(X, y, w) {
    ### this function implement he direct numerical intergartion of the Bayes factor in Liang et. al. (2008), which use
    ### Zellner-Siow prior.  The Bayes factor calculated is based on the Null model.  - X design matrix - y vector of the
    ### response variable - w weights from voom

    X = as.matrix(X)
    y = as.vector(y)
    n = length(y)
    p <- dim(X)[2] - 1

    if (dim(X)[1] != n) {
        print("err:in function postmodelprob, dimension of X and y does not match!")
        return(-1e+05)
    }

    mylm = lm(y ~ X, weights = w)
    rsquared = summary(mylm)$r.squared

    exp.h <- function(g) exp((n - p - 1)/2 * log(1 + g) + (1 - n)/2 * log(1 + (1 - rsquared) * g) + log(sqrt(n/2)/sqrt(pi)) -
        1.5 * log(g) - n/2/g)

    result = integrate(exp.h, lower = 0, upper = Inf)$value
    return(result)
}
