Bayesfactor.Laplace <- function(X, y, w) {

    ### this function implement the Laplace approximation of the Bayes factor in Liang et.  al. (2008), which use
    ### Zellner-Siow prior.  The Bayes factor calculated is based on the Null model.

    ## X design matrix y vector of the response variable set parameter values
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

    ## find the mode
    myfunction <- function(g) -(1 - rsquared) * (p + 3) * g^3 + (n - p - 4 - 2 * (1 - rsquared)) * g^2 + (n * (2 - rsquared) -
        3) * g + n
    myroot = uniroot(myfunction, interval = c(0, 1e+10))
    ghat = myroot$root

    ## define the function h(g)
    h <- function(g) (n - p - 1)/2 * log(1 + g) - (n - 1)/2 * log(1 + (1 - rsquared) * g) + 0.5 * log(n/2) - log(sqrt(pi)) -
        1.5 * log(g) - n/2/g

    ## calculate 2nd derivative
    ddh = 0.5 * (n - 1) * (1 - rsquared)^2/(1 + ghat * (1 - rsquared))^2 - 0.5 * (n - p - 1)/(1 + ghat)^2 + 1.5/ghat^2 - n/ghat^3

    if (ddh > 0) {
        ## numerical correction.
        ddh = 0
    }

    logBF = log(sqrt(2 * pi)) - log(sqrt(-ddh)) + h(ghat)

    return(logBF)
}

