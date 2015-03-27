# Perform vectorized logistic regression across all columns of a data-frame

vecLogi <- function(y,X0.char){
  try(silent = TRUE, {

    X <- eval(parse(text=X0.char))
    XtY = crossprod(X, y) # XtY
    Bxy = solve(crossprod(X), XtY) # Use OLS approximation to reduce required iterations
    p <- as.vector(X%*%Bxy) # rough fitted probabilities

    # Refine OLS fit of probabilities with 2 iterations of maximum likelihood weighting
    # Iteration 1:
    w <- p * (1 - p) # weights for iteration 1
    z <- log(abs(p / (1 - p))) + (y - p) / (p * (1 - p)) # calculate deviance
    xtw <- t(X * w) # apply weighting
    U1 <- xtw %*% z # adjust using likelihood based estimates
    U2 <- solve(xtw %*% X, U1) # iteration 1 estimates
    o <- as.matrix(X %*% U2) # iteration 1 predicted odds
    p2 <- as.vector(1 / (1+exp(-o))) # odds to probs
    # Iteration 2:
    wi <- p2 * (1 - p2)
    zi <- log(abs(p2 / (1 - p2))) + (y - p2) / (p2 * (1 - p2))
    xtwi <- t(X * wi)
    M1 <- xtwi %*% zi
    XtwX <- xtwi%*%X
    # Results:
    M2 <- solve(XtwX) # predictor covariance matrix
    ors <- solve(XtwX, M1) # odds ratios
    ses <- sqrt(diag(M2)) # standard errors
    dn <- sum((y-mean(y))^2) # null model deviance
    ndf <- n-1 # null df
    da <- sum((y-(exp(X %*% ors)/(1+exp(X %*% ors))))^2) # full model deviance
    adf <- n-k # model df
    tvals = t(ors)/ses # coefficient t-stats
    mstats <- cbind(c(da,dn),c(adf,ndf)) # model df and deviance
    pvals <- 2*pnorm(tvals)
    fx <- list(as.numeric(ors),as.numeric(ses),as.numeric(tvals),as.numeric(pvals))#,mstats)
    names(fx) <- c("OddsRatio","StandardError","Tstat","Pvalue")#,"ModelDeviance")
    #
    return(fx)})
}
