make_corr <- function(x1,rho){ # output vector that is correlated with input vector to the degree specified by "rho"
n     <- length(x1)            # length of vector                 # desired correlation = cos(angle)
theta <- acos(rho)             # corresponding angle        # fixed given data
x2    <- rnorm(n, mean=mean(x1), sd=sd(x1))      # new random data
X     <- cbind(x1, x2)         # matrix
Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)

Id   <- diag(n)                               # identity matrix
Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1

x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector
return(x)
} # end function "make_corr" 
