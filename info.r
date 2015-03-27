# The central function is "info", dimx" is required for "info" to function

dimx <- function (dd) 
{
    if (is.null(dim(dd))) 
        length(dd)
    else dim(dd)
}

info <- function (x) 
{
    cat("MODE:          ")
    cat(mode(x))
    cat("\n")
    cat("CLASS:         ")
    cat(class(x))
    cat("\n")
    cat("DIM or LENGTH: ")
    cat(dimx(x))
    cat("\n")
    cat("NAMES:         ")
    if (length(names(x)) > 6) 
        cat(names(x)[1:6], " ...")
    else cat(names(x))
    cat("\n")
}


