
  
  X0.char <- paste("X0 <- cbind(1,",paste(IVs,collapse=","),")",collapse=",") # specify predictors
  
  fastSNPlm <- function(SNP){         
    X <- eval(parse(text=X0.char))
    XtY = crossprod(X, y) # XtY
    Bxy = solve(crossprod(X), XtY) 
    res <- y-as.vector(X%*%Bxy) 
    ser = sqrt(diag(1/(n-k) * as.numeric(t(res)%*%res) * solve(t(X)%*%X)))
    tvals = t(Bxy)/ser
    pvals <- 2*pnorm(abs(tvals),lower.tail=FALSE)
    fx <- cbind(as.numeric(Bxy),as.numeric(ser),as.numeric(tvals),as.numeric(pvals))
    return(fx) # t(rbind(t(U2),StdErr,tval,pval))
  }
  
  gwas.list  <- mclapply(1:22,function(i){
    options(warn=-1)
    gwas.chrom <- eval(parse(text=paste("apply.gdsn(index.gdsn('chr",i,".gds','genotype'),margin=2,FUN=fastSNPlm)",sep="")))
    return(gwas.chrom)
  },mc.cores=cores)
  
