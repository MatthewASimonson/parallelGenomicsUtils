##
# SNP Associations Analysis Utilities
# By: Matthew A. Simonson


########## GWAS.frame(plkfile,phe,mod,cores)
# plkfile = character vector with name of binary plink files before extension .bim,...
# phe = phenotype file with missing values set to NA
# mod = character vector specifying linear or logistic regression model to be tested, example: c("BP ~ SNP + BMI + SNP*BMI")
# cores = number of processors to use

GWAS.frame <- function(plkfile,phe,mod,pheno.type=c("case.control","quantitative"),cores){
  require("OmicKriging",quietly = TRUE)
  require("parallel",quietly = TRUE)


  # function called for glm models (when case/control phenotype), executed during Map/Reduce:
  myglm <- function(i){
    SNP <- genotypes[[i]]
    data <- gwas.covars
    data$SNP <- SNP
    GenoMod <- glm(as.formula(mod),data=data,family=binomial())
    #return the vector
    cbind(matrix(t(summary(GenoMod)$coefficients),nrow=1,byrow=TRUE),matrix(unlist(t(anova(GenoMod,test="Chisq"))),nrow=1,byrow=TRUE))
  }

  # function called for lm models (when quantitative phenotype), executed during Map/Reduce:
  mylm <- function(i){
    SNP <- genotypes[[i]]
    data <- gwas.covars
    data$SNP <- SNP
    GenoMod <- lm(as.formula(mod),data=data)
    #return the vector
    cbind(matrix(t(summary(GenoMod)$coefficients),nrow=1,byrow=TRUE),matrix(unlist(t(anova(GenoMod))),nrow=1,byrow=TRUE))
  }

  # Main body of function:
  bim <- read.table(paste("",plkfile,".bim",sep=""),colClasses=c("numeric","character","numeric","numeric","character","character"))
  gdsFile <- "temp"
  convert_genotype_data(plkfile,gdsFile) # convert SNP data to an extern .gds format file for efficient storage and access
  gds <- openfn.gds(gdsFile) # open a connection to the .gds format file and store the location reference in the variable 'gds'
  snp <- bim$V2[1] # only read in first SNP column to check formatting
  g1 <- snpgdsGetGeno(gds,snp.id=snp, sample.id=phe$IID,verbose=TRUE)
  SNP.matrix <- cbind.data.frame(phe[,1],g1) # merge SNPs with ID info
  names(SNP.matrix) <- c(c("IID"))
  gwas.covars <-  merge(phe,SNP.matrix,by='IID')# Exclude subjects that do not have relevant phenotype and covariate information
  subj.ids <- gwas.covars$IID
  j <- as.list(1:22)
  snp.ids <- lapply(j,function(x){
    chrom.index <- bim$V1==x
    return(bim[chrom.index,2])
  })
  gwas.list <- list()

  if(pheno.type=="case.control"){
    for(i in 1:22){
      print(paste("Running logistic regression models on chromosome ",i,"",sep=""))
      geno.chrom <-  snpgdsGetGeno(gds,sample.id=subj.ids,snp.id=snp.ids[[i]],verbose=TRUE)
      genotypes <- as.list(as.data.frame(geno.chrom))
      Out <- mclapply(1:length(genotypes), myglm,mc.cores=cores) # returns list
      fitted <- data.frame(do.call(rbind,Out)) # make list a matrix

      # assign names to rows
      row.names(fitted) <- snp.ids[[i]]
      print(paste("Chromosome ",i," complete",sep=""))
      gwas.list[[i]] <- fitted
    }
  } else if (pheno.type=="quantitative"){
    for(i in 1:22){
      print(paste("Running linear regression models on chromosome ",i,"",sep=""))
      geno.chrom <-  snpgdsGetGeno(gds,sample.id=subj.ids,snp.id=snp.ids[[i]],verbose=TRUE)
      genotypes <- as.list(as.data.frame(geno.chrom))

      # Map and Reduce using mclapply:
      Out <- mclapply(1:length(genotypes), mylm,mc.cores=cores) # returns list
      fitted <- data.frame(do.call(rbind,Out)) # make list a matrix

      # assign names to rows
      row.names(fitted) <- snp.ids[[i]]
      print(paste("Chromosome ",i," complete",sep=""))
      gwas.list[[i]] <- fitted
    }
  }
  gwas.frame <- data.frame(do.call(rbind,gwas.list)) # convert list to single frame
  # assign names to columns
  SNP <- genotypes[[1]]
  dat.tmp <- cbind.data.frame(gwas.covars[,1:5],SNP)
  GenoMod.temp <- glm(as.formula(mod),data=dat.tmp,family=binomial())
  glmLab <- t(outer(row.names(summary(GenoMod.temp)$coefficients),row.names(t(summary(GenoMod.temp)$coefficients)),paste,sep="_"))
  lrtLab <- t(outer(row.names(anova(GenoMod.temp)),row.names(t(anova(GenoMod.temp,test="Chisq"))),paste,sep="_"))
  colnam <- c(as.vector(glmLab),as.vector(lrtLab))
  names(gwas.frame) <- colnam
  return(gwas.frame)
}

########## GPB.frame()



########## countCPU()

countCPU <- function(silent=FALSE){
  #----------------------------------------
  # INPUTS:
  # silent=TRUE/FALSE to print message that lists the CPU count, default is FALSE. When set to TRUE, 'countCPU' instead returns an integer count
  #
  # OUTPUT: printed message stating the number of CPUs on the current machine, or an integer count of the number of CPUs available
  #----------------------------------------
  time.stamp <- substr(Sys.time(),start=12,stop=27) # have time-stamp identifier for all temp files to prevent read/write mix-up
  system(paste("less /proc/cpuinfo | grep 'cpu cores' | head -n 1 > ",time.stamp,".cpu.tmp",sep=""))
  cpus <- read.table(paste("",time.stamp,".cpu.tmp",sep=""),header=FALSE)$V4
  system(paste("rm ",time.stamp,".cpu.tmp",sep=""))
  if(silent!=TRUE){
    print(paste("There are currently ",cpus," CPUs available",sep=""))
  }else{
    return(cpus)
  }
}

########## lit(x)

lit <-  function(x){
  require("NCBI2R",quietly = TRUE)
  tot <- GetPublishedGWAS(term=x)
  a <- c(which(names(tot)=="DiseaseTrait")) #
  b <- c(which(names(tot)=="Mapped_gene"))
  c <- c(which(names(tot)=="StrongestSNPRiskAllele"))
  d <- c(which(names(tot)=="RiskAlleleFrequency"))
  e <- c(which(names(tot)=="pValue"))
  f <- c(which(names(tot)=="ORorbeta"))
  g <- c(which(names(tot)=="X95CItext"))
  h <- c(which(names(tot)=="PUBMEDID"))
  y <- tot[tot$SNPs==x,c(a,b,c,d,e,f,g,h)]
  return(y)
}

########## readChroms(dataset,chrom)

readChroms <- function(dataset,chrom){
  require("snpStats",quietly = TRUE)

  bimt <- read.table(paste("",dataset,".bim",sep=""),header=FALSE)
  chrom.snps <- which(bimt[,1]==chrom)
  fam <- paste("",dataset,".fam",sep="")
  bim <- paste("",dataset,".bim",sep="")
  bed <- paste("",dataset,".bed",sep="")
  sample <- read.plink(bed, bim, fam,select.snps=chrom.snps)
  print("Reading Plink files complete!")
  print(sample$genotypes)
  return(sample)
}

########## readTraits(sample,cvFile,ID.col,trait.cols)

readTraits <- function(sample,cvFile,ID.col,trait.cols){
  covars.raw <- read.table("gain_BMI_covar.txt",header=TRUE)
  covars <- covars.raw[,trait.cols]
  row.names(covars) <- as.character(as.matrix(covars.raw[,ID.col]))

  phen <- sample$fam[,c(2,5,6)] # <--- Read in sex and phenotype information from Plink format .fam file
  inds <- as.data.frame(covars.raw[,ID.col])
  names(inds) <- c('member')
  phe <- merge(phen,inds,by="member") # <--- subject ID column to merge files and also exclude subjects without relevant variables
  subject.support <- cbind.data.frame(covars,phe[,2:ncol(phe)])
  return(subject.support)
}

########## genoSummary(sample)

genoSummary<- function(sample){
  require("snpStats",quietly = TRUE)
  snps.10 <- sample$genotypes
  sample.sum <- summary(snps.10)
  snpsum <- col.summary(snps.10)
  gsum <- list()
  gsum[[1]] <- sample.sum
  gsum[[2]] <- snpsum
  names(gsum) <- c("Across_SNP_summary","Individual_SNP_summaries")
  return(gsum)
}


#' Create a manhattan plot of genetic analysis p-values
#'
#' Create a manhattan plot of genetic analysis p-values
#' @param ifile the name of the output file of the genetic analysis
#' @param title title of plot
#' @param ... arguments to go into gap::mhtplot
#' @export
#' @import gap
manhattan<-function(ifile,title="Manhattan Plot",...){
  gwas<-read.table(ifile,header=T)
  d<-gwas[complete.cases(gwas),c("CHR", "BP","P")]
  ord <- order(d$CHR,d$BP)
  d <- d[ord,]
  top=ceiling(-log10(min(d[,"P"])))
  colors <- c(rep(c("blue","red"),15),"red")
  oldpar<-par()
  par(cex=0.6)
  mhtplot(d,control=mht.control(colors=colors,gap=1000),pch=19,srt=0,ylim=c(0,top),...)
  axis(2,cex.axis=2)
  title(title)
  #sline<--log10(3.60036E-05)
  #gline<--log10(1.8E-06)
  #abline(h=sline,col="blue")
  #abline(h=gline,col="green")
  abline(h=0)
}

#' Summarize top snp information
#'
#' Summarize top snp information. Can get top n snps, a proportion of top snps, or all snps with p-values
#' less than a specified cutoff. Each option will write to a unique file
#'@param ifile String corresponding to the name of the results file you want to get the top snps from
#'@param ofile String corresponding to the name of the resulting files (by default same as ifile)
#'@param topn Specify how many snps you want to retreive
#'@param topprop specifiy what proportion of snps you want to retreive
#'@param topcut specify the p-value cutoff you want to retreive
#'@param wd directory that all files are located in. Is by default set to the current R working directory

#'@export
topsnps<-function(ifile,ofile=NULL,topn=0,topprop=0,topcut=0,wd=""){
  if(topprop>1 | topprop<0) stop("topprop must be between 0 and 1")
  if(is.null(ofile)) ofile=ifile
  data=read.table(ifile,header=T,stringsAsFactors=F)
  ord=order(data[,"P"],decreasing=F)
  temp=which(colnames(data)=="P")
  data=data[ord,c(setdiff(1:ncol(data),temp),temp)]
  if(topn!=0){
    write.csv(data[1:topn,],paste0(wd,ofile,"_ntopsnps.csv"),row.names=F,quote=F)
  }
  if(topprop!=0){
    write.csv(data[1:ceiling(nrow(data)*topprop),],paste0(wd,ofile,"_proptopsnps.csv"),row.names=F,quote=F)
  }
  if(topcut!=0){
    write.csv(data[which(data[,"P"]<=topcut),],paste0(wd,ofile,"_cutofftopsnps.csv"),row.names=F,quote=F)
  }
}

#' Inserts a matched column containing rs format SNP identifiers into a data-frame that contains Affymetrix SNP identifiers
#' @param data-frame with column labeled "SNP" containing Affymetrix SNP identifiers
#' @export returns data-frame sorted by chromosome and physical position with input data and an additional column containing rs format SNP identifiers in column labeled "rs" matched to relevant Affymetrix identifiers
affy2rsID <- function(data){
  load("Affymetrix6.0Key.RData")
  r.mat <- cbind.data.frame(affy.key$SNP,affy.key$rsID)
  names(r.mat) <- c("SNP","rs")
  mdat <- merge(data,r.mat,by="SNP")
  done <- mdat[order(mdat$CHR,mdat$POS),]
  return(done)
}

#' Create a qq plot of genetic analysis p-values
#'
#' Create a qq plot of genetic analysis p-values
#' @param ifile the name of the output file of the genetic analysis
#' @param title title of plot
#' @param ... arguments to go into plot
#' @export
qq<-function(results,title="qq plot of p-values",...){
  obs=sort(results[,"P"],decreasing=F)
  ept=c(1:length(obs))/(length(obs))
  plot(-log10(ept),-log10(obs),col=4,xlab="Expected -log10(pvalue)",ylab="Observed -log10(pvalue)",main=title,...)
  abline(0,1,col="red")
}


#' This function takes file paths for PLINK binary files and loads them into
#' a Genomic Data Structure (GDS) stored on disk. This is a wrapper function
#' aroud utilities provided by the gdsfmt and SNPRelate packages.
convert_genotype_data <- function(bFile, gdsFile) {
  # Load dependencies:
  require("gdsfmt")
  require("SNPRelate")
  '%&%' <- function(a, b) paste(a, b, sep="")
  bedFile <- bFile %&% '.bed'
  bimFile <- bFile %&% '.bim'
  famFile <- bFile %&% '.fam'
  snpgdsBED2GDS(bed.fn=bedFile, fam.fn=famFile, bim.fn=bimFile, out.gdsfn=gdsFile)
}
