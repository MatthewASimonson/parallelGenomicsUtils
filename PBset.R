##########Gene Set Parametric bootstrap
# Summary: Functions to efficiently reduce a set of SNPs that are located within a set of genes into a single measure of association using parametric bootstrap
# Method:
# 1). Map SNPs to genes
# 2). Specify relevant covariates and null/alternative models. An empirical estimate of model significance is generated for each SNP by projecting
#     the respective models and recording the residual variance from a model comparison.

#     Implementation applies resampling with replacement so that the same individuals are used within each iteration across all SNPs in the gene set.
#     Any between SNP correlation effects (linkage-disequilibrium) are controlled for by holding individuals constant within each bootstrap iteration,
#     although between SNP correlation exists within each iteration, all iterations are independent.
#
#   2.a.) Note:  method is robust for evaluating non-additive interaction effects (UNLIKE PERMUTATION TESTS see Bůžková et. al. 2010 Annals of Human Genetics)
# 3.) Use PBSet() function to apply map/reduce of linear or logistic models within bootstap framework
#



#' Map SNPs to a list of genes
#'
#' Return a List with each element containing a character vector of the SNPs within gene ranges
#' @param name of Plink format binary files, exclude extension
#' @param data-frame containing phenotype information; with columns "IID" containing subject identifiers, and all other columns containing other predictors to be include in regression model; names of columsn should match model predictor names
#' @param character vector containing regression model to be examined
#' @param
#' @export data frame containing GWAS results; each row contains parameter summaries normally returned by functions lm() or glm() in addition to model comparison summaries returned by function anova()
#' @import OmicKriging
#' @import parallel

snps.to.genes <- function(snp.frame,plkfile,gene.info,distance){
  print(paste("Reading ",plkfile,".bim...",sep=""))
  bim <- read.table(paste("",plkfile,".bim",sep=""),colClasses=c("numeric","character","numeric","numeric","character","character"))
  names(snp.frame) <- c("V2")
  print("Extracting position information from relevant SNPs...")
  snp.info <- merge(snp.frame,bim,by="V2")
  names(snp.info) <- c("Name","Chr","CM","Position")
  j=1
  gene.snps=c()
  print("Mapping SNPs to list of genes...")
  for(j in 1:nrow(gene.info)){
    snps.inds=snp.info$Chr==gene.info$Chr[j] & snp.info$Position>=(gene.info[j,"Start"]-distance) & snp.info$Position<=(gene.info[j,"End"]+distance)
    gene.snps=c(gene.snps,list(snp.info$Name[which(snps.inds=="TRUE")]))
  }
  rm(j)

  names(gene.snps)=gene.info$Name
  gsf.list <- lapply(1:length(gene.snps),function(i){
    gene.snps <- gene.snps # create reference to object inside scope of function so not re-creating external on each iteration
    gene.snp.fr <- as.data.frame(as.character(unlist(gene.snps[[i]])))
    gene.snp.fr$gene <- rep(names(gene.snps[i]),nrow(gene.snp.fr))
    names(gene.snp.fr) <- c("Name","GENE")
    gene.snp.temp <- merge(gene.snp.fr,snp.info[,c(1,2,4)],by="Name")
    o.index <- order(gene.snp.temp$Chr,gene.snp.temp$Position)
    gene.snp.fr <- gene.snp.temp[o.index,]
    return(gene.snp.fr)
  })

  gene.snp.fr <- data.frame(do.call(rbind,gsf.list)) # make list a matrix
  return(gene.snp.fr)
} # End function


############
# Read in relevant phenotype, covariate, and subject identifier data
############

covars <- read.table("gain_BMI_covar.txt",header=TRUE)[,2:4] # example covariate/phenotype file with columns containing BMI and ancestry PCA

names(covars) <- c("IID","BMI","ANC") # <--- Any covariates can be used as long as subject ID column name matches across all files for easy merging

phen <- read.table("EA.GAIN.filtered.fam",header=FALSE)[,c(2,5,6)] # <--- Read in sex and phenotype information from Plink format .fam file
names(phen) <- c("IID","SEX","BP")
phe <- merge(phen,covars,by="IID") # <--- subject ID column to merge files and also exclude subjects without relevant variables
phe$BP[phe$BP==-9] <- NA # assign missing phenotype values of NA (Plink has these coded -9 by default)
phe$BP[phe$BP==1] <- 0 #
phe$BP[phe$BP==2] <- 1 #

# NOTE: model should be specified with predictors in the order most easily suited for model comparisons using anova()
# The predictor of greatest interest should probably be the last term

pheno.type <-c("case.control")
boots <- 1000
cores <- 24
gene.plkfile <- "WNT.snps"# file name that comes before binary plink extensions (exclude extensions)

snp.frame <- read.table("WNT.snps",header=FALSE)
gene.info <- read.table("WNT.genes.txt",header=FALSE)
names(gene.info) <- c("Chr","Start","End","Name")
distance = 20000 # 20KB up and downstream


altmod <- c("BP ~  ANC + SEX + BMI + SNP + SNP*ANC + SNP*SEX + BMI*ANC + BMI*SEX + BMI*SNP")
nullmod <- c("BP ~  ANC + SEX + BMI + SNP + SNP*ANC + SNP*SEX + BMI*ANC + BMI*SEX")

PBSet <- function(plkfile,snp.frame,gene.info,distance,phe,altmod,nullmod,boots,cores){
  require("OmicKriging",quietly = TRUE)
  require("parallel",quietly = TRUE)

  # Main body of function:
  bim <- read.table(paste("",plkfile,".bim",sep=""),colClasses=c("numeric","character","numeric","numeric","character","character"))
  snps.ids <- bim$V2
  gdsFile <- "temp"
  convert_genotype_data(plkfile,gdsFile) # convert SNP data to an extern .gds format file for efficient storage and access
  gds <- openfn.gds(gdsFile) # open a connection to the .gds format file and store the location reference in the variable 'gds'
  g <- snpgdsGetGeno(gds, sample.id=phe$IID,verbose=TRUE)
  SNP.matrix <- cbind.data.frame(phe[,1],g) # merge SNPs with ID info
  names(SNP.matrix) <- c("IID",snps.ids)
  geno.path <-  merge(phe,SNP.matrix,by='IID')# Exclude subjects that do not have relevant phenotype and covariate information
  path.covars <- geno.path[,1:(ncol(geno.path)-length(snps.ids))]
  path.snps <-geno.path[,((ncol(path.covars)+1):(ncol(geno.path)-ncol(path.covars)))]
  path.snps.list <- as.list(path.snps)
  snp.id.list <- as.list(names(path.snps))
  boot.index <- lapply(1:boots,function(b){return(sample(1:nrow(geno.path),replace=TRUE))})# generate list with index for each bootstrap iteration

  # Perform multi-threaded association analysis and parametric bootstrap (Map phase of map/reduce):
      pboot.tests <- mclapply(1:length(snp.id.list),function(i){
      snp.id.list <- snp.id.list
      path.snps.list <- path.snps.list
      path.covars <- path.covars
      SNP <- path.snps.list[[i]]
      data <- path.covars
      data$SNP <- SNP
      modA <- lm(altmod,data=data)
      modN <- lm(nullmod,data=data)
      lr.p<-anova(modA,modN)$"Pr(>F)"[2] # observed p-value for given SNP

      # apply parametric bootstrap to same re-sampling index across all SNPs
      boot.data <- data
      obs.probs <-predict(object = modN, type = "response") # projected probability of response=1 estimated from null model of observed data
      SNP.p.boot.list <- lapply(1:length(boot.index),function(j){
        temp.data <- boot.data # creating new instance of object inside function scope REALLY speeds it up, prevents mememory access bottlenecks
        temp.index <- boot.index
        probs.boot <- obs.probs[temp.index[[j]]]# probability of response=1 estimated from null model of given bootstrap sample
        ystar.boot <- rbinom(length(probs.boot),1,probs.boot) # predicted response from bootstrap iteration under the null
        temp.data$BP <- ystar.boot
        modA.boot <- lm(altmod,data=temp.data)
        modN.boot <- lm(nullmod,data=temp.data)
        lr.p.boot <-anova(modA.boot,modN.boot)$"Pr(>F)"[2] # observed p-value for given SNP
        return(lr.p.boot)
        })


      SNP.p.boot <- unlist(SNP.p.boot.list)
      pboot.compare <- list()
      observed <- lr.p
      null.pboot <- SNP.p.boot

      pboot.compare[[1]] <- observed
      pboot.compare[[2]] <- null.pboot
      pboot.compare[[3]] <- snp.id.list[[i]]
      return(pboot.compare)
    },mc.cores=cores) # End parametric bootstrap call

    # SAVE OUTPUT OF TESTS!!!
    save(pboot.tests,file="pboot.tests.RData") # With 24 cores and ~1500 samples, parametric bootstrap runs at 100 SNPs per minute; for Wnt ~ 40 minutes
    #

    # Extract and organize elements of results list (Reduce results)
    obs.pvals <- lapply(1:length(pboot.tests),function(i){
      pboot.tests <- pboot.tests
      return(pboot.tests[[i]][[1]])})

    obs.SNP <- lapply(1:length(pboot.tests),function(i){
      pboot.tests <- pboot.tests
      return(pboot.tests[[i]][[3]])})

    null.pvals <- lapply(1:length(pboot.tests),function(i){
      pboot.tests <- pboot.tests
      return(pboot.tests[[i]][[2]])})

    o.pframe <- data.frame(do.call(rbind,obs.pvals)) # make list a matrix
    names(o.pframe) <- "Observed P-value"
    o.pframe$SNP <- as.character(unlist(obs.SNP))

    # SNP LEVEL EMPIRICAL P-VALUES:
    # calculate empirical p-values for each SNP by examining the percent of null p-values that are more significant than the observed WITHIN EACH SNP ACROSS ITERATIONS (iteration ordering doesn't matter)
    SNP.emp <- unlist(lapply(1:nrow(o.pframe),function(i){
      o.pframe <- o.pframe
      null.pvals <- null.pvals

      emp <- signif(mean(null.pvals[[i]]<o.pframe[i,1]),digits=4)
      if(emp==0){
        quants <- qnorm(null.pvals[[i]])
        mq <- mean(quants)
        sq <- sd(quants)
        oq <- qnorm(o.pframe[i,1])
        qa<- (oq-mq)/sq # normal approximation of extreme p-value
        emp <- pnorm(qa)
        print(i)
      }
      if(emp==1){
        emp <- (boots/(boots+1))
      }
      emp <- emp-0.00001
      return(emp)
    }))

    emp.fr <- as.data.frame(SNP.emp)
    emp.fr$SNP <- o.pframe$SNP

    SNP.pval.fr <- merge(o.pframe,emp.fr,by="SNP")
    names(SNP.pval.fr) <- c("SNP","ObsP","EmpP")

    # Read in and format SNP data from genes within pathway:
    gene.map <- snps.to.genes(snp.frame,plkfile,gene.info,distance)
    names(gene.map) <- c("SNP","GENE","CHR","POS")
    snp.pv.map <- merge(gene.map,SNP.pval.fr,by="SNP")
    oi <- order(snp.pv.map$CHR,snp.pv.map$POS)
    snp.map <- snp.pv.map[oi,]
    snp.map <- snp.map[!duplicated(snp.map$SNP),]
    snp.map$EmpLogP <- -log10(snp.map$EmpP)

    # Add rs identifier column
    final.gm <- affy2rsID(snp.map)
    SNP.results <- final.gm # THIS CONTAINS FINAL SNP ASSOCIATIONS
    save(final.gm, file="pboot.pvals.RData")

















