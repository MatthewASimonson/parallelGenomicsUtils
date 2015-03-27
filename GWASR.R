#####
# Run GWAS Models in EA GAIN Bipolar

setwd("/data5/bsi/psych/s105109.goa/Bipolar_WNT_Pathway/Analysis/EA/GWAS")

# Read in relevant phenotype, covariate, and subject identifier data
covars <- read.table("gain_BMI_covar.txt",header=TRUE)[,2:4] # example covariate/phenotype file with columns containing BMI and ancestry PCA

names(covars) <- c("IID","BMI","ANC") # <--- Any covariates can be used as long as subject ID column name matches across all files for easy merging

phen <- read.table("EA.GAIN.filtered.fam",header=FALSE)[,c(2,5,6)] # <--- Read in sex and phenotype information from Plink format .fam file
names(phen) <- c("IID","SEX","BP")
phe <- merge(phen,covars,by="IID") # <--- subject ID column to merge files and also exclude subjects without relevant variables
phe$BP[phe$BP==-9] <- NA # assign missing phenotype values of NA (Plink has these coded -9 by default)
phe$BP[phe$BP==1] <- 0 #
phe$BP[phe$BP==2] <- 1 # 


plkfile <- "EA.GAIN.filtered"# file name that comes before binary plink extensions (exclude extensions)
cores <-24

altmod <- c("BP ~  ANC + SEX + BMI + SNP + SNP*ANC + SNP*SEX + BMI*ANC + BMI*SEX + BMI*SNP") # Interaction SNP effect 
nullmod <- c("BP ~  ANC + SEX + BMI + SNP + BMI*SNP") # effect of not controlling for interaction covariates

gwas.full <- GWASlog(plkfile,phe,altmod,nullmod,cores)
save(gwas.full, file = "gwas.full.Rdata")

altmod.b <- c("BP ~  ANC + SEX + BMI + SNP") # SNP main effect controlling for BMI
nullmod.b <- c("BP ~  ANC + SEX + SNP") #  effect of not controlling for BMI

gwas.bmi <-  GWASlog(plkfile,phe,altmod.b,nullmod.b,cores)
save(gwas.bmi, file = "gwas.reduced.RData")
#
altmod.c <- c("BP ~  ANC + SEX + SNP") # SNP main effect
nullmod.c <- c("BP ~ SNP") # effect of not controlling for Ancestry and sex

gwas.add <-  GWASlog(plkfile,phe,altmod.c,nullmod.c,cores)


######################################
# Create separate R script for below

#####
# Run GWAS Models in AA GAIN Bipolar

setwd("/data5/bsi/psych/s105109.goa/Bipolar_WNT_Pathway/Analysis/AA/GWAS")


# Read in relevant phenotype, covariate, and subject identifier data
covars <- read.table("AA_GAIN_pheno.txt",header=TRUE)[,c(2,4,5)] # example covariate/phenotype file with columns containing BMI and ancestry PCA

names(covars) <- c("IID","BMI","ANC") # <--- Any covariates can be used as long as subject ID column name matches across all files for easy merging

phen <- read.table("AA.GAIN.filtered.fam",header=FALSE)[,c(2,5,6)] # <--- Read in sex and phenotype information from Plink format .fam file
names(phen) <- c("IID","SEX","BP")
phe <- merge(phen,covars,by="IID") # <--- subject ID column to merge files and also exclude subjects without relevant variables
phe$BP[phe$BP==-9] <- NA # assign missing phenotype values of NA (Plink has these coded -9 by default)
phe$BP[phe$BP==1] <- 0 #
phe$BP[phe$BP==2] <- 1 # 


plkfile <- "AA.GAIN.filtered"# file name that comes before binary plink extensions (exclude extensions)

altmod.intcv <- c("BP ~  ANC + SEX + BMI + SNP + SNP*ANC + SNP*SEX + BMI*ANC + BMI*SEX + BMI*SNP")
nullmod.intcv <- c("BP ~  ANC + SEX + BMI + SNP + SNP*ANC + SNP*SEX + BMI*ANC + BMI*SEX")

pheno.type <-c("case.control")
cores <- 24
altmod <- altmod.intcv
nullmod <- nullmod.intcv

gwas.mod.intxcv <- GWAS.frame(plkfile,phe,altmod,nullmod,pheno.type,cores)
save(gwas.mod.intxcv, file = "AAgwas.mod.intxcv.RData")
#
altmod.int <- c("BP ~  ANC + SEX + BMI + SNP + BMI*SNP")
nullmod.int <- c("BP ~  ANC + SEX + BMI + SNP")
altmod <- altmod.int
nullmod <- nullmod.int

gwas.mod.int <- GWAS.frame(plkfile,phe,altmod,nullmod,pheno.type,cores)
save(gwas.mod.intx, file = "AAgwas.mod.intx.RData")
#