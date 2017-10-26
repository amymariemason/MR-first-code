######
# Author: Amy Mason ; based heavily on code by James Staley  
# Date: Oct 2017
# Goal: Perform Mendelian randomization  
# Inputs: list of variants, data to search for those variants
# Outputs: log file reporting missing variants, dataframe containing the subset of larger data matching the variants given 
######

#####################################################
##### Set-up #####
#####################################################

rm(list=ls())
setwd("C://Users/am2609/Dropbox (Personal)/lpamaster/")
library(MendelianRandomization)

data_rho_master <- read.table("lpa/data/LPA_master_dataset_pcs_EUwinsor_withoutM.txt", header=T, sep="\t", colClasses="character")

ivw_corr <- function(bx, by, byse, rho){
  Omega <- byse%o%byse*rho
  beta <- solve(t(bx)%*%solve(Omega)%*%bx)*t(bx)%*%solve(Omega)%*%by
  se <- sqrt(solve(t(bx)%*%solve(Omega)%*%bx))
  df <- length(bx) - 1
  lci <- beta - qt(df=df, 0.975)*se
  uci <- beta + qt(df=df, 0.975)*se
  p <- 2*(1-pt(abs(beta/se),df))
  results <- data.frame(beta,se,lci,uci,p)
  return(results)
}

## External

lpa <- read.table("lpa/data/LPA_Variants_EUwinsor_withoutM.txt", sep="\t", header=T, colClasses="character")
# lpa <- read.table("lpa/data/LPA_Variants_EUnowinsor_withoutM.txt", sep="\t", header=T, colClasses="character")
# lpa <- read.table("lpa/data/LPA_Variants_EUquantile_withoutM.txt", sep="\t", header=T, colClasses="character")
# lpa <- read.table("lpa/data/LPA_Variants_EUrobust_withoutM.txt", sep="\t", header=T, colClasses="character")
fsteps <- read.table("lpa/data/fstep_snps_0.4_EUwinsor.txt", sep="\t", header=T, colClasses="character")
# fsteps = list(snp = c("rs10455872_g", "rs3798220_t"))
lpa <- lpa[(lpa$variantID %in% fsteps$snp),]
lpa$chr.pos <- paste0("chr", lpa$chr, ":", lpa$pos)
lpa$snp <- lpa$chr.pos; lpa$a1 <- lpa$allele1; lpa$a2 <- lpa$allele2
lpa <- lpa[, c("variantID", "snp", "chr.pos", "chr", "pos" , "a1", "a2", "beta", "se")]
#CHECK WITH STEVE: Assumption is a1 is the "effect allele", and that is the same as "alternate allele" in my data



load('//me-filer1/home$/am2609/My Documents/Programs/Amy 1/Output/20002_1067subset.Rda')
outcome<-output
# list of fields 
# : ref allele: alternate allele
#variant (hg19) [CHROM:POS:REF:ALT]	Variant position as [Chromosome : hg19 Position : reference allele : alternate allele]
#rsid	SNP rsID as provided by UK Biobank
#nCompleteSamples	Number of samples analyzed with non-missing phenotypes
#AC	Dosage allele count from all samples with non-missing phenotypes
#ytx	Dosage alternate allele count (ytx means "y * x" where y=phenotype and x=alternate allele dosage) 
#beta	linear regression beta coefficient
#se	linear regression standard error
#pval	linear regression p-value

# extract alleles
seperates<- read.table(text = as.character(outcome$variant), sep = ":")
names(seperates)<- c("chr","pos","a2","a1")
seperates$chr<-as.character(seperates$chr)
seperates$pos<-as.character(seperates$pos)
seperates$a1<-as.character(seperates$a1)
seperates$a2<-as.character(seperates$a2)
outcome2<-cbind(outcome, seperates)

# attempt to merge with lpa dataset to highlight any alleles in wrong direction
stopifnot(nrow(merge(outcome2, lpa[,c("chr", "pos", "a1","a2")])) != nrow(outcome2))


chd <- read.table("lpa/program_files/mr/phenoscanner/lpa_PhenoScanner_GWAS_winsor.tsv", header=T, sep="\t", stringsAsFactors=F)
# chd <- read.table("lpa/program_files/mr/phenoscanner/clark_phenoscan.csv", header=T, sep=",", stringsAsFactors=F)
chd <- chd[(chd$PMID=="26343387" & chd$Trait=="Coronary artery disease"),]
lpa <- lpa[(lpa$chr.pos %in% chd$Pos..hg19.),]
data_rho <- data_rho_master[(data_rho_master$prev_chd==0), ]
data_rho <- data_rho[, match(lpa$variantID, names(data_rho))]
data_rho <- as.matrix(data_rho)
class(data_rho) <- "numeric"
rho <- cor(data_rho, use="complete.obs")
 

stop() 
 
 
# Set-up
bx_effect <- lpa$a1
by_effect <- chd$Effect.Allele
bx <- as.numeric(lpa$beta)
by <- as.numeric(chd$Beta)
by = ifelse(bx_effect == by_effect, by, -by)
byse <- as.numeric(chd$SE)
rho <- as.matrix(rho)
# Analysis
model_external <- ivw_corr(bx,by,byse,rho)
beta_external = model_external$beta
se_external = model_external$se

bxse <- as.numeric(lpa$se)
mr_external = mr_ivw(mr_input(bx, bxse, by, byse, corr=rho))

mr_external30 = mr_ivw(mr_input(bx[which(lpa$variantID%in%fsteps$snp[1:30])],
                              bxse[which(lpa$variantID%in%fsteps$snp[1:30])],
                                by[which(lpa$variantID%in%fsteps$snp[1:30])], 
                              byse[which(lpa$variantID%in%fsteps$snp[1:30])],
  corr=rho[which(lpa$variantID%in%fsteps$snp[1:30]),which(lpa$variantID%in%fsteps$snp[1:30])])) 

mr_external20 = mr_ivw(mr_input(bx[which(lpa$variantID%in%fsteps$snp[1:20])],
                              bxse[which(lpa$variantID%in%fsteps$snp[1:20])],
                                by[which(lpa$variantID%in%fsteps$snp[1:20])], 
                              byse[which(lpa$variantID%in%fsteps$snp[1:20])],
  corr=rho[which(lpa$variantID%in%fsteps$snp[1:20]),which(lpa$variantID%in%fsteps$snp[1:20])])) 

mr_external10 = mr_ivw(mr_input(bx[which(lpa$variantID%in%fsteps$snp[1:10])],
                              bxse[which(lpa$variantID%in%fsteps$snp[1:10])],
                                by[which(lpa$variantID%in%fsteps$snp[1:10])], 
                              byse[which(lpa$variantID%in%fsteps$snp[1:10])],
  corr=rho[which(lpa$variantID%in%fsteps$snp[1:10]),which(lpa$variantID%in%fsteps$snp[1:10])])) 

mr_external5 = mr_ivw(mr_input(bx[which(lpa$variantID%in%fsteps$snp[1:5])],
                              bxse[which(lpa$variantID%in%fsteps$snp[1:5])],
                                by[which(lpa$variantID%in%fsteps$snp[1:5])], 
                              byse[which(lpa$variantID%in%fsteps$snp[1:5])],
  corr=rho[which(lpa$variantID%in%fsteps$snp[1:5]),which(lpa$variantID%in%fsteps$snp[1:5])])) 


## Internal

lpa <- read.table("lpa/data/LPA_Variants_EUwinsor_withoutM.txt", sep="\t", header=T, colClasses="character")
fsteps <- read.table("lpa/data/fstep_snps_0.4_EUwinsor.txt", sep="\t", header=T, colClasses="character")
# fsteps = list(snp = c("rs10455872_g", "rs3798220_t"), snps = c("rs10455872", "rs3798220"))
lpa <- lpa[(lpa$variantID %in% fsteps$snp),]
lpa$chr.pos <- paste0("chr", lpa$chr, ":", lpa$pos)
lpa$snp <- lpa$chr.pos; lpa$a1 <- lpa$allele1; lpa$a2 <- lpa$allele2
lpa <- lpa[, c("variantID", "snp", "chr.pos", "chr", "pos" , "a1", "a2", "beta", "se")]
chd_internal <- read.table("lpa/data/LPA_Variants_CHD_EUwinsor_withoutM.txt", sep="\t", header=T, colClasses="character")
chd_internal <- chd_internal[(chd_internal$variantID %in% fsteps$snp),]
data_rho <- data_rho_master[(data_rho_master$prev_chd==0), ]
data_rho <- data_rho[, match(lpa$variantID, names(data_rho))]
data_rho <- as.matrix(data_rho)
class(data_rho) <- "numeric"
rho <- cor(data_rho, use="complete.obs")

# Set-up
bx <- as.numeric(lpa$beta)
by <- as.numeric(chd_internal$beta)
byse <- as.numeric(chd_internal$se)
rho <- as.matrix(rho)

model_internal <- ivw_corr(bx,by,byse,rho)
beta_internal = model_internal$beta
se_internal = model_internal$se

bxse <- as.numeric(lpa$se)
mr_internal = mr_ivw(mr_input(bx, bxse, by, byse, corr=rho))

mr_internal30 = mr_ivw(mr_input(bx[which(lpa$variantID%in%fsteps$snp[1:30])],
                              bxse[which(lpa$variantID%in%fsteps$snp[1:30])],
                                by[which(lpa$variantID%in%fsteps$snp[1:30])], 
                              byse[which(lpa$variantID%in%fsteps$snp[1:30])],
  corr=rho[which(lpa$variantID%in%fsteps$snp[1:30]),which(lpa$variantID%in%fsteps$snp[1:30])])) 

mr_internal20 = mr_ivw(mr_input(bx[which(lpa$variantID%in%fsteps$snp[1:20])],
                              bxse[which(lpa$variantID%in%fsteps$snp[1:20])],
                                by[which(lpa$variantID%in%fsteps$snp[1:20])], 
                              byse[which(lpa$variantID%in%fsteps$snp[1:20])],
  corr=rho[which(lpa$variantID%in%fsteps$snp[1:20]),which(lpa$variantID%in%fsteps$snp[1:20])])) 

mr_internal10 = mr_ivw(mr_input(bx[which(lpa$variantID%in%fsteps$snp[1:10])],
                              bxse[which(lpa$variantID%in%fsteps$snp[1:10])],
                                by[which(lpa$variantID%in%fsteps$snp[1:10])], 
                              byse[which(lpa$variantID%in%fsteps$snp[1:10])],
  corr=rho[which(lpa$variantID%in%fsteps$snp[1:10]),which(lpa$variantID%in%fsteps$snp[1:10])])) 

mr_internal5 = mr_ivw(mr_input(bx[which(lpa$variantID%in%fsteps$snp[1:5])],
                              bxse[which(lpa$variantID%in%fsteps$snp[1:5])],
                                by[which(lpa$variantID%in%fsteps$snp[1:5])], 
                              byse[which(lpa$variantID%in%fsteps$snp[1:5])],
  corr=rho[which(lpa$variantID%in%fsteps$snp[1:5]),which(lpa$variantID%in%fsteps$snp[1:5])])) 

###

100-100*exp((mr_internal$Estimate)*-5)
100-100*exp((mr_internal$Estimate)*-10)
100-100*exp((mr_internal$Estimate)*-20)
100-100*exp((mr_internal$Estimate)*-30)
100-100*exp((mr_internal$Estimate)*-50)
100-100*exp((mr_internal$Estimate)*-80)
100-100*exp((mr_internal$Estimate)*-100)
100-100*exp((mr_internal$Estimate)*-120)

100-100*exp((mr_internal$Estimate-1.96*mr_internal$StdError)*-5)
100-100*exp((mr_internal$Estimate-1.96*mr_internal$StdError)*-10)
100-100*exp((mr_internal$Estimate-1.96*mr_internal$StdError)*-20)
100-100*exp((mr_internal$Estimate-1.96*mr_internal$StdError)*-30)
100-100*exp((mr_internal$Estimate-1.96*mr_internal$StdError)*-50)
100-100*exp((mr_internal$Estimate-1.96*mr_internal$StdError)*-80)
100-100*exp((mr_internal$Estimate-1.96*mr_internal$StdError)*-100)
100-100*exp((mr_internal$Estimate-1.96*mr_internal$StdError)*-120)

100-100*exp((mr_internal$Estimate+1.96*mr_internal$StdError)*-5)
100-100*exp((mr_internal$Estimate+1.96*mr_internal$StdError)*-10)
100-100*exp((mr_internal$Estimate+1.96*mr_internal$StdError)*-20)
100-100*exp((mr_internal$Estimate+1.96*mr_internal$StdError)*-30)
100-100*exp((mr_internal$Estimate+1.96*mr_internal$StdError)*-50)
100-100*exp((mr_internal$Estimate+1.96*mr_internal$StdError)*-80)
100-100*exp((mr_internal$Estimate+1.96*mr_internal$StdError)*-100)
100-100*exp((mr_internal$Estimate+1.96*mr_internal$StdError)*-120)

###

exp(-10*mr_internal$Estimate)
exp(-10*mr_internal30$Estimate)
exp(-10*mr_internal20$Estimate)
exp(-10*mr_internal10$Estimate)
exp(-10*mr_internal5$Estimate)

exp(-10*(mr_internal$Estimate+1.96*mr_internal$StdError))
exp(-10*(mr_internal30$Estimate+1.96*mr_internal30$StdError))
exp(-10*(mr_internal20$Estimate+1.96*mr_internal20$StdError))
exp(-10*(mr_internal10$Estimate+1.96*mr_internal10$StdError))
exp(-10*(mr_internal5$Estimate +1.96*mr_internal5$StdError))

exp(-10*(mr_internal$Estimate  -1.96*mr_internal$StdError))
exp(-10*(mr_internal30$Estimate-1.96*mr_internal30$StdError))
exp(-10*(mr_internal20$Estimate-1.96*mr_internal20$StdError))
exp(-10*(mr_internal10$Estimate-1.96*mr_internal10$StdError))
exp(-10*(mr_internal5$Estimate -1.96*mr_internal5$StdError))

###

exp(-10*mr_external$Estimate)
exp(-10*mr_external30$Estimate)
exp(-10*mr_external20$Estimate)
exp(-10*mr_external10$Estimate)
exp(-10*mr_external5$Estimate)

exp(-10*(mr_external$Estimate+1.96*mr_external$StdError))
exp(-10*(mr_external30$Estimate+1.96*mr_external30$StdError))
exp(-10*(mr_external20$Estimate+1.96*mr_external20$StdError))
exp(-10*(mr_external10$Estimate+1.96*mr_external10$StdError))
exp(-10*(mr_external5$Estimate +1.96*mr_external5$StdError))

exp(-10*(mr_external$Estimate  -1.96*mr_external$StdError))
exp(-10*(mr_external30$Estimate-1.96*mr_external30$StdError))
exp(-10*(mr_external20$Estimate-1.96*mr_external20$StdError))
exp(-10*(mr_external10$Estimate-1.96*mr_external10$StdError))
exp(-10*(mr_external5$Estimate -1.96*mr_external5$StdError))


#################


## External CARDIoGRAM (2011)

lpa <- read.table("lpa/data/LPA_Variants_EUwinsor_withoutM.txt", sep="\t", header=T, colClasses="character")
# lpa <- read.table("lpa/data/LPA_Variants_EUnowinsor.txt", sep="\t", header=T, colClasses="character")
# lpa <- read.table("lpa/data/LPA_Variants_EUquantile.txt", sep="\t", header=T, colClasses="character")
fsteps <- read.table("lpa/data/fstep_snps_0.4_EUwinsor.txt", sep="\t", header=T, colClasses="character")
lpa <- lpa[(lpa$variantID %in% fsteps$snp),]
lpa$chr.pos <- paste0("chr", lpa$chr, ":", lpa$pos)
lpa$snp <- lpa$chr.pos; lpa$a1 <- lpa$allele1; lpa$a2 <- lpa$allele2
lpa <- lpa[, c("variantID", "snp", "chr.pos", "chr", "pos" , "a1", "a2", "beta", "se")]
chd <- read.table("lpa/program_files/mr/phenoscanner/lpa_PhenoScanner_GWAS_winsor.tsv", header=T, sep="\t", stringsAsFactors=F)
chd <- chd[(chd$PMID=="21378990" & chd$Trait=="Coronary artery disease"),]
lpa <- lpa[(lpa$chr.pos %in% chd$Pos..hg19.),]
data_rho <- data_rho_master[(data_rho_master$prev_chd==0), ]
data_rho <- data_rho[, match(lpa$variantID, names(data_rho))]
data_rho <- as.matrix(data_rho)
class(data_rho) <- "numeric"
rho <- cor(data_rho, use="complete.obs")
 
# Set-up
bx_effect <- lpa$a1
by_effect <- chd$Effect.Allele
bx <- as.numeric(lpa$beta)
by <- as.numeric(chd$Beta)
by = ifelse(bx_effect == by_effect, by, -by)
byse <- as.numeric(chd$SE)
rho <- as.matrix(rho)
bxse <- as.numeric(lpa$se)
mr_external_2011 = mr_ivw(mr_input(bx, bxse, by, byse, corr=rho))


exp((mr_internal$Estimate)*10)
exp((mr_internal$Estimate-1.96*mr_internal$StdError)*10)
exp((mr_internal$Estimate+1.96*mr_internal$StdError)*10)

exp((mr_external$Estimate)*10)
exp((mr_external$Estimate-1.96*mr_external$StdError)*10)
exp((mr_external$Estimate+1.96*mr_external$StdError)*10)

exp((mr_external_2011$Estimate)*100)
exp((mr_external_2011$Estimate-1.96*mr_external_2011$StdError)*100)
exp((mr_external_2011$Estimate+1.96*mr_external_2011$StdError)*100)

###

exp((mr_external30$Estimate)*10)
exp((mr_external30$Estimate-1.96*mr_external30$StdError)*10)
exp((mr_external30$Estimate+1.96*mr_external30$StdError)*10)

exp((mr_external20$Estimate)*10)
exp((mr_external20$Estimate-1.96*mr_external20$StdError)*10)
exp((mr_external20$Estimate+1.96*mr_external20$StdError)*10)

exp((mr_external10$Estimate)*10)
exp((mr_external10$Estimate-1.96*mr_external10$StdError)*10)
exp((mr_external10$Estimate+1.96*mr_external10$StdError)*10)

exp((mr_external5$Estimate)*10)
exp((mr_external5$Estimate-1.96*mr_external5$StdError)*10)
exp((mr_external5$Estimate+1.96*mr_external5$StdError)*10)

exp((mr_internal30$Estimate)*10)
exp((mr_internal30$Estimate-1.96*mr_internal30$StdError)*10)
exp((mr_internal30$Estimate+1.96*mr_internal30$StdError)*10)

exp((mr_internal20$Estimate)*10)
exp((mr_internal20$Estimate-1.96*mr_internal20$StdError)*10)
exp((mr_internal20$Estimate+1.96*mr_internal20$StdError)*10)

exp((mr_internal10$Estimate)*10)
exp((mr_internal10$Estimate-1.96*mr_internal10$StdError)*10)
exp((mr_internal10$Estimate+1.96*mr_internal10$StdError)*10)

exp((mr_internal5$Estimate)*10)
exp((mr_internal5$Estimate-1.96*mr_internal5$StdError)*10)
exp((mr_internal5$Estimate+1.96*mr_internal5$StdError)*10)