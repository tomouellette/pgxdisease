# ALLOPURINOL RESPONSE
######################################################################################################################################################################
######################################################################################################################################################################
# Given the Wen Allopurinol response GWAS and Kottgen Serum Urate GWAS - perform GPA joint/multiple association analysis
# at rs2231142 --- 4:89052323:G:T

# Allopurinol_PGxAllopurinolResponse.tsv = Allopurinol response GWAS
# Allopurinol_SerumUrate-GCST001791.csv = Kottgen Serum Urate

#Allopurinol response can be downloaded from https://www.pgrn.org/riken-gwas-statistics.html
#Serum urate summary statistics can be downloaded at the EBI-GWAS catalog under the ID GCST001791

setwd("~/Desktop/NATUREPGX/Data/APR")

library(data.table)
library(dplyr)

# Read in AINOD Data

APR_response <- fread("Allopurinol_PGxAllopurinolResponse.tsv")
APR_EBI_urate <- fread("Allopurinol_SerumUrate-GCST001791.csv")

# Select relevant columns for curation

APR_response <- APR_response[, c("rsid", "pval", "beta", "se", "Effect_Allele(A1)", "Other_Allele(A2)")]
APR_EBI_urate <- APR_EBI_urate[, c("MarkerName", "p_gc", "beta", "se", "A1", "A2")]

colnames(APR_response) <- c("rsid", "pval", "beta", "se", "effectAllele", "otherAllele")
colnames(APR_EBI_urate) <- c("rsid", "pval", "beta", "se", "effectAllele", "otherAllele")

#Analyze GPA
library(GPA)
C <- merge(APR_response, APR_EBI_urate, by = c("rsid"))
pmat <- C[,c("rsid", "pval.x", "pval.y")]
colnames(pmat) <- c("rsid", "P_response", "P_urate")
library(tibble)
pmat <- column_to_rownames(pmat, "rsid")
pmat <- as.matrix(pmat)
pmat <- na.omit(pmat)

fit.GPA <- GPA(pmat[,1:2])
fit.GPA.H0 <- GPA(pmat[,1:2], pleiotropyH0 = TRUE)

test_p <- pTest(fit.GPA, fit.GPA.H0)
fdrs <- fdr(fit.GPA, pattern="11")
tablez <- data.table(rsid = rownames(pmat), fdr = fdrs)
