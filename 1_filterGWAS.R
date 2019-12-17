# A LOOP TO FILTER GWAS DATASETS BASED ON A GIVEN P-VALUE THRESHOLD
# =================================================================
# 1. Set working directory to the folder containing your GWAS text files
# 2. Set the output directory for the filtered files on line 17

setwd('YOUR_DIRECTORY_TO_GWAS_SUMMARY_STATISTICS')
Sys.setenv('R_MAX_VSIZE'=32000000000)
library(data.table)
library(dplyr)

pb = txtProgressBar(min = 0, max = length(list.files()), initial = 0) 

for (i in 1:length(list.files())) {
  setTxtProgressBar(pb,i)
  gwas <- data.table::fread(list.files()[i])
  gwas <- gwas %>% filter(pval < 1e-05)
  write.table(gwas, file = paste0('YOUR_OUTPUT_DIRECTORY/', list.files()[i], sep=""), row.names=FALSE, quote= FALSE, sep = "\t") 
  rm(gwas)
  gwas <- 1
}