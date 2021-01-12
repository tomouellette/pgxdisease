# A SCRIPT TO IDENTIFY WHICH SNPS ARE ASSOCIATED WITH WHICH TRAITS
# ================================================================
# NOTE: THIS SCRIPT CAN BE USED DIRECTLY ON FILTERED OR UNFILTERED GWAS DATASETS BY MODIFYING P-VALUE THRESHOLD ON LINE 46
# ---------------------------------------------------------------------------------------------------------------------
#### This code is separated into 3 sections.
## Section 1. Compiling Variant Counts across all GWAS Summary Statistics
## Section 2. Filtering for Shared Variants (Count >= 2)
## Section 3. Add Complete Column Names and Write Table for Shared Variants

### Command Line Arguments
#args[1] = directory where GWAS summary statistics are located with quotes
#args[2] = name of output directory with quotes
#args[2] = name of output file with quotes

#Example command line: Rscript 2/_identifySharedVariants.R "~/Desktop/GWAS_data" "~/Desktop/Folder" "Shared_Variants.tsv"

#################### LOAD LIBRARIES ####################

library(dplyr)
library(data.table)
library(tidyverse)

#################### SECTION 1 ####################

#Argument for GWAS Summary Statistics Directory
args <- commandArgs(trailing = TRUE)

setwd(args[1])

#For loop to read in all GWAS Summary Statistics

pb = txtProgressBar(min = 0, max = length(list.files()), initial = 0) 

for (i in 1:length(list.files())) {
  setTxtProgressBar(pb,i)
  temp <- list.files()[i]
  assign(temp, filter(data.table::fread(list.files()[i]))) 
}

#For loop to build list of filtered variants
#Caution: for eval(parse) to work correctly, ensure GWAS text file names have no special characters e.g no _, -, etc.

variant_list <- list()
library(dplyr)
for (i in 1:length(list.files())) {
  variants <- eval(parse(text = list.files()[i]))
  variants <- variants %>% filter(PVAL < 1e-05)
  variants <- as.list(variants$VARIANT)
  variant_list[[i]] <- variants
}

#Tabulate all shared variants
table <- addmargins(
  table(gene=unlist(variant_list), vec=rep(list.files(), times=sapply(variant_list,length))),
  2,
  list(Count=function(x) sum(x[x>0])))

#################### SECTION 2 ####################

##Convert table from matrix type to data frame
t2 <- as.data.frame.matrix(table)

##Move rownames to a column
t2 <- rownames_to_column(t2)

##Filter for counts greater or equal to 2
t2 <- t2 %>% filter(Count >= 2)

#################### SECTION 3 ####################

###Add complete ICD10 column names

###Writing a table with second argument as name
setwd(args[2])
write.table(t2, file = args[3], sep = "\t", quote = FALSE, row.names = FALSE)

















