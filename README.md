# <span style="color:red">Code will be added following manuscript publication/initial review.</span>

# Integrating disease and drug-related phenotypes for improved identification of pharmacogenomic variants

## Data

The required data to replicate statistically significant findings for visualization and colocalization analysis can be retrieved from:

**Pharmacogenomics Research Network**

```
wget https://www.pgrn.org/uploads/1/0/7/8/107807723/wenc_eur_sumstats_a2.txt_1.zip 
wget https://www.pgrn.org/uploads/1/0/7/8/107807723/wheeler_et_al_ccr_2017_sum_stats.txt_5.zip 
wget https://www.pgrn.org/uploads/1/0/7/8/107807723/meta_w_h_b_n1.tbl_1.zip
wget https://www.pgrn.org/uploads/1/0/7/8/107807723/adv_adno_assoc.assoc.txt.zip 
```

[ Files in order of (1) allopurinol response, (2) cisplatin ototoxicity, (3) anti-hypertensive induced new-onset diabetes, and (4) celecoxib prevention of colorectal adenoma ]

**UK Biobank**

```
wget https://www.dropbox.com/s/k5o6xn6uw1hvwpm/K50.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O K50.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://www.dropbox.com/s/x060zthr9agkp2h/J33.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O J33.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://www.dropbox.com/s/p697in77wsizlvz/R31.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O R31.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://www.dropbox.com/s/7l26kh3kfduu7gh/M10.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O M10.gwas.imputed_v3.both_sexes.tsv.bgz
```

[ Files in order of (1) Crohn's disease (K50), (2) nasal polyps (J33), (3) unspecified hematuria (R31), and (4) ]

**EBI-GWAS Catalog**

```
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/KottgenA_23263486_GCST001791/GUGC_MetaAnalysis_Results_UA.csv.zip
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/KottgenA_23263486_GCST001790/GUGC_MetaAnalysis_Results_Gout.csv.zip
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/deLangeKM_28067908_GCST004131/ibd_build37_59957_20161107.txt.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/deLangeKM_28067908_GCST004132/cd_build37_40266_20161107.txt.gz
```

[ Files in order of (1) serum urate _Kottgen_, (2) gout _Kottgen_, (3) inflammatory bowel disease _de Lange_, and (4) Crohn's disease _de Lange_ ]


## Visualizing regional association plots and performing colocalization with MultiZoom.R

**Step One** 

- All GWAS Summary statistics must have SNP, PVAL, CHR, POS, and N columns. 

**Step Two** 

- An information file with columns FILENAME, STUDY, FILTER, CHR, POS_MIN, POS_MAX, TYPE, SDY, N_CASE, S
  - FILENAME is the name of the GWAS summary statistic file
  - STUDY is a simpler name for the study
  - FILTER is used to group specific GWAS together for analysis - base this off shared SNPS from comparison of all GWAS
  - CHR is the chromosome you want to analyze colocalization on
  - POS_MIN and POS_MAX are the range you want to examine colocalization on - we recommend choosing +/-100kb the top shared SNP between a group of traits
  - TYPE *must* be either "quant" or "cc" to indicate study type; required for colocalization
  - SDY is the population standard deviation of the quantitative trait; NA if trait == "cc" or if unavailable
  - N_CASE is the number of cases for binary/case-control traits; NA if trait == "quant"
  - S is the proportion of samples that are cases (i.e. N_CASES / Total N); required if type =="cc" else put NA for quant

**Step Three**

- Run MultiZoom on the command-line the following way

```
Rscript MultiZoom.R [~/Directory/of/GWAS/Data] [~/Link/to/INFO/File] [hg19 or hg38] [~/Output/Directory]

Example:

Rscript MultiZoom.R "~/Desktop" "~/Desktop/INFOFILE.tsv" hg19 "~/Desktop/OutputFolder"
```

Having trouble? Please contact t.ouellette@mail.utoronto.ca
