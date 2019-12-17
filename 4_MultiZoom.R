#########################################################
# MULTI-ZOOM: Integrated Visualization and Colocalization
# - Requires pre-formatted GWAS data described below
# - Requires information file described below
# - Pre-subsetting at +/- 550kb around causative SNP will enable faster run times
#################################################################################

# PRE-FORMATTED GWAS DATA
# -----------------------
# Minimum GWAS Summary Statistic Labels Required: [ CHR, POS, PVAL, MAF, N ]
# Recommended GWAS Summary Statistic Labels: [ CHR, POS, PVAL, MAF, N, BETA, SE ]
# Note: N can be the same for all SNPs if only study size is available

#################################################################################

# INFORMATION FILE
# -----------------------
# Requires text file (e.g. tsv, txt, etc) with following columns:
# FILENAME, STUDY, FILTER, CHR, POS_MIN, POS_MAX, TYPE, SDY, N_CASE, S

# Notes
# -----
# FILENAME is the actual name of GWAS text file in directory containing GWAS
# STUDY is a name you would like in visualization
# FILTER can be used to produce multiple multizoom, each FILTER group will be compared
# CHR is the chromosome you would like to evaluate for colocalization
# POS_MIN and POS_MAX are the maximum and minimum positions you would like to analyze for colocalization
# TYPE requires either an input of cc for case-control or quant for QTL study
# SDY is needed if TYPE == quant and is the population standard deviation for that trait (if unknown enter NA)
# N_CASE or S is needed for TYPE == cc; if N varies across each SNP in GWAS file then provide N_CASE
# - if N doesn't vary across SNPs in gwas then S = N_CASE / (N_CASE + N_CONTROL)

#################################################################################

#Usage: Rscript "~/Directory/of/GWAS/files" "~/InformationFile.textfile" [genome build: hg19 or hg38] "~/Directory/For/Plots/Tables/"
#Example: Rscript "~/Desktop/GWAS" "~/Documents/INFO.tsv" hg19 "~/Desktop/Output_Folder/"

#################################################################################

args <- commandArgs(trailingOnly = TRUE)
print(args[1])
print(args[2])
print(args[3])
# 0. Load frequently used libraries and information file
# ------------------------------------------------------

library(data.table)
library(dplyr)

INFO <- data.table::fread(args[2])

# 1. Use the RACER package to obtain gene track data from hg19 or hg38
# --------------------------------------------------------------------

library(RACER)
if (args[3] == "hg19") {
  utils::data(hg19)
  genes <- subset(data.table::data.table(hg19), gene_type == "protein_coding")
} else if (args[3] == "hg38") {
  utils::data(hg18)
  genes <- subset(data.table::data.table(hg38), gene_type == "protein_coding")
} else {
  stop("The build selection is incorrect or missing. Add hg19 or hg38 to command and re-run.")
}

genes <- dplyr::select(genes, chromosome, gene_start, gene_end, gene_name)
colnames(genes) <- c("CHR", "START", "END", "GENE")

# 2. Build an association data frame using INFORMATION FILE
# ---------------------------------------------------------

setwd(args[1])

association_data <- data.table()
for (i in 1:length(INFO$FILENAME)) {
  
  #Read in GWAS file
  A <- data.table::fread(INFO$FILENAME[i])
  
  #Get position information from INFO file
  GET_INFO <- INFO[i,]
  
  #Subset GWAS file based on GET_INFO
  
  A <- A %>% filter(CHR == GET_INFO$CHR)
  A <- A %>% subset(POS > (GET_INFO$POS_MIN-400000) & POS < (GET_INFO$POS_MAX+400000))
  
  #Restructure to make large association data frame with all groups
  CHR <- as.numeric(A$CHR)
  POS <- as.numeric(A$POS)
  PVAL <- as.numeric(A$PVAL)
  MAF <- as.numeric(A$MAF)
  N <- as.numeric(A$N)
  GROUP <- rep(GET_INFO$FILENAME, dim(A)[1])
  STUDY <- rep(GET_INFO$STUDY, dim(A)[1])
  FILTER <- rep(GET_INFO$FILTER, dim(A)[1])
  if ("BETA" %in% colnames(A)) {
    BETA <- A$BETA
  } else {
    BETA <- rep(NA, dim(A)[1])
  }
  if ("SE" %in% colnames(A)) {
    SE <- A$SE
  } else {
    SE <- rep(NA, dim(A)[1])
  }
  
  column_bind <- cbind(CHR, POS, PVAL, BETA, SE, N, MAF, GROUP, STUDY, FILTER)
  association_data <- rbind(association_data, column_bind)
}

# 3A. Build MultiZoom function for layered regional association plots
# -------------------------------------------------------------------

############## MULTI-ZOOM FUNCTION #############
################################################
multi_zoom <- function(association_data, 
                       gene_position_data, 
                       chromosome, 
                       plot_start, 
                       plot_end, 
                       range, 
                       ymin, 
                       ymax, 
                       colors) {
  
  # DATA format: 
  # -- association data must contain CHR, POS, P, GROUP columns
  # -- gene_position data must contain CHR, START, END, GENE columns where START and END are base pair positions
  
  ### 1. Load and filter gene position data
  ### =====================================
  ### =====================================
  ### =====================================
  
  library(dplyr)
  library(data.table)
  
  # Read data
  gene_data <- gene_position_data
  gene_data <- dplyr::select(gene_data, CHR, START, END, GENE)
  gene_data$CHR <- as.numeric(as.character(gene_data$CHR))
  gene_data$START <- as.numeric(gene_data$START)
  gene_data$END <- as.numeric(gene_data$END)
  
  # Subset data
  gene_data <- gene_data %>% 
    filter(CHR == chromosome) %>%
    filter(START >= plot_start - range) %>%
    filter(END <= plot_end + range)
  
  # Build pseudo-coordinates for gene blocks using factors
  
  gene_data <- melt(gene_data, id.vars = c("CHR", "GENE"))
  gene_data$pseudo_y_coordinate <- as.numeric(as.factor(gene_data$GENE))
  
  # Generate label coordinates for the plot where, value = x coordinate, pseudo_y_coordinate = y coordinate
  label_genes <- subset(gene_data, gene_data$variable == "END")
  
  ### 2. Load and filter association data
  ### ===================================
  ### ===================================
  ### ===================================
  
  #library <dplyr>
  #library <data.table>
  
  # Read data
  pvalues <- association_data
  pvalues$POS <- as.numeric(pvalues$POS)
  pvalues$NEGLOG10P <- -log10(as.numeric(pvalues$PVAL))
  pvalues$CHR <- as.numeric(as.character(pvalues$CHR))
  
  # Subset data
  pvalues <- pvalues %>% 
    filter(CHR == chromosome) %>%
    filter(POS > plot_start - range) %>%
    filter(POS < plot_end + range)
  pvalues <- data.table(pvalues)
  
  ### 3. Make top hits for shape annotation
  ### ===================================
  ### ===================================
  ### ===================================
  
  tophits <- pvalues[order(pvalues$PVAL),]
  tophits <- distinct(tophits, GROUP, .keep_all = TRUE)
  pvalues$tophit <- "not"
  for (i in 1:dim(tophits)[1]) {
    index <- which(pvalues$POS == tophits$POS[i] & 
                     pvalues$NEGLOG10P == -log10(as.numeric(tophits$PVAL[i])) & 
                     pvalues$GROUP == tophits$GROUP[i])
    pvalues$tophit[index] <- "tophit"
  }
  pvalues <- pvalues[order(pvalues$tophit),]
  
  
  ### 4. Build plots
  ### ===================================
  ### ===================================
  ### ===================================
  
  
  if (dim(gene_data)[1] >= 1) {
        library(ggplot2)
        library(wesanderson)
        gene_plot <- ggplot(gene_data, aes(x = value, y = pseudo_y_coordinate)) +
          geom_line(aes(group = GENE), size = 2) + 
          theme_bw() +
          geom_text(data = label_genes, aes(x = value, y = pseudo_y_coordinate, label = GENE), 
                    hjust = -0.1,
                    vjust = 0.5, 
                    size = 2.5) +
          theme(axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                panel.grid.major.y = element_blank(),
                axis.ticks.y = element_blank(),
                plot.margin = margin(l = 25, r = 0)) + 
          xlab(paste0("Chromosome ", chromosome, " Position (bp)")) +
          coord_cartesian(xlim = c(plot_start - range, plot_end + range), ylim = c(0,(max(as.numeric(gene_data$pseudo_y_coordinate))+1)))
        
        pvalue_plot <- ggplot(pvalues, aes(x = POS, y = NEGLOG10P)) +
          geom_point(aes(color=GROUP,
                         shape = tophit,
                         size = tophit,
                         alpha = tophit)) + 
          theme_bw() + 
          xlab("") +
          ylab("-log10(pval)") +
          scale_color_manual(values = colors) +
          scale_size_manual(values = c(1.5, 8)) +
          scale_shape_manual(values = c("circle", "diamond")) +
          scale_alpha_manual(values = c(0.4, 0.9)) +
          coord_cartesian(xlim = c(plot_start - range, plot_end + range), ylim = c(ymin, ymax)) +
          guides(color = guide_legend("MultiZoom")) +
          guides(shape = FALSE, size = FALSE, fill = FALSE, alpha = FALSE) +
          theme(axis.text.x = element_blank(),
                panel.grid = element_blank(),
                plot.margin = margin(b = -14, l = 0, r = 0, t=14),
                legend.position = "bottom")
        
        ### 4. Combine plots for final visual
        ### ===================================
        ### ===================================
        ### ===================================
        
        library(ggpubr)
        ggarrange(pvalue_plot, gene_plot, heights = c(3,1), nrow = 2, ncol = 1, common.legend = TRUE, legend = "bottom", vjust = 0)
  
    
   
  } else {
    message("There were no genes located within your specified range. Plotting only layered regional association plot.")
    pvalue_plot <- ggplot(pvalues, aes(x = POS, y = NEGLOG10P)) +
      geom_point(aes(color=GROUP,
                     shape = tophit,
                     size = tophit,
                     alpha = tophit)) + 
      theme_bw() + 
      xlab("") +
      ylab("-log10(pval)") +
      scale_color_manual(values = colors) +
      scale_size_manual(values = c(1.5, 8)) +
      scale_shape_manual(values = c("circle", "diamond")) +
      scale_alpha_manual(values = c(0.4, 0.9)) +
      coord_cartesian(xlim = c(plot_start - range, plot_end + range), ylim = c(ymin, ymax)) +
      guides(color = guide_legend("MultiZoom")) +
      guides(shape = FALSE, size = FALSE, fill = TRUE, alpha = FALSE) +
      theme(legend.position = "bottom")+ 
      xlab(paste0("Chromosome ", chromosome, " Position (bp)"))
  
  }
}

# 4. Generate layered regional associations plots
# -----------------------------------------------

plotting_groups <- distinct(INFO, FILTER, .keep_all = TRUE)

for (i in 1:dim(plotting_groups)[1]) {
  
  #Make distinct plotting group
  one_plotting_group <- subset(INFO, FILTER == plotting_groups$FILTER[i])
  
  #Filter association data base on specific filter
  one_association_group <- subset(association_data, FILTER == plotting_groups$FILTER[i])
  one_association_group$CHR <- as.numeric(one_association_group$CHR)
  one_association_group$POS <- as.numeric(one_association_group$POS)
  one_association_group$PVAL <- as.numeric(one_association_group$PVAL)
  one_association_group$BETA <- as.numeric(one_association_group$BETA)
  one_association_group$SE <- as.numeric(one_association_group$SE)
  one_association_group$MAF <- as.numeric(one_association_group$MAF)
  
  #Generate plot
  library(wesanderson)
  layered_regional <- multi_zoom(association_data = one_association_group, 
                                 gene_position_data = genes, 
                                 chromosome = one_association_group$CHR[1], 
                                 plot_start = as.numeric(one_plotting_group$POS_MIN[1]), 
                                 plot_end = as.numeric(one_plotting_group$POS_MAX[1]), 
                                 range = 400000, 
                                 ymin = 0, 
                                 ymax = max(-log10(na.omit(one_association_group$PVAL))), 
                                 colors = c(wes_palette("Darjeeling1", length(subset(INFO, FILTER == plotting_groups$FILTER[i])$FILTER))))
  
  plot_name <- as.character(paste(args[4],"/Plot_for_FILTER_", one_plotting_group$FILTER[1], ".pdf", sep = ""))
  ggsave(plot_name, layered_regional, width = 12, height = 8, device = "pdf")
}

# 5. Generate colocalization results for each FILTER group

library(coloc)

complete_colocalization_table <- data.table::data.table()

for (i in 1:dim(plotting_groups)[1]) {
  
  one_plotting_group <- subset(INFO, FILTER == plotting_groups$FILTER[i])
  one_association_group <- subset(association_data, GROUP %in% one_plotting_group$FILENAME)
  one_association_group$CHR <- as.numeric(one_association_group$CHR)
  one_association_group$POS <- as.numeric(one_association_group$POS)
  one_association_group$PVAL <- as.numeric(one_association_group$PVAL)
  one_association_group$BETA <- as.numeric(one_association_group$BETA)
  one_association_group$SE <- as.numeric(one_association_group$SE)
  one_association_group$MAF <- as.numeric(one_association_group$MAF)
  one_association_group$N <- as.numeric(one_association_group$N)
  
  distinct_groups <- distinct(one_association_group, GROUP)$GROUP
  combos <- data.table::data.table(t(combn(distinct_groups, 2)))
  
  inner_colocalization_table <- data.table::data.table()
  for (a in 1:dim(combos)[1]) {
    
    #Organize and get required data for A
    A <- one_association_group[GROUP == combos[a, 1]]
    A$TYPE <- as.character(INFO[FILENAME == combos[a, 1]]$TYPE)
    A$SDY <- as.numeric(INFO[FILENAME == combos[a, 1]]$SDY)
    A$N_CASE <- as.numeric(INFO[FILENAME == combos[a, 1]]$N_CASE)
    A$S <- as.numeric(INFO[FILENAME == combos[a, 1]]$S)
    
    A <- A %>% subset(CHR == one_plotting_group$CHR[1] & POS > one_plotting_group$POS_MIN[1] & POS < one_plotting_group$POS_MAX[1])
    
    B <- one_association_group[GROUP == combos[a, 2]]
    B$TYPE <- as.character(INFO[FILENAME == combos[a, 2]]$TYPE)
    B$SDY <- as.numeric(INFO[FILENAME == combos[a, 2]]$SDY)
    B$N_CASE <- as.numeric(INFO[FILENAME == combos[a, 2]]$N_CASE)
    B$S <- as.numeric(INFO[FILENAME == combos[a, 2]]$S)
    
    B <- B %>% subset(CHR == one_plotting_group$CHR[1] & POS > one_plotting_group$POS_MIN[1] & POS < one_plotting_group$POS_MAX[1])
    
    #Merge to get common intersect of data
    
    C <- merge(A, B, by = "POS")
    
    #Build lists for colocalization for A
    if (A$TYPE[1] == "cc" | A$TYPE[1] == "quant") {
      if (A$TYPE[1] == "quant") {
        dataset_A <- list(MAF = C$MAF.x,
                          N = C$N.x,
                          pvalues = C$PVAL.x, 
                          type = "quant")
      } else {
        dataset_A <- list(MAF = C$MAF.x,
                          N = C$N.x,
                          pvalues = C$PVAL.x, 
                          type = "cc")
      }
      
      if (is.numeric(C$BETA.x) == TRUE) {
        dataset_A$beta <- C$BETA.x
      } else {
        0
      }
      
      if (is.numeric(C$SE.x) == TRUE) {
        dataset_A$varbeta <- C$SE.x^2
      } else {
        0
      }
      
      if (is.numeric(C$SDY.x) == TRUE) {
        dataset_A$sdY <- C$SDY.x
      } else {
        dataset_A$sdY <- "NA"
      }
      
      if (dim(distinct(C, N.x))[1] > 1  & is.numeric(C$N_CASE.x) == TRUE) {
        dataset_A$s <- C$N_CASE.x/C$N.x
      } else {
        dataset_A$s <- C$S.x
      }
    } else {
      stop("The type in the information file is neither quant nor cc. Please fix before proceeding.")
    } 
    
    # Build list for colocalization for B
    
    if (B$TYPE[1] == "cc" | B$TYPE[1] == "quant") {
      if (B$TYPE[1] == "quant") {
        dataset_B <- list(MAF = C$MAF.y,
                          N = C$N.y,
                          pvalues = C$PVAL.y, 
                          type = "quant")
      } else {
        dataset_B <- list(MAF = C$MAF.y,
                          N = C$N.y,
                          pvalues = C$PVAL.y, 
                          type = "cc")
      }
      
      if (is.numeric(C$BETA.y) == TRUE) {
        dataset_B$beta <- C$BETA.y
      } else {
        0
      }
      
      if (is.numeric(C$SE.y) == TRUE) {
        dataset_B$varbeta <- C$SE.y^2
      } else {
        0
      }
      
      if (is.numeric(C$SDY.y) == TRUE) {
        dataset_B$sdY <- C$SDY.y
      } else {
        dataset_B$sdY <- "NA"
      }
      
      if (dim(distinct(C, N.y))[1] > 1  & is.numeric(C$N_CASE.y) == TRUE) {
        dataset_B$s <- C$N_CASE.y/C$N.y
      } else {
        dataset_B$s <- C$S.y
      }
    } else {
      stop("The type in the information file is neither quant nor cc. Please fix before proceeding.")
    } 
    
    #Perform colocalization
    
    colocalization <- list(coloc.abf(dataset_A, dataset_B))
    coloc_test <- colocalization[[1]]$summary
    
    dt <- data.table::data.table(TRAIT_A = combos[a, 1], 
                                 TRAIT_B = combos[a, 2], 
                                 N_SNPS_TESTED = coloc_test[1][[1]],
                                 PP.H0 = coloc_test[2][[1]],
                                 PP.H4 = coloc_test[6][[1]],
                                 CHR = C$CHR.x[1],
                                 POS_MIN = min(as.numeric(C$POS)),
                                 POS_MAX = max(as.numeric(C$POS)))
    
    inner_colocalization_table <- rbind(inner_colocalization_table, dt)
  }
  
  complete_colocalization_table <- rbind(complete_colocalization_table, inner_colocalization_table)
}

# Write colocalization results into table

write.table(complete_colocalization_table, 
            file = paste(args[4], "/colocalization_results.tsv", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
