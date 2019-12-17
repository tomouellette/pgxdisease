
##### Effect Size Comparisons for Allopurinol Response, Gout, and Serum Urate Concentrations
############################################################################################
# rs10011796 (FDR = 0.00048), 
# rs3114020 (FDR = 0.0057), 
# rs2231142 (FDR = 0.0063), 
# rs21999936 (FDR = 0.0096),
# rs4148155 

setwd("~/Desktop/NATUREPGX/DATA/APR")
library(dplyr)
library(data.table)

# Standardize columns names
gout <- data.table::fread(list.files()[1])
urate <- data.table::fread(list.files()[4])
apr <- data.table::fread(list.files()[3])
apr <- dplyr::select(apr, rsid, N, 'Effect_Allele(A1)', 'Other_Allele(A2)', beta, se, pval)
colnames(apr) <- colnames(gout)

# Standardize UKBB
ukbb <- data.table::fread(list.files()[2])
ukbb <- ukbb %>% filter(low_confidence_variant == FALSE)
variants <- data.table::fread("variants.tsv")
variants <- variants[,c("variant", "rsid")]
ukbb <- merge(ukbb, variants, by = "variant")
ukbb <- dplyr::select(ukbb, rsid, n_complete_samples, minor_allele, minor_AF, beta, se, pval)
colnames(ukbb) <- colnames(gout)

# Convert case-controls effect sizes to odds ratios
gout$OR <- exp(gout$beta)
gout$upper95 <- exp(gout$beta + gout$se)
gout$lower95 <- exp(gout$beta - gout$se)

ukbb$OR <- exp(ukbb$beta)
ukbb$upper95 <- exp(ukbb$beta + ukbb$se)
ukbb$lower95 <- exp(ukbb$beta - ukbb$se)

#Normalize quantitative effects based on min-max feature scaling
# Equation: X' = a + (X - min(x))*(b - a) / (max(x) - min(x)) where [a, b] is the new range of the distribution

urate$normB <- urate$beta/max(urate$beta)
urate$normSE <- urate$se/max(urate$beta)
urate$upperSE <- urate$normB + urate$normSE
urate$lowerSE <- urate$normB - urate$normSE

apr$normB <- apr$beta/max(apr$beta)
apr$normSE <- apr$se/max(apr$beta)
apr$upperSE <- apr$normB + apr$normSE
apr$lowerSE <- apr$normB - apr$normSE

# Filter for rsids
#markers <- c("rs10011796", "rs3114020", "rs2231142", "rs2199936", "rs4148155")
markers <- c("rs10011796", "rs2231142")
gout <- gout %>% filter(MarkerName %in% markers)
urate <- urate %>% filter(MarkerName %in% markers)
apr <- apr %>% filter(MarkerName %in% markers)
ukbb <- ukbb %>% filter(MarkerName %in% markers)

gout <- gout[order(gout$MarkerName),]
urate <- urate[order(urate$MarkerName),]
apr <- apr[order(apr$MarkerName),]
ukbb <- ukbb[order(ukbb$MarkerName),]

# Format for plotting

library(ggplot2)

urate_p <- data.table::melt(urate, id.vars = c("MarkerName", "n_total", "A1", "A2", "p_gc", "beta", "se", "normSE"))
apr_p <- data.table::melt(apr, id.vars = c("MarkerName", "n_total", "A1", "A2", "p_gc", "beta", "se", "normSE"))

library(wesanderson)
urate_p$group <- "Serum Urate (Kottgen et al. 2012)"
apr_p$group <- "Allopurinol Response"
boxing <- rbind(urate_p, apr_p)
A <- ggplot(boxing, aes(x = MarkerName, y = value)) + 
  geom_boxplot(aes(fill=group)) +
  theme_light()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 11))+
  scale_fill_manual(values = wes_palette("Darjeeling1")[c(1, 5)]) +
  labs(x = "", y = "Normalized Effect Size (+/- StdErr) \n") +
  ylim(0, 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted")

gout_p <- data.table::melt(gout, id.vars = c("MarkerName", "n_total", "A1", "A2", "p_gc", "beta", "se"))
gout_p$group <- "Gout (Kottgen et al. 2012)"

B <- ggplot(gout_p, aes(x = MarkerName, y = value)) + 
  geom_boxplot(aes(fill=group)) +
  theme_light()+
  scale_fill_manual(values = wes_palette("Darjeeling1")[c(2,3)]) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 11)) +
  xlab("") +
  ylab("Odds Ratio (+/- StdErr)\n") +
  scale_y_continuous(position = "right", limits = c(0, 2)) +
  geom_hline(yintercept = 0, linetype = "dotted")

library(ggpubr)
ggarrange(A, B, legend = "bottom",  widths=c(2, 1))

##### Effect Size Comparisons for AINOD, Crohn's, and IBD
############################################################################################

setwd("~/Desktop/NATUREPGX/DATA/AINOD")

crohns <- data.table::fread("AINOD_crohns_GCST004132.txt")
ibd <- data.table::fread("AINOD_ibd_GCST004131.txt")
ainod <- data.table::fread("AINOD_MetaDiabetesAntihypertensives.tsv")

crohns <- crohns[,c("MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P.value")]
ibd <- ibd[,c("MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P.value")]
ainod <- ainod[,c("MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P-value")]

# Format crohns and ibd to merge with variants
crohns$MarkerName <- gsub("_", ":", crohns$MarkerName)
ibd$MarkerName <- gsub("_", ":", ibd$MarkerName)

colnames(crohns)[1] <- "variant"
colnames(ibd)[1] <- "variant"

crohns <- merge(crohns, variants, by = "variant")
ibd <- merge(ibd, variants, by = "variant")
colnames(crohns)[7] <- "MarkerName"
colnames(ibd)[7] <- "MarkerName"

crohns <- dplyr::select(crohns, MarkerName, Allele1, Allele2, Effect, StdErr, P.value)
ibd <- dplyr::select(ibd, MarkerName, Allele1, Allele2, Effect, StdErr, P.value)

# Match ainod column names to crohns and ibd

colnames(ainod)[6] <- "P.value"

# Select for rsids of interest
## SLC9A4: rs17027255, rs17027258, rs56043441 and IL18RAP: rs11123929, rs56044378, rs60038017

markers <- c("rs17027255", "rs17027258", "rs56043441", "rs11123929", "rs56044378", "rs60038017")
crohns <- crohns[MarkerName %in% markers]
ibd <- ibd[MarkerName %in% markers]
ainod <- ainod[MarkerName %in% markers]

crohns <- crohns[order(crohns$MarkerName),]
ibd <- ibd[order(ibd$MarkerName),]
ainod <- ainod[order(ainod$MarkerName),]

# Add group labels

crohns$group <- "Crohn's (de Lange et al. 2017)"
ibd$group <- "Inflammatory Bowel Disease (de Lange et al. 2017)"
ainod$group <- "Antihypertensive-Induced New Onset Diabetes"

ainod_effects <- rbind(crohns, ibd, ainod)
ainod_effects$OR <- ifelse(ainod_effects$group == "Antihypertensive-Induced New Onset Diabetes", 1/(exp(as.numeric(ainod_effects$Effect))), exp(as.numeric(ainod_effects$Effect))) 
ainod_effects$upperSE <- ifelse(ainod_effects$group == "Antihypertensive-Induced New Onset Diabetes", 1/exp(as.numeric(ainod_effects$Effect+ainod_effects$StdErr)), exp(as.numeric(ainod_effects$Effect+ainod_effects$StdErr)))
ainod_effects$lowerSE <- ifelse(ainod_effects$group == "Antihypertensive-Induced New Onset Diabetes", 1/exp(as.numeric(ainod_effects$Effect-ainod_effects$StdErr)), exp(as.numeric(ainod_effects$Effect-ainod_effects$StdErr)))
ainod_effects$direction <- ifelse(ainod_effects$OR >= 1, as.numeric(0.1), as.numeric(-0.1))
ainod_effects$significance <- ifelse(ainod_effects$P.value >= 1e-08, as.numeric(-0.1), as.numeric(0.1))
ainod_effects$gene <- ifelse(ainod_effects$MarkerName %in% c("rs17027255", "rs17027258", "rs56043441"), "SLC9A4", "IL18RAP")


SLC9A4 <- ggplot(ainod_effects[ainod_effects$gene == "SLC9A4"], aes(x = MarkerName, y = direction)) +
  geom_hline(yintercept = c(-0.1, 0.1), alpha = 0.1)+
  geom_boxplot(aes(color = group), size = 4, ) +
  coord_flip()+
  theme_light()+
  scale_color_manual(values = wes_palette("Darjeeling1")[c(2,3, 5)]) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(color = "white"),
        axis.title.x = element_text(color = "white")) +
  geom_hline(yintercept = 0, linetype = "dotted")+
  ylim(-0.15, 0.15) +
  labs(x = "SLC9A4\n", y = "\nEffect Direction") +
  scale_x_discrete(position = "left")

IL18RAP <- ggplot(ainod_effects[ainod_effects$gene == "IL18RAP"], aes(x = MarkerName, y = direction)) +
  geom_hline(yintercept = c(-0.1, 0.1), alpha = 0.1) +
  geom_boxplot(aes(color = group), size = 4) +
  coord_flip()+
  theme_light()+
  scale_color_manual(values = wes_palette("Darjeeling1")[c(2,3, 5)]) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(color = "white"),
        axis.title.x = element_text(color = "white")) +
  geom_hline(yintercept = 0, linetype = "dotted")+
  ylim(-0.15, 0.15) +
  labs(x = "IL18RAP\n", y = "\nEffect Direction")+
  scale_x_discrete(position = "left")

ggarrange(SLC9A4, IL18RAP, common.legend = TRUE, nrow = 2, ncol = 1, align = "v")
