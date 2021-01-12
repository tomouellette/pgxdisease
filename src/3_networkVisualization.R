# NETWORK PLOT SCRIPT USING GGNET2
# ================================
# NOTE: THE DATA FILE MUST BE IN TIDY FORMAT WITH THE COLULMS CONTAINING SNP, ASSOCIATED PHENOTYPE, PHENOTYPE CLASS/ABBREVIATION
# -------------------------------------------------------------------------------------------------------------------------------
#Network plot
# args[1] is working directory for tidy data with quotes
# args[2] is the output directory with quotes
# args[3] is name of network plot output with quotes

library(dplyr)
library(readr)
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(tidyselect)

#Read in tidy data
args <- commandArgs(trailing = TRUE)
phe <- data.table::fread(args[1])
phe$abbreviation <- sapply(strsplit(phe$specific_phenotype, " "), '[', 1)
phe$abbreviation <- as.character(phe$abbreviation)

#Build Trait <-> Trait Connections
select <- phe %>% select(phenotype, abbreviation, variant) %>%
  inner_join(., select(., abbreviation, variant), by = "variant") %>%
  rename(sp1 = abbreviation.x, sp2 = abbreviation.y) %>%
  filter(sp1 != sp2)

#Build Edgelist
edgelist <- data.frame(select(select, sp1, sp2), check.names = FALSE)

#Remove reversed duplicates
edgelist <- edgelist[!duplicated(apply(edgelist,1,function(x) paste(sort(x),collapse='_'))),]

#Initialize network
net = network(edgelist, directed = FALSE)

#Match primary phenotype
x = data.frame(abbreviation = network.vertex.names(net))
x = merge(x, distinct(phe, abbreviation, .keep_all = TRUE), by="abbreviation")
x = x$phenotype
net %v% "ICD10" = as.character(x)

#Set edge attribute
set.edge.attribute(net, 
                   "color", 
                   ifelse(net %e% "weights" > 1, 
                          "black", 
                          "grey75"))
# TO RANDOMLY SAMPLE COLOURS USE THE CODE BELOW (To sample 20 of available R colors used code that is blocked out)
#y = sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)], 20)
#Alternatively, check out http://phrogz.net/css/distinct-colors.html - algorithmically generates
#visually distinct array of colours

y = c(rgb(191,0,102, max=255), 
      rgb(242,182,230, max=255), 
      rgb(255,0,238, max=255), 
      rgb(77,57,75, max=255), 
      rgb(86,48,191, max=255), 
      rgb(19,27,77, max=255), 
      rgb(89,125,179, max=255), 
      rgb(96,172,191, max=255), 
      rgb(64,255,242, max=255), 
      rgb(19,77,73, max=255), 
      rgb(163,217,170, max=255), 
      rgb(29,217,0, max=255), 
      rgb(97,102,26, max=255), 
      rgb(255,242,64, max=255), 
      rgb(229,122,0, max=255), 
      rgb(191,147,96, max=255), 
      rgb(76,10,0, max=255), 
      rgb(255,145,128, max=255), 
      rgb(255,0,0, max=255), 
      rgb(204,0,0, max=255))

names(y) = levels(factor(x))

networkplot <- ggnet2(net,
       color = "ICD10",
       palette = y,
       alpha = 0.8, 
       edge.alpha = 0.5,
       edge.color = c("color", "black"),
       label.size = 5,
       legend.size = 6,
       max_size = 10,
       label.alpha = 0.6,
       layout.par = list(cell.jitter = 0.75)
       ) + 
      geom_text(
        label = network.vertex.names(net), 
        #position=position_jitter(width=0.01,height=0.01), 
        fontface = "bold",
        alpha = 1)

setwd(args[2])
ggsave(filename = args[3], plot = networkplot, height = 14, width = 20, units = "in")