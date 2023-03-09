###load libraries panels tabs####
library(phyloseq) #
library(DESeq2) #
library(kableExtra) #
library(genefilter) #
library(microbiome) #
library(ggplot2)  #
library(vegan) #
library(ggpubr) #
library(ggplot2) #
library(plyr)
library(multcompView)

#####Load libraries correlations####
library(ggplot2)
library(phyloseq)
library(stringr)
library(dplyr)

#####Load libraries LR####
library(tidyverse)
library(ape)
library(ggpubr)

####Load libraries categorical histo####
library(kableExtra)
library(plyr)
#library("fantaxtic")
library(data.table)

####load libraries Ratios-histo boxplots####
library(multcompView)
library(reshape)

### OTRICAS
library(plotly)
library(plyr)
library(flexdashboard)
library(shiny)
library(DT)
library(stringr)
library(rstatix)
#library(microViz)
library(lsr)
library(readxl)
library(openxlsx)

### summaries ###
library(psych)

#EAFIT surveillance tables (same for all reports)
species_taxonomy_info <-read.csv("/Users/sebastianbedoyamazo/Documents/siwa_git/Reports/Version1.0/Input_data/species_metabolic_effects.csv",check.names = FALSE, sep=";")
genera_taxonomy_info <- read.csv( "/Users/sebastianbedoyamazo/Documents/siwa_git/Reports/Version1.0/Input_data/genus_metabolic_effects.csv", check.names = FALSE, sep=";")
broad_taxonomy_info <- read.csv("/Users/sebastianbedoyamazo/Documents/siwa_git/Reports/Version1.0/Input_data/broad_groups_metabolic_effects.csv", check.names = FALSE, sep=";")


#open phyloseq object for extract performance data
input_dir = "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/inputdata/"

ODLEPobj <- readRDS(paste0(input_dir, "PhyloseqObject_perf.rds"))


#performances summary
metadata <- as.data.frame(sample_data(ODLEPobj))
meta_exp <- metadata
colnames(meta_exp)


#BOXPLOT by projectid
ggplot(metadata, aes(x = projectid, y = BW42, fill = projectid)) +
  geom_violin()  #geom_boxplot()+  #or geom_violin() +
  #geom_jitter(shape=16, position=position_jitter())
  
ggplot(metadata, aes(x = projectid, y = FCR0_42, fill = projectid)) +
  geom_violin()  #geom_boxplot()+  #or geom_violin() 

ggplot(metadata, aes(x = projectid, y = FCR36_42, fill = projectid)) +
  geom_violin()  #geom_boxplot()+  #or geom_violin() 

ggplot(metadata, aes(x = projectid, y = BWG0_42, fill = projectid)) +
  geom_violin()  #geom_boxplot()+  #or geom_violin() 

ggplot(metadata, aes(x = projectid, y = BWG36_42, fill = projectid)) +
  geom_violin()  #geom_boxplot()+  #or geom_violin() 




#summary by projectid
meta_exp <- meta_exp[, c("projectid", "BW42", "FCR0_42", "FCR36_42", "BWG0_42", "BWG36_42")]
summary_meta_exp <- describe.by(meta_exp, group = meta_exp$projectid)

#tables by projectid
E267 <- summary_meta_exp$E267
E267 <- E267[, c("n", "mean", "sd", "min", "max")]
E267 <- E267[-which(rownames(E267) == "projectid*"),]

E271 <- summary_meta_exp$E271
E271 <- E271[, c("n", "mean", "sd", "min", "max")]
E271 <- E271[-which(rownames(E271) == "projectid*"),]

E345 <- summary_meta_exp$E345
E345 <- E345[, c("n", "mean", "sd", "min", "max")]
E345 <- E345[-which(rownames(E345) == "projectid*"),]


E347 <- summary_meta_exp$E347
E347 <- E347[, c("n", "mean", "sd", "min", "max")]
E347 <- E347[-which(rownames(E347) == "projectid*"),]

options(digits = 2)



#relation between Genera Score and Performance correlation
genera_taxonomy_info <- genera_taxonomy_info[, c("Genus", "Probiotic_potential")]
genera_taxonomy_info <- genera_taxonomy_info %>% rename(Taxa = Genus)

genus_list <- as.list(genera_taxonomy_info$Taxa)
analysis_probiotic_genus <- subset(corr_analysis, Taxa %in% genus_list)


merged_genus_corr <- merge(corr_analysis_probiotic, genera_taxonomy_info, by = "Taxa")


merged_genus_corr_filter <-  merged_genus_corr[merged_genus_corr$AdjPvalue < 0.05, ]

merged_genus_corr_filter <- split(merged_genus_corr_filter, f = ifelse(merged_genus_corr_filter$Correlation >= 0, "positive_cor", "negative_cor"))
merged_genus_corr_filter$negative_cor
merged_genus_corr_filter$positive_cor
merged_genus_corr_filter_neg <- merged_genus_corr_filter$negative_cor
merged_genus_corr_filter_pos <- merged_genus_corr_filter$positive_cor




#relation between Species Score and Performance correlation
species_taxonomy_info <- species_taxonomy_info[, c("Species", "Probiotic_potential")]
species_taxonomy_info <- species_taxonomy_info %>% rename(Taxa = Species)

species_list_corr <- as.list(corr_analysis$Taxa)


species_list_tax <- as.list(species_taxonomy_info$Taxa)
#revisar algunas species 
corr_analysis %>% filter(str_detect(Taxa, "sp."))

#cambiar " " por "_"
species_list_tax <- lapply(species_list_tax, function(x) gsub(" ", "_", x))



analysis_probiotic_species <- subset(corr_analysis, Taxa %in% species_list_tax)
as.list(unique(analysis_probiotic_species$Taxa))


species_taxonomy_info <- species_taxonomy_info %>% mutate(Taxa = gsub(" ", "_", Taxa))
merged_species_corr <- merge(analysis_probiotic_species, species_taxonomy_info, by = "Taxa")
as.list(unique(merged_species_corr$Taxa))

merged_species_corr_filter <-  merged_species_corr[merged_species_corr$AdjPvalue < 0.20, ]


merged_species_corr_filter <- split(merged_species_corr_filter, f = ifelse(merged_species_corr_filter$Correlation >= 0, "positive_cor", "negative_cor"))
merged_species_corr_filter$negative_cor
merged_species_corr_filter$positive_cor



merged_species_corr_filter_neg <- merged_species_corr_filter$negative_cor[, c("Taxa", "Env", "Correlation", "AdjPvalue", "projectid","Probiotic_potential")]
merged_species_corr_filter_pos <- merged_species_corr_filter$positive_cor[, c("Taxa", "Env", "Correlation", "AdjPvalue", "projectid","Probiotic_potential")]




#filter by Species in corr df
corr_analysis_species <- filter(corr_analysis, Taxa_aggre == "species")
corr_analysis_species <-corr_analysis_species[corr_analysis_species$AdjPvalue < 0.05, ]
corr_analysis_species_list <- as.list(unique(corr_analysis_species$Taxa))



