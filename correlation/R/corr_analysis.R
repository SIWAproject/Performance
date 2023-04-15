library(phyloseq)
library(dplyr)

source("/Users/sebastianbedoyamazo/Documents/siwa_git/Methods-review/functions.R")
phyl_dir <- "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/performance_df/pq_obj/"

#EAFIT surveillance tables (same for all reports)
species_taxonomy_info <-read.csv("/Users/sebastianbedoyamazo/Documents/siwa_git/Reports/Version1.0/Input_data/species_metabolic_effects.csv",check.names = FALSE, sep=";")
genera_taxonomy_info <- read.csv( "/Users/sebastianbedoyamazo/Documents/siwa_git/Reports/Version1.0/Input_data/genus_metabolic_effects.csv", check.names = FALSE, sep=";")

#put "species_taxonomy_info" or "genera_taxonomy_info" based on analysis
kgbg = genera_taxonomy_info
kgbsp = species_taxonomy_info


# correlation tables
Broiler_g <- readRDS(paste0(phyl_dir, paste0("CorrTableGenus_broiler.rds")))
Broiler_sp <- readRDS(paste0(phyl_dir, paste0("CorrTableSp_broiler.rds")))

# correlation tables by taxa
Broiler_gg <- readRDS(paste0(phyl_dir, paste0("CorrTableGnGn_broiler.rds")))
Broiler_spsp <- readRDS(paste0(phyl_dir, paste0("CorrTableSpSp_broiler.rds")))




colnames(kgbsp)
kgbsp <- kgbsp[, c("Species", "Probiotic_potential")]
colnames(kgbsp)[colnames(kgbsp)=="Species"] <- "Taxa"
kgbsp$Taxa <- gsub(" ", "_", kgbsp$Taxa)



#merge corr with KGB 
Broiler_sp <- merge(Broiler_sp, kgbsp, by = "Taxa", all.x = TRUE)



colnames(kgbg)
kgbg <- kgbg[, c("Genus", "Probiotic_potential")]
colnames(kgbg)[colnames(kgbg)=="Genus"] <- "Taxa"
kgbg$Taxa <- gsub(" ", "_", kgbg$Taxa)

#merge corr with KGB 
Broiler_g <- merge(Broiler_g, kgbg, by = "Taxa", all.x = TRUE)



#filter by p-value
Broiler_sp_filter <-  Broiler_sp[Broiler_sp$AdjPvalue < 0.05, ]
unique(Broiler_sp_filter$Taxa)

Broiler_g_filter <-  Broiler_g[Broiler_g$AdjPvalue < 0.05, ]
unique(Broiler_g_filter$Taxa)



# #filter by repeated name
# kgb_corr_filter <- kgb_corr_filter %>%
#   group_by(Taxa) %>%
#   filter(n() > 1)
# unique(kgb_corr_filter$Taxa)

broiler_score <- Broiler_g_filter
  
broiler_score <- broiler_score %>% 
  mutate(Score = case_when(
    grepl("FCR", Env) & Correlation > 0 ~ "negative",
    grepl("FCR", Env) & Correlation <= 0 ~ "positive",
    grepl("BW", Env) & Correlation > 0 ~ "positive",
    grepl("BW", Env) & Correlation <= 0 ~ "negative",
    TRUE ~ "negative"
))


#ages list
colnames(broiler_score)
list_taxas <- unique(broiler_score$Taxa)
list_samples <- unique(broiler_score$Type)

#empty list to store the correlation results for each age group
taxa_list <- list()

# Iterate over each age group and perform the same operations
for (i in list_taxas) {
  for (j in list_samples) {
    # Subconjunto del dataframe que cumple las condiciones
    subset_taxa <- subset(broiler_score, Taxa == i & Type == j)
    
    # Verificar si hay filas en el subconjunto antes de agregar la nueva columna
    if (nrow(subset_taxa) > 0) {
      subset_taxa$Potential <- ifelse(length(unique(subset_taxa$Score)) == 1, "consistent", "no consistent")
      taxa_list <- rbind(taxa_list, subset_taxa)
    }
  }
}




taxa_list_sp = taxa_list  
taxa_list_genus = taxa_list 

taxa_list_sp$Taxonomy = 'Species'  
taxa_list_genus$Taxonomy = 'Genus' 


unique(taxa_list_sp$ProjectID)
unique(taxa_list_genus$ProjectID)
unique(taxa_list_genus$Taxa)
unique(taxa_list_sp$Taxa)
list_taxa <- rbind(taxa_list_sp, taxa_list_genus)


f = paste0(phyl_dir, paste0("CorrTableFilter_", exp,".rds"))
f
saveRDS(list_taxa, file=f)


# 
# list_positive_genus = subset(taxa_list_genus, Score == "positive" & Potential == "consistent")
# list_negative_genus = subset(taxa_list_genus, Score == "negative" & Potential == "consistent")
# list_positive_sp = subset(taxa_list_sp, Score == "positive" & Potential == "consistent")
# list_negative_sp = subset(taxa_list_genus, Score == "negative" & Potential == "consistent")
# 
# 
# list_positive <- rbind(list_positive_genus, list_positive_sp)
# list_negative <- rbind(list_negative_genus, list_negative_sp)
# 
# unique_list_positive <- unique(list_positive$Taxa)
# unique_list_negative <- unique(list_negative$Taxa)






# mask boolean for filtrar los valores de df1
mask <- Broiler_gg$Taxa %in% Broiler_g$Taxa & Broiler_gg$Env %in% Broiler_g$Env & Broiler_gg$Type %in% Broiler_g$Type & Broiler_gg$Age %in% Broiler_g$Age & Broiler_gg$ProjectID %in% Broiler_g$ProjectID

# Filtramos los valores de df1 utilizando la máscara booleana y la función subset()
df_filtered <- subset(Broiler_gg, mask)

