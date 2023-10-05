## script para hacer una tabla con todas las correlations a nivel de OTU

source("/Users/sebastianbedoyamazo/Documents/siwa_git/Methods-review/functions.R")
source("/Users/sebastianbedoyamazo/Documents/siwa_git/Methods-review/correlations.R")

input_dir <- "/Users/sebastianbedoyamazo/Documents/siwa_git/Pendientes/"
exp <- "Broiler"
ODLEPobj <- readRDS(paste0(input_dir, paste0("PhyloseqObject_", exp,".rds")))

meta_all <- as.data.frame(as.matrix(sample_data(ODLEPobj)))
tax_all <- as.data.frame(as.matrix(tax_table(ODLEPobj)))
phyloseq_rel <-microbiome::transform(ODLEPobj, "compositional") 


corr_df_OTU <- data_frame()
project_ids <- c("E321", "E325", "E326", "E335", "E345", "E347", "E267", "E271", "Fieldale")
# i <- "E345"
for (i in project_ids) {

  
  
  phy_exp <- subset_samples(phyloseq_rel, ProjectID== i)
  phy_exp <- prune_taxa(taxa_sums(phy_exp) > 0, phy_exp)
  meta_exp <- as.data.frame(as.matrix(sample_data(phy_exp)))
  otu_table <- as.matrix(as.data.frame(otu_table(phy_exp)))
  otu_table <- t(otu_table) 
  dim(otu_table)
  ### filtro: que tenga al menos dos OTUs diferentes de cero para que Rmult no falle
  otu_table_filtered <- otu_table[, colSums(otu_table > 0) >= 2]
  otu_table_filtered <- otu_table_filtered[rowSums(otu_table_filtered) != 0, ]
  dim(otu_table_filtered)
  otu_table_filtered_transformed <- cmultRepl(otu_table_filtered)
  otu_table_filtered_transformed_log <- clr(otu_table_filtered_transformed)
  otu_table_filtered_transformed_log <- as.matrix(otu_table_filtered_transformed_log)
  ### coger solo muestras de contexto 
  ages <- unique(meta_exp$Age)
  locations <- unique(meta_exp$SampleLocation)
  # j <- "42"
  # k <- "C"
  for (j in ages){
    print(j)
    for (k in locations){
      print(paste0(j, "-", k))
      meta_ctx <-
        meta_exp[(meta_exp$Age == j & meta_exp$SampleLocation == k),]
      meta_ctx <- meta_ctx[meta_ctx$SampleID %in% rownames(otu_table_filtered),]
      samples_list <- meta_ctx$SampleID
      
      sel_vars <- c( "BW", "BWG", "BWGbefore", "FCR", "FCRbefore")
      # Dejar solo las columnas que no contienen NaNs
      col_nan <- colnames(meta_ctx)[colSums(is.na(meta_ctx)) > 0]
      print(col_nan)
      sel_vars <- sel_vars[!sel_vars %in% col_nan]
      print(paste0("Columnas con todos los valores -->", toString(sel_vars)))
      
      data_perf <- dplyr::select(meta_ctx, sel_vars)
      X <- otu_table_filtered_transformed_log[rownames(otu_table_filtered_transformed_log) %in% samples_list,]
      Y <- data_perf
      Y[, sel_vars] <- sapply(Y[, sel_vars], as.numeric)
      groups_ <- meta_ctx$SampleLocation
      corr_ctx <- correlation_siwa(X, Y, groups_, method = "pearson", adj_method = "BH")
      corr_ctx$Age <- j
      corr_ctx$ProjectID <- i
      corr_df_OTU <- rbind(corr_df_OTU, corr_ctx)
      rm(meta_ctx, X, Y, corr_ctx, samples_list, data_perf)
    }
  }
  rm(phy_exp, meta_exp, otu_table, otu_table_filtered, otu_table_filtered_transformed, otu_table_filtered_transformed_log)
  
}




# FILTER BY COLUMN (SIGNIFICANT) ----- 
column <- "Pvalue"
corr_df_OTU_filter <-  corr_df_OTU[corr_df_OTU[column] < 0.05, ]

broiler_score_OTU <- corr_df_OTU_filter %>% 
  mutate(Corr.relation = case_when(
    grepl("FCR", Env) & Correlation > 0 ~ "negative",
    grepl("FCR", Env) & Correlation <= 0 ~ "positive",
    grepl("BW", Env) & Correlation > 0 ~ "positive",
    grepl("BW", Env) & Correlation <= 0 ~ "negative",
  ))

colnames(broiler_score_OTU)
broiler_score_OTU$Context <- paste0(broiler_score_OTU$ProjectID, "-", broiler_score_OTU$Age, "-", broiler_score_OTU$Type) 
broiler_score_OTU$Context <- as.factor(broiler_score_OTU$Context)

# CREATE CONSISTENCY COLUMN ----- 
df_contexts <- broiler_score_OTU %>% group_by(Context)

list_df_contexts <- group_split(df_contexts)
contexts_order <- levels(df_contexts$Context)
names(list_df_contexts) <- contexts_order

### loop sobre los contextos 
corr_table_OTU <- data.frame()
for (df_con in list_df_contexts) {
  print(as.character(df_con$Context[[1]]))
  list_taxas <- unique(df_con$Taxa)
  print(length(list_taxas))
  for (b in list_taxas) {
    df_exp_bacteria <- df_con[df_con$Taxa == b,]
    df_exp_bacteria$ContCons <-
      ifelse(length(unique(df_exp_bacteria$Corr.relation)) == 1,  "C", "NC")
    corr_table_OTU <- rbind(corr_table_OTU, df_exp_bacteria)
  }
}

corr_dir <- "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/correlation/results/"
corr_table_OTU <- readRDS(paste0(corr_dir, "OtuCorr_120723.rds"))


tax_all$OTU <- rownames(tax_all)
taxa_all_otu <- dplyr::select(tax_all, OTU, Species, Genus)
corr_table_OTU <- merge(corr_table_OTU, taxa_all_otu, by.x = "Taxa", by.y = "OTU", all.x = TRUE)


### guardar tabla por genero 
saveRDS(corr_table_OTU, 
        "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/correlation/results/OtuCorr_120723.rds")

corr_table_OTU <- readRDS("/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/correlation/results/OtuCorr_120723.rds")
colnames(corr_table_OTU)
context_UCG_014 <- c("E267-28-C", "E271-28-C", "E326-21-C", "E335-42-C", "E347-42-C")
columna = "AdjPvalue"
corr_table_OTU_filter <-  corr_table_OTU[corr_table_OTU[columna] < 0.05, ]

corr_table_OTU_014 <- subset(corr_table_OTU_filter, Genus == "Clostridia_UCG-014")
length(unique(corr_table_OTU_014$Taxa))


corr_table_OTU_014_ctx <- subset(corr_table_OTU_014, Genus == "Clostridia_UCG-014" & Context %in% context_UCG_014)
length(unique(corr_table_OTU_014_ctx$Taxa))


taxa_score <- list()
list_taxas <- unique(corr_table_OTU_014_ctx$Taxa)


for (i in list_taxas) {
  df_taxa <- subset(corr_table_OTU_014_ctx, Taxa == i)
  list_samples <- unique(df_taxa$Type)
  for (j in list_samples) {
    df_taxa_loc <- subset(df_taxa, Type == j)
    unique_values <- unique(df_taxa_loc$Corr.relation)
    if (length(unique_values) != 1) {
      # Hay más de un valor único, continuar con la siguiente iteración del bucle
      next
    }
    score_sign <- ifelse(unique_values == "negative", -1, 1)
    Projects <- length(unique(df_taxa_loc$ProjectID))
    Context <- length(unique(df_taxa_loc$Context)) 
    Score <- Context * Projects * score_sign
    
    df_temp <- data.frame(Taxa = i,
                          SampleLocation = j,
                          Contexts = Context,
                          Projects = Projects,
                          Score = Score,
                          stringsAsFactors = FALSE)
    
    taxa_score <- rbind(taxa_score, df_temp)
  }
}


UCG014_corr.relation <- taxa_score
length(unique(UCG014_corr.relation$Taxa))
UCG014_corr.relation_positive <- subset(UCG014_corr.relation, SampleLocation == "C" & Score > 0)
length(unique(UCG014_corr.relation_positive$Taxa))


list_UCG014_otus <- UCG014_corr.relation$Taxa
list_UCG014_otus_positive <- UCG014_corr.relation_positive$Taxa

sum(UCG014_corr.relation_C$Score < 0)
sum(UCG014_corr.relation_C$Score > 0)


otu_table <- as.matrix(otu_table(phyloseq_rel))
otu_table <- t(otu_table) 
otu_table <- as.data.frame(otu_table)
colnames(otu_table)

meta_all <- as.data.frame(as.matrix(sample_data(ODLEPobj)))

meta_all$Context <- paste0(meta_all$ProjectID, "-", meta_all$Age, "-", meta_all$SampleLocation)
meta_all_OTU_014_ctx <- subset(meta_all, Context %in% context_UCG_014)
samples_ctx <- meta_all_OTU_014_ctx$SampleID


data_OTUs <- dplyr::select(otu_table, list_UCG014_otus)

data_OTUs_samples <- data_OTUs[rownames(data_OTUs) %in% samples_ctx,]
sum_OTUs <- as.data.frame(colSums(data_OTUs_samples))
sum_OTUs$OTUs <- rownames(sum_OTUs)
sum_OTUs$corr.relation <- ifelse(sum_OTUs$OTUs %in% list_UCG014_otus_positive, "positive", "negative")
colnames(sum_OTUs) <- c("Total.samples.abun", "Corr.relation", "OTUs")



prevalence <- as.data.frame(colSums(data_OTUs_samples > 0) / nrow(data_OTUs_samples))
prevalence$OTUs <- rownames(prevalence)
prevalence$corr.relation <- ifelse(prevalence$OTUs %in% list_UCG014_otus_positive, "positive", "negative")
colnames(prevalence) <- c("prevalence", "OTUs", "Corr.relation")



select(data_OTUs_samples, c931df285c2e86cf166b513d458bb187)

##### RF39

context_UCG_RF39 <- c("E271-28-C", "E321-21-C", "E326-21-C", "E335-42-C", "E347-42-C")
corr_table_OTU_RF39 <- subset(corr_table_OTU, Genus == "RF39")
length(unique(corr_table_OTU_RF39$Taxa))


corr_table_OTU_RF39_ctx <- subset(corr_table_OTU_RF39, Genus == "RF39" & Context %in% context_UCG_RF39)
length(unique(corr_table_OTU_RF39_ctx$Taxa))


taxa_score <- list()
list_taxas <- unique(corr_table_OTU_RF39_ctx$Taxa)

i <- "RF39"
for (i in list_taxas) {
  df_taxa <- subset(corr_table_OTU_RF39_ctx, Taxa == i)
  list_samples <- unique(df_taxa$Type)
  for (j in list_samples) {
    df_taxa_loc <- subset(df_taxa, Type == j)
    unique_values <- unique(df_taxa_loc$Corr.relation)
    if (length(unique_values) != 1) {
      # Hay más de un valor único, continuar con la siguiente iteración del bucle
      next
    }
    score_sign <- ifelse(unique_values == "negative", -1, 1)
    Projects <- length(unique(df_taxa_loc$ProjectID))
    Context <- length(unique(df_taxa_loc$Context)) 
    Score <- Context * Projects * score_sign
    
    df_temp <- data.frame(Taxa = i,
                          SampleLocation = j,
                          Contexts = Context,
                          Projects = Projects,
                          Score = Score,
                          stringsAsFactors = FALSE)
    
    taxa_score <- rbind(taxa_score, df_temp)
  }
}


RF39_corr.relation <- taxa_score
length(unique(RF39_corr.relation$Taxa))
RF39_corr.relation_positive <- subset(RF39_corr.relation, SampleLocation == "C" & Score > 0)
length(unique(UCGRF39_corr.relation_positive$Taxa))


list_RF39_otus <- RF39_corr.relation$Taxa
list_RF39_otus_positive <- RF39_corr.relation_positive$Taxa

sum(RF39_corr.relation$Score < 0)
sum(RF39_corr.relation$Score > 0)


otu_table <- as.matrix(otu_table(phyloseq_rel))
otu_table <- t(otu_table) 
otu_table <- as.data.frame(otu_table)
colnames(otu_table)

meta_all <- as.data.frame(as.matrix(sample_data(ODLEPobj)))

meta_all$Context <- paste0(meta_all$ProjectID, "-", meta_all$Age, "-", meta_all$SampleLocation)
meta_all_OTU_RF39_ctx <- subset(meta_all, Context %in% context_UCG_RF39)
samples_ctx <- meta_all_OTU_RF39_ctx$SampleID


data_OTUs <- dplyr::select(otu_table, list_RF39_otus)

data_OTUs_samples <- data_OTUs[rownames(data_OTUs) %in% samples_ctx,]
sum_OTUs <- as.data.frame(colSums(data_OTUs_samples))
sum_OTUs$OTUs <- rownames(sum_OTUs)
sum_OTUs$corr.relation <- ifelse(sum_OTUs$OTUs %in% list_UCGRF39_otus_positive, "positive", "negative")
colnames(sum_OTUs) <- c("Total.samples.abun", "Corr.relation", "OTUs")



prevalence <- as.data.frame(colSums(data_OTUs_samples > 0) / nrow(data_OTUs_samples))
prevalence$OTUs <- rownames(prevalence)
prevalence$corr.relation <- ifelse(prevalence$OTUs %in% list_UCGRF39_otus_positive, "positive", "negative")
colnames(prevalence) <- c("prevalence", "OTUs", "Corr.relation")



select(data_OTUs_samples, c931df285c2e86cf166b513d458bb187)
