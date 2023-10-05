source("/Users/sebastianbedoyamazo/Documents/siwa_git/Methods-review/functions.R")
source("/Users/sebastianbedoyamazo/Documents/siwa_git/Methods-review/correlations.R")

input_dir <- "/Users/sebastianbedoyamazo/Documents/siwa_git/Pendientes/"
exp <- "Broiler"
ODLEPobj <- readRDS(paste0(input_dir, paste0("PhyloseqObject_", exp,".rds")))

meta_all <- as.data.frame(as.matrix(sample_data(ODLEPobj)))
tax_all <- as.data.frame(as.matrix(tax_table(ODLEPobj)))
subset(tax_all, Genus == "Clostridia_UCG")
phyloseq_rel <-microbiome::transform(ODLEPobj, "compositional") 


filter_otus <- function(matrix, th=0.90){
  ### matrix: FEATURES as columns
  rowz<-apply(matrix == 0, 2, sum)/dim(matrix)[1] #sum of zeros over columns (OTUS)
  # dejar otus frecuentes, quitando los que están llenos de CEROS
  p <- which(rowz> th) ## p SON LOS QUE SALEN
  matrix[, -p]
}

corr_df_G <- data_frame()
project_ids <- c("E321", "E325", "E326", "E335", "E345", "E347", "E267", "E271", "Fieldale")
#i <- "E345"
for (i in project_ids) {
  print(i)
  phy_exp <- subset_samples(phyloseq_rel, ProjectID== i)
  phy_exp <- prune_taxa(taxa_sums(phy_exp) > 0, phy_exp)
  # prune_samples(sample_sums(phy_exp) > 0, phy_exp)
  meta_exp <- as.data.frame(as.matrix(sample_data(phy_exp)))
  dfagg_exp <- aggregate_taxa_siwa(phy_exp, "Genus") 
  otu_table <- as.matrix(dfagg_exp)
  otu_table <- t(otu_table) 
  dim(otu_table)
  ### filtro: que tenga al menos dos OTUs diferentes de cero para que Rmult no falle
  otu_table_filtered <- otu_table[, colSums(otu_table > 0) >= 2]
  otu_table_filtered <- filter_otus(otu_table)
  otu_table_filtered <- otu_table_filtered[rowSums(otu_table_filtered) != 0, ]
  dim(otu_table_filtered)
  otu_table_filtered_transformed <- cmultRepl(otu_table_filtered)
  otu_table_filtered_transformed_log <- clr(otu_table_filtered_transformed)
  otu_table_filtered_transformed_log <- as.data.frame(otu_table_filtered_transformed_log)
  ### coger solo muestras de contexto 
  ages <- unique(meta_exp$Age)
  locations <- unique(meta_exp$SampleLocation)
  #j <- "42"
  #k <- "C"
  for (j in ages){
    print(j)
    for (k in locations){
      print(paste0(j, "-", k))
      meta_ctx <-
        meta_exp[(meta_exp$Age == j & meta_exp$SampleLocation == k),]
      meta_ctx <- meta_ctx[meta_ctx$SampleID %in% rownames(otu_table_filtered),]
      samples_list <- meta_ctx$SampleID
      sel_vars <-c( "BW", "BWG", "BWGbefore", "FCR", "FCRbefore")
      # Dejar solo las columnas que no contienen NaNs
      col_nan <- colnames(meta_ctx)[colSums(is.na(meta_ctx)) > 0]
      sel_vars <- sel_vars[!sel_vars %in% col_nan]
      
      print(paste0("Columnas con todos los valores -->", toString(sel_vars)))
      data_perf <- as.data.frame(dplyr::select(meta_ctx, sel_vars)) 
      Y <- data_perf
      Y[, sel_vars] <- sapply(Y[, sel_vars], as.numeric)
      X <-
        otu_table_filtered_transformed_log[rownames(otu_table_filtered_transformed_log) %in% samples_list,]
      Y <-
        otu_table_filtered_transformed_log[rownames(otu_table_filtered_transformed_log) %in% samples_list,]
      groups_ <- meta_ctx$SampleLocation
      #correlation siwa for each ProjectID, Age and SampleLocation
      corr_ctx <- correlation_siwa(X, Y, groups_, method = "pearson", adj_method = "BH")
      # Store the correlation results for the current age group in the corr_list
      corr_ctx$Age <- j
      corr_ctx$ProjectID <- i
      corr_df_G <- rbind(corr_df_G, corr_ctx)
      rm(meta_ctx, X, Y, corr_ctx, samples_list, data_perf)
    }
  }
  print("------------------------")
  rm(phy_exp, meta_exp,dfagg_exp,otu_table, otu_table_filtered, otu_table_filtered_transformed, otu_table_filtered_transformed_log)
}


# FILTER BY COLUMN (SIGNIFICANT) ----- 
columna <- "Pvalue"
corr_df_G_filter <-  corr_df_G[corr_df_G[columna] < 0.05, ]

broiler_score_G <- corr_df_G_filter %>% 
  mutate(Corr.relation = case_when(
    grepl("FCR", Env) & Correlation > 0 ~ "negative",
    grepl("FCR", Env) & Correlation <= 0 ~ "positive",
    grepl("BW", Env) & Correlation > 0 ~ "positive",
    grepl("BW", Env) & Correlation <= 0 ~ "negative",
  ))

broiler_score_G$Context <- paste0(broiler_score_G$ProjectID, "-", broiler_score_G$Age, "-", broiler_score_G$Type) 
broiler_score_G$Context <- as.factor(broiler_score_G$Context)

# CREATE CONSISTENCY COLUMN ----- 
df_contexts <- broiler_score_G %>% group_by(Context)

list_df_contexts <- group_split(df_contexts)
contexts_order <- levels(df_contexts$Context)
names(list_df_contexts) <- contexts_order

### loop sobre los contextos para ver consistencia (esta casi que siempre debe cumplirse) ---- 
corr_table_G <- data.frame()
for (df_con in list_df_contexts) {
  print(as.character(df_con$Context[[1]]))
  list_taxas <- unique(df_con$Taxa)
  print(length(list_taxas)) 
  for (b in list_taxas) {
    df_exp_bacteria <- df_con[df_con$Taxa == b,]
    df_exp_bacteria$ContCons <-
      ifelse(length(unique(df_exp_bacteria$Corr.relation)) == 1,  "C", "NC")
    corr_table_G <- rbind(corr_table_G, df_exp_bacteria)
  }
}


### guardar tabla por genero 
saveRDS(corr_df_G_filter, 
        "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/correlation/results/GeneraVsGeneraCorr_280723.rds")


corr_dir <- "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/correlation/results/"
corr_table_G <- readRDS(paste0(corr_dir, "GeneraCorr_120723.rds"))
colnames(corr_table_G)


### loop sobre EXP+LOC para ver consistencia (esta es menos obvia) ---- 

corr_table_G$ExpLoc <- paste0(corr_table_G$ProjectID, "-", corr_table_G$Type) 
corr_table_G$ExpLoc <- as.factor(corr_table_G$ExpLoc)
df_contexts <- corr_table_G %>% group_by(ExpLoc)

list_df_contexts <- group_split(df_contexts)
contexts_order <- levels(df_contexts$ExpLoc)
names(list_df_contexts) <- contexts_order

corr_table_G2 <- data.frame()
for (df_con in list_df_contexts) {
  print(as.character(df_con$Context[[1]]))
  list_taxas <- unique(df_con$Taxa)
  print(length(list_taxas)) 
  for (b in list_taxas) {
    df_exp_bacteria <- df_con[df_con$Taxa == b,]
    df_exp_bacteria$ExpLocCons <-
      ifelse(length(unique(df_exp_bacteria$Corr.relation)) == 1,  "C", "NC")
    corr_table_G2 <- rbind(corr_table_G2, df_exp_bacteria)
  }
}


corr_dir <- ("/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/correlation/results/GeneraCorr_120723.rds")
table(corr_table_G2$ExpLocCons)
saveRDS(corr_table_G2, 
        "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/correlation/results/GeneraCorr_120723.rds")

corr_dir <- "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/correlation/results/"
corr_table_G <- readRDS(paste0(corr_dir, "GeneraCorr_120723.rds"))


## PRUEBAS EXTRAS ---- 
taxa_score <- list()
colnames(corr_table_G)
list_taxas <- unique(corr_table_G$Taxa)

i <- "RF39"
for (i in list_taxas) {
  df_taxa <- subset(corr_table_G, Taxa == i)
  list_samples <- unique(df_taxa$Type)
  for (j in list_samples) {
    df_taxa_loc <- subset(df_taxa, Type == j)
    unique_values <- unique(df_taxa_loc$ExpLocCons)
    if (length(unique_values) != 1) {
      # Hay más de un valor único, continuar con la siguiente iteración del bucle
      next
    }
      #score_sign <- ifelse(unique_values == "negative", -1, 1)
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

taxa_score_corr.relation 
taxa_score_corr.relation_C <- subset(taxa_score_corr.relation, SampleLocation == "C")
length(unique(taxa_score_corr.relation_C$Taxa))
taxa_score_corr.relation_I <- subset(taxa_score_corr.relation, SampleLocation == "I")
length(unique(taxa_score_corr.relation_I$Taxa))
taxa_score_corr.relation_F <- subset(taxa_score_corr.relation, SampleLocation == "F")
length(unique(taxa_score_corr.relation_F$Taxa))


taxa_score_ContCons 
taxa_score_ContCons_C <- subset(taxa_score_ContCons, SampleLocation == "C")
length(unique(taxa_score_ContCons_C$Taxa))
taxa_score_ContCons_I <- subset(taxa_score_ContCons, SampleLocation == "I")
length(unique(taxa_score_ContCons_I$Taxa))
taxa_score_ContCons_F <- subset(taxa_score_ContCons, SampleLocation == "F")
length(unique(taxa_score_ContCons_F$Taxa))


taxa_score_ExpLocCons 
taxa_score_ExpLocCons_C <- subset(taxa_score_ExpLocCons, SampleLocation == "C")
length(unique(taxa_score_ExpLocCons_C$Taxa))
taxa_score_ExpLocCons_I <- subset(taxa_score_ExpLocCons, SampleLocation == "I")
length(unique(taxa_score_ExpLocCons_I$Taxa))
taxa_score_ExpLocCons_F <- subset(taxa_score_ExpLocCons, SampleLocation == "F")
length(unique(taxa_score_ExpLocCons_F$Taxa))


corr_table_G_filter <-
  subset(corr_table_G, Taxa == "RF39" & Type == "C")
unique(corr_table_G_filter$Context)




for (i in list_taxas) {
  df_taxa <- subset(corr_table_G, Taxa == i)
  list_samples <- unique(df_taxa$Type)
  for (j in list_samples) {
    df_taxa_loc <- subset(df_taxa, Type == j)
    if (nrow(df_taxa_loc) == 0) {
      # No hay filas para el valor de j, continuar con la siguiente iteración del bucle
      next
    }
    if (all(unique(df_taxa_loc$Corr.relation == 1))) {
      score_sign <- ifelse(all(unique(df_taxa_loc$Corr.relation) == "negative"), -1, 1)
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
}


