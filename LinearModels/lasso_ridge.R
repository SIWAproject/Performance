library(phyloseq)
library(zCompositions)
library(compositions)
library(microbiomeMarker)
library(glmnet)
library(psych)#summaries


source("/Users/sebastianbedoyamazo/Documents/siwa_git/Methods-review/functions.R")
source("/Users/sebastianbedoyamazo/Documents/siwa_git/Methods-review/correlations.R")

input_dir <- "/Users/sebastianbedoyamazo/Documents/siwa_git/Pendientes/"
exp <- "Broiler"
ODLEPobj <- readRDS(paste0(input_dir, paste0("PhyloseqObject_", exp,".rds")))

meta_all <- as.data.frame(as.matrix(sample_data(ODLEPobj)))
tax_all <- as.data.frame(as.matrix(tax_table(ODLEPobj)))
phyloseq_rel <-microbiome::transform(ODLEPobj, "compositional") 

meta_all$ProjectID_Age  <- paste0(meta_all$ProjectID, "_", meta_all$Age)

length(unique(meta_all$ProjectID_Age))

filter_otus <- function(matrix, th=0.90){
  ### matrix: FEATURES as columns
  rowz<-apply(matrix == 0, 2, sum)/dim(matrix)[1] #sum of zeros over columns (OTUS)
  # dejar otus frecuentes, quitando los que están llenos de CEROS
  p <- which(rowz> th) ## p SON LOS QUE SALEN
  matrix[, -p]
}

#loop por ProjectID para transformaciones y almacenar en lista cada otu_table
#se debe correr por location para garantizar que se compartan la mayor cantidad 
# de otus entre PojectsID
Sample.location <- "F"
otu_table_project <- list()
project_ids <- c("E321", "E325", "E335", "E345", "E347", "E267", "E271", "Fieldale")
for (i in project_ids) {
  print(i)
  phy_exp <- subset_samples(phyloseq_rel, ProjectID== i)
  phy_exp <- subset_samples(phy_exp, SampleLocation == Sample.location)
  phy_exp <- prune_taxa(taxa_sums(phy_exp) > 0, phy_exp)
  # prune_samples(sample_sums(phy_exp) > 0, phy_exp)
  dfagg_exp <- aggregate_taxa_siwa(phy_exp, "Genus") 
  otu_table <- as.matrix(dfagg_exp)
  otu_table <- t(otu_table) 
  dim(otu_table)
  ### filtro: que tenga al menos dos OTUs diferentes de cero para que Rmult no falle
  otu_table_filtered <- otu_table[, colSums(otu_table > 0) >= 2]
  dim(otu_table_filtered)
  otu_table_filtered <- filter_otus(otu_table)
  dim(otu_table_filtered)
  otu_table_filtered <- otu_table_filtered[rowSums(otu_table_filtered) != 0, ]
  dim(otu_table_filtered)
  otu_table_filtered_transformed <- cmultRepl(otu_table_filtered)
  dim(otu_table_filtered_transformed)
  otu_table_filtered_transformed_log <- clr(otu_table_filtered_transformed)
  otu_table_filtered_transformed_log <- as.data.frame(otu_table_filtered_transformed_log)
  otu_table_project[[i]] <- otu_table_filtered_transformed_log
  
}

# a = otu_table_project$E267 #126 61
# colnames(a)
# a = otu_table_project$E271 #118 69
# colnames(a)
# a = otu_table_project$E321 #99 67
# colnames(a)
# a = otu_table_project$E325 #36 82
# colnames(a)
# a = otu_table_project$E326 #33 29
# colnames(a)
# a = otu_table_project$E335 #83 72
# colnames(a)
# a = otu_table_project$E345 #60 75
# colnames(a)
# a = otu_table_project$E347 #57 59
# colnames(a)
# a = otu_table_project$Fieldale #94 101
# colnames(a)


# Loop para convertir los row names en una columna en cada dataframe
for (i in 1:length(otu_table_project)) {
  otu_table_project[[i]]$SampleID <- rownames(otu_table_project[[i]])
  #rownames(otu_table_project[[i]]) <- NULL
}


# Obtener las columnas comunes entre todos los dataframes en la lista
common_columns <- Reduce(intersect, lapply(otu_table_project, colnames))

# Inicializar el dataframe resultante con el primer dataframe de la lista
otu_table_all_cecum <- otu_table_project[[1]]
otu_table_all_cecum <- dplyr::select(otu_table_all_cecum, common_columns)

# Realizar el merge con los demás dataframes en la lista
for (i in 2:length(otu_table_project)) {
  otu_table_project_filter <- dplyr::select(otu_table_project[[i]], common_columns)
  otu_table_all_cecum <- rbind(otu_table_all_cecum, otu_table_project_filter)
}




#filtrar meta_all
colnames(meta_all)
meta_all[, 10:14] <- lapply(meta_all[, 10:14], as.numeric)
str(meta_all)

meta_all_cecum <- meta_all[meta_all$SampleID %in% rownames(otu_table_all_cecum),]
sel_vars <-c( "BW", "BWG", "BWGbefore", "FCR", "FCRbefore")
meta_all_cecum <- dplyr::select(meta_all_cecum, sel_vars, ProjectID, SampleLocation, Age)
XY <- merge(otu_table_all_cecum, meta_all_cecum, by = "row.names", all = TRUE)
row.names(XY) <- XY$Row.names
XY$Row.names <- NULL
XY$SampleID <- NULL


################################################################################
# Función para normalizar columnas específicas en función de columna categórica
normalize_columns <- function(data, categorical_col, cols_norm) {
  normalized_data <- data %>%
    group_by({{categorical_col}}) %>%
    mutate(across({{cols_norm}}, scale)) %>%
    ungroup()
  return(normalized_data)
}

colnames(XY)
# Columnas numéricas que se normalizarán
columns_to_normalize <- colnames(XY)[1:17]

# Normalizar las columnas var1 y var2 en función de la columna categoria_col
normalized_df <- normalize_columns(XY, ProjectID, columns_to_normalize)
describeBy(normalized_df, group = normalized_df$ProjectID)

################################################################################
colnames(XY)
lasso_table_list <- list()
for (i in sel_vars) {
  #complete cases para cada performance variable
  XY_complcases <- normalized_df[complete.cases(normalized_df[i]), ]
  
  #X y Y para linear model
  colnames(XY_complcases)
  Y <- XY_complcases[[i]]
  X <- as.matrix(XY_complcases[, 1:12])
  
  # Entrenar un modelo de Regresión Lasso
  lasso_model <- glmnet(X, Y, alpha = 1)  # alpha = 1 para Lasso
  
  
  n_lambda <- cv.glmnet(X, Y, alpha = 1)
  best_lambda <- n_lambda$lambda.min
  
  
  # Calcula los coeficientes para el valor óptimo de lambda
  best_model <- glmnet(X, Y, alpha = 1, lambda = best_lambda)
  lasso_coeff <- as.data.frame(as.matrix(coef(best_model)))
  #lasso_coeff$Taxa <- row.names(lasso_coeff)
  # Renombrar las columnas
  colnames(lasso_coeff)[1] <- as.character(i) 
  class(lasso_coeff)
  
  # Almacena el dataframe en la lista
  lasso_table_list[[i]] <- lasso_coeff
}



# Combina todos los dataframes en uno solo
lasso_table <- do.call(cbind, lasso_table_list)
lasso_table$Taxa <- row.names(lasso_table)


long_lasso_table <- gather(lasso_table, PerfVariables, Slopes, -Taxa)
long_lasso_table$SampleLocation <- Sample.location


################################################################################
#consistency LASSO COEFF
# FILTER BY COLUMN (SIGNIFICANT) ----- 
long_lasso_table_score <- long_lasso_table %>% 
  mutate(Coeff.relation = case_when(
    grepl("FCR", PerfVariables) & Slopes > 0 ~ "negative",
    grepl("FCR", PerfVariables) & Slopes <= 0 ~ "positive",
    grepl("BW", PerfVariables) & Slopes > 0 ~ "positive",
    grepl("BW", PerfVariables) & Slopes <= 0 ~ "negative",
  ))


# CREATE CONSISTENCY COLUMN ----- 
df_contexts <- long_lasso_table_score %>% group_by(Taxa)

list_df_contexts <- group_split(df_contexts)
contexts_order <- levels(df_contexts$Taxa)
names(list_df_contexts) <- contexts_order

### loop sobre los contextos para ver consistencia (esta casi que siempre debe cumplirse) ---- 
lasso_consist <- data.frame()
for (df_con in list_df_contexts) {
  print(as.character(df_con$SampleLocation[[1]]))
  list_taxas <- unique(df_con$Taxa)
  print(length(list_taxas)) 
  for (b in list_taxas) {
    df_exp_bacteria <- df_con[df_con$Taxa == b,]
    df_exp_bacteria$ContCons <-
      ifelse(length(unique(df_exp_bacteria$Coeff.relation)) == 1,  "C", "NC")
    lasso_consist <- rbind(lasso_consist, df_exp_bacteria)
  }
}




unique(lasso_consist$Taxa)


THRESHOLD = 0.75
taxa_score <- data.frame()
list_taxas <- unique(lasso_consist$Taxa)
for (tax in list_taxas) {
  #tax = "Bacteroides"
  df_taxa <- subset(lasso_consist, Taxa == tax)
  max_freq <- max(table(df_taxa$Coeff.relation))
  max_prop <- max_freq/dim(df_taxa)[1]
  #prop.table(table(df_taxa$Corr.relation))
  if (max_prop < THRESHOLD){next}
  ## quito las inconsistencias del score
  score_sign <- names(table(df_taxa$Coeff.relation)[table(df_taxa$Coeff.relation) == max_freq])
  df_taxa <- df_taxa[df_taxa$Coeff.relation == score_sign,] 
  PerfVariables <- paste0(unique(df_taxa$PerfVariables), collapse = ", ")
  ProjectIDs <- paste0(unique(XY_complcases$ProjectID), collapse = ", ")
  PropConsistency <- max_prop
  df_row <- data.frame(Taxa = tax,
                       ProjectID = ProjectIDs,
                       SampleLocation = Sample.location,
                       PerfVariables = PerfVariables,
                       PropConsistency= max_prop,
                       Sign = score_sign,
                       stringsAsFactors = FALSE)
  taxa_score <- rbind(taxa_score, df_row)
}
# unique(taxa_score$Taxa) 
# 
# saveRDS(taxa_score,
#         "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/linear_mixed_models/codes and tables/consistencytable_globalcor_genus_cecum.rds")

################################################################################
colnames(normalized_df)
ridge_table_list <- list()
for (i in sel_vars) {
  #complete cases para cada performance variable
  XY_complcases <- normalized_df[complete.cases(normalized_df[i]), ]
  
  #X y Y para linear model
  colnames(XY_complcases)
  Y <- XY_complcases[[i]]
  X <- as.matrix(XY_complcases[, 1:12])
  
  # Entrenar un modelo de Regresión ridge
  ridge_model <- glmnet(X, Y, alpha = 0)  # alpha = 1 para ridge
  
  
  n_lambda <- cv.glmnet(X, Y, alpha = 0)
  best_lambda <- n_lambda$lambda.min
  
  
  # Calcula los coeficientes para el valor óptimo de lambda
  best_model <- glmnet(X, Y, alpha = 0, lambda = best_lambda)
  ridge_coeff <- as.data.frame(as.matrix(coef(best_model)))
  #ridge_coeff$Taxa <- row.names(ridge_coeff)
  # Renombrar las columnas
  colnames(ridge_coeff)[1] <- as.character(i) 
  class(ridge_coeff)
  
  # Almacena el dataframe en la lista
  ridge_table_list[[i]] <- ridge_coeff
}



# Combina todos los dataframes en uno solo
ridge_table <- do.call(cbind, ridge_table_list)
ridge_table$Taxa <- row.names(ridge_table)


long_ridge_table <- gather(ridge_table, PerfVariables, Slopes, -Taxa)
long_ridge_table$SampleLocation <- Sample.location


################################################################################
#consistency ridge COEFF
# FILTER BY COLUMN (SIGNIFICANT) ----- 
long_ridge_table_score <- long_ridge_table %>% 
  mutate(Coeff.relation = case_when(
    grepl("FCR", PerfVariables) & Slopes > 0 ~ "negative",
    grepl("FCR", PerfVariables) & Slopes <= 0 ~ "positive",
    grepl("BW", PerfVariables) & Slopes > 0 ~ "positive",
    grepl("BW", PerfVariables) & Slopes <= 0 ~ "negative",
  ))


# CREATE CONSISTENCY COLUMN ----- 
df_contexts <- long_ridge_table_score %>% group_by(Taxa)

list_df_contexts <- group_split(df_contexts)
contexts_order <- levels(df_contexts$Taxa)
names(list_df_contexts) <- contexts_order

### loop sobre los contextos para ver consistencia (esta casi que siempre debe cumplirse) ---- 
ridge_consist <- data.frame()
for (df_con in list_df_contexts) {
  print(as.character(df_con$SampleLocation[[1]]))
  list_taxas <- unique(df_con$Taxa)
  print(length(list_taxas)) 
  for (b in list_taxas) {
    df_exp_bacteria <- df_con[df_con$Taxa == b,]
    df_exp_bacteria$ContCons <-
      ifelse(length(unique(df_exp_bacteria$Coeff.relation)) == 1,  "C", "NC")
    ridge_consist <- rbind(ridge_consist, df_exp_bacteria)
  }
}




unique(ridge_consist$Taxa)


THRESHOLD = 0.75
taxa_score <- data.frame()
list_taxas <- unique(ridge_consist$Taxa)
for (tax in list_taxas) {
  #tax = "Bacteroides"
  df_taxa <- subset(ridge_consist, Taxa == tax)
  max_freq <- max(table(df_taxa$Coeff.relation))
  max_prop <- max_freq/dim(df_taxa)[1]
  #prop.table(table(df_taxa$Corr.relation))
  if (max_prop < THRESHOLD){next}
  ## quito las inconsistencias del score
  score_sign <- names(table(df_taxa$Coeff.relation)[table(df_taxa$Coeff.relation) == max_freq])
  df_taxa <- df_taxa[df_taxa$Coeff.relation == score_sign,] 
  PerfVariables <- paste0(unique(df_taxa$PerfVariables), collapse = ", ")
  ProjectIDs <- paste0(unique(XY_complcases$ProjectID), collapse = ", ")
  PropConsistency <- max_prop
  df_row <- data.frame(Taxa = tax,
                       ProjectID = ProjectIDs,
                       SampleLocation = Sample.location,
                       PerfVariables = PerfVariables,
                       PropConsistency= max_prop,
                       Sign = score_sign,
                       stringsAsFactors = FALSE)
  taxa_score <- rbind(taxa_score, df_row)
}
# unique(taxa_score$Taxa) 
# 
# saveRDS(taxa_score,
#         "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/linear_mixed_models/codes and tables/consistencytable_globalcor_genus_cecum.rds")

