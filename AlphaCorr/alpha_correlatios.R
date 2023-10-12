source("/Users/sebastianbedoyamazo/Documents/siwa_git/Methods-review/functions.R")
source("/Users/sebastianbedoyamazo/Documents/siwa_git/Methods-review/correlations.R")

folder = "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/correlation/alpha/"
metadata_alpha <- read.csv(file = paste0(folder, "metadata_alpha.csv"),check.names = FALSE, sep = ",")



metadata_alpha$SampleLocation <- gsub("cecum", "C", metadata_alpha$SampleLocation)
metadata_alpha$SampleLocation <- gsub("ileum", "I", metadata_alpha$SampleLocation)
metadata_alpha$SampleLocation <- gsub("feces", "F", metadata_alpha$SampleLocation)

metadata_alpha$Context <- paste0(metadata_alpha$ProjectID, "-", metadata_alpha$Age, "-", metadata_alpha$SampleLocation) 
metadata_alpha$Context <- as.factor(metadata_alpha$Context)


colnames(metadata_alpha)

ggplot(metadata_alpha, aes(x = ProjectID, y = alphaobserved, fill=ProjectID)) +
  geom_violin() +
  facet_wrap(~SampleLocation) + 
  xlab("") +
  ylab("alphaobserved") +
  ggtitle("Distribution of alphaobserved across Siwa projects") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1),
        text = element_text(size = 15), axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 15))+
  geom_jitter(color="black", size=0.8, alpha=0.9)
  
  


metadata_alpha_cecum <- subset(metadata_alpha, SampleLocation == "C")
library(ggstatsplot)
ggbetweenstats(
  data  = metadata_alpha_cecum,
  x     = ProjectID,
  y     = alphashannon,
  title = "Distribution of Alphashannon across Cecum Siwa-projects"
)



metadata_alpha_cecum <- subset(metadata_alpha, SampleLocation == "C")
metadata_alpha_ileum <- subset(metadata_alpha, SampleLocation == "I")
metadata_alpha_feces <- subset(metadata_alpha, SampleLocation == "F")


ggplot(metadata_alpha_cecum, aes(x = alphashannon)) +
  geom_histogram(binwidth = 0.2, fill = "#5cb08e", color = "black") +
  labs(title = "Raw Data - Cecum") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(metadata_alpha_cecum, aes(x = alphashannon)) +
  geom_histogram(binwidth = 0.2, fill = "#5cb08e", color = "black") +
  labs(title = "Raw Data by ProjectID - Cecum") +
  facet_wrap(~ ProjectID) +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(metadata_alpha_cecum, aes(x = alphashannon, y = BW, col = Context)) +
  geom_point() +  # Agregar puntos de dispersión
  geom_smooth(method = "lm", se = FALSE) +  # Agregar la línea de tendencia lineal
  labs(title = "Raw Data")  +
  theme(plot.title = element_text(hjust = 0.5))


################################################################################
# Función para normalizar columnas específicas en función de columna categórica
normalize_columns <- function(data, categorical_col, cols_norm) {
  normalized_data <- data %>%
    group_by({{categorical_col}}) %>%
    mutate(across({{cols_norm}}, scale)) %>%
    ungroup()
  return(normalized_data)
}

colnames(metadata_alpha)
# Columnas numéricas que se normalizarán
columns_to_normalize <- c("BW", "BWGbefore", "BWG", "FCRbefore", "FCR", "alphashannon", "alphaobserved")

# Normalizar las columnas var1 y var2 en función de la columna categoria_col
normalized_df <- as.data.frame(normalize_columns(metadata_alpha_cecum, ProjectID, columns_to_normalize))
row.names(normalized_df) <- normalized_df$SampleID
colnames(normalized_df)



ggplot(normalized_df, aes(x = alphashannon)) +
  geom_histogram(binwidth = 0.2, fill = "#5cb08e", color = "black") +
  labs(title = "Scaled Data - Cecum") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(normalized_df, aes(x = alphashannon)) +
  geom_histogram(binwidth = 0.2, fill = "#5cb08e", color = "black") +
  labs(title = "Scaled Data by ProjectID - Cecum") +
  facet_wrap(~ ProjectID) +
  theme(plot.title = element_text(hjust = 0.5))


ggplot(normalized_df, aes(x = alphashannon, y = BW, col = Context)) +
  geom_point() +  # Agregar puntos de dispersión
  geom_smooth(method = "lm", se = FALSE) +  # Agregar la línea de tendencia lineal
  labs(title = "Scaled data")  +
  theme(plot.title = element_text(hjust = 0.5))

################################################################################



project_ids <- unique(metadata_alpha$ProjectID)
corr_df <- data_frame()
for (i in project_ids) {
  metadata_alpha_pj <- subset(metadata_alpha, ProjectID == i)
  ages <- unique(metadata_alpha_pj$Age)
  locations <- unique(metadata_alpha_pj$SampleLocation)
  #j <- "42"
  #k <- "C"
  for (j in ages){
    print(j)
    for (k in locations){
      print(paste0(j, "-", k))
      meta_ctx <-
        metadata_alpha_pj[(metadata_alpha_pj$Age == j & metadata_alpha_pj$SampleLocation == k),]
      sel_vars <-c( "BW", "BWG", "BWGbefore", "FCR", "FCRbefore")
      # Dejar solo las columnas que no contienen NaNs
      col_nan <- colnames(meta_ctx)[colSums(is.na(meta_ctx)) > 0]
      sel_vars <- sel_vars[!sel_vars %in% col_nan]
      
      print(paste0("Columnas con todos los valores -->", toString(sel_vars)))
      data_perf <- as.data.frame(dplyr::select(meta_ctx, sel_vars)) 
      Y <- data_perf
      Y[, sel_vars] <- sapply(Y[, sel_vars], as.numeric)
      
      alpha_vars <-c( "alphashannon", "alphaobserved")
      data_alpha <- as.data.frame(dplyr::select(meta_ctx, alpha_vars))
      X <- data_alpha
      groups_ <- meta_ctx$SampleLocation
      #correlation siwa for each ProjectID, Age and SampleLocation
      corr_ctx <- correlation_siwa(X, Y, groups_, method = "pearson", adj_method = "BH")
      # Store the correlation results for the current age group in the corr_list
      corr_ctx$Age <- j
      corr_ctx$ProjectID <- i
      corr_df <- rbind(corr_df, corr_ctx)
      rm(meta_ctx, X, Y, corr_ctx, samples_list, data_perf)
    }
  }
  print("------------------------")
  rm(metadata_alpha_pj, meta_ctx, data_alpha, data_perf)
}




################################################################################
#consistency correlations
# FILTER BY COLUMN (SIGNIFICANT) ----- 
columna <- "Pvalue"
corr_all_filter <-  corr_df[corr_df[columna] < 0.05, ]

corr_all_filter_score <- corr_all_filter %>% 
  mutate(Corr.relation = case_when(
    grepl("FCR", Env) & Correlation > 0 ~ "negative",
    grepl("FCR", Env) & Correlation <= 0 ~ "positive",
    grepl("BW", Env) & Correlation > 0 ~ "positive",
    grepl("BW", Env) & Correlation <= 0 ~ "negative",
  ))


colnames(corr_all_filter_score)

# CREATE CONSISTENCY COLUMN ----- 
df_contexts <- corr_all_filter_score %>% group_by(Taxa)

list_df_contexts <- group_split(df_contexts)
contexts_order <- levels(df_contexts$Taxa)
names(list_df_contexts) <- contexts_order

### loop sobre los contextos para ver consistencia (esta casi que siempre debe cumplirse) ---- 
corr_consist <- data.frame()
for (df_con in list_df_contexts) {
  print(as.character(df_con$Type[[1]]))
  list_taxas <- unique(df_con$Taxa)
  print(length(list_taxas)) 
  for (b in list_taxas) {
    df_exp_bacteria <- df_con[df_con$Taxa == b,]
    df_exp_bacteria$ContCons <-
      ifelse(length(unique(df_exp_bacteria$Corr.relation)) == 1,  "C", "NC")
    corr_consist <- rbind(corr_consist, df_exp_bacteria)
  }
}


saveRDS(corr_consist,
        "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/AlphaCorr/significant_alphacorrelations.rds")


unique(corr_consist$Taxa)

colnames(corr_consist)
THRESHOLD = 0.75
taxa_score <- data.frame()
list_taxas <- unique(corr_consist$Taxa)
list_locations <- unique(corr_consist$Type)
for (tax in list_taxas) {
  for (loc in list_locations) {
  #tax = "Bacteroides"
  df_taxa <- subset(corr_consist, Taxa == tax & Type == loc)
  max_freq <- max(table(df_taxa$Corr.relation))
  max_prop <- max_freq/dim(df_taxa)[1]
  #prop.table(table(df_taxa$Corr.relation))
  if (max_prop < THRESHOLD){next}
  ## quito las inconsistencias del score
  score_sign <- names(table(df_taxa$Corr.relation)[table(df_taxa$Corr.relation) == max_freq])
  df_taxa <- df_taxa[df_taxa$Corr.relation == score_sign,] 
  PerfVariables <- paste0(unique(df_taxa$Env), collapse = ", ")
  ProjectIDs <- paste0(unique(df_taxa$ProjectID), collapse = ", ")
  PropConsistency <- max_prop
  df_row <- data.frame(Taxa = tax,
                       ProjectID = ProjectIDs,
                       SampleLocation = loc,
                       PerfVariables = PerfVariables,
                       PropConsistency= max_prop,
                       Sign = score_sign,
                       stringsAsFactors = FALSE)
  taxa_score <- rbind(taxa_score, df_row)
  }
}


unique(taxa_score$Taxa)


saveRDS(taxa_score,
        "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/AlphaCorr/consistency_alphacorrelations.rds")
