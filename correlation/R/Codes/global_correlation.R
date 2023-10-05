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
Sample.location <- "C"
otu_table_project <- list()
project_ids <- c("E321", "E325", "E326", "E335", "E345", "E347", "E267", "E271", "Fieldale")
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

a = otu_table_project$E267 
a = otu_table_project$E271 
a = otu_table_project$E321 
a = otu_table_project$E325 
a = otu_table_project$E326 
a = otu_table_project$E335 
a = otu_table_project$E345 
a = otu_table_project$E347 
a = otu_table_project$Fieldale 

# Loop para convertir los row names en una columna en cada dataframe
for (i in 1:length(otu_table_project)) {
  otu_table_project[[i]]$SampleID <- rownames(otu_table_project[[i]])
  rownames(otu_table_project[[i]]) <- NULL
}

# Inicializar el DataFrame resultado con el primer DataFrame de la lista
otu_table_all_cecum <- otu_table_project[[1]]
dim(otu_table_project[[1]])

# Loop para realizar full join con cada dataframe en otu_table_project
for (i in 2:length(otu_table_project)) {
  common_columns <- intersect(colnames(otu_table_all_cecum), colnames(otu_table_project[[i]]))
  otu_table_all_cecum <- otu_table_all_cecum %>%
    full_join(otu_table_project[[i]], by = common_columns)
}
dim(otu_table_all_cecum)

otu_table_all_cecum$SampleID
row.names(otu_table_all_cecum) <-  otu_table_all_cecum$SampleID
otu_table_all_cecum$SampleID <- NULL

#filtrar meta_all
colnames(meta_all)
meta_all[, 10:14] <- lapply(meta_all[, 10:14], as.numeric)
str(meta_all)

meta_all_cecum <- meta_all[meta_all$SampleID %in% rownames(otu_table_all_cecum),]
sel_vars <-c( "BW", "BWG", "BWGbefore", "FCR", "FCRbefore")
meta_all_cecum <- dplyr::select(meta_all_cecum, sel_vars, ProjectID, SampleLocation)
XY <- merge(otu_table_all_cecum, meta_all_cecum, by = "row.names", all = TRUE)
row.names(XY) <- XY$Row.names
XY$Row.names <- NULL
colnames(XY)

ggplot(XY, aes(x = Bacteroides)) +
  geom_histogram(binwidth = 0.2, fill = "#5cb08e", color = "black") +
  labs(title = "Raw Data") +
  theme(plot.title = element_text(hjust = 0.5))


ggplot(XY, aes(x = Bacteroides)) +
  geom_histogram(binwidth = 0.2, fill = "#5cb08e", color = "black") +
  labs(title = "Raw Data by ProjectID") +
  facet_wrap(~ ProjectID) +
  theme(plot.title = element_text(hjust = 0.5))



# saveRDS(XY,
#         "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/xy_clr_byproject_feces.rds")

################################################################################
#correlation
# Función para normalizar columnas específicas en función de columna categórica
XY$SampleID <- rownames(XY)
normalize_columns <- function(data, categorical_col, cols_norm) {
  normalized_data <- data %>%
    group_by({{categorical_col}}) %>%
    mutate(across({{cols_norm}}, scale)) %>%
    ungroup()
  return(normalized_data)
}

colnames(XY)
# Columnas numéricas que se normalizarán
columns_to_normalize <- column_names <- colnames(XY)[1:125]

# Normalizar las columnas var1 y var2 en función de la columna categoria_col
normalized_df <- as.data.frame(normalize_columns(XY, ProjectID, columns_to_normalize))
row.names(normalized_df) <- normalized_df$SampleID
colnames(normalized_df)

ggplot(normalized_df, aes(x = Bacteroides)) +
  geom_histogram(binwidth = 0.2, fill = "#5cb08e", color = "black") +
  labs(title = "Scaled Data") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(normalized_df, aes(x = Bacteroides)) +
  geom_histogram(binwidth = 0.2, fill = "#5cb08e", color = "black") +
  labs(title = "Scaled Data by ProjectID") +
  facet_wrap(~ ProjectID) +
  theme(plot.title = element_text(hjust = 0.5))




summary_table <- describe(normalized_df)


Xvariables <- as.data.frame(normalized_df[, 1:140])

Yvariables <- as.data.frame(normalized_df[, 141:145])
str(Yvariables)

x_colnames <- colnames(Xvariables)

df_all_cor <- data.frame()
method <- "pearson"
adj_method <- "BH"
for (i in x_colnames) {
  for (j in sel_vars) {
    X1 <- dplyr::select(Xvariables, i)
    Y1 <- dplyr::select(Yvariables, j)
    XY1 <- cbind(X1, Y1)
    XY1.completecases <- XY1[complete.cases(XY1), ]
    if (nrow(XY1.completecases) == 0) {
      cat("there is not data for", paste0(i, "_", j), "combination", "\n")
      next  
    }
    
    Xcor <- dplyr::select(XY1.completecases, i)
    Ycor <- dplyr::select(XY1.completecases, j)
    corr.pval <- cor.test(Xcor[complete.cases(Ycor),], Ycor[complete.cases(Ycor),], method = method)
    Sample_list <- row.names(XY1.completecases)
    ProjectIDs = unique(subset(normalized_df, SampleID %in% Sample_list)$ProjectID)
    df_cor <- data.frame(Taxa = i,
                         SampleLocation = Sample.location,
                         ProjectID = paste0(ProjectIDs, collapse = ", "),
                         PerfVariables = j,
                         Observations = nrow(XY1.completecases),
                         Correlation = corr.pval$estimate,
                         Pvalue =  corr.pval$p.value,
                         #AdjPvalue = p.adjust(Pvalue, method = adj_method), 
                         stringsAsFactors = FALSE)
    df_all_cor <- rbind(df_all_cor, df_cor)
  }
}


################################################################################
#consistency correlations
# FILTER BY COLUMN (SIGNIFICANT) ----- 
columna <- "Pvalue"
corr_all_filter <-  df_all_cor[df_all_cor[columna] < 0.05, ]

corr_all_filter_score <- corr_all_filter %>% 
  mutate(Corr.relation = case_when(
    grepl("FCR", PerfVariables) & Correlation > 0 ~ "negative",
    grepl("FCR", PerfVariables) & Correlation <= 0 ~ "positive",
    grepl("BW", PerfVariables) & Correlation > 0 ~ "positive",
    grepl("BW", PerfVariables) & Correlation <= 0 ~ "negative",
  ))


# CREATE CONSISTENCY COLUMN ----- 
df_contexts <- corr_all_filter_score %>% group_by(Taxa)

list_df_contexts <- group_split(df_contexts)
contexts_order <- levels(df_contexts$Taxa)
names(list_df_contexts) <- contexts_order

### loop sobre los contextos para ver consistencia (esta casi que siempre debe cumplirse) ---- 
corr_consist <- data.frame()
for (df_con in list_df_contexts) {
  print(as.character(df_con$SampleLocation[[1]]))
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
        "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/linear_mixed_models/codes and tables/globalcorrelations_feces.rds")


unique(corr_consist$Taxa)


THRESHOLD = 0.75
taxa_score <- data.frame()
list_taxas <- unique(corr_consist$Taxa)
  for (tax in list_taxas) {
    #tax = "Bacteroides"
    df_taxa <- subset(corr_consist, Taxa == tax)
    max_freq <- max(table(df_taxa$Corr.relation))
    max_prop <- max_freq/dim(df_taxa)[1]
    #prop.table(table(df_taxa$Corr.relation))
    if (max_prop < THRESHOLD){next}
    ## quito las inconsistencias del score
    score_sign <- names(table(df_taxa$Corr.relation)[table(df_taxa$Corr.relation) == max_freq])
    df_taxa <- df_taxa[df_taxa$Corr.relation == score_sign,] 
    PerfVariables <- paste0(unique(df_taxa$PerfVariables), collapse = ", ")
    ProjectIDs <- paste0(unique(df_taxa$ProjectID), collapse = ", ")
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
unique(taxa_score$Taxa) 

saveRDS(taxa_score,
        "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/linear_mixed_models/codes and tables/consistencytable_globalcor_genus_ileum.rds")

input_dir <- "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/linear_mixed_models/codes and tables/"
GlobalCorrResults <- readRDS(paste0(input_dir, "consistencytable_globalcor_genus_cecum.rds"))
