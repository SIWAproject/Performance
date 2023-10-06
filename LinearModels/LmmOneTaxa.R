library(phyloseq)
library(zCompositions)
library(compositions)
library(microbiomeMarker)
library(glmnet)
library(psych)#summaries
library(lme4)

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
Sample.location <- "F"
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

# a = otu_table_project$E267 
# a = otu_table_project$E271 
# a = otu_table_project$E321 
# a = otu_table_project$E325 
# a = otu_table_project$E326 
# a = otu_table_project$E335 
# a = otu_table_project$E345 
# a = otu_table_project$E347 
# a = otu_table_project$Fieldale 

# Loop para convertir los row names en una columna en cada dataframe
for (i in 1:length(otu_table_project)) {
  otu_table_project[[i]]$SampleID <- rownames(otu_table_project[[i]])
  #rownames(otu_table_project[[i]]) <- NULL
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
meta_all_cecum <- dplyr::select(meta_all_cecum, sel_vars, ProjectID, SampleLocation, Age)
XY <- merge(otu_table_all_cecum, meta_all_cecum, by = "row.names", all = TRUE)
row.names(XY) <- XY$Row.names
XY$Row.names <- NULL

saveRDS(XY,
        "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/ILUE_03102023.rds")


################################################################################
#LINEAR MIXED MODEL
XY$ProjectID_Age <- paste0(XY$ProjectID, "_", XY$Age)
XY$SampleID <- row.names(XY)

colnames(XY)

list_otus <- colnames(XY)[1:120]

new_otus_names <- as.data.frame(paste("Genus", 1:120, sep = "_"))

taxa_names <- cbind(list_otus, new_otus_names)
colnames(taxa_names) <- c("Genus", "NewGenus")
list_taxa_reg <- taxa_names$NewGenus

colnames(XY)[1:120] <- list_taxa_reg

sel_vars <- colnames(XY)[121:125]
sample_locations <- unique(XY$SampleLocation)



# Definir la función para calcular RMSE
calculate_rmse <- function(observed, predicted) {
  sqrt(mean((observed - predicted)^2))
}


results <- data.frame()
for (i in list_taxa_reg) {
  for (j in sel_vars) {
    for (k in sample_locations) {
      XY_filter <- subset(XY, SampleLocation == k) 
      XY_filter <- XY_filter[complete.cases(XY_filter[, c(i, j)]), ]
      if (nrow(XY_filter) == 0) {
        cat("XY_filter has no rows for", paste0(i, "_", j), "\n")
        next  # Saltar a la siguiente iteración del bucle
      }
      
      if (length(unique(XY_filter$ProjectID_Age)) == 1) {
        cat("There is just one category in the random variable for", paste0(i, "_", j), "\n")
        next  # Saltar a la siguiente iteración del bucle
      }
      equation <- as.formula(paste(j, "~", i, "+ (1 | ProjectID_Age)"))
      model <- lmer(equation, data = XY_filter, na.action = na.omit)
      sum_mod <- summary(model)
      
      #unique projects and n observations
      Sample_list <- row.names(XY_filter)
      ProjectIDs = unique(subset(XY_filter, SampleID %in% Sample_list)$ProjectID)
      
      mod_coeff <- sum_mod$coefficients
      mod_coeff <- as.data.frame(mod_coeff)
      mod_coeff$bic <- BIC(model)
      mod_coeff$NewGenus <- i
      mod_coeff$Observations <- nrow(XY_filter)
      mod_coeff$perf_var <- j
      mod_coeff$ProjectID <- paste0(ProjectIDs, collapse = ", ")
      mod_coeff$SampleLocation <- k
      mod_coeff$variables <- rownames(mod_coeff)
      
      # Calcular RMSE
      observed <- XY_filter[, j]  # Los valores observados
      predicted <- fitted(model)   # Los valores predichos por el modelo
      rmse <- calculate_rmse(observed, predicted)
      mod_coeff$rmse <- rmse
      
      
      results <- rbind(results, mod_coeff)
      
    }
  }
}

colnames(results)
colnames(taxa_names)
results <- merge(results, taxa_names, by.x = "NewGenus", by.y = "NewGenus")
results$RandomVar <- "ProjectID_Age"

input_dir <- "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/linear_mixed_models/"
lmm_results <- readRDS(paste0(input_dir, "LinearReg_feces_results_agefix.rds"))


################################################################################
#consistency coefficients
# FILTER BY COLUMN (SIGNIFICANT) ----- 
colnames(lmm_results)
columna <- "Pr(>|t|)"
lmm_results_filter <-  lmm_results[lmm_results[columna] < 0.05, ]
lmm_results_filter <-  subset(lmm_results_filter, variables != "(Intercept)")
lmm_results_filter <- subset(lmm_results_filter, !grepl("Age", variables))

lmm_results_score <- lmm_results_filter %>% 
  mutate(Coef.relation = case_when(
    grepl("FCR", perf_var) & Estimate > 0 ~ "negative",
    grepl("FCR", perf_var) & Estimate <= 0 ~ "positive",
    grepl("BW", perf_var) & Estimate > 0 ~ "positive",
    grepl("BW", perf_var) & Estimate <= 0 ~ "negative",
  ))


# CREATE CONSISTENCY COLUMN ----- 
df_contexts <- lmm_results_score %>% group_by(Genus)
unique(lmm_results_score$Genus)

list_df_contexts <- group_split(df_contexts)
contexts_order <- levels(df_contexts$Genus)
names(list_df_contexts) <- contexts_order

### loop sobre los contextos para ver consistencia (esta casi que siempre debe cumplirse) ---- 
reg_consist <- data.frame()
for (df_con in list_df_contexts) {
  print(as.character(df_con$Genus[[1]]))
  list_taxas <- unique(df_con$Genus)
  print(length(list_taxas)) 
  for (b in list_taxas) {
    df_exp_bacteria <- df_con[df_con$Genus == b,]
    df_exp_bacteria$ContCons <-
      ifelse(length(unique(df_exp_bacteria$Coef.relation)) == 1,  "C", "NC")
    reg_consist <- rbind(reg_consist, df_exp_bacteria)
  }
}

# saveRDS(reg_consist,
#         "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/linear_mixed_models/codes and tables/CoeffRegresion_projectIDAge_ileum.rds")



unique(reg_consist$Genus)


THRESHOLD = 0.75
taxa_score <- data.frame()
list_taxas <- unique(reg_consist$Genus)
for (tax in list_taxas) {
  #tax = "Bacteroides"
  df_taxa <- subset(reg_consist, Genus == tax)
  max_freq <- max(table(df_taxa$Coef.relation))
  max_prop <- max_freq/dim(df_taxa)[1]
  #prop.table(table(df_taxa$Corr.relation))
  if (max_prop < THRESHOLD){next}
  ## quito las inconsistencias del score
  score_sign <- names(table(df_taxa$Coef.relation)[table(df_taxa$Coef.relation) == max_freq])
  df_taxa <- df_taxa[df_taxa$Coef.relation == score_sign,] 
  PerfVariables <- paste0(unique(df_taxa$perf_var), collapse = ", ")
  ProjectIDs <- paste0(unique(df_taxa$ProjectID), collapse = ", ")
  PropConsistency <- max_prop
  df_row <- data.frame(Taxa = tax,
                       ProjectID = ProjectIDs,
                       SampleLocation = unique(df_taxa$SampleLocation),
                       PerfVariables = PerfVariables,
                       PropConsistency= max_prop,
                       Sign = score_sign,
                       stringsAsFactors = FALSE)
  taxa_score <- rbind(taxa_score, df_row)
}
unique(taxa_score$Taxa) 


# Función para eliminar duplicados y mantener valores únicos en columna de ProjectID
eliminar_duplicados <- function(string) {
  unique_strings <- unique(unlist(strsplit(string, ", ")))
  return(paste(unique_strings, collapse = ", "))
}

# Aplicar la función a cada fila de la columna PRojectID
taxa_score$ProjectID <- apply(taxa_score, 1, function(row) {
  eliminar_duplicados(row["ProjectID"])
})


saveRDS(taxa_score,
        "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/linear_mixed_models/codes and tables/consistencytable_lmm_onetaxa_genus_feces.rds")
