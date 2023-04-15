library(phyloseq)
library(zCompositions)
library(compositions)

source("/Users/sebastianbedoyamazo/Documents/siwa_git/Methods-review/functions.R")
source("/Users/sebastianbedoyamazo/Documents/siwa_git/Methods-review/correlations.R")


phyl_dir <- "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/performance_df/pq_obj/"

exp <- "broiler"

ODLEPobj <- readRDS(paste0(phyl_dir, paste0("PhyloseqObject_", exp,".rds")))


pseq.rel <- microbiome::transform(ODLEPobj, "compositional")
colnames(tax_table(pseq.rel))
dfagg_taxa <- aggregate_taxa_siwa(pseq.rel, "Genus") 
class(dfagg_taxa)
dim(dfagg_taxa) ### TAXAS x SAMPLES


meta_exp <- as.data.frame(as.matrix(sample_data(pseq.rel)))
colnames(meta_exp)


# Create a vector with uniques ProjectIDs
project_ids <- unique(meta_exp$ProjectID)
# Create a vector with uniques Ages
ages <- unique(meta_exp$Age)
# Create a vector with uniques SampleLocations
sample_locations <- unique(meta_exp$SampleLocation)


filter_otus <- function(matrix, th=0.8){
  ### matrix: FEATURES as columns
  rowz<-apply(matrix == 0, 2, sum)/dim(matrix)[1] #sum of zeros over columns (OTUS)
  # dejar otus frecuentes, quitando los que están llenos de CEROS
  p <- which(rowz> th) ## p SON LOS QUE SALEN
  matrix[, -p]
}

# Iterate over the unique values of ProjectID, Age and SampleLocation

#empty list to store the correlation results for ProjectID, Ages, and SampleLocation
corr_df <- data_frame()
for (i in project_ids) {
  print(paste0("Running project = ", i))
  df_exp <- meta_exp[meta_exp$ProjectID == i, ]
  for (j in ages) {
    print(paste0("Running Age = ", j))
    if (!j %in% unique(df_exp$Age)){
      print(paste0("no tengo Age ", j))
    } else {
      df_exp_age <- df_exp[df_exp$Age == j, ]
      for (k in sample_locations) {
        print(paste0("Location = ", k))
        # Filtrar el dataframe por las condiciones
        meta_exp_filter <- df_exp_age[df_exp_age$SampleLocation == k, ]
        otu_table <- as.matrix(dfagg_taxa)
        otu_table <- t(otu_table) 
        otu_table <- otu_table[rownames(meta_exp_filter), ]
        otu_table_filtered <- filter_otus(otu_table)
        otu_table_filtered_transformed <- cmultRepl(otu_table_filtered)
        otu_table_filtered_transformed_log <- clr(otu_table_filtered_transformed)
        sel_vars <-
          c( "BW", "BWG", "BWGbefore", "FCR", "FCRbefore")
        # Identificar las columnas que no contienen NaNs en ninguna fila
        print("Columnas con NAN")
        col_nan <- colnames(meta_exp_filter)[colSums(is.na(meta_exp_filter)) > 0]
        print(col_nan)
        # Seleccionar solo las columnas que no tienen NaNs
        sel_vars <- sel_vars[!sel_vars %in% col_nan]
        meta_exp_filtered <- dplyr::select(meta_exp_filter, sel_vars) 
        X <- otu_table_filtered_transformed_log
        dim(X) 
        otus_list <- colnames(X)
        X <- data.frame(X[, colnames(X) %in% otus_list], check.names = FALSE)
        Y <- meta_exp_filtered
        for (col in sel_vars){
          Y[[col]] <- as.numeric(Y[[col]])
        }
        groups_ <- meta_exp_filter$SampleLocation
        #correlation siwa for each ProjectID, Age and SampleLocation
        corr_exp <- correlation_siwa(X, X, groups_, method = "pearson")
        # Store the correlation results for the current age group in the corr_list
        print(dim(corr_exp))
        corr_exp$Age <- j
        corr_exp$ProjectID <- i
        corr_df <- rbind(corr_df, corr_exp)
      }

    }
  } 
}


f = paste0(phyl_dir, paste0("CorrTableGnGn_", exp,".rds"))
f
saveRDS(corr_df, file=f)






colnames(meta_exp)
meta_exp_filter <- subset(meta_exp, Age == "42" & SampleLocation == "cecum" & ProjectID == "E267")


otu_table <- as.matrix(dfagg_taxa)
otu_table <- t(otu_table) 
otu_table <- otu_table[rownames(meta_exp_filter), ]
otu_table_filtered <- filter_otus(otu_table)
otu_table_filtered_transformed <- cmultRepl(otu_table_filtered)
otu_table_filtered_transformed_log <- clr(otu_table_filtered_transformed)

sel_vars <-
  c( "BW", "BWG", "BWGbefore", "FCR", "FCRbefore")

col_nan <- colnames(meta_exp_filter)[colSums(is.na(meta_exp_filter)) > 0]
print(col_nan)
# Seleccionar solo las columnas que no tienen NaNs
sel_vars <- sel_vars[!sel_vars %in% col_nan]
meta_exp_filtered <- dplyr::select(meta_exp_filter, sel_vars) 
X <- otu_table_filtered_transformed_log
dim(X) 
otus_list <- colnames(X)
X <- data.frame(X[, colnames(X) %in% otus_list], check.names = FALSE)
Y <- meta_exp_filtered
for (col in sel_vars){
  Y[[col]] <- as.numeric(Y[[col]])
}



dfmerge <- merge(Y, X, by = "row.names")
str(dfmerge)

cor(dfmerge$BWGbefore, dfmerge$Helicobacter_pullorum)



