library(psych)

source("/Users/sebastianbedoyamazo/Documents/siwa_git/Methods-review/functions.R")
source("/Users/sebastianbedoyamazo/Documents/siwa_git/Methods-review/correlations.R")


phyl_dir <- "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/performance_df/pq_obj/"

exp <- "broiler"

ODLEPobj <- readRDS(paste0(phyl_dir, paste0("PhyloseqObject_", exp,".rds")))

meta_exp <- as.data.frame(as.matrix(sample_data(pseq.rel)))
colnames(meta_exp)

# Create a vector with uniques ProjectIDs
project_ids <- unique(meta_exp$ProjectID)
# Create a vector with uniques Ages
ages <- unique(meta_exp$Age)


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
      
      sel_vars <-
        c( "BW", "BWG", "BWGbefore", "FCR", "FCRbefore")
      
      sel_vars <- sel_vars[sel_vars %in% colnames(df_exp_age)] ### por si no existen la columna
      meta_exp_filtered<-df_exp_age[, sel_vars]
      
      ## ENTRADAS summary
      Y <- meta_exp_filtered
      for (col in sel_vars){
        Y[[col]] <- as.numeric(Y[[col]])
      }
      summ_exp <- describe(Y , na.rm = TRUE)
      summ_exp$Age <- j
      summ_exp$ProjectID <- i
      
      corr_df <- rbind(corr_df, summ_exp)
      
    }
  }
}





meta_exp_filter <- subset(meta_exp, SampleLocation == "cecum")

sel_vars <-
  c( "BW", "BWG", "BWGbefore", "FCR", "FCRbefore")

sel_vars <- sel_vars[sel_vars %in% colnames(meta_exp_filter)] ### por si no existen la columna
meta_exp_filtered<-meta_exp_filter[, sel_vars]

## ENTRADAS summary
Y <- meta_exp_filtered
for (col in sel_vars){
  Y[[col]] <- as.numeric(Y[[col]])
}

exp_summ = describe.by(Y, group = meta_exp_filter$Age, na.rm = T)

age_14_summary <- exp_summ$'14' %>%
  as.data.frame() %>%
  select(n, mean, sd, min, max)%>%
  round(digits = 2)


age_21_summary <- exp_summ$'21' %>%
  as.data.frame() %>%
  select(n, mean, sd, min, max)%>%
  round(digits = 2)


age_28_summary <- exp_summ$'28' %>%
  as.data.frame() %>%
  select(n, mean, sd, min, max)%>%
  round(digits = 2)


age_42_summary <- exp_summ$'42' %>%
  as.data.frame() %>%
  select(n, mean, sd, min, max)%>%
  round(digits = 2)
