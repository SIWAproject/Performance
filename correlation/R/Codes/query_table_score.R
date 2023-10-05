corr_dir <- "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/correlation/results/"
corr_table_G <- readRDS(paste0(corr_dir, "GeneraCorr_120723.rds"))


x = corr_table_G[corr_table_G$Taxa == "[Eubacterium]_brachy_group", ]
### 1. filtro solo los significativos de acuerdo a la columns que escoja 
colnames(corr_table_G)
column = "Pvalue" #AdjPvalue - más estricto

corr_table_G_filter <-  corr_table_G[corr_table_G[column] < 0.05, ]
unique(corr_table_G_filter$Taxa) #159 TAXAS


unique(corr_table_G_filter$Taxa) #156 TAXAS

sample_locations <- c("C", "I", "F")
THRESHOLD = 0.75
taxa_score <- data.frame()
for (loc in sample_locations) {
  print(loc)
  ## SACAR INCONSISTENTES DENTRO DE LOS MISMOS CONTEXTOS
  ### quitar las BACTERIAS que tengan inconsistencias dentro del contexto
  corr_table_G_filter_loc <- corr_table_G_filter[corr_table_G_filter$Type == loc,]
  corr_table_G_filter_loc <-
    corr_table_G_filter_loc[!corr_table_G_filter_loc$Taxa %in% 
                              corr_table_G_filter_loc[corr_table_G_filter_loc$ContCons == "NC",]$Taxa,]
  df_loc <- corr_table_G_filter_loc %>%
    group_by(Taxa, Context) %>%
    summarize(PerfVariables = paste(Env, collapse = ", "),
              ProjectID = first(ProjectID),
              Corr.relation = first(Corr.relation)) %>%
    ungroup()
  list_taxas <- unique(df_loc$Taxa)
  for (tax in list_taxas) {
    #tax = "Bacteroides"
    df_taxa <- subset(df_loc, Taxa == tax)
    max_freq <- max(table(df_taxa$Corr.relation))
    max_prop <- max_freq/dim(df_taxa)[1]
    #prop.table(table(df_taxa$Corr.relation))
    if (max_prop < THRESHOLD){next}
    ## quito las inconsistencias del score
    score_sign <- names(table(df_taxa$Corr.relation)[table(df_taxa$Corr.relation) == max_freq])
    df_taxa <- df_taxa[df_taxa$Corr.relation == score_sign,] 
    Projects <- length(unique(df_taxa$ProjectID))
    Context <- length(unique(df_taxa$Context)) 
    ContextsID <- paste0(unique(df_taxa$Context), collapse = ", ")
    PerfVariables <- paste0(unique(df_taxa$PerfVariables), collapse = ", ")
    PropConsistency <- max_prop
    df_row <- data.frame(Taxa = tax,
                         SampleLocation = loc,
                         PerfVariables = PerfVariables,
                         Contexts = Context,
                         ContextsID = ContextsID,
                         Projects = Projects,
                         PropConsistency= max_prop,
                         Sign = score_sign,
                         stringsAsFactors = FALSE)
    taxa_score <- rbind(taxa_score, df_row)
  }
}
unique(taxa_score$Taxa) ## 152 taxas


# taxa_score --> 267 rows

taxa_score <- subset(taxa_score, Taxa != "Incertae_Sedis")
saveRDS(taxa_score, "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/correlation/results/query_table_score_190723.rds")


contexts <- unlist(str_split(taxa_score[taxa_score$SampleLocation == "F",]$ContextsID, ","))
contexts <- str_replace(contexts, " ", "")
length(unique(contexts))


# https://www.notion.so/GENERA-PERFORMANCE-545f63fd2c8b429193bb40ebac4c9a75



