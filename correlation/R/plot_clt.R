library(phyloseq)
library(zCompositions)
library(compositions)
library(gridExtra)
library(ggpmisc)


source("/Users/sebastianbedoyamazo/Documents/siwa_git/Methods-review/functions.R")
source("/Users/sebastianbedoyamazo/Documents/siwa_git/Methods-review/correlations.R")


phyl_dir <- "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/performance_df/pq_obj/"

exp <- "broiler"
ODLEPobj <- readRDS(paste0(phyl_dir, paste0("PhyloseqObject_", exp,".rds")))

cor_filter_table <- readRDS(paste0(phyl_dir, paste0("CorrTableFilter_broiler.rds")))
colnames(cor_filter_table)

corr_genus <- subset(cor_filter_table, Taxonomy == "Genus" & Score == "positive")
unique(corr_genus$Taxa)

corr_specie <- subset(cor_filter_table, Taxonomy == "Species" & Score == "positive")
unique(corr_specie$Taxa)


pseq.rel <- microbiome::transform(ODLEPobj, "compositional")
colnames(tax_table(pseq.rel))
dfagg_taxa <- aggregate_taxa_siwa(pseq.rel, "Species") 
class(dfagg_taxa)
dim(dfagg_taxa) ### TAXAS x SAMPLES


meta_exp <- as.data.frame(as.matrix(sample_data(pseq.rel)))
colnames(meta_exp)

filter_otus <- function(matrix, th=0.8){
  ### matrix: FEATURES as columns
  rowz<-apply(matrix == 0, 2, sum)/dim(matrix)[1] #sum of zeros over columns (OTUS)
  # dejar otus frecuentes, quitando los que están llenos de CEROS
  p <- which(rowz> th) ## p SON LOS QUE SALEN
  matrix[, -p]
}


ages = "28"
locations = "ileum"
projects = "E271"

taxa = "Lactobacillus_coleohominis"
perf_var = "BWGbefore"
perf_var <- c("BW", "BWGbefore")


meta_exp_filter <- subset(meta_exp, Age == ages & SampleLocation == locations & ProjectID == projects)
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


### revisar con transformación
a <- X[taxa]
class(a)
a <- a %>% rename(clr_t = taxa)

b <- Y[perf_var]
class(b)

otu_table_filtered <- as.data.frame(otu_table_filtered)
c <- otu_table_filtered[taxa]
c <- c %>% rename(rawdata = taxa)

d <- otu_table_filtered_transformed[taxa]
d <- d %>% rename(cmultRepl_t = taxa)


df_plot <- cbind(b,c,d,a)
cor(b,a)
cor(b,c)

# convert to long format
df_long <- df_plot %>%
  pivot_longer(cols = -perf_var, names_to = c("data_transf"))


ggplot(df_long, aes(x = value, y = eval(as.symbol(perf_var)))) +
  geom_point()+ geom_smooth(method="lm", se = FALSE)+
  facet_wrap(~ data_transf, scales = "free") + #"free" for autoscale
  stat_poly_line(formula = y ~ x) +
  stat_poly_eq(mapping = use_label(c("R2", "eq")))  +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Relation y ~ x for Different Data Transformation") +
  labs(y = perf_var, x = paste0(taxa, "(Relative Abundance)")) 










