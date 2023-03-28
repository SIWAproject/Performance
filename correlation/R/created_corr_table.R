library(phyloseq)
library(zCompositions)
library(compositions)

source("/Users/sebastianbedoyamazo/Documents/siwa_git/Methods-review/functions.R")
source("/Users/sebastianbedoyamazo/Documents/siwa_git/Methods-review/correlations.R")


phyl_dir <- "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/performance_df/pq_obj/"

exp <- "E271"

ODLEPobj <- readRDS(paste0(phyl_dir, paste0("PhyloseqObject_", exp,".rds")))
### filtrar singletons: otus solo en 1 muestra 

pseq.rel <- microbiome::transform(ODLEPobj, "compositional")
colnames(tax_table(pseq.rel))
pseq_genus <- aggregate_taxa_siwa(pseq.rel, "species")
class(pseq_genus)

otu_table <- as.matrix(pseq_genus)
otu_table <- t(otu_table)
rowz<-apply(otu_table == 0, 2, sum)/dim(otu_table)[1] #sum over columns (OTUS)
# dejar otus frecuentes, quitando los que están llenos de CEROS
p <- which(rowz>0.90)
otu_table_filtered <-otu_table[, -p]
colnames(otu_table_filtered)
### transformation 
otu_table_filtered_transformed <- cmultRepl(otu_table_filtered)
otu_table_filtered_transformed_log <- clr(otu_table_filtered_transformed)


meta_exp <- as.data.frame(as.matrix(sample_data(pseq.rel)))
meta_exp <- meta_exp[rownames(otu_table_filtered_transformed_log), ]
sel_vars <-
  c( "BW42","FCR0_42","FCR36_42","BWG0_42","BWG36_42"
  )
sel_vars <- sel_vars[sel_vars %in% colnames(meta_exp)] ### por si no existen la columna
meta_exp_filtered<-meta_exp[, sel_vars]

## ENTRADAS CORRELATION 
X <- otu_table_filtered_transformed_log
#X <- X[, order(colSums(X), decreasing = TRUE)] organizar por abundancia, pero no tiene sentido en log
dim(X) 
otus_list <- colnames(X)
X <- data.frame(X[, colnames(X) %in% otus_list], check.names = FALSE)
Y <- meta_exp_filtered
groups_ <- meta_exp$SampleLocation
for (col in sel_vars){
  Y[[col]] <- as.numeric(Y[[col]])
}

corr_exp <- correlation_siwa(X, Y, groups_,  method="pearson")
corr_exp$projectid <- exp
unique(corr_exp$Taxa)
subset(corr_exp, Taxa == "Alistipes_finegoldii")


f = paste0(phyl_dir, paste0("CorrTableSp_", exp,".rds"))
f
saveRDS(corr_exp, file=f)




# # list_corr_tables_genus <- list() SOLO CORRER PARA EL PRIMER EXP 
# list_corr_tables_genus[[exp]] <- corr_exp_genus
# 
# my_data <- readRDS(f)



