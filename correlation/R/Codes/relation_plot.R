#---El objetivo de este cod es graficar la relación entre bacterias y variables
#--- de performance.
#--- pasos:
# 1. Llamar phyloseq (versión de interes)
# 2. Llamar metadatos con performance as.data.frame(as.matrix(sample_data(ODLEPobj)))
# 3. 

library(phyloseq)
library(zCompositions)
library(compositions)
library(microbiomeMarker)
library(glmnet)
library(psych)#summaries

cvi_colours = list(
  cvi_siwa = c("#03343a", "#4e8e74","#f99b35",  "#e5c217",  
                        "#075b44", "#f9b870", "#f7e76d", 
                        "#017fb1", "#5cb08e" , "#fcd8b6", "#fcf5cd", "#ABF4D4",
                        "#8CDBF4","#F7927F", "red"),
                        
  
  samplelocation_colors = c( "#f99b35",  "#017fb1"),
  age_colors=c("#6732a8"),
  treatment_colors = c( "#4e8e74", "#e5c217", "#f99b35", "#017fb1", "#F7927F"), 
  groups=c("#4e8e74", "#035060", "#f99b35", "#BC8808")
)

cvi_palettes = function(name, n, all_palettes = cvi_colours, type = c("discrete", "continuous")) {
  palette = all_palettes[[name]]
  if (missing(n)) {
    n = length(palette)
  }
  type = match.arg(type)
  out = switch(type,continuous = grDevices::colorRampPalette(palette)(n),discrete = palette[1:n]
  )
  structure(out, name = name, class = "palette")
}

scale_color_cvi_d = function(name) {
  ggplot2::scale_colour_manual(values = cvi_palettes(name, type = "discrete"))
}
scale_fill_cvi_d = function(name) {
  ggplot2::scale_fill_manual(values = cvi_palettes(name,type = "discrete"))
}




source("/Users/sebastianbedoyamazo/Documents/siwa_git/Methods-review/functions.R")
source("/Users/sebastianbedoyamazo/Documents/siwa_git/Methods-review/correlations.R")

file_path <- "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/performance_tables/excel/DE_1_2_6.xlsx"

DE_context <- openxlsx::read.xlsx(file_path, sheet = "cons_table", colNames = TRUE)

input_dir <- "/Users/sebastianbedoyamazo/Documents/siwa_git/Pendientes/"
exp <- "Broiler"
ODLEPobj <- readRDS(paste0(input_dir, paste0("PhyloseqObject_", exp,".rds")))
meta_all <- as.data.frame(as.matrix(sample_data(ODLEPobj)))


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
  phy_exp <- subset_samples(phy_exp, SampleLocation== Sample.location)
  phy_exp <- prune_taxa(taxa_sums(phy_exp) > 0, phy_exp)
  dfagg_exp <- aggregate_taxa_siwa(phy_exp, "Genus") 
  otu_table <- as.matrix(as.data.frame(dfagg_exp))
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


# Loop para convertir los row names en una columna en cada dataframe
for (i in 1:length(otu_table_project)) {
  otu_table_project[[i]]$SampleID <- rownames(otu_table_project[[i]])
  rownames(otu_table_project[[i]]) <- NULL
}

# Inicializar el DataFrame resultado con el primer DataFrame de la lista
otu_table_all <- otu_table_project[[1]]
dim(otu_table_project[[1]])

# Loop para realizar full join con cada dataframe en otu_table_project
for (i in 2:length(otu_table_project)) {
  common_columns <- intersect(colnames(otu_table_all), colnames(otu_table_project[[i]]))
  otu_table_all <- otu_table_all %>%
    full_join(otu_table_project[[i]], by = common_columns)
}


dim(otu_table_all)

otu_table_all$SampleID
row.names(otu_table_all) <-  otu_table_all$SampleID
otu_table_all$SampleID <- NULL

#filtrar meta_all
colnames(meta_all)
meta_all[, 10:14] <- lapply(meta_all[, 10:14], as.numeric)
str(meta_all)

meta_all_location <- meta_all[meta_all$SampleID %in% rownames(otu_table_all),]
sel_vars <-c( "BW", "BWG", "BWGbefore", "FCR", "FCRbefore")
meta_all_location <- dplyr::select(meta_all_location, sel_vars, ProjectID, SampleLocation, Age)
XY <- merge(otu_table_all, meta_all_location, by = "row.names", all = TRUE)
row.names(XY) <- XY$Row.names
XY$Row.names <- NULL
colnames(XY)

XY$Context <- paste0(XY$ProjectID, "-", XY$Age, "-", XY$SampleLocation)



perf_var <- "BW"
bacteria <- "RF39"
col_by <- "ProjectID"
Context_list <- unlist(strsplit(subset(DE_context, Genus == bacteria & SampleLocation == Sample.location)$ContextIDs_ContextCorr, ", "))


#plot de complete data
ggplot(XY, aes(x = eval(as.symbol(bacteria)), y = eval(as.symbol(perf_var)))) +
  geom_point()+ geom_smooth(method="lm", se = FALSE)+
  #facet_wrap(~ data_transf, scales = "free") + #"free" for autoscale
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y = paste0(perf_var), x = paste0(bacteria, " (RA, clr transformation)", fill = paste0(col_by))) 


#plot complete data col by categorical variable
ggplot(XY, aes(x = eval(as.symbol(bacteria)), y = eval(as.symbol(perf_var)), col = eval(as.symbol(col_by)))) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    y = paste0(perf_var),
    x = paste0(bacteria, " (RA, CLR transformation)"),
    col = paste0(col_by)
  ) 

col_by <- "Context"
#plot context data col by categorical variable
ggplot(subset(XY, Context %in% Context_list), aes(x = eval(as.symbol(bacteria)), y = eval(as.symbol(perf_var)), col = eval(as.symbol(col_by)))) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    y = paste0(perf_var),
    x = paste0(bacteria, " (RA, CLR transformation)"),
    col = paste0(col_by)
  ) 


