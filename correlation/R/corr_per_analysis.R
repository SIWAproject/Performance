###load libraries panels tabs####
library(phyloseq) #
library(DESeq2) #
library(kableExtra) #
library(genefilter) #
library(microbiome) #
library(ggplot2)  #
library(vegan) #
library(ggpubr) #
library(ggplot2) #
library(plyr)
library(multcompView)
library(RColorBrewer)
library(plotly)

#####Load libraries correlations####
library(ggplot2)
library(phyloseq)
library(stringr)
library(dplyr)

#####Load libraries LR####
library(tidyverse)
library(ape)
library(ggpubr)

####Load libraries categorical histo####
library(kableExtra)
library(plyr)
#library("fantaxtic")
library(data.table)

####load libraries Ratios-histo boxplots####
library(multcompView)
library(reshape)

### OTRICAS
library(plotly)
library(plyr)
library(flexdashboard)
library(shiny)
library(DT)
library(stringr)
library(rstatix)
library(patchwork)
library(readxl)
library(openxlsx)


#library(microViz)
library(lsr)

### summaries ###
library(psych)


#EAFIT surveillance tables (same for all reports)
species_taxonomy_info <-read.csv("/Users/sebastianbedoyamazo/Documents/siwa_git/Reports/Version1.0/Input_data/species_metabolic_effects.csv",check.names = FALSE, sep=";")
genera_taxonomy_info <- read.csv( "/Users/sebastianbedoyamazo/Documents/siwa_git/Reports/Version1.0/Input_data/genus_metabolic_effects.csv", check.names = FALSE, sep=";")


#open phyloseq object for extract performance data
input_dir = "/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/correlation/inputdata/"

ODLEPobj <- readRDS(paste0(input_dir, "PhyloseqObject_perf.rds"))

corr_analysis = read_excel("/Users/sebastianbedoyamazo/Documents/siwa_git/Performance/correlation/inputdata/corr_analysis.xlsx")
#performances summary


#plot analysis
df_species_clr <-read.csv("/Users/sebastianbedoyamazo/Library/CloudStorage/GoogleDrive-sebastian.bedoya@premexcorp.com/.shortcut-targets-by-id/1j_mPOvs6DPPqRINFw74_ARcBZ9HwkQpE/Performance exploration/corr_results/species_otu_table_filtered_transformed_log.csv",check.names = FALSE, sep=",", row.names = 1)
df_genus_clr <-read.csv("/Users/sebastianbedoyamazo/Library/CloudStorage/GoogleDrive-sebastian.bedoya@premexcorp.com/.shortcut-targets-by-id/1j_mPOvs6DPPqRINFw74_ARcBZ9HwkQpE/Performance exploration/corr_results/genus_otu_table_filtered_transformed_log.csv",check.names = FALSE, sep=",", row.names = 1)


#metadata as.data.frame
metadata <- as.data.frame(sample_data(ODLEPobj))
meta_exp <- metadata
colnames(meta_exp)


#anova and tukey for performance and projects
TukeyHSD(aov(meta_exp$BW42 ~ meta_exp$projectid))
TukeyHSD(aov(meta_exp$BWG0_42 ~ meta_exp$projectid))
TukeyHSD(aov(meta_exp$BWG36_42 ~ meta_exp$projectid))
TukeyHSD(aov(meta_exp$FCR0_42 ~ meta_exp$projectid))
TukeyHSD(aov(meta_exp$FCR36_42 ~ meta_exp$projectid))



#boxplot perfomance ~ projects
my_comparisons <- list( c("E345", "E347"), c("E347", "E267"), c("E267", "E271"),
                        c("E345", "E267"),  c("E347", "E271"))
p1 <- ggboxplot(meta_exp, x = "projectid", y = "BW42",
                fill = "projectid")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "anova")+     # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.") 
#geom_jitter(position = position_nudge(x = -0.5))




my_comparisons <- list( c("E345", "E347"), c("E347", "E267"), c("E267", "E271"),
                        c("E345", "E267"),  c("E347", "E271"))
p2 <- ggboxplot(meta_exp, x = "projectid", y = "BWG0_42",
                fill = "projectid")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "anova")+     # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.") 
#geom_jitter(position = position_nudge(x = -0.5))


my_comparisons <- list( c("E345", "E347"), c("E347", "E267"), c("E267", "E271"),
                        c("E345", "E267"),  c("E347", "E271"))
p3 <- ggboxplot(meta_exp, x = "projectid", y = "BWG36_42",
                fill = "projectid")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "anova")+     # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.") 
#geom_jitter(position = position_nudge(x = -0.5))


my_comparisons <- list( c("E345", "E267"), c("E267", "E271"), c("E345", "E271"))
p4 <- ggboxplot(meta_exp, x = "projectid", y = "FCR0_42",
                fill = "projectid")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "anova")+     # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.") 
#geom_jitter(position = position_nudge(x = -0.5))



my_comparisons <- list( c("E345", "E347"), c("E347", "E267"), c("E267", "E271"),
                        c("E345", "E267"),  c("E347", "E271"))
p5 <- ggboxplot(meta_exp, x = "projectid", y = "FCR36_42",
                fill = "projectid")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "anova")+     # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.") 
#geom_jitter(position = position_nudge(x = -0.5))


# Unir los gráficos con patchwork
plot_unique <- p1 + p2 + p3 + p4 + p5
plot_unique <- plot_unique + plot_layout(ncol=3, nrow = 2)

# Imprimir el gráfico unido
plot_unique



#summary by projectid
meta_exp_sum <- meta_exp[, c("projectid", "BW42", "FCR0_42", "FCR36_42", "BWG0_42", "BWG36_42")]
summary_meta_exp <- describe.by(meta_exp_sum, group = meta_exp_sum$projectid)

#tables by projectid
E267 <- summary_meta_exp$E267
E267[-which(rownames(E267) == "projectid*"),][, c("n", "mean", "sd", "min", "max")]

E271 <- summary_meta_exp$E271
E271[-which(rownames(E271) == "projectid*"),][, c("n", "mean", "sd", "min", "max")]

E345 <- summary_meta_exp$E345
E345[-which(rownames(E345) == "projectid*"),][, c("n", "mean", "sd", "min", "max")]


E347 <- summary_meta_exp$E347
E347[-which(rownames(E347) == "projectid*"),][, c("n", "mean", "sd", "min", "max")]





#relation between Genera Score and Performance correlation
#select Genus and score, rename and filter genus in correlation table
library(dplyr)
genera_taxonomy_info_sc <- genera_taxonomy_info[, c("Genus", "Probiotic_potential")]
genera_taxonomy_info_sc <- genera_taxonomy_info_sc %>% rename(Taxa = Genus)

corr_analysis_info_genus <- merge(corr_analysis, genera_taxonomy_info_sc, by = "Taxa", all.x = TRUE)

#complete genus correlation analysis
corr_analysis_info_genus <- filter(corr_analysis_info_genus, Taxa_aggre == "genus")
x = corr_analysis_info_genus[complete.cases(corr_analysis_info_genus), ]
x =x[x$Probiotic_potential < -3,]


#filter by p-value
genus_corr_analysis_filter <-  genus_corr_analysis[genus_corr_analysis$AdjPvalue < 0.05, ]
dim(genus_corr_analysis_filter)
#filter by repeated name
genus_corr_analysis_filter_rep <- genus_corr_analysis_filter %>%
  group_by(Taxa) %>%
  filter(n() > 1)

#### así se le pega una columna con el conteo de repeticiones
# genus_corr_analysis_filter %>%
#   group_by(Taxa) %>%
#   mutate(n = n())

unique(genus_corr_analysis_filter_rep$Taxa)
genus_corr_analysis_filter


# #split positive and negative correlation
# merged_genus_corr_filter <-
#   split(
#     merged_genus_corr_filter,
#     f = ifelse(
#       merged_genus_corr_filter$Correlation >= 0,
#       "positive_cor",
#       "negative_cor"
#     )
#   )




#relation between Species Score and Performance correlation
species_taxonomy_info_sc <- species_taxonomy_info[, c("Species", "Probiotic_potential")]
species_taxonomy_info_sc <- species_taxonomy_info_sc %>% rename(Taxa = Species)

#complete genus correlation analysis
corr_analysis_info_sp <- filter(corr_analysis, Taxa_aggre == "species")

#cambiar " " por "_"
species_taxonomy_info_sc$Taxa <- gsub(" ", "_", species_taxonomy_info_sc$Taxa)


colnames(corr_analysis_info_sp)
colnames(species_taxonomy_info_sc)
dim(corr_analysis_info_sp)
dim(species_taxonomy_info_sc)

sp_corr_analysis <- merge(corr_analysis_info_sp, species_taxonomy_info_sc, by = "Taxa", all.x = TRUE)

#filter by p-value
sp_corr_analysis_filter <-  sp_corr_analysis[sp_corr_analysis$AdjPvalue < 0.10, ]
unique(sp_corr_analysis_filter$Taxa)


#filter by repeated name
sp_corr_analysis_filter_rep <- sp_corr_analysis_filter %>%
  group_by(Taxa) %>%
  filter(n() > 1)





#some graphs
X <- df_species_clr
Y <- meta_exp

bacteria <- "Alistipes_finegoldii"
variable <- "BW42"


#scatterplot 

# Select columns from dataframes

a <- X[, bacteria, drop = F]
b <- Y[, c( variable, "SampleLocation", "projectid"), drop = F]

# merge with rownames
a <- as.data.frame(a)
b <- as.data.frame(b)
b$SampleLocation <- as.factor(b$SampleLocation)
b$projectid <- as.factor(b$projectid)
df_merge <- merge(a, b, by = "row.names")
colnames(df_merge) 
unique(df_merge$projectid)

ggplot(df_merge, aes(x = eval(as.symbol(bacteria)), y = eval(as.symbol(variable)), color = projectid)) +
  facet_wrap(~ SampleLocation) +
  geom_smooth(method='lm', se = FALSE)+
  labs(x='Relative Abundace (CLR)', y='BWG (36-42 days)', title= "Lactobacillus ingluviei") +
  theme(plot.title = element_text(hjust=0.5, size=15, face='bold')) +
  geom_point() +
  labs(color = "Project ID")



df_filtered <- df_merge %>% 
  filter(SampleLocation == "cecum" & projectid == "E271")

ggplot(df_filtered, aes(x = eval(as.symbol(f)), y = eval(as.symbol(p)), color = SampleLocation)) +
  facet_wrap(~ projectid) +
  geom_smooth(method='lm', se = FALSE)+
  labs(x='Relative Abundace (CLR)', y='BW (42 days)', title= "Alistipes_finegoldii") +
  theme(plot.title = element_text(hjust=0.5, size=15, face='bold')) +
  geom_point() +
  labs(color = "Project ID")



```{r correlation method 1, include=FALSE}
df <- filter(corr_analysis_filter, Taxa_aggre == "species")


Env_labeller <- function(variable,value){
  return(sel_vars_label[as.character(value),"Trans"])
}

df$label <- df$Taxa

method="pearson"

```


```{r crear plot 1, include=TRUE,echo=FALSE, warning=FALSE, fig.dim = c(10, 6)}

p <- ggplot(aes(x = Type, y = Taxa, fill = Correlation), data = df)

p <-
  p + geom_tile() + scale_fill_gradient2(low = "#2C7BB6", mid = "white", high =
                                           "#D7191C")
                                         p <-
                                           p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
                                         
                                         p <-
                                           p + geom_text(aes(label = Significance), color = "black", size = 3) + labs(y =
                                                                                                                        NULL, x = NULL, fill = method)
                                         p <-
                                           p + facet_grid(
                                             . ~ Env,
                                             drop = TRUE,
                                             scale = "free",
                                             space = "free_x",
                                             labeller = Env_labeller
                                           )+ theme(text=element_text(size=9))
                                         
                                         p
                                         ```





