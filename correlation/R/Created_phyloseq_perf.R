# Performance phyloseq object


#####Load libraries#####
library(phyloseq)
library(tibble)
library(dplyr)

#####Open data#####
input_dir<- "/Users/sebastianbedoyamazo/Documents/siwa_git/performance/performance_df/pq_obj/"
taxonomy <- read.table(paste0(input_dir,"taxonomy_broiler.csv"), header=T, sep=",", check.names = F)
features <- read.table(paste0(input_dir,"otu_table_broiler.csv"), header=T, sep=",", check.names = F)
metadata <- read.table(paste0(input_dir,"metadata_broiler.csv"), header=T, sep=",", check.names = F) 

#####Organize data#####
#check repeated columns
colnames(metadata)
#metadata <- metadata[ , -8]
#metadata <- metadata %>% rename(SampleID = sampleid, SampleLocation = samplelocation, Treatment = Trt)
# taxonomy <- taxonomy %>% rename(OTU = otu)
# features <- features %>% rename(OTU = otu)

# #delete letter X from colnames
# colnames(features) <- gsub("X", "", colnames(features))

#Organize metadata
metadata <- metadata[metadata$SampleID %in% colnames(features),]
metadata$SampleLocation <- as.factor(metadata$SampleLocation)
metadata$Treatment <- as.factor(metadata$Trt)

#Make OTUs row names in otutable and taxo 
rownames(metadata) <- metadata$SampleID


# metadata <- metadata %>% 
#   tibble::column_to_rownames(var = "SampleID") %>% 
#   tibble::add_column(SampleID = rownames(.), .before = 1)

# features <- features %>% 
#   tibble::column_to_rownames(var = "OTU") %>% 
#   tibble::add_column(OTU = rownames(.), .before = 1)


features <- features %>% column_to_rownames(var = "OTU")
# dim(features)


# taxonomy <- taxonomy %>% 
#   tibble::column_to_rownames(var = "OTU") %>% 
#   tibble::add_column(otu = rownames(.), .before = 1)

taxonomy <- taxonomy %>% column_to_rownames(var = "OTU")


taxonomy <- taxonomy[order(row.names(taxonomy)), ]
features <- features[order(row.names(features)), ]
metadata <- metadata[order(metadata$SampleID), ]
idx <- order(metadata$SampleID)
features <- features[, idx]



dim(features)
featuresM <- as(features, "matrix")
class(featuresM) <- "numeric"
taxa_matrix<-as.matrix(taxonomy)
colnames(taxonomy)

#Check dims:
dim(features) #lines X, cols Y
dim(metadata) #lines Y, cols Z
dim(taxonomy) #lines X, cols W


#####Build Phyloseq Object#####
ASVobj = otu_table(featuresM, taxa_are_rows = TRUE)  
TAXobj = tax_table(taxa_matrix)                        
SAMobj = sample_data(metadata)
ODLEPobj = phyloseq(ASVobj, TAXobj, SAMobj)     

#####Review Phyloseq object#####
rank_names(ODLEPobj)
sample_variables(ODLEPobj)
ntaxa(ODLEPobj)


#####Save phyloseq object#####
exp <-"broiler"
f = paste0(input_dir, paste0("PhyloseqObject_", exp,".rds"))
f
saveRDS(ODLEPobj, file=f)

