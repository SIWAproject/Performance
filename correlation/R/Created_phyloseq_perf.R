# Performance phyloseq object


#####Load libraries#####
library(phyloseq)
library(tibble)

#####Open data#####
input_dir<- "/Users/sebastianbedoyamazo/Documents/siwa/siwa2022/performance/performance_df/pq_obj/"
taxonomy <- read.table(paste0(input_dir,"taxonomy_E347.csv"), header=T, sep=",", check.names = F)
features <- read.table(paste0(input_dir,"otu_table_E347.csv"), header=T, sep=",", check.names = F)
metadata <- read.table(paste0(input_dir,"metadata_E347.csv"), header=T, sep=",", check.names = F) 

#####Organize data#####
metadata <- metadata %>% rename(SampleID = sampleid, SampleLocation = samplelocation, Treatment = Trt)
taxonomy <- taxonomy %>% rename(otu = OTU)
features <- features %>% rename(otu = OTU)

#delete letter X from colnames
colnames(features) <- gsub("X", "", colnames(features))

#Organize metadata
metadata <- metadata[metadata$SampleID %in% colnames(features),]
metadata$SampleLocation <- as.factor(metadata$SampleLocation)
metadata$Treatment <- as.factor(metadata$Treatment)

#Make OTUs row names in otutable and taxo 
metadata <- metadata %>% 
  tibble::column_to_rownames(var = "SampleID") %>% 
  tibble::add_column(SampleID = rownames(.), .before = 1)

# features <- features %>% 
#   tibble::column_to_rownames(var = "OTU") %>% 
#   tibble::add_column(OTU = rownames(.), .before = 1)


features <- features %>% column_to_rownames(var = "otu")
# dim(features)


taxonomy <- taxonomy %>% 
  tibble::column_to_rownames(var = "otu") %>% 
  tibble::add_column(otu = rownames(.), .before = 1)


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
saveRDS(ODLEPobj, paste0(input_dir,"PhyloseqObject_E347.rds"))
