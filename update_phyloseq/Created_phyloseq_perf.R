# Performance phyloseq object


#####Load libraries#####
library(phyloseq)
library(tibble)
library(dplyr)

#####Open data#####
input_dir<- "/Users/sebastianbedoyamazo/Documents/siwa_git/swine_data/para_merge/"
taxonomy <- read.table(paste0(input_dir,"taxonomy_table_E350.csv"), header=T, sep=",", check.names = F)
features <- read.table(paste0(input_dir,"otu_table_E350.csv"), header=T, sep=",", check.names = F)
metadata <- read.table(paste0(input_dir,"metadata_E350_pen_perf.csv"), header=T, sep=",", check.names = F) 

####subset metadata ····
colnames(metadata)
length(unique(metadata$sampleid))
metadata <- dplyr::select(metadata, sampleid, samplelocation, alphashannon,
                          alphaobserved, ProjectID, Age, treatment, Trt,
                          animaltype, animalnumber, Pen, House, Sex,
                          BW, ADGbefore, ADFIbefore, FCRbefore)

metadata <- metadata %>% rename(SampleID = sampleid, SampleLocation = samplelocation, 
                                TreatmentNumber = Trt, AlphaShannon = alphashannon,
                                AlphaObserved = alphaobserved, Treatment = treatment,
                                AnimalType = animaltype, AnimalNumber = animalnumber)




#####Organize data#####
#check repeated columns
colnames(metadata)
unique(metadata$ProjectID)
#Organize metadata
metadata <- metadata[metadata$SampleID %in% colnames(features),]
metadata$SampleLocation <- as.factor(metadata$SampleLocation)
metadata$Treatment <- as.factor(metadata$Treatment)

#Make OTUs row names in otutable and taxo 
rownames(metadata) <- metadata$SampleID


features <- features %>% column_to_rownames(var = "OTU")
features <- features[, intersect(rownames(metadata), colnames(features))]


taxonomy <- taxonomy %>% column_to_rownames(var = "OTU")
taxonomy <- taxonomy[rownames(taxonomy) %in% rownames(features), ]


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
exp <-"E350"
f = paste0(input_dir, paste0("PhySwine_", exp,".rds"))
f
saveRDS(ODLEPobj, file=f)


###merge phyloseq ------
#llamar phy actual
exp <- "E350"
PhySwine_E350 <- readRDS(paste0(input_dir, paste0("PhySwine_", exp,".rds")))

#llamar phy ultima actualizacion
exp <- "14112023"
PhySwine_14112023 <- readRDS(paste0(input_dir, paste0("PhySwine_", exp,".rds")))


#hacer merge de phy actual con phy ultima actualizacion
physeq_merged <- merge_phyloseq(physeq_merged, PhySwine_E350)


#Guardar nueva version del phy, siempre poner la fecha de actualizacion en el nombre
exp <-"14112023"
f = paste0(input_dir, paste0("PhySwine_", exp,".rds"))
f
saveRDS(physeq_merged, file=f)
