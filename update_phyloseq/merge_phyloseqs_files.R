#####Load libraries#####
library(phyloseq)
library(tibble)
library(dplyr)

###merge phyloseq ------
### Hacer merge de los phyoseq con data de performance, los pasos son:
# 1. Llamar phyloseq por experimento
# 2. Hacer merge
# 3. Guardar phyloseq


# 1. Directorio de phyloseq
input_dir <- "/Users/sebastianbedoyamazo/Documents/siwa_git/swine_data/phyloseq/phy_files"
files <- list.files(path = input_dir, full.names = TRUE)  # Obtener rutas completas


phy_list <- list()

# 2. Loop para cargar y renombrar los phyloseq presentes en la carpeta
for (i in 1:length(files)) {
  nombre_archivo <- paste0("Phy", i)
  
  # Leer y asignar el archivo como un objeto en el entorno de trabajo
  assign(nombre_archivo, readRDS(files[i])) 
  
  # Imprimir un mensaje para verificar que se ha cargado el archivo
  cat("Cargado phyloseq:", nombre_archivo, "\n")
  phy_list <- c(phy_list, nombre_archivo)
}

# 3. Iterar a través de la lista de los phyloseq para hacer merge 
# sobre el Phy1, osea el primer phyloseq file que llamo el loop anterior
for (i in phy_list[-1]) {
  # Acceder al DataFrame por su nombre

  Phy_n <- get(i)

  # Combinar el DataFrame actual con el resultado final utilizando rbind
  Phy1 <- merge_phyloseq(Phy1, Phy_n)
}



# 4. Guardar nuevo complete phyloseq con fecga de creación

phy_dir <- "/Users/sebastianbedoyamazo/Documents/siwa_git/swine_data/phyloseq/phy_files/"
date <-"14112023"
f = paste0(phy_dir, paste0("PhySwine_", date,".rds"))
f
saveRDS(Phy1, file=f)
