####################################################################################################
# SPTIAL OBJECT CREATION WITH MANUAL AREA ANNOTATION
#####################################################################################################
rm(list=ls())

library(Seurat)
library(dplyr)
library(ggplot2)

#Variables---------------------------------
DIR_ROOT <- file.path(getwd())
DIR_DATA <- file.path(DIR_ROOT, "data/data/")

#Data--------------------------------------

samples <- dir(path = DIR_DATA)
samples <- samples[! samples %in% c("HGP study 1 metadata 20230509.xlsx","zip")]

# Create each individual Seurat object for every sample
lista <- c(samples)
names(lista) <- samples
prueba <-c()

for (i in lista){
  a <- Load10X_Spatial(data.dir = paste0(DIR_DATA, i),
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial",
                       slice = i,
                       filter.matrix = TRUE)
  
  ######add area data
  area_a <- read.csv(paste0(DIR_DATA, i, "/area.csv"))
  area_a <- as.vector(area_a$area)
  a@meta.data["area"] <- as.factor(area_a)
  
  ######add replicate data
  replicate_a <- read.csv(paste0(DIR_DATA, i, "/replicate.csv"))
  replicate_a <- as.vector(replicate_a$replicate)
  a@meta.data["replicate"] <- as.factor(replicate_a)
  
  # Check if additional file exists (sample info)
  extra_file <- paste0(DIR_DATA, i, "/sample.csv")
  if (file.exists(extra_file)) {
    # Read and add the extra file data
    sample_a <- read.csv(extra_file)
    sample_a <- as.vector(sample_a$sample)
    a@meta.data["sample"] <- as.factor(sample_a)
  }
  
  ######subset data
  Seurat::Idents(object = a) <- a@meta.data[["area"]]
  a <- subset(x =a, idents = c("HZ", "PZ", "RZ", "SOC"))
  a@meta.data[["area"]] <- a@active.ident
  
  Seurat::Idents(object = a) <- a@meta.data[["replicate"]]
  ######add object to list
  prueba[[length(prueba) + 1]] <- a

}
names(prueba) <- samples

###################################################OPTION1####################################
####Separate data in replicates
## separate list
list2env(prueba,envir=.GlobalEnv)

#separate in replicates
`V11B18-363_B1`@meta.data[["orig.ident"]] <- "DQ9500"
DQ9500_1 <- subset(x = `V11B18-363_B1` , idents = c("one"))
DQ9500_1@meta.data[["id"]] <- "DQ9500_1"
DQ9500_2 <- subset(x = `V11B18-363_B1` , idents = c("two"))
DQ9500_2@meta.data[["id"]] <- "DQ9500_2"
`V11B04-032_C1`@meta.data[["orig.ident"]] <- "CC5217"
CC5217_1 <- subset(x = `V11B04-032_C1` , idents = c("one"))
CC5217_1@meta.data[["id"]] <- "CC5217_1"
CC5217_2 <- subset(x = `V11B04-032_C1` , idents = c("two"))
CC5217_2@meta.data[["id"]] <- "CC5217_2"
`V11B04-032_D1`@meta.data[["orig.ident"]] <- "CC5222"
CC5222_1 <- subset(x = `V11B04-032_D1` , idents = c("one"))
CC5222_1@meta.data[["id"]] <- "CC5222_1"
CC5222_2 <- subset(x = `V11B04-032_D1` , idents = c("two"))
CC5222_2@meta.data[["id"]] <- "CC5222_2"
`V11B04-129_C1`@meta.data[["orig.ident"]] <- "CC5205"
CC5205_1 <- subset(x = `V11B04-129_C1` , idents = c("one"))
CC5205_1@meta.data[["id"]] <- "CC5205_1"
CC5205_2 <- subset(x = `V11B04-129_C1` , idents = c("two"))
CC5205_2@meta.data[["id"]] <- "CC5205_2"

#sample A is two samples
Seurat::Idents(object = `V42L22-361_A1`) <- `V42L22-361_A1`@meta.data[["sample"]]
DQ9502 <- subset(x = `V42L22-361_A1` , idents = c("two"))
DQ9502@meta.data[["orig.ident"]] <- "DQ9502"
CC5218 <- subset(x = `V42L22-361_A1` , idents = c("one"))
CC5218@meta.data[["orig.ident"]] <- "CC5218"
Seurat::Idents(object = DQ9502) <- DQ9502@meta.data[["replicate"]]
DQ9502_1 <- subset(x = DQ9502 , idents = c("one"))
DQ9502_1@meta.data[["id"]] <- "DQ9502_1"
DQ9502_2 <- subset(x = DQ9502 , idents = c("two"))
DQ9502_2@meta.data[["id"]] <- "DQ9502_2"
Seurat::Idents(object = CC5218) <- CC5218@meta.data[["replicate"]]
CC5218_1 <- subset(x = CC5218 , idents = c("one"))
CC5218_1@meta.data[["id"]] <- "CC5218_1"
CC5218_2 <- subset(x = CC5218 , idents = c("two"))
CC5218_2@meta.data[["id"]] <- "CC5218_2"

#sample D is two samples
Seurat::Idents(object = `V42L22-361_D1`) <- `V42L22-361_D1`@meta.data[["sample"]]
DQ9511 <- subset(x = `V42L22-361_D1` , idents = c("two"))
DQ9511@meta.data[["orig.ident"]] <- "DQ9511"
DQ9512 <- subset(x = `V42L22-361_D1` , idents = c("one"))
DQ9512@meta.data[["orig.ident"]] <- "DQ9512"
Seurat::Idents(object = DQ9511) <- DQ9511@meta.data[["replicate"]]
DQ9511_1 <- subset(x = DQ9511 , idents = c("one"))
DQ9511_1@meta.data[["id"]] <- "DQ9511_1"
DQ9511_2 <- subset(x = DQ9511 , idents = c("two"))
DQ9511_2@meta.data[["id"]] <- "DQ9511_2"
DQ9511_3 <- subset(x = DQ9511 , idents = c("two"))
DQ9511_3@meta.data[["id"]] <- "DQ9511_3"
Seurat::Idents(object = DQ9512) <- DQ9512@meta.data[["replicate"]]
DQ9512_1 <- subset(x = DQ9512 , idents = c("one"))
DQ9512_1@meta.data[["id"]] <- "DQ9512_1"
DQ9512_2 <- subset(x = DQ9512 , idents = c("two"))
DQ9512_2@meta.data[["id"]] <- "DQ9512_2"
  
#############
prueba <- c(DQ9500_1, DQ9500_2, CC5217_1, CC5217_2, CC5222_1, CC5222_2,
            CC5205_1 ,CC5205_2 , DQ9502_1, DQ9502_2, CC5218_1, CC5218_2,
            DQ9511_1, DQ9511_2, DQ9511_3, DQ9512_1, DQ9512_2)
names_indiv <- c("DQ9500_1", "DQ9500_2", "CC5217_1", "CC5217_2", "CC5222_1", "CC5222_2",
                 "CC5205_1" ,"CC5205_2" , "DQ9502_1", "DQ9502_2", "CC5218_1", "CC5218_2",
                 "DQ9511_1", "DQ9511_2", "DQ9511_3", "DQ9512_1", "DQ9512_2")
names(prueba) <- names_indiv

####Mix each samples same labeled areas together
x <- prueba
objects <- c()
for (i in 1:length(x)){
  ####analysis
  object <- x[[i]]
  object <- SCTransform(object, assay = "Spatial", verbose = FALSE)
  ######add object to list
  objects[[length(objects) + 1]] <- object
}

names(objects) <- names_indiv
list2env(objects,envir=.GlobalEnv)


##Merge them
combined <- merge(DQ9500_1, y = c(DQ9500_2, CC5217_1, CC5217_2, CC5222_1, CC5222_2,
                                  CC5205_1 ,CC5205_2 , DQ9502_1, DQ9502_2, CC5218_1, CC5218_2,
                                  DQ9511_1, DQ9511_2, DQ9511_3, DQ9512_1, DQ9512_2), 
                  add.cell.ids = c("DQ9500_1", "DQ9500_2", "CC5217_1", "CC5217_2", "CC5222_1", "CC5222_2",
                                   "CC5205_1" ,"CC5205_2" , "DQ9502_1", "DQ9502_2", "CC5218_1", "CC5218_2",
                                   "DQ9511_1", "DQ9511_2", "DQ9511_3", "DQ9512_1", "DQ9512_2"), project = "Bone")

saveRDS(combined, "./objects/sp/combined_nofilter_discovery.rds")

############PREPARE PSEUDO DATA#########
## prepare the data
Seurat::Idents(object = combined) <- combined@meta.data[["area"]]
cts <- AggregateExpression(combined, 
                           group.by = c("area", "id"),
                           assays = 'Spatial',
                           slot = "counts",
                           return.seurat = FALSE)

cts <- cts$Spatial
head(cts)
df <- as.data.frame(cts)

write.csv(cts, "./count_replicates_discovery_spatial.csv", row.names=TRUE)

##Create metadata info
# First, we will get column names from your count_table
column_names <- colnames(df)

# Now, we will split column names into areas, sampleIDs, and replicates
metadata_list <- strsplit(column_names, "_")

# Now we create the dataframe from the list. 
# Note that the structure of column names has to be consistent for this to work.
metadata <- data.frame(
  Area = sapply(metadata_list, `[`, 1),
  SampleID = sapply(metadata_list, `[`, 2),
  Replicate = sapply(metadata_list, `[`, 3),
  row.names = column_names
)

write.csv(metadata, "./meta_discovery.csv", row.names=TRUE)


