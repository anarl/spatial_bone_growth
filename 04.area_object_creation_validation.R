####################################################################################################
# SPTIAL OBJECT CREATION WITH MANUAL AREA ANNOTATION
#####################################################################################################
rm(list=ls())

set.seed(123)
library(Seurat)
library(dplyr)
library(ggplot2)

  #Variables---------------------------------
DIR_ROOT <- file.path(getwd())
DIR_DATA <- file.path(DIR_ROOT, "Validation RRST/Validation RRST/")

  #Data--------------------------------------

samples <- dir(path = DIR_DATA)
samples <- samples[! samples %in% c("HGP study 1 metadata 20230509.xlsx","zip")]

  # Create each individual Seurat object for every sample
lista <- c(samples)
names(lista) <- c("B032", "A032","A120","B120","A129","B129","D129","A363", "D288")
prueba <-list()

  for (i in 1:length(lista)){
    a <- Load10X_Spatial(data.dir = paste0(DIR_DATA,"/", lista[i]),
                         filename = "filtered_feature_bc_matrix.h5",
                         assay = "Spatial",
                         slice = names(lista)[i],
                         image = Read10X_Image(image.dir = paste0(DIR_DATA,"/", lista[i], "/spatial"),
                                               image.name = "tissue_hires_image.png"),
                         filter.matrix = TRUE)
    
    ######add area data
    area_a <- read.csv(list.files(path = paste0(DIR_DATA,"/", lista[i]), pattern = "area", full.names = T))
    area_a <- as.vector(area_a$area)
    a@meta.data["area"] <- as.factor(area_a)
    
    ######add replicate data
    replicate_a <- read.csv(paste0(DIR_DATA,"/", lista[i], "/replicate.csv"))
    replicate_a <- as.vector(replicate_a$replicate)
    a@meta.data["replicate"] <- as.factor(replicate_a)
    
    ######subset data
    Seurat::Idents(object = a) <- a@meta.data[["area"]]
    a <- subset(x =a, idents = c("HZ","PZ","RZ", "SOC"))
    a@meta.data[["area"]] <- a@active.ident
    
    prueba[[i]] <- a

  }
i
###### for sample 8
a <- subset(x =a, idents = c("RZ", "SOC"))
a@meta.data[["area"]] <- a@active.ident

prueba[[i]] <- a
i=9

a <- Load10X_Spatial(data.dir = paste0(DIR_DATA,"/", lista[i]),
                     filename = "filtered_feature_bc_matrix.h5",
                     assay = "Spatial",
                     slice = names(lista)[i],
                     image = Read10X_Image(image.dir = paste0(DIR_DATA,"/", lista[i], "/spatial"),
                                           image.name = "tissue_hires_image.png"),
                     filter.matrix = TRUE)

######add area data
area_a <- read.csv(list.files(path = paste0(DIR_DATA,"/", lista[i]), pattern = "area", full.names = T))
area_a <- as.vector(area_a$area)
a@meta.data["area"] <- as.factor(area_a)

######add replicate data
replicate_a <- read.csv(paste0(DIR_DATA,"/", lista[i], "/replicate.csv"))
replicate_a <- as.vector(replicate_a$Sample)
a@meta.data["replicate"] <- as.factor(replicate_a)

######subset data
Seurat::Idents(object = a) <- a@meta.data[["area"]]
a <- subset(x =a, idents = c("HZ","PZ","RZ", "SOC"))
a@meta.data[["area"]] <- a@active.ident

prueba[[i]] <- a


names(prueba) <- names(lista)

  ###################################################OPTION1####################################
####Separate data in replicates
## separate list
list2env(prueba,envir=.GlobalEnv)

  #separate in replicates
  Seurat::Idents(object = B032) <- B032@meta.data[["replicate"]]
  B032@meta.data[["orig.ident"]] <- "DQ9500rep2"
  DQ9500rep2_1 <- subset(x = B032 , idents = c("one"))
  DQ9500rep2_1@meta.data[["id"]] <- "DQ9500rep2_1"
  DQ9500rep2_2 <- subset(x = B032 , idents = c("two"))
  DQ9500rep2_2@meta.data[["id"]] <- "DQ9500rep2_2"

  Seurat::Idents(object = A120) <- A120@meta.data[["replicate"]]
  A120@meta.data[["orig.ident"]] <- "CC5218rep2"
  CC5218rep2_1 <- subset(x = A120 , idents = c("one"))
  CC5218rep2_1@meta.data[["id"]] <- "CC5218rep2_1"
  
  Seurat::Idents(object = D129) <- D129@meta.data[["replicate"]]
  D129@meta.data[["orig.ident"]] <- "DQ9500rep4"
  DQ9500rep4_1 <- subset(x = D129 , idents = c("one"))
  DQ9500rep4_1@meta.data[["id"]] <- "DQ9500rep4_1"
  DQ9500rep4_2 <- subset(x = D129 , idents = c("two"))
  DQ9500rep4_2@meta.data[["id"]] <- "DQ9500rep4_2"
  
  Seurat::Idents(object = B120) <- B120@meta.data[["replicate"]]
  B120@meta.data[["orig.ident"]] <- "CC5218rep3"
  CC5218rep3_1 <- subset(x = B120 , idents = c("one"))
  CC5218rep3_1@meta.data[["id"]] <- "CC5218rep3_1"
  
  Seurat::Idents(object = A032) <- A032@meta.data[["replicate"]]
  A032@meta.data[["orig.ident"]] <- "CC5205rep2"
  CC5205rep2_1 <- subset(x = A032 , idents = c("one"))
  CC5205rep2_1@meta.data[["id"]] <- "CC5205rep2_1"
  CC5205rep2_2 <- subset(x = A032 , idents = c("two"))
  CC5205rep2_2@meta.data[["id"]] <- "CC5205rep2_2"
  CC5205rep2_3 <- subset(x = A032 , idents = c("three"))
  CC5205rep2_3@meta.data[["id"]] <- "CC5205rep2_3"
  
  Seurat::Idents(object = A129) <- A129@meta.data[["replicate"]]
  A129@meta.data[["orig.ident"]] <- "CC5215"
  CC5215_1 <- subset(x = A129 , idents = c("one"))
  CC5215_1@meta.data[["id"]] <- "CC5215_1"
  CC5215_2 <- subset(x = A129 , idents = c("two"))
  CC5215_2@meta.data[["id"]] <- "CC5215_2"
  
  Seurat::Idents(object = A363) <- A363@meta.data[["replicate"]]
  A363@meta.data[["orig.ident"]] <- "DQ9502rep2"
  DQ9502rep2_1 <- subset(x = A363 , idents = c("one"))
  DQ9502rep2_1@meta.data[["id"]] <- "DQ9502rep2_1"
  
  Seurat::Idents(object = B129) <- B129@meta.data[["replicate"]]
  B129@meta.data[["orig.ident"]] <- "CC5216"
  CC5216_1 <- subset(x = B129 , idents = c("one"))
  CC5216_1@meta.data[["id"]] <- "CC5216_1"
  CC5216_2 <- subset(x = B129 , idents = c("two"))
  CC5216_2@meta.data[["id"]] <- "CC5216_2"
  
  Seurat::Idents(object = D288) <- D288@meta.data[["replicate"]]
  CC5222_1 <- subset(x = D288 , idents = c("CC5222_rep1"))
  CC5222_1@meta.data[["orig.ident"]] <- "CC5222"
  CC5222_1@meta.data[["replicate"]] <- "one"
  CC5222_1@meta.data[["id"]] <- "CC5222cyt_1"
  CC5222_2 <- subset(x = D288 , idents = c("CC5222_rep2"))
  CC5222_2@meta.data[["orig.ident"]] <- "CC5222"
  CC5222_2@meta.data[["replicate"]] <- "two"
  CC5222_2@meta.data[["id"]] <- "CC52222cyt_1"
  DQ9500_1 <- subset(x = D288 , idents = c("DQ9500_rep1"))
  DQ9500_1@meta.data[["orig.ident"]] <- "DQ9500"
  DQ9500_1@meta.data[["replicate"]] <- "one"
  DQ9500_1@meta.data[["id"]] <- "DQ9500cyt_1"
  DQ9500_2 <- subset(x = D288 , idents = c("DQ9500_rep2"))
  DQ9500_2@meta.data[["orig.ident"]] <- "DQ9500"
  DQ9500_2@meta.data[["replicate"]] <- "two"
  DQ9500_2@meta.data[["id"]] <- "DQ9500cyt_2"
  

#############
prueba <- c(DQ9500rep2_1, DQ9500rep2_2, DQ9500rep4_1, DQ9500rep4_2,DQ9502rep2_1, 
            CC5205rep2_1, CC5205rep2_2,CC5205rep2_3,
            CC5215_1,CC5215_2,CC5216_1,CC5216_2,CC5218rep2_1,CC5218rep3_1)
names_indiv <- c("DQ9500rep2_1", "DQ9500rep2_2", "DQ9500rep4_1", "DQ9500rep4_2","DQ9502rep2_1", 
                 "CC5205rep2_1", "CC5205rep2_2","CC5205rep2_3",
                 "CC5215_1","CC5215_2","CC5216_1","CC5216_2","CC5218rep2_1","CC5218rep3_1")
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
  combined <- merge(DQ9500rep2_1, y = c(DQ9500rep2_2, DQ9500rep4_1, DQ9500rep4_2,DQ9502rep2_1, 
                                        CC5205rep2_1, CC5205rep2_2,CC5205rep2_3,
                                        CC5215_1,CC5215_2,CC5216_1,CC5216_2,CC5218rep2_1,CC5218rep3_1), 
                    add.cell.ids = c("DQ9500rep2_1", "DQ9500rep2_2", "DQ9500rep4_1", "DQ9500rep4_2","DQ9502rep2_1", 
                                     "CC5205rep2_1", "CC5205rep2_2","CC5205rep2_3",
                                     "CC5215_1","CC5215_2","CC5216_1","CC5216_2","CC5218rep2_1","CC5218rep3_1"), project = "Bone")

  saveRDS(combined, "combined_nofilter_validation_RRST.rds")

  
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

  write.csv(cts, "./count_replicates_validation_spatial_RRST.csv", row.names=TRUE)

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

  write.csv(metadata, "meta_validation_RRST.csv", row.names=TRUE)
  
  
  
  
  
  
  

  #############
  prueba <- c(DQ9500_1,DQ9500_2,CC5222_1,CC5222_2)
  names_indiv <- c("DQ9500cyt_1","DQ9500cyt_2","CC5222cyt_1","CC5222cyt_2")
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
  combined <- merge(DQ9500_1, y = c(DQ9500_2,CC5222_1,CC5222_2), 
                    add.cell.ids = c("DQ9500cyt_1","DQ9500cyt_2","CC5222cyt_1","CC5222cyt_2"), project = "Bone")
  
  saveRDS(combined, "combined_nofilter_validation_Cyt.rds")
  
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
  
  write.csv(cts, "./count_replicates_validation_spatial_Cyt.csv", row.names=TRUE)
  
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
  
  write.csv(metadata, "meta_validation_Cyt.csv", row.names=TRUE)
  
  
  