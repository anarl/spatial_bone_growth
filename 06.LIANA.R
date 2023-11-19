#### LIANA on Spatial

rm(list=ls())

setwd(workdir)
set.seed(123)
library(tidyverse)
library(magrittr)
library(liana)
library(Seurat)

combined <- readRDS("combined_nofilter_discovery.rds")
DefaultAssay(combined) <- "SCT"
metaData <- read.table("meta_discovery_2.csv", sep = ",", header = T, row.names = 1)
head(metaData)

seurat_RRST <- subset(combined, subset = orig.ident %in% unique(metaData$SampleID[which(metaData$Platform == "RRST")]))


liana_RRST <- liana_wrap(seurat_RRST, idents_col = "area")
liana_RRST <- liana_RRST %>% liana_aggregate()


seurat_CytAssist <- subset(combined, subset = orig.ident %in% unique(metaData$SampleID[which(metaData$Platform == "CytAssist")]))


liana_CytAssist <- liana_wrap(seurat_CytAssist, idents_col = "area")
liana_CytAssist <- liana_CytAssist %>% liana_aggregate()


liana_RRST %>% liana_dotplot(source_groups = c("HZ","PZ","RZ","SOC"), target_groups = c("HZ","PZ","RZ","SOC"), ntop = 10)
liana_CytAssist %>% liana_dotplot(source_groups = c("HZ","PZ","RZ","SOC"), target_groups = c("HZ","PZ","RZ","SOC"), ntop = 10)


liana_RRST_trunc <- liana_RRST %>% filter(aggregate_rank <= 0.05) # note that these pvals are already corrected
liana_CytAssist_trunc <- liana_CytAssist %>% filter(aggregate_rank <= 0.05) # note that these pvals are already corrected


liana_RRST_trunc %>% liana_dotplot(source_groups = c("HZ","PZ","RZ","SOC"),
                                   target_groups = c("HZ","PZ","RZ","SOC"), ntop = 30) 
liana_CytAssist_trunc %>% liana_dotplot(source_groups = c("HZ","PZ","RZ","SOC"), 
                                        target_groups = c("HZ","PZ","RZ","SOC"), ntop = 40) 



liana_RRST$name <- paste(liana_RRST$source, liana_RRST$target, 
                         liana_RRST$ligand.complex, liana_RRST$receptor.complex, sep = "_")
liana_CytAssist$name <- paste(liana_CytAssist$source, liana_CytAssist$target, 
                              liana_CytAssist$ligand.complex, liana_CytAssist$receptor.complex, sep = "_")




liana_total <- merge(liana_RRST, liana_CytAssist, by="name", all = T)

liana_both <- subset(liana_total, aggregate_rank.x <= 0.05 & aggregate_rank.y <= 0.05)



### LIANA analysis removing SOC
seurat_RRST_noSOC <- subset(seurat_RRST, subset = area %in% c("HZ","PZ","RZ"))
table(seurat_RRST_noSOC$area)
liana_RRST_noSOC <- liana_wrap(seurat_RRST_noSOC, idents_col = "area")
liana_RRST_noSOC <- liana_RRST_noSOC %>% liana_aggregate()

seurat_CytAssist_noSOC <- subset(seurat_CytAssist, subset = area %in% c("HZ","PZ","RZ"))
table(seurat_CytAssist_noSOC$area)
liana_CytAssist_noSOC <- liana_wrap(seurat_CytAssist_noSOC, idents_col = "area")
liana_CytAssist_noSOC <- liana_CytAssist_noSOC %>% liana_aggregate()


liana_RRST_noSOC %>% liana_dotplot(source_groups = c("HZ","PZ","RZ"), 
                             target_groups = c("HZ","PZ","RZ"), ntop = 10)
liana_CytAssist_noSOC %>% liana_dotplot(source_groups = c("HZ","PZ","RZ"), 
                                  target_groups = c("HZ","PZ","RZ"), ntop = 10)


liana_RRST_noSOC_trunc <- liana_RRST_noSOC %>% filter(aggregate_rank <= 0.05) # note that these pvals are already corrected
liana_CytAssist_noSOC_trunc <- liana_CytAssist_noSOC %>% filter(aggregate_rank <= 0.05) # note that these pvals are already corrected

pdf("LIANA_RRST_noSOC.pdf")
liana_RRST_noSOC_trunc %>% liana_dotplot(source_groups = c("HZ","PZ","RZ","SOC"),
                                   target_groups = c("HZ","PZ","RZ","SOC"), ntop = 30) 
dev.off()
pdf("LIANA_CytAssit_noSOC.pdf")
liana_CytAssist_noSOC_trunc %>% liana_dotplot(source_groups = c("HZ","PZ","RZ","SOC"), 
                                        target_groups = c("HZ","PZ","RZ","SOC"), ntop = 40) 
dev.off()



#### LIANA only in SOC and RZ

seurat_RRST_RZSOC <- subset(seurat_RRST, subset = area %in% c("RZ", "SOC"))
table(seurat_RRST_RZSOC$area)
liana_RRST_RZSOC <- liana_wrap(seurat_RRST_RZSOC, idents_col = "area")
liana_RRST_RZSOC <- liana_RRST_RZSOC %>% liana_aggregate()

seurat_CytAssist_RZSOC <- subset(seurat_CytAssist, subset = area %in% c("RZ","SOC"))
table(seurat_CytAssist_RZSOC$area)
liana_CytAssist_RZSOC <- liana_wrap(seurat_CytAssist_RZSOC, idents_col = "area")
liana_CytAssist_RZSOC <- liana_CytAssist_RZSOC %>% liana_aggregate()

pdf("LIANA_RRST_top10RZSOC.pdf")
liana_RRST_RZSOC %>% liana_dotplot(source_groups = c("RZ","SOC"), 
                                   target_groups = c("RZ","SOC"), ntop = 10)
dev.off()
pdf("LIANA_CytAssist_top10RZSOC.pdf")
liana_CytAssist_RZSOC %>% liana_dotplot(source_groups = c("RZ","SOC"), 
                                  target_groups = c("RZ","SOC"), ntop = 10)
dev.off()

liana_RRST_RZSOC_trunc <- liana_RRST_RZSOC %>% filter(aggregate_rank <= 0.05) # note that these pvals are already corrected
liana_CytAssist_RZSOC_trunc <- liana_CytAssist_RZSOC %>% filter(aggregate_rank <= 0.05) # note that these pvals are already corrected

pdf("LIANA_RRST_RZSOC.pdf")
liana_RRST_RZSOC_trunc %>% liana_dotplot(source_groups = c("RZ","SOC"),
                                         target_groups = c("RZ","SOC"), ntop = 30) 
dev.off()
pdf("LIANA_CytAssit_RZSOC.pdf")
liana_CytAssist_RZSOC_trunc %>% liana_dotplot(source_groups = c("HZ","PZ","RZ","SOC"), 
                                              target_groups = c("HZ","PZ","RZ","SOC"), ntop = 40) 
dev.off()


#### LIANA only in RZ and PZ

seurat_RRST_PZRZ <- subset(seurat_RRST, subset = area %in% c("RZ", "PZ"))
table(seurat_RRST_PZRZ$area)
liana_RRST_PZRZ <- liana_wrap(seurat_RRST_PZRZ, idents_col = "area")
liana_RRST_PZRZ <- liana_RRST_PZRZ %>% liana_aggregate()

seurat_CytAssist_PZRZ <- subset(seurat_CytAssist, subset = area %in% c("RZ","PZ"))
table(seurat_CytAssist_PZRZ$area)
liana_CytAssist_PZRZ <- liana_wrap(seurat_CytAssist_PZRZ, idents_col = "area")
liana_CytAssist_PZRZ <- liana_CytAssist_PZRZ %>% liana_aggregate()

pdf("LIANA_RRST_top10PZRZ.pdf")
liana_RRST_PZRZ %>% liana_dotplot(source_groups = c("PZ","RZ"), 
                                   target_groups = c("PZ","RZ"), ntop = 10)
dev.off()
pdf("LIANA_CytAssist_top10PZRZ.pdf")
liana_CytAssist_PZRZ %>% liana_dotplot(source_groups = c("PZ","RZ"), 
                                  target_groups = c("PZ","RZ"), ntop = 10)
dev.off()

liana_RRST_PZRZ_trunc <- liana_RRST_PZRZ %>% filter(aggregate_rank <= 0.05) # note that these pvals are already corrected
liana_CytAssist_PZRZ_trunc <- liana_CytAssist_PZRZ %>% filter(aggregate_rank <= 0.05) # note that these pvals are already corrected


#### LIANA only in HZ and PZ

seurat_RRST_HZPZ <- subset(seurat_RRST, subset = area %in% c("HZ", "PZ"))
table(seurat_RRST_HZPZ$area)
liana_RRST_HZPZ <- liana_wrap(seurat_RRST_HZPZ, idents_col = "area")
liana_RRST_HZPZ <- liana_RRST_HZPZ %>% liana_aggregate()

seurat_CytAssist_HZPZ <- subset(seurat_CytAssist, subset = area %in% c("HZ","PZ"))
table(seurat_CytAssist_HZPZ$area)
liana_CytAssist_HZPZ <- liana_wrap(seurat_CytAssist_HZPZ, idents_col = "area")
liana_CytAssist_HZPZ <- liana_CytAssist_HZPZ %>% liana_aggregate()

pdf("LIANA_RRST_top10HZPZ.pdf")
liana_RRST_HZPZ %>% liana_dotplot(source_groups = c("HZ","PZ"), 
                                   target_groups = c("HZ","PZ"), ntop = 10)
dev.off()
pdf("LIANA_CytAssist_top10HZPZ.pdf")
liana_CytAssist_HZPZ %>% liana_dotplot(source_groups = c("HZ","PZ"), 
                                  target_groups = c("HZ","PZ"), ntop = 10)
dev.off()

liana_RRST_HZPZ_trunc <- liana_RRST_HZPZ %>% filter(aggregate_rank <= 0.05) # note that these pvals are already corrected
liana_CytAssist_HZPZ_trunc <- liana_CytAssist_HZPZ %>% filter(aggregate_rank <= 0.05) # note that these pvals are already corrected

pdf("LIANA_RRST_HZPZ.pdf")
liana_RRST_HZPZ_trunc %>% liana_dotplot(source_groups = c("HZ","PZ"),
                                         target_groups = c("HZ","PZ"), ntop = 30) 
dev.off()
pdf("LIANA_CytAssit_HZPZ.pdf")
liana_CytAssist_HZPZ_trunc %>% liana_dotplot(source_groups = c("HZ","PZ"), 
                                              target_groups = c("HZ","PZ"), ntop = 40) 
dev.off()



### Check spots that are marginal SOC
head(seurat_RRST@meta.data)
barcodes <- data.frame(Row.names = row.names(seurat_RRST@meta.data),
                       Barcode = unlist(lapply(strsplit(row.names(seurat_RRST@meta.data), split = "_"),
                                               function(x) as.character(x)[3])),
                       names = seurat_RRST$orig.ident)
head(barcodes)
area <- list()
area$CC5217 <- read.csv("Discovery/Discovery/V11B04-032_C1/C032_marginalSOC.csv")
area$CC5222 <- read.csv("Discovery/Discovery/V11B04-032_D1/D032_marginalSOC.csv")
area$CC5205 <- read.csv("Discovery/Discovery/V11B04-129_C1/C129_marginalSOC.csv")
area$DQ9500 <- read.csv("Discovery/Discovery/V11B18-363_B1/B363_marginalSOC.csv")

area_new <- list()
for (i in 1:length(area)) {
  bar_sample <- subset(barcodes, names %in% names(area)[i])
  area_new[[i]] <- merge(area[[i]], bar_sample, by="Barcode")
  rownames(area_new[[i]]) <- area_new[[i]]$Row.names
  colnames(area_new[[i]])[2] <- "NewArea"
}

area_new <- rbind(area_new[[1]], area_new[[2]], area_new[[3]], area_new[[4]])
seurat_RRST <- AddMetaData(seurat_RRST, area_new)

seurat_RRST_marginalSOC <- subset(seurat_RRST, subset = NewArea %in% c("RZ","marginalSOC"))

table(seurat_RRST_marginalSOC$NewArea)
liana_RRST_RZmarginalSOC <- liana_wrap(seurat_RRST_marginalSOC, idents_col = "NewArea")
liana_RRST_RZmarginalSOC <- liana_RRST_RZmarginalSOC %>% liana_aggregate()


pdf("LIANA_RRST_top10RZmarginalSOC.pdf", width = 10)
liana_RRST_RZmarginalSOC %>% liana_dotplot(source_groups = c("RZ", "marginalSOC"), 
                                  target_groups = c("RZ","marginalSOC"), ntop = 10)
dev.off()

liana_RRST_RZmarginalSOC_trunc <- liana_RRST_RZmarginalSOC %>% filter(aggregate_rank <= 0.05) # note that these pvals are already corrected

pdf("LIANA_RRST_RZmarginalSOC.pdf")
liana_RRST_RZmarginalSOC_trunc %>% liana_dotplot(source_groups = c("RZ","marginalSOC"),
                                        target_groups = c("RZ","marginalSOC"), ntop = 30) 
dev.off()


head(seurat_CytAssist@meta.data)
barcodes <- data.frame(Row.names = row.names(seurat_CytAssist@meta.data),
                       Barcode = unlist(lapply(strsplit(row.names(seurat_CytAssist@meta.data), split = "_"),
                                               function(x) as.character(x)[3])),
                       names = seurat_CytAssist$orig.ident)
head(barcodes)
table(metaData$Slide, metaData$SampleID)
area <- list()
area$A361 <- read.csv("Discovery/Discovery/V42L22-361_A1/A361_marginalSOC.csv")
area$D361 <- read.csv("Discovery/Discovery/V42L22-361_D1/D361_marginalSOC.csv")

bar_sample <- subset(barcodes, names %in% c("CC5218","DQ9502"))
area_new$A361 <- merge(area$A361, bar_sample, by="Barcode")
rownames(area_new$A361) <- area_new$A361$Row.names
colnames(area_new$A361)[2] <- "NewArea"

bar_sample <- subset(barcodes, names %in% c("DQ9511","DQ9512"))
area_new$D361 <- merge(area$D361, bar_sample, by="Barcode")
rownames(area_new$D361) <- area_new$D361$Row.names
colnames(area_new$D361)[2] <- "NewArea"

area_new <- rbind(area_new[[1]], area_new[[2]])
seurat_CytAssist <- AddMetaData(seurat_CytAssist, area_new)
head(seurat_CytAssist@meta.data)

seurat_CytAssist_marginalSOC <- subset(seurat_CytAssist, subset = NewArea %in% c("RZ","marginalSOC"))

table(seurat_CytAssist_marginalSOC$NewArea)
liana_CytAssist_RZmarginalSOC <- liana_wrap(seurat_CytAssist_marginalSOC, idents_col = "NewArea")
liana_CytAssist_RZmarginalSOC <- liana_CytAssist_RZmarginalSOC %>% liana_aggregate()


pdf("LIANA_CytAssist_top10RZmarginalSOC.pdf", width = 10)
liana_CytAssist_RZmarginalSOC %>% liana_dotplot(source_groups = c("RZ", "marginalSOC"), 
                                           target_groups = c("RZ","marginalSOC"), ntop = 10)
dev.off()

liana_CytAssist_RZmarginalSOC_trunc <- liana_CytAssist_RZmarginalSOC %>% filter(aggregate_rank <= 0.05) # note that these pvals are already corrected

pdf("LIANA_CytAssist_RZmarginalSOC.pdf", width = 10)
liana_CytAssist_RZmarginalSOC_trunc %>% liana_dotplot(source_groups = c("RZ","marginalSOC"),
                                                 target_groups = c("RZ","marginalSOC"), ntop = 30) 
dev.off()



### MarginalSOC vs marginalRZ
barcodes <- data.frame(Row.names = row.names(seurat_RRST@meta.data),
                       Barcode = unlist(lapply(strsplit(row.names(seurat_RRST@meta.data), split = "_"),
                                               function(x) as.character(x)[3])),
                       names = seurat_RRST$orig.ident)
head(barcodes)
area <- list()
area$CC5217 <- read.csv("Discovery/Discovery/V11B04-032_C1/C032_marginalSOC_marginalRZ.csv")
area$CC5222 <- read.csv("Discovery/Discovery/V11B04-032_D1/D032_marginalSOC_marginalRZ.csv")
area$CC5205 <- read.csv("Discovery/Discovery/V11B04-129_C1/C129_marginalSOC_marginalRZ.csv")
area$DQ9500 <- read.csv("Discovery/Discovery/V11B18-363_B1/B363_marginalSOC_marginalRZ.csv")

area_new <- list()
for (i in 1:length(area)) {
  bar_sample <- subset(barcodes, names %in% names(area)[i])
  area_new[[i]] <- merge(area[[i]], bar_sample, by="Barcode")
  rownames(area_new[[i]]) <- area_new[[i]]$Row.names
  colnames(area_new[[i]])[2] <- "NewArea"
}

area_new <- rbind(area_new[[1]], area_new[[2]], area_new[[3]], area_new[[4]])
seurat_RRST <- AddMetaData(seurat_RRST, area_new)

seurat_RRST_marginalSOCmarginalRZ <- subset(seurat_RRST, subset = NewArea %in% c("marginalRZ","marginalSOC"))

table(seurat_RRST_marginalSOCmarginalRZ$NewArea)
liana_RRST_marginalRZmarginalSOC <- liana_wrap(seurat_RRST_marginalSOCmarginalRZ, idents_col = "NewArea")
liana_RRST_marginalRZmarginalSOC <- liana_RRST_marginalRZmarginalSOC %>% liana_aggregate()


pdf("LIANA_RRST_top10marginalRZmarginalSOC.pdf")
liana_RRST_marginalRZmarginalSOC %>% liana_dotplot(source_groups = c("marginalRZ", "marginalSOC"), 
                                           target_groups = c("marginalRZ","marginalSOC"), ntop = 10)
dev.off()

liana_RRST_marginalRZmarginalSOC_trunc <- liana_RRST_marginalRZmarginalSOC %>% filter(aggregate_rank <= 0.05) # note that these pvals are already corrected

pdf("LIANA_RRST_marginalRZmarginalSOC.pdf")
liana_RRST_marginalRZmarginalSOC_trunc %>% liana_dotplot(source_groups = c("RZ","marginalSOC"),
                                                 target_groups = c("RZ","marginalSOC"), ntop = 30) 
dev.off()


barcodes <- data.frame(Row.names = row.names(seurat_CytAssist@meta.data),
                       Barcode = unlist(lapply(strsplit(row.names(seurat_CytAssist@meta.data), split = "_"),
                                               function(x) as.character(x)[3])),
                       names = seurat_CytAssist$orig.ident)
head(barcodes)
area <- list()
area$A361 <- read.csv("Discovery/Discovery/V42L22-361_A1/A361_marginalSOC_marginalRZ.csv")
area$D361 <- read.csv("Discovery/Discovery/V42L22-361_D1/D361_marginalSOC_marginalRZ.csv")

bar_sample <- subset(barcodes, names %in% c("CC5218","DQ9502"))
area_new <- list()
area_new$A361 <- merge(area$A361, bar_sample, by="Barcode")
rownames(area_new$A361) <- area_new$A361$Row.names
colnames(area_new$A361)[2] <- "NewArea"

bar_sample <- subset(barcodes, names %in% c("DQ9511","DQ9512"))
area_new$D361 <- merge(area$D361, bar_sample, by="Barcode")
rownames(area_new$D361) <- area_new$D361$Row.names
colnames(area_new$D361)[2] <- "NewArea"

area_new <- rbind(area_new[[1]], area_new[[2]])
seurat_CytAssist <- AddMetaData(seurat_CytAssist, area_new)

seurat_CytAssist_marginalSOCmarginalRZ <- subset(seurat_CytAssist, subset = NewArea %in% c("marginalRZ","marginalSOC"))

table(seurat_CytAssist_marginalSOCmarginalRZ$NewArea)
liana_CytAssist_marginalRZmarginalSOC <- liana_wrap(seurat_CytAssist_marginalSOCmarginalRZ, idents_col = "NewArea")
liana_CytAssist_marginalRZmarginalSOC <- liana_CytAssist_marginalRZmarginalSOC %>% liana_aggregate()


pdf("LIANA_CytAssist_top10marginalRZmarginalSOC.pdf")
liana_CytAssist_marginalRZmarginalSOC %>% liana_dotplot(source_groups = c("marginalRZ", "marginalSOC"), 
                                                   target_groups = c("marginalRZ","marginalSOC"), ntop = 10)
dev.off()

liana_CytAssist_marginalRZmarginalSOC_trunc <- liana_CytAssist_marginalRZmarginalSOC %>% filter(aggregate_rank <= 0.05) # note that these pvals are already corrected

pdf("LIANA_CytAssist_marginalRZmarginalSOC.pdf")
liana_CytAssist_marginalRZmarginalSOC_trunc %>% liana_dotplot(source_groups = c("RZ","marginalSOC"),
                                                         target_groups = c("RZ","marginalSOC"), ntop = 30) 
dev.off()

save.image("LIANA.RData")
load("LIANA.RData")
