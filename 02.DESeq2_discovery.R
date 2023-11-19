####################################################################################################
# DIFFERENTIAL EXPRESSION ON PSEUDOBULK DATA USING DESeq2
#####################################################################################################
rm(list=ls())

set.seed(123)

## Read files with counts and metadata
countTable <- read.table("count_replicates_discovery_spatial.csv", sep = ",", header = T, row.names = 1)
head(countTable)
metaData <- read.table("meta_discovery.csv", sep = ",", header = T, row.names = 1)
head(metaData)

## Load packages needed
suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(rtracklayer)
  library(VennDiagram)
  library(ggpubr)
  library(ComplexHeatmap)
  library(magrittr)
  library(DOSE)
  library(ggupset)
  library(cowplot)
  library(biomaRt)
  library(sva)
  library(enrichplot)
})


## split for RRST samples
metaData_RRST <- subset(metaData, Platform == "RRST")
countTable_RRST <- countTable[,rownames(metaData_RRST)]
dds_RRST <- DESeqDataSetFromMatrix(countData = countTable_RRST, colData = metaData_RRST, design = ~Area)
## Remove samples with less than 10 reads
keep <- rowSums(counts(dds_RRST)) >= 10
dds_RRST <- dds_RRST[keep,]

## Batch correction with Combat
corrected_RRST <- ComBat_seq(counts=counts(dds_RRST), batch = dds_RRST$SampleID, group = dds_RRST$Area)
ddsCor_RRST <- DESeqDataSetFromMatrix(countData = corrected_RRST, colData = colData(dds_RRST), design = ~ Area)
ddsCor_RRST <- DESeq(ddsCor_RRST)


## Estracting transformed values
rld_RRST <- rlog(ddsCor_RRST, blind=FALSE)
vsd_RRST <- varianceStabilizingTransformation(ddsCor_RRST, blind = FALSE)
dim(assay(rld_RRST))

## Heatmap of the sample matrix
samplesDists_RRST <- dist(t(assay(rld_RRST)))
sampleDistMatrix_RRST <- as.matrix(samplesDists_RRST)
row.names(sampleDistMatrix_RRST) <- colnames(countTable_RRST)
colnames(sampleDistMatrix_RRST) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Greens")))(255)

## Plot PCA
data_RRST <- plotPCA(rld_RRST, intgroup = "Area", returnData=TRUE)
data.frame(data_RRST$name)
percentVar_RRST <- round(100*attr(data_RRST, "percentVar"))


## Extract the results for each contrast and obtain only significant genes
HZ_PZ_RRST <- results(ddsCor_RRST, contrast = c("Area", "HZ", "PZ"))
HZ_RZ_RRST <- results(ddsCor_RRST, contrast = c("Area", "HZ", "RZ"))
HZ_SOC_RRST <- results(ddsCor_RRST, contrast = c("Area", "HZ", "SOC"))
PZ_RZ_RRST <- results(ddsCor_RRST, contrast = c("Area", "PZ", "RZ"))
PZ_SOC_RRST <- results(ddsCor_RRST, contrast = c("Area", "PZ", "SOC"))
RZ_SOC_RRST <- results(ddsCor_RRST, contrast = c("Area", "RZ", "SOC"))

HZ_PZ_RRST_sig <- subset(HZ_PZ_RRST, padj < 0.05) ## dim(HZ_PZ_RRST_sig) = 37 6
HZ_RZ_RRST_sig <- subset(HZ_RZ_RRST, padj < 0.05) ## dim(HZ_RZ_RRST_sig) = 210 6
HZ_SOC_RRST_sig <- subset(HZ_SOC_RRST, padj < 0.05) ## dim(HZ_SOC_RRST_sig) = 524 6
PZ_RZ_RRST_sig <- subset(PZ_RZ_RRST, padj < 0.05) ## dim(PZ_RZ_RRST_sig) = 209 6
PZ_SOC_RRST_sig <- subset(PZ_SOC_RRST, padj < 0.05) ## dim(PZ_SOC_RRST_sig) = 540 6
RZ_SOC_RRST_sig <- subset(RZ_SOC_RRST, padj < 0.05) ## dim(RZ_SOC_RRST_sig) = 61 6

### Remane colnmaes
colnames(HZ_PZ_RRST) <- c("baseMean_HZvsPZ_RRST","log2FC_HZvsPZ_RRST","lfcSE_HZvsPZ_RRST",
                          "stat_HZvsPZ_RRST","pvalue_HZvsPZ_RRST","padj_HZvsPZ_RRST")
colnames(HZ_RZ_RRST) <- c("baseMean_HZvsRZ_RRST","log2FC_HZvsRZ_RRST","lfcSE_HZvsRZ_RRST",
                          "stat_HZvsRZ_RRST","pvalue_HZvsRZ_RRST","padj_HZvsRZ_RRST")
colnames(HZ_SOC_RRST) <- c("baseMean_HZvsSOC_RRST","log2FC_HZvsSOC_RRST","lfcSE_HZvsSOC_RRST",
                          "stat_HZvsSOC_RRST","pvalue_HZvsSOC_RRST","padj_HZvsSOC_RRST")
colnames(PZ_RZ_RRST) <- c("baseMean_PZvsRZ_RRST","log2FC_PZvsRZ_RRST","lfcSE_PZvsRZ_RRST",
                          "stat_PZvsRZ_RRST","pvalue_PZvsRZ_RRST","padj_PZvsRZ_RRST")
colnames(PZ_SOC_RRST) <- c("baseMean_PZvsSOC_RRST","log2FC_PZvsSOC_RRST","lfcSE_PZvsSOC_RRST",
                          "stat_PZvsSOC_RRST","pvalue_PZvsSOC_RRST","padj_PZvsSOC_RRST")
colnames(RZ_SOC_RRST) <- c("baseMean_RZvsSOC_RRST","log2FC_RZvsSOC_RRST","lfcSE_RZvsSOC_RRST",
                          "stat_RZvsSOC_RRST","pvalue_RZvsSOC_RRST","padj_RZvsSOC_RRST")

## Combine results tables
Results1_RRST <- merge(as.data.frame(HZ_PZ_RRST)[,c(2,5,6)], as.data.frame(HZ_RZ_RRST)[,c(2,5,6)], by="row.names")
Results2_RRST <- merge(as.data.frame(HZ_SOC_RRST)[,c(2,5,6)], as.data.frame(PZ_RZ_RRST)[,c(2,5,6)], by="row.names")
Results3_RRST <- merge(as.data.frame(PZ_SOC_RRST)[,c(2,5,6)], as.data.frame(RZ_SOC_RRST)[,c(2,5,6)], by="row.names")
Results4_RRST <- merge(Results1_RRST, Results2_RRST, by="Row.names")
Results_RRST <- merge(Results4_RRST, Results3_RRST, by="Row.names")
rm(Results1_RRST,Results2_RRST,Results3_RRST,Results4_RRST)
colnames(Results_RRST)[1] <- "Genes"

### Create table with normalized count
normalized_counts_RRST <- counts(ddsCor_RRST, normalize = TRUE)
normalized_counts_RRST <- as.data.frame(normalized_counts_RRST)
pheatmap_counts_RRST <- normalized_counts_RRST[which(rownames(normalized_counts_RRST) %in% 
                                                       c(rownames(HZ_PZ_RRST_sig),rownames(HZ_RZ_RRST_sig),
                                                        rownames(HZ_SOC_RRST_sig),rownames(PZ_RZ_RRST_sig),
                                                        rownames(PZ_SOC_RRST_sig),rownames(RZ_SOC_RRST_sig))),] 
dim(pheatmap_counts_RRST)
## Total 767 genes DE in the different pairwise comparisons



### Cytoassist

metaData_Cyt <- subset(metaData, Platform == "CytAssist")
countTable_Cyt <- countTable[,rownames(metaData_Cyt)]
dds_Cyt <- DESeqDataSetFromMatrix(countData = countTable_Cyt, colData = metaData_Cyt, design = ~Area)
keep <- rowSums(counts(dds_Cyt)) >= 10
dds_Cyt <- dds_Cyt[keep,]

corrected_Cyt <- ComBat_seq(counts=counts(dds_Cyt), batch = dds_Cyt$SampleID, group = dds_Cyt$Area)
ddsCor_Cyt <- DESeqDataSetFromMatrix(countData = corrected_Cyt, colData = colData(dds_Cyt), design = ~ Area)
ddsCor_Cyt <- DESeq(ddsCor_Cyt)


## Estracting transformed values
rld_Cyt <- rlog(ddsCor_Cyt, blind=FALSE)
vsd_Cyt <- varianceStabilizingTransformation(ddsCor_Cyt, blind = FALSE)
dim(assay(rld_Cyt))

## Heatmap of the sample matrix
samplesDists_Cyt <- dist(t(assay(rld_Cyt)))
sampleDistMatrix_Cyt <- as.matrix(samplesDists_Cyt)
row.names(sampleDistMatrix_Cyt) <- colnames(countTable_Cyt)
colnames(sampleDistMatrix_Cyt) <- NULL

## Plot PCA
data_Cyt <- plotPCA(rld_Cyt, intgroup = "Area", returnData=TRUE)
data.frame(data_Cyt$name)
percentVar_Cyt <- round(100*attr(data_Cyt, "percentVar"))


## Extract the results for each contrast and obtain only significant genes
HZ_PZ_Cyt <- results(ddsCor_Cyt, contrast = c("Area", "HZ", "PZ"))
HZ_RZ_Cyt <- results(ddsCor_Cyt, contrast = c("Area", "HZ", "RZ"))
HZ_SOC_Cyt <- results(ddsCor_Cyt, contrast = c("Area", "HZ", "SOC"))
PZ_RZ_Cyt <- results(ddsCor_Cyt, contrast = c("Area", "PZ", "RZ"))
PZ_SOC_Cyt <- results(ddsCor_Cyt, contrast = c("Area", "PZ", "SOC"))
RZ_SOC_Cyt <- results(ddsCor_Cyt, contrast = c("Area", "RZ", "SOC"))

HZ_PZ_Cyt_sig <- subset(HZ_PZ_Cyt, padj < 0.05) ## dim(HZ_PZ_Cyt_sig) = 54 6
HZ_RZ_Cyt_sig <- subset(HZ_RZ_Cyt, padj < 0.05) ## dim(HZ_RZ_Cyt_sig) = 252 6
HZ_SOC_Cyt_sig <- subset(HZ_SOC_Cyt, padj < 0.05) ## dim(HZ_SOC_Cyt_sig) = 648 6
PZ_RZ_Cyt_sig <- subset(PZ_RZ_Cyt, padj < 0.05) ## dim(PZ_RZ_Cyt_sig) = 554 6
PZ_SOC_Cyt_sig <- subset(PZ_SOC_Cyt, padj < 0.05) ## dim(PZ_SOC_Cyt_sig) = 1397 6
RZ_SOC_Cyt_sig <- subset(RZ_SOC_Cyt, padj < 0.05) ## dim(RZ_SOC_Cyt_sig) = 171 6

### Remane colnmaes
colnames(HZ_PZ_Cyt) <- c("baseMean_HZvsPZ_Cyt","log2FC_HZvsPZ_Cyt","lfcSE_HZvsPZ_Cyt",
                          "stat_HZvsPZ_Cyt","pvalue_HZvsPZ_Cyt","padj_HZvsPZ_Cyt")
colnames(HZ_RZ_Cyt) <- c("baseMean_HZvsRZ_Cyt","log2FC_HZvsRZ_Cyt","lfcSE_HZvsRZ_Cyt",
                          "stat_HZvsRZ_Cyt","pvalue_HZvsRZ_Cyt","padj_HZvsRZ_Cyt")
colnames(HZ_SOC_Cyt) <- c("baseMean_HZvsSOC_Cyt","log2FC_HZvsSOC_Cyt","lfcSE_HZvsSOC_Cyt",
                           "stat_HZvsSOC_Cyt","pvalue_HZvsSOC_Cyt","padj_HZvsSOC_Cyt")
colnames(PZ_RZ_Cyt) <- c("baseMean_PZvsRZ_Cyt","log2FC_PZvsRZ_Cyt","lfcSE_PZvsRZ_Cyt",
                          "stat_PZvsRZ_Cyt","pvalue_PZvsRZ_Cyt","padj_PZvsRZ_Cyt")
colnames(PZ_SOC_Cyt) <- c("baseMean_PZvsSOC_Cyt","log2FC_PZvsSOC_Cyt","lfcSE_PZvsSOC_Cyt",
                           "stat_PZvsSOC_Cyt","pvalue_PZvsSOC_Cyt","padj_PZvsSOC_Cyt")
colnames(RZ_SOC_Cyt) <- c("baseMean_RZvsSOC_Cyt","log2FC_RZvsSOC_Cyt","lfcSE_RZvsSOC_Cyt",
                           "stat_RZvsSOC_Cyt","pvalue_RZvsSOC_Cyt","padj_RZvsSOC_Cyt")

## Combine results tables
Results1_Cyt <- merge(as.data.frame(HZ_PZ_Cyt)[,c(2,5,6)], as.data.frame(HZ_RZ_Cyt)[,c(2,5,6)], by="row.names")
Results2_Cyt <- merge(as.data.frame(HZ_SOC_Cyt)[,c(2,5,6)], as.data.frame(PZ_RZ_Cyt)[,c(2,5,6)], by="row.names")
Results3_Cyt <- merge(as.data.frame(PZ_SOC_Cyt)[,c(2,5,6)], as.data.frame(RZ_SOC_Cyt)[,c(2,5,6)], by="row.names")
Results4_Cyt <- merge(Results1_Cyt, Results2_Cyt, by="Row.names")
Results_Cyt <- merge(Results4_Cyt, Results3_Cyt, by="Row.names")
rm(Results1_Cyt,Results2_Cyt,Results3_Cyt,Results4_Cyt)
colnames(Results_Cyt)[1] <- "Genes"

### Create table with normalized count
normalized_counts_Cyt <- counts(ddsCor_Cyt, normalize = TRUE)
normalized_counts_Cyt <- as.data.frame(normalized_counts_Cyt)
pheatmap_counts_Cyt <- normalized_counts_Cyt[which(rownames(normalized_counts_Cyt) %in% 
                                                       c(rownames(HZ_PZ_Cyt_sig),rownames(HZ_RZ_Cyt_sig),
                                                         rownames(HZ_SOC_Cyt_sig),rownames(PZ_RZ_Cyt_sig),
                                                         rownames(PZ_SOC_Cyt_sig),rownames(RZ_SOC_Cyt_sig))),] 
dim(pheatmap_counts_Cyt)
## Total 1553 genes DE in the different pairwise comparisons



### Megre RRST and CytAssist


Results <- merge(Results_RRST, Results_Cyt, by = "Genes")



save.image("DiffExpression_boneGrowth.RData")

