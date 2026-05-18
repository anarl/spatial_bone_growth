####################################################################################################
# VALIDATION OF MARKERS IN THE VALIDATION COHORT
#####################################################################################################
rm(list=ls())
set.seed(123)

## Read files with counts and metadata
countTable <- read.table("count_replicates_validation_spatial.csv", sep = ",", header = T, row.names = 1)
head(countTable)
metaData <- read.table("meta_validation.csv", sep = ",", header = T, row.names = 1)
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
  library(Seurat)
  library(ggupset)
  library(cowplot)
  library(biomaRt)
  library(sva)
  library(ComplexHeatmap)
  library(enrichplot)
})


## split for RRST samples
metaData_RRST <- subset(metaData, Platform == "RRST")
countTable_RRST <- countTable[,rownames(metaData_RRST)]
countTable_RRST <- countTable_RRST[,colSums(countTable_RRST) > 1000]
metaData_RRST <- metaData_RRST[rownames(metaData_RRST) %in% colnames(countTable_RRST),]
dds_RRST <- DESeqDataSetFromMatrix(countData = countTable_RRST, colData = metaData_RRST, design = ~Area)
## Remove samples with less than 10 reads
keep <- rowSums(counts(dds_RRST)) >= 10
dds_RRST <- dds_RRST[keep,]

## Batch correction with Combat
corrected_RRST <- ComBat_seq(counts=counts(dds_RRST), batch = dds_RRST$SampleID, group = dds_RRST$Area)
ddsCor_RRST <- DESeqDataSetFromMatrix(countData = corrected_RRST, colData = colData(dds_RRST), design = ~ Area+SampleID)
ddsCor_RRST <- DESeq(ddsCor_RRST)

### Selected marker genes
SelectedMarkers <- read.csv("SelectedMarkers.csv")
head(SelectedMarkers)

dds_SelMar_RRST <- ddsCor_RRST[rownames(ddsCor_RRST) %in% SelectedMarkers$Genes]

## Estracting transformed values
rld_RRST <- rlog(dds_SelMar_RRST, blind=FALSE)
vsd_RRST <- varianceStabilizingTransformation(dds_SelMar_RRST, blind = FALSE)
dim(assay(rld_RRST))


## Plot PCA
data_RRST <- plotPCA(rld_RRST, intgroup = c("Area","SampleID"), returnData=TRUE)
data.frame(data_RRST$name)
percentVar_RRST <- round(100*attr(data_RRST, "percentVar"))



## Extract the results for each contrast and obtain only significant genes
HZ_PZ_RRST <- results(dds_SelMar_RRST, contrast = c("Area", "HZ", "PZ"))
HZ_RZ_RRST <- results(dds_SelMar_RRST, contrast = c("Area", "HZ", "RZ"))
HZ_SOC_RRST <- results(dds_SelMar_RRST, contrast = c("Area", "HZ", "SOC"))
PZ_RZ_RRST <- results(dds_SelMar_RRST, contrast = c("Area", "PZ", "RZ"))
PZ_SOC_RRST <- results(dds_SelMar_RRST, contrast = c("Area", "PZ", "SOC"))
RZ_SOC_RRST <- results(dds_SelMar_RRST, contrast = c("Area", "RZ", "SOC"))

HZ_PZ_RRST_sig <- subset(HZ_PZ_RRST, padj < 0.05) ## dim(HZ_PZ_RRST_sig) = 5 6
HZ_RZ_RRST_sig <- subset(HZ_RZ_RRST, padj < 0.05) ## dim(HZ_RZ_RRST_sig) = 14 6
HZ_SOC_RRST_sig <- subset(HZ_SOC_RRST, padj < 0.05) ## dim(HZ_SOC_RRST_sig) = 25 6
PZ_RZ_RRST_sig <- subset(PZ_RZ_RRST, padj < 0.05) ## dim(PZ_RZ_RRST_sig) = 6 6
PZ_SOC_RRST_sig <- subset(PZ_SOC_RRST, padj < 0.05) ## dim(PZ_SOC_RRST_sig) = 40 6
RZ_SOC_RRST_sig <- subset(RZ_SOC_RRST, padj < 0.05) ## dim(RZ_SOC_RRST_sig) = 15 6



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



### Plot heatmap with markers
rowAnn <- SelectedMarkers[which(SelectedMarkers$Genes %in% rownames(vsd_RRST)), c("Genes","Area")]
rownames(rowAnn) <- rowAnn$Genes
RRST_scale = t(apply(assay(rld_RRST), 1, scale))
colnames(RRST_scale) <- colnames(assay(vsd_RRST))
HZ_finalMarkers <- sort(rowAnn$Genes[rowAnn$Area == "HZ"])
PZ_finalMarkers <- sort(rowAnn$Genes[rowAnn$Area == "PZ"])
RZ_finalMarkers <- sort(rowAnn$Genes[rowAnn$Area == "RZ"])
SOC_finalMarkers <- sort(rowAnn$Genes[rowAnn$Area == "SOC"])
rowAnn <- rowAnn[c(SOC_finalMarkers, RZ_finalMarkers, PZ_finalMarkers,HZ_finalMarkers),]
rowAnn$Genes <- NULL
pdf("NewPlots/HeatmapValidation_SelectedMarkers.pdf",  width = 10, height = 17)
Heatmap(RRST_scale[rownames(rowAnn),c(33:43,20:32,9:19,1:8)], 
        name = "RRST",cluster_columns = F, width = unit(50, "mm"), cluster_rows = F,
        col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
        row_split = rowAnn, 
        column_title = "RRST",
        show_column_names = F, row_names_gp = gpar(fontsize=5),
        row_title_gp = gpar(fontsize=0),
        left_annotation = rowAnnotation(df=rowAnn, 
                                        col= list(Area = c(SOC="blue3",RZ="red3",PZ="yellow2",HZ="green4"))),
        top_annotation = HeatmapAnnotation(df = metaData_RRST[c(33:43,20:32,9:19,1:8),c("Area","SampleID","Slide")],
                                           col = list(Area = c(HZ="green4",PZ="yellow2",RZ="red3",SOC="blue3"),
                                                      SampleID = c(CC5205="purple",CC5217="turquoise",CC5222="orange",
                                                                   DQ9500="green", CC5215="cadetblue1",CC5216="steelblue1",
                                                                   CC5218="maroon2",DQ9502="palevioletred"),
                                                      Slide = c(V11B04_032_C1="black",V11B04_032_D1="brown",
                                                                V11B04_129_C1="grey",V11B18_363_B1="darkblue",
                                                                V11B04_032_A1="bisque3",V11B04_129_A1="darkgoldenrod",
                                                                V11B04_129_B1="sienna3",V11B04_120_A1="sandybrown",
                                                                V11B18_363_A1="azure4",V11B04_032_B1="snow2",
                                                                V11B04_129_D1="wheat2",V11B04_120_B1="cornsilk")),
                                           annotation_name_side = "left")) 
dev.off()





Mahtab <- c("SFRP5","CLU","COL2A1","COL9A1","ACAN","MATN3","COL10A1","ALPL","COL1A1","CXCL12") #"UCMA","MKI67",
Mahtab %in% rownames(ddsCor_RRST)
rld_RRST_total <- rlog(ddsCor_RRST)
RRST_scale_total = t(apply(assay(rld_RRST_total), 1, scale))
colnames(RRST_scale_total) <- colnames(assay(rld_RRST_total))

pdf("NewPlots/BoneCartilageMarkers_Validation.pdf",  width = 7, height = 7)
Heatmap(RRST_scale_total[Mahtab,c(33:43,20:32,9:19,1:8)], 
        name = "RRST",cluster_columns = F, width = unit(50, "mm"), cluster_rows = F,
        col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
        show_column_names = F, row_names_gp = gpar(fontsize=8),
        row_title_gp = gpar(fontsize=0),
        top_annotation = HeatmapAnnotation(df = metaData_RRST[c(33:43,20:32,9:19,1:8),c("Area","SampleID","Slide")],
                                           col = list(Area = c(HZ="green4",PZ="blue3",RZ="red3",SOC="yellow2"),
                                                      SampleID = c(CC5205="purple",CC5217="turquoise",CC5222="orange",
                                                                   DQ9500="green", CC5215="cadetblue1",CC5216="steelblue1",
                                                                   CC5218="maroon2",DQ9502="palevioletred"),
                                                      Slide = c(V11B04_032_C1="black",V11B04_032_D1="brown",
                                                                V11B04_129_C1="grey",V11B18_363_B1="darkblue",
                                                                V11B04_032_A1="bisque3",V11B04_129_A1="darkgoldenrod",
                                                                V11B04_129_B1="sienna3",V11B04_120_A1="sandybrown",
                                                                V11B18_363_A1="azure4",V11B04_032_B1="snow2",
                                                                V11B04_129_D1="wheat2",V11B04_120_B1="cornsilk")),
                                           annotation_name_side = "left")) 
dev.off()





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



## Plot PCA
data_Cyt <- plotPCA(rld_Cyt, intgroup = "Area", returnData=TRUE)
data.frame(data_Cyt$name)
percentVar_Cyt <- round(100*attr(data_Cyt, "percentVar"))
pdf("PCA.Cyt.pdf")
ggplot(data_Cyt, aes(PC1, PC2, color= Area)) + 
  scale_color_manual(values=c("green4", "blue3", "red3", "yellow2")) +
  #geom_text(data = data, aes(x = PC1+1, y = PC2+2,label=name), cex =4) +
  geom_point(cex=3) + cowplot::theme_cowplot() +
  xlab(paste0("PC1: ", percentVar_Cyt[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_Cyt[2], "% variance"))
dev.off()

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


### Plot heatmap with markers
rowAnn <- SelectedMarkers[which(SelectedMarkers$Genes %in% rownames(vsd_Cyt)), c("Genes","Area")]
rownames(rowAnn) <- rowAnn$Genes
Cyt_scale = t(apply(assay(rld_Cyt), 1, scale))
colnames(Cyt_scale) <- colnames(assay(vsd_Cyt))
HZ_finalMarkers <- sort(rowAnn$Genes[rowAnn$Area == "HZ"])
PZ_finalMarkers <- sort(rowAnn$Genes[rowAnn$Area == "PZ"])
RZ_finalMarkers <- sort(rowAnn$Genes[rowAnn$Area == "RZ"])
SOC_finalMarkers <- sort(rowAnn$Genes[rowAnn$Area == "SOC"])
rowAnn <- rowAnn[c(SOC_finalMarkers, RZ_finalMarkers, PZ_finalMarkers,HZ_finalMarkers),]
rowAnn$Genes <- NULL
pdf("NewPlots/HeatmapValidation_SelectedMarkers.pdf",  width = 10, height = 17)
Heatmap(Cyt_scale[rownames(rowAnn),c(13:16,9:12,5:8,1:4)], 
        name = "Cyt",cluster_columns = F, width = unit(50, "mm"), cluster_rows = F,
        col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
        row_split = rowAnn, 
        column_title = "Cyt",
        show_column_names = F, row_names_gp = gpar(fontsize=5),
        row_title_gp = gpar(fontsize=0),
        left_annotation = rowAnnotation(df=rowAnn, 
                                        col= list(Area = c(SOC="blue3",RZ="red3",PZ="yellow2",HZ="green4"))),
        top_annotation = HeatmapAnnotation(df = metaData_Cyt[c(13:16,9:12,5:8,1:4),c("Area","SampleID","Slide")],
                                           col = list(Area = c(HZ="green4",PZ="yellow2",RZ="red3",SOC="blue3"),
                                                      SampleID = c(CC5205="purple",CC5217="turquoise",CC5222="orange",
                                                                   DQ9500="green", CC5215="cadetblue1",CC5216="steelblue1",
                                                                   CC5218="maroon2",DQ9502="palevioletred")),
                                           annotation_name_side = "left")) 
dev.off()

Mahtab <- c("SFRP5","CLU","COL2A1","COL9A1","ACAN","MATN3","COL10A1","ALPL","COL1A1","CXCL12") #"UCMA","MKI67",
Mahtab %in% rownames(ddsCor_Cyt)
rld_Cyt_total <- rlog(ddsCor_Cyt)
Cyt_scale_total = t(apply(assay(rld_Cyt_total), 1, scale))
colnames(Cyt_scale_total) <- colnames(assay(rld_Cyt_total))

pdf("NewPlots/BoneCartilageMarkers_Validation.pdf",  width = 7, height = 7)
Heatmap(Cyt_scale_total[Mahtab,c(13:16,9:12,5:8,1:4)], 
        name = "Cyt",cluster_columns = F, width = unit(50, "mm"), cluster_rows = F,
        col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
        show_column_names = F, row_names_gp = gpar(fontsize=8),
        row_title_gp = gpar(fontsize=0),
        top_annotation = HeatmapAnnotation(df = metaData_Cyt[c(13:16,9:12,5:8,1:4),c("Area","SampleID","Slide")],
                                           col = list(Area = c(HZ="green4",PZ="blue3",RZ="red3",SOC="yellow2"),
                                                      SampleID = c(CC5205="purple",CC5217="turquoise",CC5222="orange",
                                                                   DQ9500="green", CC5215="cadetblue1",CC5216="steelblue1",
                                                                   CC5218="maroon2",DQ9502="palevioletred")),
                                           annotation_name_side = "left")) 
dev.off()


save.image("DESeq2_validation")
