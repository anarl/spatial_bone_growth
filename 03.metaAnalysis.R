####################################################################################################
# METANALYSIS OF RRST AND CYTASSIST SAMPLES
#####################################################################################################
rm(list=ls())

set.seed(123)

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
  library(ComplexHeatmap)
  library(ggupset)
  library(cowplot)
  library(biomaRt)
  library(sva)
  library(enrichplot)
  library(metaRNASeq)
})

load("DiffExpression_boneGrowth.RData")


##  Merge pairwise results from RRST and CytAssist

HZ_PZ <- merge(as.data.frame(HZ_PZ_RRST), as.data.frame(HZ_PZ_Cyt), by="row.names")
HZ_PZ_fisher <- subset(HZ_PZ, (HZ_PZ$log2FC_HZvsPZ_RRST > 0 & HZ_PZ$log2FC_HZvsPZ_Cyt > 0) |
                         (HZ_PZ$log2FC_HZvsPZ_RRST < 0 & HZ_PZ$log2FC_HZvsPZ_Cyt < 0))
HZ_PZ_fisher <- subset(HZ_PZ_fisher, HZ_PZ_fisher$pvalue_HZvsPZ_RRST < 0.2 | HZ_PZ_fisher$pvalue_HZvsPZ_Cyt < 0.2)
fishcomb_HZvsPZ <- fishercomb(list(RRST=HZ_PZ_fisher$pvalue_HZvsPZ_RRST, CytAssist=HZ_PZ_fisher$pvalue_HZvsPZ_Cyt), BHth = 0.05)
lapply(fishcomb_HZvsPZ, head)
lapply(fishcomb_HZvsPZ, length)

HZ_PZ_fisher$fisherTest <- fishcomb_HZvsPZ$TestStatistic
HZ_PZ_fisher$fisherPval <- fishcomb_HZvsPZ$rawpval
HZ_PZ_fisher$fisherPadj <- fishcomb_HZvsPZ$adjpval

HZ_PZ_fisher <- HZ_PZ_fisher[,-c(2:13)]
HZ_PZ_fisher <- merge(HZ_PZ, HZ_PZ_fisher, by="Row.names", all=T)

######

HZ_RZ <- merge(as.data.frame(HZ_RZ_RRST), as.data.frame(HZ_RZ_Cyt), by="row.names")
HZ_RZ_fisher <- subset(HZ_RZ, (HZ_RZ$log2FC_HZvsRZ_RRST > 0 & HZ_RZ$log2FC_HZvsRZ_Cyt > 0) |
                         (HZ_RZ$log2FC_HZvsRZ_RRST < 0 & HZ_RZ$log2FC_HZvsRZ_Cyt < 0))
HZ_RZ_fisher <- subset(HZ_RZ_fisher, HZ_RZ_fisher$pvalue_HZvsRZ_RRST < 0.2 | HZ_RZ_fisher$pvalue_HZvsRZ_Cyt < 0.2)

fishcomb_HZvsRZ <- fishercomb(list(RRST=HZ_RZ_fisher$pvalue_HZvsRZ_RRST, CytAssist=HZ_RZ_fisher$pvalue_HZvsRZ_Cyt), BHth = 0.05)
lapply(fishcomb_HZvsRZ, head)
lapply(fishcomb_HZvsRZ, length)

HZ_RZ_fisher$fisherTest <- fishcomb_HZvsRZ$TestStatistic
HZ_RZ_fisher$fisherPval <- fishcomb_HZvsRZ$rawpval
HZ_RZ_fisher$fisherPadj <- fishcomb_HZvsRZ$adjpval

HZ_RZ_fisher <- HZ_RZ_fisher[,-c(2:13)]
HZ_RZ_fisher <- merge(HZ_RZ, HZ_RZ_fisher, by="Row.names", all=T)

####


HZ_SOC <- merge(as.data.frame(HZ_SOC_RRST), as.data.frame(HZ_SOC_Cyt), by="row.names")
HZ_SOC_fisher <- subset(HZ_SOC, (HZ_SOC$log2FC_HZvsSOC_RRST > 0 & HZ_SOC$log2FC_HZvsSOC_Cyt > 0) |
                         (HZ_SOC$log2FC_HZvsSOC_RRST < 0 & HZ_SOC$log2FC_HZvsSOC_Cyt < 0))
HZ_SOC_fisher <- subset(HZ_SOC_fisher, HZ_SOC_fisher$pvalue_HZvsSOC_RRST < 0.2 | HZ_SOC_fisher$pvalue_HZvsSOC_Cyt < 0.2)

fishcomb_HZvsSOC <- fishercomb(list(RRST=HZ_SOC_fisher$pvalue_HZvsSOC_RRST, CytAssist=HZ_SOC_fisher$pvalue_HZvsSOC_Cyt), BHth = 0.05)
lapply(fishcomb_HZvsSOC, head)
lapply(fishcomb_HZvsSOC, length)

HZ_SOC_fisher$fisherTest <- fishcomb_HZvsSOC$TestStatistic
HZ_SOC_fisher$fisherPval <- fishcomb_HZvsSOC$rawpval
HZ_SOC_fisher$fisherPadj <- fishcomb_HZvsSOC$adjpval

HZ_SOC_fisher <- HZ_SOC_fisher[,-c(2:13)]
HZ_SOC_fisher <- merge(HZ_SOC, HZ_SOC_fisher, by="Row.names", all=T)

#####

PZ_RZ <- merge(as.data.frame(PZ_RZ_RRST), as.data.frame(PZ_RZ_Cyt), by="row.names")
PZ_RZ_fisher <- subset(PZ_RZ, (PZ_RZ$log2FC_PZvsRZ_RRST > 0 & PZ_RZ$log2FC_PZvsRZ_Cyt > 0) |
                         (PZ_RZ$log2FC_PZvsRZ_RRST < 0 & PZ_RZ$log2FC_PZvsRZ_Cyt < 0))
PZ_RZ_fisher <- subset(PZ_RZ_fisher, PZ_RZ_fisher$pvalue_PZvsRZ_RRST < 0.2 | PZ_RZ_fisher$pvalue_PZvsRZ_Cyt < 0.2)

fishcomb_PZvsRZ <- fishercomb(list(RRST=PZ_RZ_fisher$pvalue_PZvsRZ_RRST, CytAssist=PZ_RZ_fisher$pvalue_PZvsRZ_Cyt), BHth = 0.05)
lapply(fishcomb_PZvsRZ, head)
lapply(fishcomb_PZvsRZ, length)

PZ_RZ_fisher$fisherTest <- fishcomb_PZvsRZ$TestStatistic
PZ_RZ_fisher$fisherPval <- fishcomb_PZvsRZ$rawpval
PZ_RZ_fisher$fisherPadj <- fishcomb_PZvsRZ$adjpval

PZ_RZ_fisher <- PZ_RZ_fisher[,-c(2:13)]
PZ_RZ_fisher <- merge(PZ_RZ, PZ_RZ_fisher, by="Row.names", all=T)

#####


PZ_SOC <- merge(as.data.frame(PZ_SOC_RRST), as.data.frame(PZ_SOC_Cyt), by="row.names")
PZ_SOC_fisher <- subset(PZ_SOC, (PZ_SOC$log2FC_PZvsSOC_RRST > 0 & PZ_SOC$log2FC_PZvsSOC_Cyt > 0) |
                          (PZ_SOC$log2FC_PZvsSOC_RRST < 0 & PZ_SOC$log2FC_PZvsSOC_Cyt < 0))
PZ_SOC_fisher <- subset(PZ_SOC_fisher, PZ_SOC_fisher$pvalue_PZvsSOC_RRST < 0.2 | PZ_SOC_fisher$pvalue_PZvsSOC_Cyt < 0.2)

fishcomb_PZvsSOC <- fishercomb(list(RRST=PZ_SOC_fisher$pvalue_PZvsSOC_RRST, CytAssist=PZ_SOC_fisher$pvalue_PZvsSOC_Cyt), BHth = 0.05)
lapply(fishcomb_PZvsSOC, head)
lapply(fishcomb_PZvsSOC, length)

PZ_SOC_fisher$fisherTest <- fishcomb_PZvsSOC$TestStatistic
PZ_SOC_fisher$fisherPval <- fishcomb_PZvsSOC$rawpval
PZ_SOC_fisher$fisherPadj <- fishcomb_PZvsSOC$adjpval

PZ_SOC_fisher <- PZ_SOC_fisher[,-c(2:13)]
PZ_SOC_fisher <- merge(PZ_SOC, PZ_SOC_fisher, by="Row.names", all=T)

####

RZ_SOC <- merge(as.data.frame(RZ_SOC_RRST), as.data.frame(RZ_SOC_Cyt), by="row.names")
RZ_SOC_fisher <- subset(RZ_SOC, (RZ_SOC$log2FC_RZvsSOC_RRST > 0 & RZ_SOC$log2FC_RZvsSOC_Cyt > 0) |
                          (RZ_SOC$log2FC_RZvsSOC_RRST < 0 & RZ_SOC$log2FC_RZvsSOC_Cyt < 0))
RZ_SOC_fisher <- subset(RZ_SOC_fisher, RZ_SOC_fisher$pvalue_RZvsSOC_RRST < 0.2 | RZ_SOC_fisher$pvalue_RZvsSOC_Cyt < 0.2)

fishcomb_RZvsSOC <- fishercomb(list(RRST=RZ_SOC_fisher$pvalue_RZvsSOC_RRST, CytAssist=RZ_SOC_fisher$pvalue_RZvsSOC_Cyt), BHth = 0.05)
lapply(fishcomb_RZvsSOC, head)
lapply(fishcomb_RZvsSOC, length)

RZ_SOC_fisher$fisherTest <- fishcomb_RZvsSOC$TestStatistic
RZ_SOC_fisher$fisherPval <- fishcomb_RZvsSOC$rawpval
RZ_SOC_fisher$fisherPadj <- fishcomb_RZvsSOC$adjpval

RZ_SOC_fisher <- RZ_SOC_fisher[,-c(2:13)]
RZ_SOC_fisher <- merge(RZ_SOC, RZ_SOC_fisher, by="Row.names", all=T)

####

DEgenes <- list(upHZvsPZ = HZ_PZ_fisher[which(HZ_PZ_fisher$log2FC_HZvsPZ_RRST > 0 & HZ_PZ_fisher$log2FC_HZvsPZ_Cyt > 0 & 
                                         HZ_PZ_fisher$fisherPadj < .1),]$Row.names,
                upHZvsRZ = HZ_RZ[which(HZ_RZ_fisher$log2FC_HZvsRZ_RRST > 0 & HZ_RZ_fisher$log2FC_HZvsRZ_Cyt > 0 & 
                                         HZ_RZ_fisher$fisherPadj < .1),]$Row.names,
                upHZvsSOC = HZ_SOC[which(HZ_SOC_fisher$log2FC_HZvsSOC_RRST > 0 & HZ_SOC_fisher$log2FC_HZvsSOC_Cyt > 0 & 
                                           HZ_SOC_fisher$fisherPadj < .1),]$Row.names,
                upPZvsRZ = PZ_RZ[which(PZ_RZ_fisher$log2FC_PZvsRZ_RRST > 0 & PZ_RZ_fisher$log2FC_PZvsRZ_Cyt > 0 & 
                                         PZ_RZ_fisher$fisherPadj < .1),]$Row.names,
                upPZvsSOC = PZ_SOC[which(PZ_SOC_fisher$log2FC_PZvsSOC_RRST > 0 & PZ_SOC_fisher$log2FC_PZvsSOC_Cyt > 0 & 
                                           PZ_SOC_fisher$fisherPadj < .1),]$Row.names,
                upRZvsSOC = RZ_SOC[which(RZ_SOC_fisher$log2FC_RZvsSOC_RRST > 0 & RZ_SOC_fisher$log2FC_RZvsSOC_Cyt > 0 & 
                                           RZ_SOC_fisher$fisherPadj < .1),]$Row.names,
                downHZvsPZ = HZ_PZ[which(HZ_PZ_fisher$log2FC_HZvsPZ_RRST < 0 & HZ_PZ_fisher$log2FC_HZvsPZ_Cyt < 0 & 
                                           HZ_PZ_fisher$fisherPadj < .1),]$Row.names,
                downHZvsRZ = HZ_RZ[which(HZ_RZ_fisher$log2FC_HZvsRZ_RRST < 0 & HZ_RZ_fisher$log2FC_HZvsRZ_Cyt < 0 & 
                                           HZ_RZ_fisher$fisherPadj < .1),]$Row.names,
                downHZvsSOC = HZ_SOC[which(HZ_SOC_fisher$log2FC_HZvsSOC_RRST < 0 & HZ_SOC_fisher$log2FC_HZvsSOC_Cyt < 0 &  
                                             HZ_SOC_fisher$fisherPadj < .1),]$Row.names,
                downPZvsRZ = PZ_RZ[which(PZ_RZ_fisher$log2FC_PZvsRZ_RRST < 0 & PZ_RZ_fisher$log2FC_PZvsRZ_Cyt < 0 & 
                                           PZ_RZ_fisher$fisherPadj < .1),]$Row.names,
                downPZvsSOC = PZ_SOC[which(PZ_SOC_fisher$log2FC_PZvsSOC_RRST < 0 & PZ_SOC_fisher$log2FC_PZvsSOC_Cyt < 0 & 
                                             PZ_SOC_fisher$fisherPadj < .1),]$Row.names,
                downRZvsSOC = RZ_SOC[which(RZ_SOC_fisher$log2FC_RZvsSOC_RRST < 0 & RZ_SOC_fisher$log2FC_RZvsSOC_Cyt < 0 & 
                                             RZ_SOC_fisher$fisherPadj < .1),]$Row.names)

## Plot the total number of genes up/down for each pairwise
DEgenes_plot <- data.frame(number = c(unlist(lapply(DEgenes[1:6], length)),
                                           -unlist(lapply(DEgenes[7:12], length))),
                                Comparison = rep(c("HZvsPZ","HZvsRZ","HZvsSOC","PZvsRZ","PZvsSOC","RZvsSOC"), times=2),
                                Reg = rep(c("up-regulated","down-regualted"), each=6))

pdf("DEgenes_metaAnalysis.pdf")
ggplot(DEgenes_plot,aes(x=number, y = Comparison, fill=Reg))+ 
  geom_bar(stat="identity", position="identity") + theme_cowplot() + 
  scale_fill_manual(values = c("red3","green3"))
dev.off()


### HZ markers : genes DE in HZ area compare to the other 3 
lt_HZ = as.data.frame(list_to_matrix(DEgenes[c(1:3,7:9)]))
HZmarkers <- rownames(lt_HZ)[rowSums(lt_HZ[1:3]) >= 3]
length(HZmarkers)

venn.diagram(x=DEgenes[1:3], main = "HZ markers up", disable.logging = T,
             filename = 'NewPlots/HZmarkers_MetaAnalysis_up_pval.png', output = T,
             category.names = c("up-reg HZvsPZ", "up-reg HZvsRZ", "up-reg HZvsSOC"),
             lwd = 2, lty = 'blank', fill = c("#B7D6A3","#F6BE98","#FFDF7F"), 
             cex = 1.6, fontface = "bold", fontfamily = "sans")



### PZ markers : genes DE in PZ area compare to the other 3 
lt_PZ = as.data.frame(list_to_matrix(DEgenes[c(7,4,5,1,10,11)]))
PZmarkers <- rownames(lt_PZ)[rowSums(lt_PZ[1:3]) >= 3]
length(PZmarkers)

## Plot Venn diagrams showing the numbers of up/down genes
venn.diagram(x=DEgenes[c(7,4,5)], main = "PZ markers up", disable.logging = T,
             filename = 'NewPlots/PZmarkers_MetaAnalysis_up_pval.png', output = T,
             category.names = c("up-reg PZvsHZ", "up-reg PZvsRZ", "up-reg PZvsSOC"),
             lwd = 2, lty = 'blank', fill = c("#B7D6A3","#F6BE98","#FFDF7F"), 
             cex = 1.6, fontface = "bold", fontfamily = "sans")


### RZ markers : genes DE in RZ area compare to the other 3 
lt_RZ = as.data.frame(list_to_matrix(DEgenes[c(8,10,6,2,4,12)]))
RZmarkers <- rownames(lt_RZ)[rowSums(lt_RZ[1:3]) >= 3]

length(RZmarkers)

## Plot Venn diagrams showing the numbers of up/down genes

venn.diagram(x=DEgenes[c(8,10,6)], main = "RZ markers up", disable.logging = T,
             filename = 'NewPlots/RZmarkers_MetaAnalysis_up_pval.png', output = T,
             category.names = c("up-reg RZvsHZ", "up-reg RZvsPZ", "up-reg RZvsSOC"),
             lwd = 2, lty = 'blank', fill = c("#B7D6A3","#F6BE98","#FFDF7F"), 
             cex = 1.6, fontface = "bold", fontfamily = "sans")


### SOC markers : genes DE in SOC area compare to the other 3 
lt_SOC = as.data.frame(list_to_matrix(DEgenes[c(9,11,12,3,5,6)]))
SOCmarkers <- rownames(lt_SOC)[rowSums(lt_SOC[1:3]) >= 3]
length(SOCmarkers)


## Plot Venn diagrams showing the numbers of up/down genes
venn.diagram(x=DEgenes[c(9,11,12)], main = "SOC markers up", disable.logging = T,
             filename = 'NewPlots/SOCmarkers_MetaAnalysis_up_pval.png', output = T,
             category.names = c("up-reg SOCvsHZ", "up-reg SOCvsPZ", "up-reg SOCvsRZ"),
             lwd = 2, lty = 'blank', fill = c("#B7D6A3","#F6BE98","#FFDF7F"), 
             cex = 1.6, fontface = "bold", fontfamily = "sans")


### FINAL MARKES SELECTED
SelectedMarkers <- Results[which(Results$Genes %in% 
                                                c(HZmarkers, PZmarkers, RZmarkers, SOCmarkers)),]
SelectedMarkers$Area <- "NA"
SelectedMarkers[which(SelectedMarkers$Genes %in% HZ_finalMarkers),"Area"] <- "HZ"
SelectedMarkers[which(SelectedMarkers$Genes %in% PZ_finalMarkers),"Area"] <- "PZ"
SelectedMarkers[which(SelectedMarkers$Genes %in% RZ_finalMarkers),"Area"] <- "RZ"
SelectedMarkers[which(SelectedMarkers$Genes %in% SOC_finalMarkers),"Area"] <- "SOC"

write.csv(SelectedMarkers_MetaAnalysis, "SelectedMarkers.csv", row.names = F)



#### Plots
rowAnn <- SelectedMarkers[which(c(SOC_finalMarkers, RZ_finalMarkers, PZ_finalMarkers[c(12,1:11)], HZ_finalMarkers ) %in% rownames(vsd_RRST)), c("Genes","Area")]
rownames(rowAnn) <- rowAnn$Genes
RRST_scale = t(apply(assay(rld_RRST), 1, scale))
colnames(RRST_scale) <- colnames(assay(vsd_RRST))
HZ_finalMarkers <- sort(rowAnn$Genes[rowAnn$Area == "HZ"])
PZ_finalMarkers <- sort(rowAnn$Genes[rowAnn$Area == "PZ"])
RZ_finalMarkers <- sort(rowAnn$Genes[rowAnn$Area == "RZ"])
SOC_finalMarkers <- sort(rowAnn$Genes[rowAnn$Area == "SOC"])
rowAnn <- rowAnn[c(SOC_finalMarkers, RZ_finalMarkers, PZ_finalMarkers,HZ_finalMarkers),]
rowAnn$Genes <- NULL
RRST_scale = t(apply(assay(vsd_RRST), 1, scale))
colnames(RRST_scale) <- colnames(assay(vsd_RRST))
Cyt_scale = t(apply(assay(vsd_Cyt), 1, scale))
colnames(Cyt_scale) <- colnames(assay(vsd_Cyt))

ht_list = Heatmap(RRST_scale[rownames(rowAnn),c(25:32,17:24,9:16,1:8)], 
                  name = "RRST",cluster_columns = F, width = unit(50, "mm"), cluster_rows = F,
                  col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
                  row_split = rowAnn, column_title = "RRST",
                  show_column_names = F, row_names_gp = gpar(fontsize=5),
                  row_title_gp = gpar(fontsize=0),
                  left_annotation = rowAnnotation(df=rowAnn, 
                                                  col= list(Area = c(SOC="blue3",RZ="red3",PZ="yellow2",HZ="green4"))),
                 top_annotation = HeatmapAnnotation(df = metaData_RRST[c(25:32,17:24,9:16,1:8),c("Area","SampleID","Slide")],
                                                     col = list(Area = c(HZ="green4",PZ="yellow2",RZ="red3",SOC="blue3"),
                                                                SampleID = c(CC5205="purple",CC5217="turquoise",
                                                                             CC5222="orange",DQ9500="green"),
                                                                Slide = c(V11B04_032_C1="black",V11B04_032_D1="brown",
                                                                          V11B04_129_C1="grey",V11B18_363_B1="darkblue")),
                                                     annotation_name_side = "left")) +
  Heatmap(Cyt_scale[rownames(rowAnn),c(26:34,17:25,8:16,1:7)], 
          name = "CytAssist", cluster_columns = F,width = unit(50, "mm"),
          row_names_gp = gpar(fontsize = 8),
          col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
          top_annotation = HeatmapAnnotation(df = metaData_Cyt[c(26:34,17:25,8:16,1:7),c("Area","SampleID","Slide")],
                                             col = list(Area = c(HZ="green4",PZ="yellow2",RZ="red3",SOC="blue3"),
                                                        SampleID = c(CC5218="purple4",DQ9511="turquoise1",
                                                                     DQ9512="orange4",DQ9502="aquamarine4"),
                                                        Slide = c(V42L22_361_A1="#6E6E6E",V42L22_361_D1="brown2")),
                                             annotation_name_side = "left", show_annotation_name = F),
          show_column_names = FALSE, row_title = NULL, show_row_dend = FALSE)
pdf("Heatmap with selected markers.pdf", width = 10, height = 5)
ht_list
dev.off()

CartMarkers <- c("SFRP5","UCMA","CLU","COL2A1","COL9A1","ACAN","MATN3","MKI67","COL10A1","ALPL","COL1A1","CXCL12")
ht_list_Cart = Heatmap(RRST_scale[CartMarkers,c(25:32,17:24,9:16,1:8)], 
                  name = "RRST",cluster_columns = F, width = unit(50, "mm"), cluster_rows = F,
                  col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
                  show_column_names = F, row_names_gp = gpar(fontsize=5),
                  row_title_gp = gpar(fontsize=0),
                  top_annotation = HeatmapAnnotation(df = metaData_RRST[c(25:32,17:24,9:16,1:8),c("Area","SampleID","Slide")],
                                                     col = list(Area = c(HZ="green4",PZ="yellow2",RZ="red3",SOC="blue3"),
                                                                SampleID = c(CC5205="purple",CC5217="turquoise",
                                                                             CC5222="orange",DQ9500="green"),
                                                                Slide = c(V11B04_032_C1="black",V11B04_032_D1="brown",
                                                                          V11B04_129_C1="grey",V11B18_363_B1="darkblue")),
                                                     annotation_name_side = "left")) +
  Heatmap(Cyt_scale[CartMarkers,c(26:34,17:25,8:16,1:7)], 
          name = "CytAssist", cluster_columns = F,width = unit(50, "mm"),row_names_gp = gpar(fontsize = 8),
          col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
          top_annotation = HeatmapAnnotation(df = metaData_Cyt[c(26:34,17:25,8:16,1:7),c("Area","SampleID","Slide")],
                                             col = list(Area = c(HZ="green4",PZ="yellow2",RZ="red3",SOC="blue3"),
                                                        SampleID = c(CC5218="purple4",DQ9511="turquoise1",
                                                                     DQ9512="orange4",DQ9502="aquamarine4"),
                                                        Slide = c(V42L22_361_A1="#6E6E6E",V42L22_361_D1="brown2")),
                                             annotation_name_side = "left", show_annotation_name = F),
          show_column_names = FALSE)
pdf("NewPlots/Heatmap with cartilage markers.pdf", width = 10, height = 7)
ht_list_Cart
dev.off()



colnames(HZ_PZ_fisher)[14:16] <- c("fisherTest_HZvsPZ", "fisherPval_HZvsPZ", "fisherPadj_HZvsPZ")
colnames(HZ_RZ_fisher)[14:16] <- c("fisherTest_HZvsRZ", "fisherPval_HZvsRZ", "fisherPadj_HZvsRZ")
colnames(HZ_SOC_fisher)[14:16] <- c("fisherTest_HZvsSOC", "fisherPval_HZvsSOC", "fisherPadj_HZvsSOC")
colnames(PZ_RZ_fisher)[14:16] <- c("fisherTest_PZvsRZ", "fisherPval_PZvsRZ", "fisherPadj_PZvsRZ")
colnames(PZ_SOC_fisher)[14:16] <- c("fisherTest_PZvsSOC", "fisherPval_PZvsSOC", "fisherPadj_PZvsSOC")
colnames(RZ_SOC_fisher)[14:16] <- c("fisherTest_RZvsSOC", "fisherPval_RZvsSOC", "fisherPadj_RZvsSOC")
Results1 <- merge(as.data.frame(HZ_PZ_fisher)[,c(1,3,6,7,9,12:16)], as.data.frame(HZ_RZ_fisher)[,c(1,3,6,7,9,12:16)], by="Row.names", all=T)
Results2 <- merge(as.data.frame(HZ_SOC_fisher)[,c(1,3,6,7,9,12:16)], as.data.frame(PZ_RZ_fisher)[,c(1,3,6,7,9,12:16)], by="Row.names", all=T)
Results3 <- merge(as.data.frame(PZ_SOC_fisher)[,c(1,3,6,7,9,12:16)], as.data.frame(RZ_SOC_fisher)[,c(1,3,6,7,9,12:16)], by="Row.names", all=T)
Results4 <- merge(Results1, Results2, by="Row.names", all=T)
Results <- merge(Results4, Results3, by="Row.names", all=T)
rm(Results1,Results2,Results3,Results4)
colnames(Results)[1] <- "Genes"
Markers_MetaAnalysis <- Results[Results$Genes %in% unique(c(HZmarkers, PZmarkers, RZmarkers, SOCmarkers)),]
write.table(Markers_MetaAnalysis, file="Markers_MetaAnalysis_pval.csv", quote = F, sep = ",", row.names = F)

save.image("metaAnalysis.RData")


