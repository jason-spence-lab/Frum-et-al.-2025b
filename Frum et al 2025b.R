library(Seurat)
library(ggplot2)
library(patchwork)

Dotplot_Zhiwei_Version <- function(seurat_object, gene_list) {
  DotPlot(seurat_object, features = gene_list,  dot.scale = 20, scale = TRUE, group.by = "name") + RotatedAxis() +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.4) +
    scale_colour_distiller(palette = "YlGnBu", direction = 1) +
    guides(size=guide_legend(title = "Percent Expressed",override.aes=list(shape=21, colour="black", fill="black"))) +
    labs(y=NULL, x= NULL) +
    guides(color = guide_colourbar(title = "Scaled Expression", ticks = TRUE, frame.colour = "black")) +
    theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())+
    theme(panel.border=element_rect(colour="black",fill=NA,size=0.8))
}

###Load Fetal Data
setwd("~/scRNAseq/Data/Cell/Daysha/Fetal/Sample_HT231_ESO/filtered_feature_bc_matrix") #file:///run/user/1002/gvfs/smb-share:server=umms-spencejr-win.turbo.storage.umich.edu,share=umms-spencejr/01_RAW_RNASEQ_AGC_SHARE/Run_2598/Sample_HT-231_esophagus/outs/web_summary.html
ht231eso.data <- Read10X(data.dir = "./") #117 day male
setwd("~/scRNAseq/Data/Cell/Daysha/Fetal/Sample_HT232_ESO/filtered_feature_bc_matrix") #file:///run/user/1002/gvfs/smb-share:server=umms-spencejr-win.turbo.storage.umich.edu,share=umms-spencejr/01_RAW_RNASEQ_AGC_SHARE/Run_2598/Sample_HT-231_esophagus/outs/web_summary.html
ht232eso.data <- Read10X(data.dir = "./") #100 day male

ht231eso <- CreateSeuratObject(counts = ht231eso.data, project = "ht231eso", min.cells = 3, min.features = 200) 
ht232eso <- CreateSeuratObject(counts = ht232eso.data, project = "ht232eso", min.cells = 3, min.features = 200)

ht231eso <- RenameCells(object = ht231eso, add.cell.id = "ht231eso")
ht232eso <- RenameCells(object = ht232eso, add.cell.id = "ht232eso")

ht231eso[["percent.mt"]] <- PercentageFeatureSet(ht231eso, pattern = "^MT-") 
ht232eso[["percent.mt"]] <- PercentageFeatureSet(ht232eso, pattern = "^MT-") 

ht231eso <- AddMetaData(ht231eso, "fetal", "age")
ht232eso <- AddMetaData(ht232eso, "fetal", "age")

###Process Without QC Filtering To Recover Multiciliated Cells
esounfiltered.fetal.combined <- merge(ht231eso, ht232eso)
esounfiltered.fetal.combined <- JoinLayers(esounfiltered.fetal.combined)
ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = esounfiltered.fetal.combined), value = TRUE);
percent.ribo <- Matrix::colSums(x = GetAssayData(object = esounfiltered.fetal.combined, slot = 'counts')[ribo.genes, ]) / Matrix::colSums(x = GetAssayData(object = esounfiltered.fetal.combined, slot = 'counts'));
esounfiltered.fetal.combined[['percent.ribo']] <- percent.ribo;


esounfiltered.fetal.combined <- SplitObject(esounfiltered.fetal.combined, split.by = "orig.ident")
for (i in seq_along(esounfiltered.fetal.combined)) {
  esounfiltered.fetal.combined[[i]] <- NormalizeData(esounfiltered.fetal.combined[[i]]) %>% FindVariableFeatures()
}
features <- SelectIntegrationFeatures(esounfiltered.fetal.combined)
for (i in seq_along(along.with = esounfiltered.fetal.combined)) {
  esounfiltered.fetal.combined[[i]] <- ScaleData(esounfiltered.fetal.combined[[i]], features = features) %>% RunPCA(features = features)
}
anchors <- FindIntegrationAnchors(esounfiltered.fetal.combined, anchor.features = features, reduction = "rpca", dims = 1:30)
esounfiltered.fetal.combined.integrated <- IntegrateData(anchors, dims = 1:30)
DefaultAssay(esounfiltered.fetal.combined.integrated) <- "integrated"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
esounfiltered.fetal.combined.integrated <- CellCycleScoring(esounfiltered.fetal.combined.integrated, s.features = s.genes, g2m.features = g2m.genes)
esounfiltered.fetal.combined.integrated <- ScaleData(esounfiltered.fetal.combined.integrated, vars.to.regress = c("S.Score", "G2M.Score"))
esounfiltered.fetal.combined.integrated <- RunPCA(esounfiltered.fetal.combined.integrated)
ElbowPlot(esounfiltered.fetal.combined.integrated, ndims = 50)
esounfiltered.fetal.combined.integrated <- RunUMAP(esounfiltered.fetal.combined.integrated, dims = 1:14, reduction.name = "umap", return.model = TRUE)
esounfiltered.fetal.combined.integrated <- FindNeighbors(esounfiltered.fetal.combined.integrated, dims = 1:14)
esounfiltered.fetal.combined.integrated <- FindClusters(esounfiltered.fetal.combined.integrated, resolution = 0.15)
DimPlot(esounfiltered.fetal.combined.integrated, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)

pdf(file.path("./", paste0("Fetal ESO No Filtering Louvain", ".pdf")), w=11, h=8.5)
DimPlot(esounfiltered.fetal.combined.integrated, reduction = "umap", label = TRUE, pt.size = 2)
dev.off()

DefaultAssay(esounfiltered.fetal.combined.integrated) <- "RNA"
esounfiltered.fetal.combined.integrated <- NormalizeData(esounfiltered.fetal.combined.integrated)
esounfiltered.fetal.combined.integrated <- JoinLayers(esounfiltered.fetal.combined.integrated)

pdf(file.path("./", paste0("Fetal ESO FOXJ1", ".pdf")), w=11, h=8.5)
FeaturePlot(esounfiltered.fetal.combined.integrated, features = c("FOXJ1"),  pt.size = 0.2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

esounfiltered.fetal.unfiltered.ciliated.cellids <- WhichCells(esounfiltered.fetal.combined.integrated, idents = 11)
esounfiltered.fetal.unfiltered.ciliated.subset <- subset(esounfiltered.fetal.combined.integrated, idents = 11) #82 cells


#####Process with Normal QC Filering Parameters and Merge in Multiciliated Cells
eso.fetal.combined <- merge(ht231eso, ht232eso)

eso.fetal.combined <- subset(eso.fetal.combined, subset = nFeature_RNA > 1500 & nFeature_RNA < 5000 & percent.mt < 10)
eso.fetal.combined.with.ciliated <- merge(eso.fetal.combined, esounfiltered.fetal.unfiltered.ciliated.subset)

eso.fetal.combined.with.ciliated <- JoinLayers(eso.fetal.combined.with.ciliated)

eso.fetal.combined.with.ciliated <- SplitObject(eso.fetal.combined.with.ciliated, split.by = "orig.ident")
for (i in seq_along(eso.fetal.combined.with.ciliated)) {
  eso.fetal.combined.with.ciliated[[i]] <- NormalizeData(eso.fetal.combined.with.ciliated[[i]]) %>% FindVariableFeatures()
}
features <- SelectIntegrationFeatures(eso.fetal.combined.with.ciliated)
for (i in seq_along(along.with = eso.fetal.combined.with.ciliated)) {
  eso.fetal.combined.with.ciliated[[i]] <- ScaleData(eso.fetal.combined.with.ciliated[[i]], features = features) %>% RunPCA(features = features)
}
anchors <- FindIntegrationAnchors(eso.fetal.combined.with.ciliated, anchor.features = features, reduction = "rpca", dims = 1:30)
eso.fetal.combined.with.ciliated.integrated <- IntegrateData(anchors, dims = 1:30)
DefaultAssay(eso.fetal.combined.with.ciliated.integrated) <- "integrated"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
eso.fetal.combined.with.ciliated.integrated <- CellCycleScoring(eso.fetal.combined.with.ciliated.integrated, s.features = s.genes, g2m.features = g2m.genes)
eso.fetal.combined.with.ciliated.integrated <- ScaleData(eso.fetal.combined.with.ciliated.integrated, vars.to.regress = c("S.Score", "G2M.Score"))
eso.fetal.combined.with.ciliated.integrated <- RunPCA(eso.fetal.combined.with.ciliated.integrated)
ElbowPlot(eso.fetal.combined.with.ciliated.integrated, ndims = 50)
eso.fetal.combined.with.ciliated.integrated <- RunUMAP(eso.fetal.combined.with.ciliated.integrated, dims = 1:14, reduction.name = "umap", return.model = TRUE)
eso.fetal.combined.with.ciliated.integrated <- FindNeighbors(eso.fetal.combined.with.ciliated.integrated, dims = 1:14)
eso.fetal.combined.with.ciliated.integrated <- FindClusters(eso.fetal.combined.with.ciliated.integrated, resolution = 0.3)
####Figure S1A
pdf(file.path("./", paste0("Figure S1A Left", ".pdf")), w=11, h=8.5)
DimPlot(eso.fetal.combined.with.ciliated.integrated, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
dev.off()
DefaultAssay(eso.fetal.combined.with.ciliated.integrated) <- "RNA"
eso.fetal.combined.with.ciliated.integrated <- NormalizeData(eso.fetal.combined.with.ciliated.integrated)
eso.fetal.combined.with.ciliated.integrated <- JoinLayers(eso.fetal.combined.with.ciliated.integrated)
pdf(file.path("./", paste0("Figure S1A Middle and Right", ".pdf")), w=22, h=8.5)
FeaturePlot(eso.fetal.combined.with.ciliated.integrated, features = c("CDH1", "VIM"),  pt.size = 0.2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

###Figure S1B
eso.fetal.combined.with.ciliated.integrated.named <- RenameIdents(eso.fetal.combined.with.ciliated.integrated, "0" = "Epibasal", "1" = "KRT4+ Luminal", "2" = "Basal", "3" = "GDF10+ Mesenchyme", "4" = "WNT2B+ Mesenchyme", "5" = "Myofibroblast", "6" = "SCGB1A1+HOPX+ Luminal", "7" = "Neurons", "8" = "Skeletal Muscle Precursors", "9" = "Schwann Cells", "10" = "Multiciliated", "11" = "T-Cells", "12" = "Pericytes", "13" = "Endothelial", "14" = "Proliferative Mesenchyme")
eso.fetal.combined.with.ciliated.integrated  <- AddMetaData(eso.fetal.combined.with.ciliated.integrated, col.name = "annotation_lvl1", Idents(eso.fetal.combined.with.ciliated.integrated.named))
eso.fetal.combined.with.ciliated.integrated$annotation_lvl1 <- factor(eso.fetal.combined.with.ciliated.integrated$annotation_lvl1, levels = c("Basal", "Epibasal",   "KRT4+ Luminal", "SCGB1A1+HOPX+ Luminal", "Multiciliated", "Pericytes", "Myofibroblast", "WNT2B+ Mesenchyme", "GDF10+ Mesenchyme","Endothelial","Skeletal Muscle Precursors",     "Neurons", "Schwann Cells",  "T-Cells", "B-Cells", "Macrophages", "Proliferative Mesenchyme"))
eso.fetal.discriminating.markers.simple <- c("TP63", "COL17A1", "LY6D", "KRT4", "SCGB1A1", "HOPX", "FOXJ1","COL1A1", "PDGFRA", "PDGFRB",  "ACTA2", "TAGLN","CXCL14" , "WNT2B", "ANGPTL1", "GDF10", "CD34", "CDH5", "PAX7", "TNNT1", "TUBB3.1", "ELAVL4", "GAP43", "SOX10", "MPZ", "CDH19", "PTPRC", "CD3D", "CD7", "MS4A1", "CD19", "CD68", "CD163", "MRC1", "MARCO", "TOP2A")
pdf(file.path("./", paste0("Figure S1B Left", ".pdf")), w=11, h=8.5)
DimPlot(eso.fetal.combined.with.ciliated.integrated, group.by = "annotation_lvl1", reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()
pdf(file.path("./", paste0("Figure S1B Right", ".pdf")), w=38, h=13)
Dotplot_Zhiwei_Version(eso.fetal.combined.with.ciliated.integrated, gene_list = eso.fetal.discriminating.markers.simple, "annotation_lvl1") + theme(axis.text = element_text(size = 30))
dev.off()


####Process Adult Data
setwd("~/scRNAseq/Data/Cell/Daysha/Adult/Sample_HT239_ESO/filtered_feature_bc_matrix") #file:///run/user/1002/gvfs/smb-share:server=umms-spencejr-win.turbo.storage.umich.edu,share=umms-spencejr/01_RAW_RNASEQ_AGC_SHARE/Run_2598/Sample_HT239_ESO/outs/web_summary.html
ht239eso.data <- Read10X(data.dir = "./")
setwd("~/scRNAseq/Data/Cell/Daysha/Adult/Sample_HT328_ESO/filtered_feature_bc_matrix") #file:///run/user/1002/gvfs/smb-share:server=umms-spencejr-win.turbo.storage.umich.edu,share=umms-spencejr/01_RAW_RNASEQ_AGC_SHARE/NovaA-147/Sample_HT328_ESO/outs/web_summary.html
ht328eso.data <- Read10X(data.dir = "./")

ht239eso <- CreateSeuratObject(counts = ht239eso.data, project = "ht239eso", min.cells = 3, min.features = 200) 
ht328eso <- CreateSeuratObject(counts = ht328eso.data, project = "ht328eso", min.cells = 3, min.features = 200)

ht239eso <- RenameCells(object = ht239eso, add.cell.id = "ht239eso")
ht328eso <- RenameCells(object = ht328eso, add.cell.id = "ht328eso")

ht239eso[["percent.mt"]] <- PercentageFeatureSet(ht239eso, pattern = "^MT-") 
ht328eso[["percent.mt"]] <- PercentageFeatureSet(ht328eso, pattern = "^MT-") 

ht239eso <- AddMetaData(ht239eso, "adult", "age")
ht328eso <- AddMetaData(ht328eso, "adult", "age")


eso.adult.combined <- merge(ht239eso, ht328eso)
eso.adult.combined <- JoinLayers(eso.adult.combined)
eso.adult.combined <- subset(eso.adult.combined, subset = nFeature_RNA > 1500 & nFeature_RNA < 7500 & percent.mt < 20)
ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = eso.adult.combined), value = TRUE);
percent.ribo <- Matrix::colSums(x = GetAssayData(object = eso.adult.combined, slot = 'counts')[ribo.genes, ]) / Matrix::colSums(x = GetAssayData(object = eso.adult.combined, slot = 'counts'));

eso.adult.combined[['percent.ribo']] <- percent.ribo;
eso.adult.combined <- SplitObject(eso.adult.combined, split.by = "orig.ident")
for (i in seq_along(eso.adult.combined)) {
  eso.adult.combined[[i]] <- NormalizeData(eso.adult.combined[[i]]) %>% FindVariableFeatures()
}
features <- SelectIntegrationFeatures(eso.adult.combined)
for (i in seq_along(along.with = eso.adult.combined)) {
  eso.adult.combined[[i]] <- ScaleData(eso.adult.combined[[i]], features = features) %>% RunPCA(features = features)
}
anchors <- FindIntegrationAnchors(eso.adult.combined, anchor.features = features, reduction = "rpca", dims = 1:30)
eso.adult.combined.integrated <- IntegrateData(anchors, dims = 1:30)
DefaultAssay(eso.adult.combined.integrated) <- "integrated"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
eso.adult.combined.integrated <- CellCycleScoring(eso.adult.combined.integrated, s.features = s.genes, g2m.features = g2m.genes)
eso.adult.combined.integrated <- ScaleData(eso.adult.combined.integrated, vars.to.regress = c("S.Score", "G2M.Score"))
eso.adult.combined.integrated <- RunPCA(eso.adult.combined.integrated)
ElbowPlot(eso.adult.combined.integrated, ndims = 50)
eso.adult.combined.integrated <- RunUMAP(eso.adult.combined.integrated, dims = 1:12, reduction.name = "umap", return.model = TRUE)
eso.adult.combined.integrated <- FindNeighbors(eso.adult.combined.integrated, dims = 1:12)
eso.adult.combined.integrated <- FindClusters(eso.adult.combined.integrated, resolution = 0.15)
DimPlot(eso.adult.combined.integrated, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)
#Figure S1C
pdf(file.path("./", paste0("Figure S1C Left Panel", ".pdf")), w=11, h=8.5)
DimPlot(eso.adult.combined.integrated, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 2)
dev.off()
pdf(file.path("./", paste0("Figure S1C Right Panel", ".pdf")), w=11, h=8.5)
DimPlot(eso.adult.combined.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 2)
dev.off()

DefaultAssay(eso.adult.combined.integrated) <- "RNA"
eso.adult.combined.integrated <- NormalizeData(eso.adult.combined.integrated)
eso.adult.combined.integrated <- JoinLayers(eso.adult.combined.integrated)

#Figure S1D
pdf(file.path("./", paste0("Figure S1D Top and Bottom", ".pdf")), w=22, h=8.5)
FeaturePlot(eso.fetal.combined.with.ciliated.integrated, features = c("CDH1", "VIM"),  pt.size = 0.2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

#Extract and Subcluster All Adult Epithelial Cells
eso.adult.epithelial <- subset(eso.adult.combined.integrated, idents = c(0, 1, 2, 3, 4))
DefaultAssay(eso.adult.epithelial) <- "integrated"
eso.adult.epithelial <- RunPCA(eso.adult.epithelial, verbose = FALSE)
ElbowPlot(eso.adult.epithelial, ndims = 50)
eso.adult.epithelial <- RunUMAP(eso.adult.epithelial, reduction = "pca", dims = 1:12, return.model = TRUE)
eso.adult.epithelial <- FindNeighbors(eso.adult.epithelial, dims = 1:12)
eso.adult.epithelial <- FindClusters(eso.adult.epithelial, resolution = 0.15) 

###Figure 1A FeaturePlots
pdf(file.path("./", paste0("Adult TP63", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.adult.epithelial, features = "TP63", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Adult COL17A1", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.adult.epithelial, features = "COL17A1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Adult LY6D", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.adult.epithelial, features =  "LY6D" , pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Adult KRT4", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.adult.epithelial, features = "KRT4", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Adult HOPX", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.adult.epithelial, features = "HOPX", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Adult CRNN", ".pdf")), w=11, h=8.5)
DimPlot(eso.adult.epithelial, features = "CRNN",  pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Adult FOXJ1", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.adult.epithelial, features = "FOXJ1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Adult PCNA", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.adult.epithelial, features = "PCNA", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

###Figure 1A Louvain
eso.adult.epithelial <- RenameIdents(object = eso.adult.epithelial, '0' = "Mid", '1' = "Basal", '2' = "Early", '3' = "Proliferative", '4' = "Late")
eso.adult.epithelial <- AddMetaData(eso.adult.epithelial, eso.adult.epithelial@active.ident, "cell_type")
pdf(file.path("./", paste0("Figure 1A Labelled Louvain", ".pdf")), w=11, h=8.5)
DimPlot(eso.adult.epithelial, group.by = "cell_type", pt.size = 2, label = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

#Extract and Subcluster All Fetal Epithelial Cells
eso.fetal.epithelial.with.ciliated <- subset(eso.fetal.combined.with.ciliated.integrated, idents = c(1, 0, 2, 6, 10))
eso.fetal.epithelial.with.ciliated.cellids <- WhichCells(eso.fetal.epithelial.with.ciliated)
DefaultAssay(eso.fetal.epithelial.with.ciliated) <- "integrated"
eso.fetal.epithelial.with.ciliated <- RunPCA(eso.fetal.epithelial.with.ciliated, verbose = FALSE)
ElbowPlot(eso.fetal.epithelial.with.ciliated, ndims = 50)
eso.fetal.epithelial.with.ciliated <- RunUMAP(eso.fetal.epithelial.with.ciliated, reduction = "pca", dims = 1:6, return.model = TRUE)
eso.fetal.epithelial.with.ciliated <- FindNeighbors(eso.fetal.epithelial.with.ciliated, dims = 1:6)
eso.fetal.epithelial.with.ciliated <- FindClusters(eso.fetal.epithelial.with.ciliated, resolution = 0.4) #use 0.4 to seperate out proliferative cluster in basal cells
DimPlot(eso.fetal.epithelial.with.ciliated, group.by = c("orig.ident", "ident"), label = TRUE, pt.size = 0.5)
eso.fetal.epithelial.with.ciliated$seurat_clusters <- factor(eso.fetal.epithelial.with.ciliated$seurat_clusters, levels = c("2",  "1", "3", "0", "5", "7", "6", "4"))

DefaultAssay(eso.fetal.epithelial.with.ciliated) <- "RNA"
eso.fetal.epithelial.with.ciliated <- NormalizeData(eso.fetal.epithelial.with.ciliated)

###Figure 1B FeaturePlots
pdf(file.path("./", paste0("Fetal TP63", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.fetal.epithelial.with.ciliated, features = "TP63", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Fetal COL17A1", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.fetal.epithelial.with.ciliated, features = "COL17A1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Fetal LY6D", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.fetal.epithelial.with.ciliated, features =  "LY6D" , pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Fetal KRT4", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.fetal.epithelial.with.ciliated, features = "KRT4", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Fetal HOPX", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.fetal.epithelial.with.ciliated, features = "HOPX", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Fetal CRNN", ".pdf")), w=11, h=8.5)
DimPlot(eso.fetal.epithelial.with.ciliated, "CRNN", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Fetal FOXJ1", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.fetal.epithelial.with.ciliated, features = "FOXJ1", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Fetal PCNA", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.fetal.epithelial.with.ciliated, features = "PCNA", pt.size = 2, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

##Figure 1B Louvain (identities aligned by hand to adult)
eso.fetal.epithelial.with.ciliated <- RenameIdents(eso.fetal.epithelial.with.ciliated, "0" = "Middle", "1" = "Early 1", "2" = "Basal", "3" = "Early 2", "4" = "Proliferative 2", "5" = "Late", "6" = "Proliferative 1", "7" = "Luminal Multiciliated")
eso.fetal.epithelial.with.ciliated <- AddMetaData(eso.fetal.epithelial.with.ciliated, col.name = "annotation_lvl3", Idents(eso.fetal.epithelial.with.ciliated))
eso.fetal.epithelial.with.ciliated <- AddMetaData(eso.fetal.epithelial.with.ciliated, col.name = "annotation_lvl1", Idents(eso.fetal.epithelial.with.ciliated))

pdf(file.path("./", paste0("Figure 1B Labelled Louvain", ".pdf")), w=11, h=8.5)
DimPlot(eso.fetal.epithelial.with.ciliated, group.by = "annotation_lvl3", pt.size = 2, label = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

###Transfer Labels from Adult to Fetal Data
common.features <- intersect(rownames(eso.adult.epithelial), rownames(eso.fetal.epithelial.with.ciliated))
length(x = common.features)
eso.adult.epithelial.anchors <- FindTransferAnchors(
  reference = eso.adult.epithelial,
  query = eso.fetal.epithelial.with.ciliated,
  features = rownames(eso.adult.epithelial),
  
)
eso.fetal.epithelial.predictions <- TransferData(anchorset = eso.adult.epithelial.anchors, refdata = eso.adult.epithelial$cell_type,  store.weights = TRUE)
eso.fetal.epithelial.with.ciliated <- AddMetaData(eso.fetal.epithelial.with.ciliated, metadata =eso.fetal.epithelial.predictions)

#Figure 1D
pdf(file.path("./", paste0("Figure 1D Louvain", ".pdf")), w=11, h=8.5)
DimPlot(eso.fetal.epithelial.predictions, group.by = "predicted.id", pt.size = 2, label = TRUE) 
dev.off()

##Compare By Age within Label Transfered Annotations between Fetal and Adult to Define Markers for Figure 1E
eso.adult.basal.subset <- subset(eso.adult.epithelial, idents = "Basal")
eso.fetal.epithelial.basal.label <- subset(eso.fetal.epithelial.with.ciliated, subset = predicted.id == "Basal")

eso.basal <- merge(eso.adult.basal.subset, eso.fetal.epithelial.basal.label)
DefaultAssay(eso.basal) <- "RNA"
eso.basal <- NormalizeData(eso.basal)
Idents(eso.basal) <- "age"
eso.basal.adultvfetal <- FindMarkers(eso.basal, ident.1 = "adult", ident.2 = "fetal", min.pct = 0.25)
write.csv(eso.basal.adultvfetal, "AdultvFetalBasal.csv")

##early cell comparison
eso.adult.early.subset <- subset(eso.adult.epithelial, idents = "Early")
eso.fetal.epithelial.early.label <- subset(eso.fetal.epithelial.with.ciliated, subset = predicted.id == "Early")
eso.early <- merge(eso.adult.early.subset, eso.fetal.epithelial.early.label)
DefaultAssay(eso.early) <- "RNA"
eso.early <- NormalizeData(eso.early)
Idents(eso.early) <- "age"
eso.early.adultvfetal <- FindMarkers(eso.early, ident.1 = "adult", ident.2 = "fetal", min.pct = 0.25)
write.csv(eso.early.adultvfetal, "AdultvFetalEarly.csv")

#middle cell comparison
eso.adult.mid.subset <- subset(eso.adult.epithelial, idents = "Mid")
eso.fetal.epithelial.mid.label <- subset(eso.fetal.epithelial.with.ciliated, subset = predicted.id == "Mid")
eso.mid <- merge(eso.adult.mid.subset, eso.fetal.epithelial.mid.label)
DefaultAssay(eso.mid) <- "RNA"
eso.mid <- NormalizeData(eso.mid)
Idents(eso.mid) <- "age"
eso.mid.adultvfetal <- FindMarkers(eso.mid, ident.1 = "adult", ident.2 = "fetal", min.pct = 0.25)
write.csv(eso.mid.adultvfetal, "AdultvFetalMid.csv")

#luminal or late cell comparison
eso.adult.late.subset <- subset(eso.adult.epithelial, idents = "Late")
eso.fetal.epithelial.late.label <- subset(eso.fetal.epithelial.with.ciliated, idents = 5)
eso.late <- merge(eso.adult.late.subset, eso.fetal.epithelial.late.label)
DefaultAssay(eso.late) <- "RNA"
eso.late <- NormalizeData(eso.late)
Idents(eso.late) <- "age"
eso.late.adultvfetal <- FindMarkers(eso.late, ident.1 = "adult", ident.2 = "fetal", min.pct = 0.25)
write.csv(eso.late.adultvfetal, "AdultvFetalLate.csv")

#Figure 1E
eso.adult.v.fetal.goi <- c("KRT14", "MAFF", "IGF2", "SCGB1A1","TP63", "COL17A1", "CAV1", "CXCL14", "GPC3", "KRT17", "LY6D", "PLAUR", "DLK1", "CLDN10", "KRT4", "KRT16", "CRABP2", "KRT7", "UPK1B", "HOPX", "CRNN", "KRT78", "TSPAN1", "CRIP1", "FOXJ1", "TPP3", "CDC20B")
eso.adult.basal.subset <- AddMetaData(eso.adult.basal.subset, "Basal-Adult", "name")
eso.adult.early.subset <- AddMetaData(eso.adult.early.subset, "Early-Adult", "name")
eso.adult.mid.subset <- AddMetaData(eso.adult.mid.subset, "Middle-Adult", "name")
eso.adult.late.subset <- AddMetaData(eso.adult.late.subset, "Luminal-Adult", "name")
eso.fetal.epithelial.basal.label <- AddMetaData(eso.fetal.epithelial.basal.label, "Basal-Fetal", "name")
eso.fetal.epithelial.early.label <- AddMetaData(eso.fetal.epithelial.early.label, "Early-Fetal", "name")
eso.fetal.epithelial.mid.label <- AddMetaData(eso.fetal.epithelial.mid.label, "Middle-Fetal", "name")
eso.fetal.epithelial.late.label <- AddMetaData(eso.fetal.epithelial.late.label, "Luminal-Fetal", "name")
eso.fetal.ciliated.subset <- subset(eso.fetal.epithelial.with.ciliated, ident = 7)
eso.fetal.ciliated.subset <- AddMetaData(eso.fetal.ciliated.subset, "Multiciliated-Fetal", "name")

eso.dotplot <- merge(eso.adult.basal.subset, c(eso.fetal.ciliated.subset, eso.fetal.epithelial.late.label,eso.fetal.epithelial.mid.label, eso.fetal.epithelial.early.label,eso.fetal.epithelial.basal.label,eso.adult.late.subset,eso.adult.mid.subset, eso.adult.early.subset))
eso.dotplot$name <- factor(eso.dotplot$name, levels = c("Multiciliated-Fetal", "Luminal-Fetal", "Middle-Fetal", "Early-Fetal", "Basal-Fetal", "Luminal-Adult", "Middle-Adult", "Early-Adult", "Basal-Adult"))

pdf(file.path("./", paste0("Figure 1E Dotplot", ".pdf")), w=20, h=7.0)
Dotplot_Zhiwei_Version(eso.dotplot, eso.adult.v.fetal.goi)
dev.off()


##Figure 2
###Extract and Subcluster Fetal Mesenchyme
eso.fetal.mesenchyme <- subset(eso.fetal.epithelial.with.ciliated, idents = c(3, 14, 4, 5))
DefaultAssay(eso.fetal.mesenchyme) <- "integrated"
eso.fetal.mesenchyme <- ScaleData(eso.fetal.mesenchyme, verbose = FALSE)
eso.fetal.mesenchyme <- RunPCA(eso.fetal.mesenchyme, verbose = FALSE)
ElbowPlot(eso.fetal.mesenchyme, ndims = 50)
eso.fetal.mesenchyme <- RunUMAP(eso.fetal.mesenchyme, reduction = "pca", dims = 1:10, return.model = TRUE)
eso.fetal.mesenchyme <- FindNeighbors(eso.fetal.mesenchyme, dims = 1:10)
eso.fetal.mesenchyme <- FindClusters(eso.fetal.mesenchyme, resolution = 0.4) 
eso.fetal.mesenchyme$seurat_clusters <- factor(eso.fetal.mesenchyme$seurat_clusters, levels = c( "4", "6","1", "0", "3", "2", "7", "5"))
eso.fetal.mesenchyme.marker <- c("PDGFRA", "GDF10", "CD34", "SPOCK3", "NRXN1", "HGF", "ETV1", "KIT", "FGL2", "KCNN3", "PEG3",  "SLC26A7", "ARHGAP15", "WNT2B", "PTGER1", "APOE", "NGFR", "PTCH2", "EBF2", "POSTN", "CXCL14","CDH6","LPPR4","VWC2","MYH11", "DES", "ACTA2", "TOP2A")
eso.fetal.mesenchyme.named <- RenameIdents(eso.fetal.mesenchyme, "4" = "APOE+ Mesenchyme 1", "6" = "KIT+ Mesenchyme 2", "1" = "WNT2B+ Mesenchyme", "0" = "KCNN3+ Mesenchyme", "3" = "NRXN1+ Mesenchyme 2", "2" = "Smooth Muscle Precursor", "7" = "Smooth Muscle", "5" = "Proliferative Mesenchyme") 
eso.fetal.mesenchyme <- AddMetaData(eso.fetal.mesenchyme, col.name = "annotation_lvl2r", Idents(eso.fetal.mesenchyme.named))
eso.fetal.mesenchyme$annotation_lvl2r <- factor(eso.fetal.mesenchyme$annotation_lvl2r, levels = c("NRXN1+ Mesenchyme", "KIT+ Mesenchyme", "KCNN3+ Mesenchyme", "WNT2B+ Mesenchyme", "APOE+ Mesenchyme", "Smooth Muscle Precursor", "Smooth Muscle", "Proliferative Mesenchyme"))

pdf(file.path("./", paste0("Figure 2A", ".pdf")), w=11, h=8.5)
DimPlot(eso.fetal.mesenchyme, reduction = "umap", group.by = "annotation_lvl2r", label = TRUE, pt.size = 2)
dev.off()

DefaultAssay(eso.fetal.mesenchyme) <- "RNA"
eso.fetal.mesenchyme <- NormalizeData(eso.fetal.mesenchyme)

pdf(file.path("./", paste0("Figure 2B.pdf")), w=19.5, h=6)
Dotplot_Zhiwei_Version(eso.fetal.mesenchyme, gene_list = eso.fetal.mesenchyme.marker, "annotation_lvl2r") + theme(axis.text.x = element_text(face = "italic", size = 16))
dev.off()

pdf(file.path("./", paste0("Fetal ESO Mesenchyme PDGFRA", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.fetal.mesenchyme, features = c("PDGFRA"), pt.size = 0.5, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Fetal ESO Mesenchyme GDF10", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.fetal.mesenchyme, features = c("GDF10"), pt.size = 0.5, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

pdf(file.path("./", paste0("Fetal ESO Mesenchyme HGF", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.fetal.mesenchyme, features = c("HGF"), pt.size = 1, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

pdf(file.path("./", paste0("Fetal ESO Mesenchyme KIT", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.fetal.mesenchyme, features = c("KIT"), pt.size = 1, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

pdf(file.path("./", paste0("Fetal ESO Mesenchyme WNT2B", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.fetal.mesenchyme, features = c("WNT2B"), pt.size = 1, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Fetal ESO Mesenchyme VWC2", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.fetal.mesenchyme, features = c("VWC2"), pt.size = 1, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("Fetal ESO Mesenchyme TOP2A", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.fetal.mesenchyme, features = c("TOP2A"), pt.size = 1, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

##Process Fetal ESO In Vitro 2D Condition
setwd("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Data/Cell/Daysha/In Vitro/DF-11431-1/filtered_feature_bc_matrix")
eso.vitro.2d.1.data <- Read10X(data.dir = "./") 
setwd("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Data/Cell/Daysha/In Vitro/DF-11290-1/filtered_feature_bc_matrix")
eso.vitro.2d.3.data <- Read10X(data.dir = "./") 

eso.vitro.2d.1 <- CreateSeuratObject(counts = eso.vitro.2d.1.data, project = "eso.vitro.2d.1", min.cells = 3, min.features = 200) 
eso.vitro.2d.3 <- CreateSeuratObject(counts = eso.vitro.2d.3.data, project = "eso.vitro.2d.3", min.cells = 3, min.features = 200)

eso.vitro.2d.1 <- RenameCells(object = eso.vitro.2d.1, add.cell.id = "eso.vitro.2d.1")
eso.vitro.2d.3 <- RenameCells(object = eso.vitro.2d.3, add.cell.id = "eso.vitro.2d.3")
eso.vitro.2d.1[["percent.mt"]] <- PercentageFeatureSet(eso.vitro.2d.1, pattern = "^MT-") 
eso.vitro.2d.3[["percent.mt"]] <- PercentageFeatureSet(eso.vitro.2d.3, pattern = "^MT-") 
eso.vitro.2d.1 <- AddMetaData(eso.vitro.2d.1, "fetal", "age")
eso.vitro.2d.3 <- AddMetaData(eso.vitro.2d.3, "fetal", "age")
eso.vitro.2d.1 <- AddMetaData(eso.vitro.2d.1, "2D", "format")
eso.vitro.2d.3 <- AddMetaData(eso.vitro.2d.3, "2D", "format")
eso.vitro.2d.1 <- AddMetaData(eso.vitro.2d.1, "651", "ht")
eso.vitro.2d.3 <- AddMetaData(eso.vitro.2d.3, "514", "ht")
eso.vitro.2d.1 <- AddMetaData(eso.vitro.2d.1, col.name = "ID", "esovitro1")
eso.vitro.2d.1 <- AddMetaData(eso.vitro.2d.1, col.name = "Type", "In Vitro")
eso.vitro.2d.3 <- AddMetaData(eso.vitro.2d.3, col.name = "ID", "esovitro3")
eso.vitro.2d.3 <- AddMetaData(eso.vitro.2d.3, col.name = "Type", "In Vitro")
eso.vitro.2d <- merge(eso.vitro.2d.1, eso.vitro.2d.3)
eso.vitro.2d <- JoinLayers(eso.vitro.2d)

ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = eso.vitro.2d), value = TRUE);
percent.ribo <- Matrix::colSums(x = GetAssayData(object = eso.vitro.2d, slot = 'counts')[ribo.genes, ]) / Matrix::colSums(x = GetAssayData(object = eso.vitro.2d, slot = 'counts'));
eso.vitro.2d[['percent.ribo']] <- percent.ribo;
eso.vitro.2d <- subset(eso.vitro.2d, subset = nFeature_RNA > 1500 & nFeature_RNA < 10000 & percent.mt < 20 & nCount_RNA > 10000)

eso.vitro.2d <- SplitObject(eso.vitro.2d, split.by = "orig.ident")
for (i in seq_along(eso.vitro.2d)) {
  eso.vitro.2d[[i]] <- NormalizeData(eso.vitro.2d[[i]]) %>% FindVariableFeatures()
}
features <- SelectIntegrationFeatures(eso.vitro.2d)
for (i in seq_along(along.with = eso.vitro.2d)) {
  eso.vitro.2d[[i]] <- ScaleData(eso.vitro.2d[[i]], features = features) %>% RunPCA(features = features)
}

anchors <- FindIntegrationAnchors(eso.vitro.2d, anchor.features = features, reduction = "rpca", dims = 1:30)
eso.vitro.2d.integrated <- IntegrateData(anchors, dims = 1:30, k.weight = 19)
DefaultAssay(eso.vitro.2d.integrated) <- "integrated"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
eso.vitro.2d.integrated <- CellCycleScoring(eso.vitro.2d.integrated, s.features = s.genes, g2m.features = g2m.genes)
eso.vitro.2d.integrated <- ScaleData(eso.vitro.2d.integrated, vars.to.regress = c("S.Score", "G2M.Score"))
eso.vitro.2d.integrated <- RunPCA(eso.vitro.2d.integrated)
ElbowPlot(eso.vitro.2d.integrated, ndims = 50)
eso.vitro.2d.integrated <- RunUMAP(eso.vitro.2d.integrated, dims = 1:14, reduction.name = "umap", return.model = TRUE)
eso.vitro.2d.integrated <- FindNeighbors(eso.vitro.2d.integrated, dims = 1:14)
eso.vitro.2d.integrated <- FindClusters(eso.vitro.2d.integrated, resolution = 0.3)

DefaultAssay(eso.vitro.2d.integrated) <- "RNA"
eso.vitro.2d.integrated <- NormalizeData(eso.vitro.2d.integrated)
eso.vitro.2d.integrated <- JoinLayers(eso.vitro.2d.integrated)

##Transfer Labels from Fetal In Vivo onto Fetal 2D In Vitro Data
common.features <- intersect(rownames(fetal.eso), rownames(eso.vitro.2d.integrated))
length(x = common.features)
fetal.eso.anchors <- FindTransferAnchors(
  reference = fetal.eso,
  query = eso.vitro.2d.integrated,
  reference.assay = "integrated",
  query.assay = "integrated",
  features = common.features,
  reference.reduction = "pca",
  k.filter = 200
)
invitro.mapped.invivo.epithelium <- MapQuery(
  anchorset = fetal.eso.anchors,
  query = eso.vitro.2d.integrated,
  reference = fetal.eso,
  refdata = "annotation_lvl3",
  reference.reduction = "pca", 
  reduction.model = "umap",
  )

eso.vitro.2d.integrated <- AddMetaData(eso.vitro.2d.integrated, col.name = "predicted.id", invitro.mapped.invivo.epithelium$predicted.id)

#Figure 3G
pdf(file.path("./", paste0("Figure 3G Louvain with In Vivo Fetal Predicted IDs", ".pdf")), w=11, h=8.5)
DimPlot(eso.vitro.2d.integrated, group.by = "predicted.id", pt.size = 2, label = TRUE) 
dev.off()

#Figure 3H
pdf(file.path("./", paste0("CDH1", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.vitro.2d.integrated, features = c("CDH1"),  pt.size = 1.0, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("TP63", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.vitro.2d.integrated, features = c("TP63"),  pt.size = 1.0, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("COL17A1", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.vitro.2d.integrated, features = c("COL17A1"),  pt.size = 1.0, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("KI67", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.vitro.2d.integrated, features = c("MKI67"),  pt.size = 1.0, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("KRT4", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.vitro.2d.integrated, features = c("KRT4"),  pt.size = 1.0, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("VIM", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.vitro.2d.integrated, features = c("VIM"),  pt.size = 1.0, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("WNT2B", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.vitro.2d.integrated, features = c("WNT2B"),  pt.size = 1.0, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("HGF", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.vitro.2d.integrated, features = c("HGF"),  pt.size = 1.0, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("PDGFRA", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.vitro.2d.integrated, features = c("PDGFRA"),  pt.size = 1.0, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()


##Figure 3I Module Scoring (In Vivo Modules on 2D In Vitro Data)
eso.fetal.combined.with.ciliated.integrated <- AddMetaData(eso.fetal.combined.with.ciliated.integrated, col.name = "annotation_lvl3", c(eso.fetal.mesenchyme$annotation_lvl2r, eso.fetal.epithelial.with.ciliated$annotation_lvl1))
Idents(eso.fetal.combined.with.ciliated.integrated) <- "annotation_lvl3"
markers <- FindAllMarkers(eso.fetal.combined.with.ciliated.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, features = trimmed.features)

##Extract Top 100 Enriched Gene per Cluster by Log2FC
top20_markers <- markers %>%
  group_by(cluster) %>%
  top_n(100, avg_log2FC)

marker_lists <- split(top20_markers$gene, top20_markers$cluster) 
DefaultAssay(eso.vitro.2d.integrated) <- "RNA"

eso.vitro.2d.integrated <- NormalizeData(eso.vitro.2d.integrated)

valid_marker_lists <- lapply(marker_lists, function(genes) {
  intersect(genes, rownames(eso.vitro.2d.integrated))
})

# Add module scores for each Cluster
for (cluster_name in names(valid_marker_lists)) {
  eso.vitro.2d.integrated <- AddModuleScore(
    object = eso.vitro.2d.integrated,
    features = list(valid_marker_lists[[cluster_name]]),
    name = paste0("ClusterScore_", cluster_name)
  )
}

#Reorder Clusters in Violin Plots
eso.vitro.2d.integrated$ordered_cluster <- factor(eso.vitro.2d.integrated$predicted.id, levels = c("Basal", "Epibasal", "KRT4+ Luminal", "SCGB1A1+HOPX+ Luminal","Multiciliated", "Pericytes", "APOE+ Mesenchyme", "Proliferative Mesenchyme"))
Idents(eso.vitro.2d.integrated) <- "ordered_cluster"

#Define Plot Categories
mesenchyme_titles <- c("NRXN1+ Mesenchyme", "KIT+ Mesenchyme", "KCNN3+ Mesenchyme", 
                       "WNT2B+ Mesenchyme", "APOE+ Mesenchyme", "Smooth Muscle Precursor", 
                       "Smooth Muscle", "Proliferative Mesenchyme")
basal_titles <- c("Basal")
epibasal_titles <- c("Epibasal")
luminal_titles <- c("KRT4+ Luminal", "SCGB1A1+HOPX+ Luminal")
multiciliated_titles <- c("Multiciliated")
other_titles <- c("Pericytes", "Endothelial", "Skeletal Muscle Precursors", "Neurons", 
                  "Schwann Cells", "T-Cells", "B-Cells", "Macrophages")

#Generate and Sort Violin Plots
feature_names <- paste0("ClusterScore_", names(valid_marker_lists), "1")
formatted_titles <- gsub("ClusterScore_|1$", "", feature_names)  # Remove extra text
formatted_titles <- paste0(formatted_titles, "\nScore")  # Add "Score" on a new line

#Create violin plots
plots <- VlnPlot(eso.vitro.2d.integrated, features = feature_names, pt.size = 0, combine = FALSE)

# Extract plot titles without "Score"
plot_titles <- gsub("\nScore", "", formatted_titles)

# Create an ordered index based on defined categories
ordered_indices <- c(
  which(plot_titles %in% mesenchyme_titles),
  which(plot_titles %in% basal_titles),
  which(plot_titles %in% epibasal_titles),
  which(plot_titles %in% luminal_titles),
  which(plot_titles %in% multiciliated_titles),
  which(plot_titles %in% other_titles)
)

# Reorder plots and titles
plots <- plots[ordered_indices]
formatted_titles <- formatted_titles[ordered_indices]

#  Adjust Each Plot
plots <- lapply(seq_along(plots), function(i) {
  plots[[i]] +
    geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA, alpha = 0.7, coef = 0) +  # Box without whiskers
    theme_minimal() +  # Clean layout
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"), 
      axis.text.y = element_text(size = 12, color = "black"), 
      axis.title.y = element_text(size = 14, color = "black"),  # No bold, slightly smaller
      axis.title.x = element_text(size = 14, color = "black", vjust = -0.1),  # Move x-axis label closer
      plot.title = element_text(size = 14, hjust = 0.5, lineheight = 1.1), # Center title & reduce size
      legend.position = "none",
      axis.line = element_line(size = 1, color = "black"),  # Darker axis lines
      axis.ticks = element_line(size = 0.8, color = "black"),    # Darker ticks
      plot.margin = margin(t = 10, b = 40)  # Extra spacing between plots
    ) +
    labs(
      x = "In Vitro Cluster",  # Label x-axis
      title = formatted_titles[i]  # Apply formatted title with "Score" on second line
    ) +
    scale_fill_manual(values = c("#00BFC4", "#DE8C00", "#00B4F0", "#7CAE00", "#00BA38", "#F8766D", "#B79F00")) +  # Seurat default colors
    ylim(NA, 1.5)  # **Set upper y-axis limit to 1.0**
})

# Arrange in Two Rows 
num_plots <- length(plots)
num_cols <- ceiling(num_plots / 2)  # Two rows

plot_layout <- wrap_plots(plots, ncol = num_cols) + 
  plot_layout(guides = "collect") + 
  plot_annotation(theme = theme(plot.margin = margin(t = 80, b = 80)))  # More space between rows

#Figure 3I
pdf(file = "./ESO_2D_Module_Scores_Top_100_from_In_Vivo_Reference_By_Cluster.pdf", width = 30, height = 12)
print(plot_layout)
dev.off()

#Figure 4
#Process 3D Esophagus In Vitro
setwd("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Data/Cell/Daysha/In Vitro/DF-11431-2/filtered_feature_bc_matrix")
eso.vitro.3d.2.data <- Read10X(data.dir = "./") 
setwd("~/University of Michigan Dropbox/Tristan Frum/0_Mac_scRNA_Seq_Processing/Data/Cell/Daysha/In Vitro/DF-11290-2/filtered_feature_bc_matrix")
eso.vitro.3d.4.data <- Read10X(data.dir = "./") 
eso.vitro.3d.2 <- CreateSeuratObject(counts = eso.vitro.3d.2.data, project = "eso.vitro.3d.2", min.cells = 3, min.features = 200)
eso.vitro.3d.4 <- CreateSeuratObject(counts = eso.vitro.3d.4.data, project = "eso.vitro.3d.4", min.cells = 3, min.features = 200)
eso.vitro.3d.2 <- RenameCells(object = eso.vitro.3d.2, add.cell.id = "eso.vitro.3d.2")
eso.vitro.3d.4 <- RenameCells(object = eso.vitro.3d.4, add.cell.id = "eso.vitro.3d.4")
eso.vitro.3d.2[["percent.mt"]] <- PercentageFeatureSet(eso.vitro.3d.2, pattern = "^MT-") 
eso.vitro.3d.4[["percent.mt"]] <- PercentageFeatureSet(eso.vitro.3d.4, pattern = "^MT-") 
eso.vitro.3d.2 <- AddMetaData(eso.vitro.3d.2, "fetal", "age")
eso.vitro.3d.4 <- AddMetaData(eso.vitro.3d.4, "fetal", "age")
eso.vitro.3d.2 <- AddMetaData(eso.vitro.3d.2, "3D", "format")
eso.vitro.3d.4 <- AddMetaData(eso.vitro.3d.4, "3D", "format")
eso.vitro.3d.2 <- AddMetaData(eso.vitro.3d.2, "651", "ht")
eso.vitro.3d.4 <- AddMetaData(eso.vitro.3d.4, "514", "ht")
eso.vitro.3d.2 <- AddMetaData(eso.vitro.3d.2, col.name = "ID", "esovitro2")
eso.vitro.3d.2 <- AddMetaData(eso.vitro.3d.2, col.name = "Type", "In Vitro")
eso.vitro.3d.4 <- AddMetaData(eso.vitro.3d.4, col.name = "ID", "esovitro4")
eso.vitro.3d.4 <- AddMetaData(eso.vitro.3d.4, col.name = "Type", "In Vitro")
so.vitro.3d <- merge(eso.vitro.3d.2, eso.vitro.3d.4)
eso.vitro.3d <- JoinLayers(eso.vitro.3d)
ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = eso.vitro.3d), value = TRUE);
percent.ribo <- Matrix::colSums(x = GetAssayData(object = eso.vitro.3d, slot = 'counts')[ribo.genes, ]) / Matrix::colSums(x = GetAssayData(object = eso.vitro.3d, slot = 'counts'));
eso.vitro.3d[['percent.ribo']] <- percent.ribo;


eso.vitro.3d <- subset(eso.vitro.3d, subset = nFeature_RNA > 1500 & nFeature_RNA < 10000 & percent.mt < 20 & nCount_RNA > 10000)
eso.vitro.3d <- SplitObject(eso.vitro.3d, split.by = "orig.ident")
for (i in seq_along(eso.vitro.3d)) {
  eso.vitro.3d[[i]] <- NormalizeData(eso.vitro.3d[[i]]) %>% FindVariableFeatures()
}
features <- SelectIntegrationFeatures(eso.vitro.3d)
for (i in seq_along(along.with = eso.vitro.3d)) {
  eso.vitro.3d[[i]] <- ScaleData(eso.vitro.3d[[i]], features = features) %>% RunPCA(features = features)
}

anchors <- FindIntegrationAnchors(eso.vitro.3d, anchor.features = features, reduction = "rpca", dims = 1:30)
eso.vitro.3d.integrated <- IntegrateData(anchors, dims = 1:30, k.weight = 19)
DefaultAssay(eso.vitro.3d.integrated) <- "integrated"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
eso.vitro.3d.integrated <- CellCycleScoring(eso.vitro.3d.integrated, s.features = s.genes, g2m.features = g2m.genes)
eso.vitro.3d.integrated <- ScaleData(eso.vitro.3d.integrated, vars.to.regress = c("S.Score", "G2M.Score"))
eso.vitro.3d.integrated <- RunPCA(eso.vitro.3d.integrated)
ElbowPlot(eso.vitro.3d.integrated, ndims = 50)
eso.vitro.3d.integrated <- RunUMAP(eso.vitro.3d.integrated, dims = 1:14, reduction.name = "umap", return.model = TRUE)
eso.vitro.3d.integrated <- FindNeighbors(eso.vitro.3d.integrated, dims = 1:14)
eso.vitro.3d.integrated <- FindClusters(eso.vitro.3d.integrated, resolution = 0.3)

DefaultAssay(eso.vitro.3d.integrated) <- "RNA"
eso.vitro.3d.integrated <- NormalizeData(eso.vitro.3d.integrated)
eso.vitro.3d.integrated <- JoinLayers(eso.vitro.3d.integrated)


##Label Transfer 
#Figure 4H
##Transfer Labels from Fetal In Vivo onto Fetal 3D In Vitro Data
common.features <- intersect(rownames(fetal.eso), rownames(eso.vitro.3d.integrated))
length(x = common.features)
fetal.eso.anchors <- FindTransferAnchors(
  reference = fetal.eso,
  query = eso.vitro.3d.integrated,
  reference.assay = "integrated",
  query.assay = "integrated",
  features = common.features,
  reference.reduction = "pca",
  k.filter = 200
)
invitro.mapped.invivo.epithelium <- MapQuery(
  anchorset = fetal.eso.anchors,
  query = eso.vitro.3d.integrated,
  reference = fetal.eso,
  refdata = "annotation_lvl3",
  reference.reduction = "pca", 
  reduction.model = "umap",
)
eso.vitro.3d.integrated <- AddMetaData(eso.vitro.3d.integrated, col.name = "predicted.id", invitro.mapped.invivo.epithelium$predicted.id)

#Figure 4H
pdf(file.path("./", paste0("Figure 4H Louvain with In Vivo Fetal Predicted IDs", ".pdf")), w=11, h=8.5)
DimPlot(eso.vitro.3d.integrated, group.by = "predicted.id", pt.size = 2, label = TRUE) 
dev.off()

#Figure 4I
pdf(file.path("./", paste0("CDH1", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.vitro.3d.integrated, features = c("CDH1"),  pt.size = 1.0, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("TP63", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.vitro.3d.integrated, features = c("TP63"),  pt.size = 1.0, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("COL17A1", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.vitro.3d.integrated, features = c("COL17A1"),  pt.size = 1.0, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("KI67", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.vitro.3d.integrated, features = c("MKI67"),  pt.size = 1.0, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("KRT4", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.vitro.3d.integrated, features = c("KRT4"),  pt.size = 1.0, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("VIM", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.vitro.3d.integrated, features = c("VIM"),  pt.size = 1.0, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("WNT2B", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.vitro.3d.integrated, features = c("WNT2B"),  pt.size = 1.0, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("HGF", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.vitro.3d.integrated, features = c("HGF"),  pt.size = 1.0, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("PDGFRA", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.vitro.3d.integrated, features = c("PDGFRA"),  pt.size = 1.0, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()
pdf(file.path("./", paste0("VWC2", ".pdf")), w=11, h=8.5)
FeaturePlot(eso.vitro.3d.integrated, features = c("VWC2"),  pt.size = 1.0, order = TRUE) & NoAxes() & NoLegend() & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58"))
dev.off()

##Figure 4J Module Scoring (In Vivo Modules on 3D In Vitro Data)
eso.fetal.combined.with.ciliated.integrated <- AddMetaData(eso.fetal.combined.with.ciliated.integrated, col.name = "annotation_lvl3", c(eso.fetal.mesenchyme$annotation_lvl2r, eso.fetal.epithelial.with.ciliated$annotation_lvl1))
Idents(eso.fetal.combined.with.ciliated.integrated) <- "annotation_lvl3"
markers <- FindAllMarkers(eso.fetal.combined.with.ciliated.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, features = trimmed.features)

##Extract Top 100 Enriched Gene per Cluster by Log2FC
top20_markers <- markers %>%
  group_by(cluster) %>%
  top_n(100, avg_log2FC)

marker_lists <- split(top20_markers$gene, top20_markers$cluster) 
DefaultAssay(eso.vitro.3d.integrated) <- "RNA"

eso.vitro.3d.integrated <- NormalizeData(eso.vitro.3d.integrated)

valid_marker_lists <- lapply(marker_lists, function(genes) {
  intersect(genes, rownames(eso.vitro.3d.integrated))
})

# Add module scores for each Cluster
for (cluster_name in names(valid_marker_lists)) {
  eso.vitro.3d.integrated <- AddModuleScore(
    object = eso.vitro.3d.integrated,
    features = list(valid_marker_lists[[cluster_name]]),
    name = paste0("ClusterScore_", cluster_name)
  )
}

#Reorder Clusters in Violin Plots
eso.vitro.3d.integrated$ordered_cluster <- factor(eso.vitro.3d.integrated$seurat_clusters, levels = c(3, 0, 1, 2))
Idents(eso.vitro.3d.integrated) <- "ordered_cluster"

#Define Plot Categories
mesenchyme_titles <- c("NRXN1+ Mesenchyme", "KIT+ Mesenchyme", "KCNN3+ Mesenchyme", 
                       "WNT2B+ Mesenchyme", "APOE+ Mesenchyme", "Smooth Muscle Precursor", 
                       "Smooth Muscle", "Proliferative Mesenchyme")
basal_titles <- c("Basal")
epibasal_titles <- c("Epibasal")
luminal_titles <- c("KRT4+ Luminal", "SCGB1A1+HOPX+ Luminal")
multiciliated_titles <- c("Multiciliated")
other_titles <- c("Pericytes", "Endothelial", "Skeletal Muscle Precursors", "Neurons", 
                  "Schwann Cells", "T-Cells", "B-Cells", "Macrophages")

#Generate and Sort Violin Plots
feature_names <- paste0("ClusterScore_", names(valid_marker_lists), "1")
formatted_titles <- gsub("ClusterScore_|1$", "", feature_names)  # Remove extra text
formatted_titles <- paste0(formatted_titles, "\nScore")  # Add "Score" on a new line

#Create violin plots
plots <- VlnPlot(eso.vitro.3d.integrated, features = feature_names, pt.size = 0, combine = FALSE)

# Extract plot titles without "Score"
plot_titles <- gsub("\nScore", "", formatted_titles)

# Create an ordered index based on defined categories
ordered_indices <- c(
  which(plot_titles %in% mesenchyme_titles),
  which(plot_titles %in% basal_titles),
  which(plot_titles %in% epibasal_titles),
  which(plot_titles %in% luminal_titles),
  which(plot_titles %in% multiciliated_titles),
  which(plot_titles %in% other_titles)
)

# Reorder plots and titles
plots <- plots[ordered_indices]
formatted_titles <- formatted_titles[ordered_indices]

#  Adjust Each Plot
plots <- lapply(seq_along(plots), function(i) {
  plots[[i]] +
    geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA, alpha = 0.7, coef = 0) +  # Box without whiskers
    theme_minimal() +  # Clean layout
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"), 
      axis.text.y = element_text(size = 12, color = "black"), 
      axis.title.y = element_text(size = 14, color = "black"),  # No bold, slightly smaller
      axis.title.x = element_text(size = 14, color = "black", vjust = -0.1),  # Move x-axis label closer
      plot.title = element_text(size = 14, hjust = 0.5, lineheight = 1.1), # Center title & reduce size
      legend.position = "none",
      axis.line = element_line(size = 1, color = "black"),  # Darker axis lines
      axis.ticks = element_line(size = 0.8, color = "black"),    # Darker ticks
      plot.margin = margin(t = 10, b = 40)  # Extra spacing between plots
    ) +
    labs(
      x = "In Vitro Cluster",  # Label x-axis
      title = formatted_titles[i]  # Apply formatted title with "Score" on second line
    ) +
    scale_fill_manual(values = c("#00BFC4", "#DE8C00", "#00B4F0", "#7CAE00", "#00BA38", "#F8766D", "#B79F00")) +  # Seurat default colors
    ylim(NA, 1.5)  # **Set upper y-axis limit to 1.0**
})

# Arrange in Two Rows 
num_plots <- length(plots)
num_cols <- ceiling(num_plots / 2)  # Two rows

plot_layout <- wrap_plots(plots, ncol = num_cols) + 
  plot_layout(guides = "collect") + 
  plot_annotation(theme = theme(plot.margin = margin(t = 80, b = 80)))  # More space between rows

#Figure 4J
pdf(file = "./ESO_3D_Module_Scores_Top_20_from_In_Vivo_Reference_By_Cluster.pdf", width = 30, height = 12)
print(plot_layout)
dev.off()