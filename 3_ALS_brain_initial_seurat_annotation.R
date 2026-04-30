options(repos = c(REPO_NAME = "https://packagemanager.posit.co/cran/__linux__/centos7/latest"))


setwd("/home/mpg02/MMNB/yigit.babal/workspace_ykb/fibroblasts/MIT_sc_dataset")
print(getwd())


source('/home/mpg08/mmnb1ngs/collab_work/workspace/ykbabal/library.R')
options(future.globals.maxSize = 3e+09)


load('whole_dataset/pmi_cutoff/workspace_checkpoint_5.RData')

seurat


seurat <- FindClusters(object = seurat, resolution = 1.5)

png(filename = "whole_dataset/pmi_cutoff/figures/umap_seurat_clusters.png", width = 10, height = 5, units = "in", res = 300)
DimPlot(seurat, reduction = "umap", group.by = "seurat_clusters", label = T, raster = F) + NoLegend()
dev.off()

png(filename = "whole_dataset/pmi_cutoff/figures/umap_condition.png", width = 10, height = 5, units = "in", res = 300)
DimPlot(seurat, reduction = "umap", group.by = "Condition", label = F, raster = F)
dev.off()

png(filename = "whole_dataset/pmi_cutoff/figures/umap_donor.png", width = 10, height = 5, units = "in", res = 300)
DimPlot(seurat, reduction = "umap", group.by = "Donor", label = F, raster = F)
dev.off()

png(filename = "whole_dataset/pmi_cutoff/figures/umap_cellytpe.png", width = 10, height = 5, units = "in", res = 300)
DimPlot(seurat, reduction = "umap", group.by = "CellType", label = T, raster = F) + NoLegend()
dev.off()

###
seurat <- FindSubCluster(seurat, cluster = "39", subcluster.name = "sub", graph.name = "SCT_snn", resolution = 0.5)

colnames(seurat@meta.data)

png(filename = "whole_dataset/pmi_cutoff/figures/umap_mtx_seurat_clusters.png", width = 10, height = 5, units = "in", res = 300)
DimPlot(seurat, reduction = "umap", group.by = "sub", label = T, raster = F) + NoLegend()
dev.off()

seurat

## cell annotation
seurat
seurat.marker <- FindAllMarkers(seurat, assay = "SCT", only.pos = T, verbose = T)
seurat.marker <- seurat.marker[seurat.marker$p_val_adj < 0.05,]
seurat.marker$pct_score <- seurat.marker$pct.1 - seurat.marker$pct.2

save.image('whole_dataset/pmi_cutoff/workspace_checkpoint_6.RData')



###
oligo.markers <- c("TF", "MOBP", "MBP", "MOG", "MAG")
astro.markers <- c("ALDH1L1","GFAP","SLC14A1","SLC1A3")
opc.marker <- c("LHFPL3","PDGFRA","VCAN")
inhibitory.markers <- c("GAD1", "GAD2")
microglia <- c("SYK", "RUNX1", "CD74", "SPP1")
endo.marker <- c("CLDN5")
fibro <- c("DCN", "COLEC12", "LAMA2", "COL1A2", "PDGFRA")
tcel.marker <- c("PTPRC", "CD247")


l2_3.markers <- c("CUX2")
l5_6.markers <- c("FOXP2", "B3GALT2", "NTNG2", "BCL11B")
l3_5.marker <- c("RORB")
l6.marker <- c("TLE4")

clusters <- c("Oligodendrocyte" = 0, "Astrocyte" = 1, "L2_3 Neuron" = 2, 
              "L2_3_Neuron" = 3, "OPC" = 4, "Inhibitory Neuron" = 5,
              "L5_6 Neuron" = 6, "L3_5 Neuron" = 7, "L4_6 Neuron" = 8,
              "L4_5 Neuron" = 9, "Microglia" = 10, "L4_6 Neuron" = 11,
              "L2_3 Neuron" = 12, "Inhibitory Neuron" = 13, "Inhibitory Neuron" = 14,
              "L4_5 Neuron" = 15, "Oligodendrocyte" = 16, "L2_3_Neuron" = 17,
              "L6 Neuron" = 18, "L6 Neuron" = 19, "Inhibitory Neuron" = 20,
              "Inhibitory Neuron" = 21, "Astrocyte" = 22, "L5_6 Neuron" = 23,
              "L5_6 Neuron" = 24, "Inhibitory Neuron" = 25, "Inhibitory Neuron" = 26,
              "L6 Neuron" = 27, "Inhibitory Neuron" = 28, "L5_6 Neuron" = 29, 
              "Inhibitory Neuron" = 30, "Endothelial cell" = 31, "Inhibitory Neuron" = 32,
              "Mural" = 33, "Inhibitory Neuron" = 34, "Fibroblast" = 35, 
              "Inhibitory Neuron" = 36, "L5_6 Neuron" = 37, "T cell" = 38)

seurat$annotated <- names(clusters[match(seurat$seurat_clusters, clusters)])
unique(seurat$annotated)

write_csv(seurat.marker, "whole_dataset/pmi_cutoff/seurat_marker_genes_total.csv")