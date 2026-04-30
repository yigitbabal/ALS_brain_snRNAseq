options(repos = c(REPO_NAME = "https://packagemanager.posit.co/cran/__linux__/centos7/latest"))


setwd("/home/mpg02/MMNB/yigit.babal/workspace_ykb/fibroblasts/MIT_sc_dataset")
print(getwd())


source('/home/mpg08/mmnb1ngs/collab_work/workspace/ykbabal/library.R')
options(future.globals.maxSize = 3e+09)

load("whole_dataset/pmi_cutoff/workspace_checkpoint_2.RData")

var.features <- VariableFeatures(seurat)


var.features <- var.features[!grepl(pattern = "^RP[SL]", x = var.features)]
var.features <- var.features[!grepl(pattern = "^HSP", x = var.features)]
var.features <- var.features[!grepl(pattern = "^LINC0", x = var.features)]
var.features <- var.features[!grepl(pattern = "^N/A", x = var.features)]


VariableFeatures(seurat) <- var.features
seurat


seurat <- RunPCA(seurat, features = VariableFeatures(seurat), npcs = 100)

png(filename = "whole_dataset/pmi_cutoff/figures/elbow_pca.png", width = 20, height = 15, units = "in", res = 300)
ElbowPlot(seurat, ndims = 100)
dev.off()

seurat$Donor <- as.factor(seurat$Donor)
Idents(seurat) <- seurat$Donor

seurat <- IntegrateLayers(
  object = seurat,
  method = RPCAIntegration,
  normalization.method = "SCT",
  verbose = T
)

save.image("whole_dataset/pmi_cutoff/workspace_checkpoint_3.RData")

seurat

seurat <- seurat %>% 
  RunUMAP(reduction = "integrated.dr", dims = 1:90) %>% 
  FindNeighbors(reduction = "integrated.dr", dims = 1:90) %>% 
  FindClusters(resolution = 0.9) %>% 
  identity()

save.image("whole_dataset/pmi_cutoff/workspace_checkpoint_4.RData")


seurat <- PrepSCTFindMarkers(seurat)

save.image("whole_dataset/pmi_cutoff/workspace_checkpoint_5.RData")


