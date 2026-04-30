options(repos = c(REPO_NAME = "https://packagemanager.posit.co/cran/__linux__/centos7/latest"))


setwd("/home/mpg02/MMNB/yigit.babal/workspace_ykb/fibroblasts/MIT_sc_dataset")
print(getwd())


source('/home/mpg08/mmnb1ngs/collab_work/workspace/ykbabal/library.R')
options(future.globals.maxSize = 3e+09)



sample.meta <- read_tsv("/home/mpg08/mmnb1ngs/collab_work/workspace/ykbabal/fibroblasts/MIT_sc_dataset/1_sample_metadata.tsv")
class(sample.meta)
sample.meta <- sample.meta %>% as.data.frame()
dim(sample.meta)
head(sample.meta)

unique(sample.meta$PMI_h)
sample.meta$PMI_h[sample.meta$PMI_h == "N/A"] <- NA
sample.meta$PMI_h <- as.numeric(sample.meta$PMI_h)
sample.meta <- sample.meta[sample.meta$Region != 'M1',]

colnames(sample.meta)
unique(sample.meta$Condition)
als.samples <- sample.meta[sample.meta$Group %in% c("SALS", "PN"), ]
als.samples <- na.omit(als.samples)

als.samples %>%
  dplyr::select(Condition, PMI_h) %>%
  table()

als.samples %>%
  dplyr::select(Condition, Sex) %>%
  table()


als.samples %>%
  dplyr::select(Condition, Donor) %>%
  group_by(Condition) %>%
  summarise(count = n_distinct(Donor))


als.samples <- als.samples[als.samples$PMI_h < 10,]

als.samples %>%
  dplyr::select(Condition, Donor) %>%
  group_by(Condition) %>%
  summarise(count = n_distinct(Donor))


raw.path <- "/home/mpg08/mmnb1ngs/collab_work/workspace/ykbabal/fibroblasts/MIT_sc_dataset/processed/"
count.list <- pblapply(als.samples$Sample_ID, function(x){
  count = ReadMtx(mtx = paste0(paste0(raw.path, x), '/counts_fil.mtx'),
                  cells = paste0(paste0(raw.path, x), '/col_metadata.tsv'),
                  features = paste0(paste0(raw.path, x), '/row_metadata.tsv'),
                  feature.column = 2, skip.cell = 1, skip.feature = 1)
  gene.anot <- read_tsv(paste0(paste0(raw.path, x), '/row_metadata.tsv'))
  gene.anot <- gene.anot[gene.anot$Biotype == "protein_coding",]
  gene.anot <- gene.anot[gene.anot$Chromosome != "MT",]
  gene.anot <- gene.anot[gene.anot$Gene != "N/A",]
  count <- count[rownames(count) %in% gene.anot$Gene,]
  return(count)
  
})

names(count.list) <- als.samples$Sample_ID



meta.list <- lapply(als.samples$Sample_ID, function(x){
  meta = read_tsv(file = paste0(paste0(raw.path, x), '/col_metadata.tsv'))
  meta = as.data.frame(meta)
  rownames(meta) = meta[,1]
  meta$PMI <- rep(sample.meta$PMI_h[sample.meta$Sample_ID == x], dim(meta)[1])
  return(meta)
})

names(meta.list) <- als.samples$Sample_ID

metadata <- Reduce(rbind, meta.list)
dim(metadata)


length(als.samples$Sample_ID)


seurat <- CreateSeuratObject(counts = count.list,
                             project = "als",
                             meta.data = metadata,
                             min.cells = 3, min.features = 200)

seurat


rm(count.list, meta.list)
gc()


write_rds(seurat, file = "whole_dataset/pmi_cutoff/seurat_als_mtx_v5_layered.rds")


seurat[["percent_ribo"]] <- PercentageFeatureSet(seurat, pattern = "^RP[SL]") # ribosomal genes
seurat[["percent_hb"]] <- PercentageFeatureSet(seurat, pattern = "^HB[^(P)]") # hemoglobin genes

png(filename = "whole_dataset/pmi_cutoff/figures/QC_before.png", width = 10, height = 40, units = "in", res = 300)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent_ribo", "percent_hb"),
        group.by = "Donor", ncol = 1, raster = F)
dev.off()

seurat <- subset(seurat, subset = nCount_RNA < 50000 & percent_ribo < 2 & percent_hb < 0.2)
seurat

png(filename = "whole_dataset/pmi_cutoff/figures/QC_after.png", width = 10, height = 40, units = "in", res = 300)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent_ribo", "percent_hb"),
        group.by = "Donor", ncol = 1, raster = F)
dev.off()


options(future.globals.maxSize = 3e+09)


save.image("whole_dataset/pmi_cutoff/workspace_checkpoint.RData")

seurat <- SCTransform(seurat, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "PMI"),
                      vst.flavor = "v2", method = "glmGamPoi", verbose = T)


save.image("whole_dataset/pmi_cutoff/workspace_checkpoint_2.RData")
