options(repos = c(REPO_NAME = "https://packagemanager.posit.co/cran/__linux__/centos7/latest"))


setwd("/home/mpg02/MMNB/yigit.babal/workspace_ykb/fibroblasts/MIT_sc_dataset")
print(getwd())


source('/home/mpg08/mmnb1ngs/collab_work/workspace/ykbabal/library.R')

seurat.list <- SplitObject(seurat, split.by = "Condition")
seurat.list
seurat.als <- seurat.list$ALS
seurat.als
seurat.pn <- seurat.list$PN
seurat.pn
rm(seurat.list)
gc()

cellchat.als <- createCellChat(object = seurat.als, group.by = "annotated")
cellchat.pn <- createCellChat(object = seurat.pn, group.by = "annotated")

CellChatDB <- CellChatDB.human

png(filename = "whole_dataset/pmi_cutoff/figures/cellchat_categories.png", width = 10, height = 5, units = "in", res = 300)
showDatabaseCategory(CellChatDB)
dev.off()

CellChatDB.use <- CellChatDB

cellchat.als@DB <- CellChatDB
cellchat.pn@DB <- CellChatDB

cellchat.als <- subsetData(cellchat.als)
cellchat.pn <- subsetData(cellchat.pn)

future::plan("multisession", workers = 10) # do parallel

cellchat.als <- identifyOverExpressedGenes(cellchat.als)
cellchat.als <- identifyOverExpressedInteractions(cellchat.als)
cellchat.als

cellchat.pn <- identifyOverExpressedGenes(cellchat.pn)
cellchat.pn <- identifyOverExpressedInteractions(cellchat.pn)

cellchat.als <- computeCommunProb(cellchat.als, population.size = TRUE)
cellchat.pn <- computeCommunProb(cellchat.pn, population.size = TRUE)

cellchat.als <- filterCommunication(cellchat.als, min.cells = 10)
cellchat.pn <- filterCommunication(cellchat.pn, min.cells = 10)

cellchat.als <- computeCommunProbPathway(cellchat.als)
cellchat.pn <- computeCommunProbPathway(cellchat.pn)

cellchat.als <- aggregateNet(cellchat.als)
cellchat.pn <- aggregateNet(cellchat.pn)

cellchat.als <- netAnalysis_computeCentrality(cellchat.als)
cellchat.pn <- netAnalysis_computeCentrality(cellchat.pn)


object.list <- list(PN = cellchat.pn, ALS = cellchat.als)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

levels(cellchat.pn@idents)
levels(cellchat.als@idents)

seurat$cell_types <- seurat$refined
DimPlot(seurat, group.by = "cell_types", label = T, repel = T) +
  DimPlot(seurat, group.by = "Condition")

gg1.1 <- netVisual_heatmap(cellchat, font.size = 15)
gg2.1 <- netVisual_heatmap(cellchat, measure = "weight", font.size = 15)

gg1.1 + gg2.1


ht1.2 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = "SPP1", title = names(object.list)[i], width = 10, height = 5, color.heatmap = "GnBu", font.size = 15)
ht2.2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = "SPP1", title = names(object.list)[i+1], width = 10, height = 5, color.heatmap = "GnBu", font.size = 15)

pdf("cellchat_pathway_heatmap_incoming.pdf", width = 10, height = 5)
draw(ht1.2 + ht2.2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1.2 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = "SPP1", title = names(object.list)[i], width = 10, height = 5, color.heatmap = "OrRd", font.size = 15)
ht2.2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = "SPP1", title = names(object.list)[i+1], width = 10, height = 5, color.heatmap = "OrRd", font.size = 15)

pdf("cellchat_pathway_heatmap_outgoing.pdf", width = 10, height = 5)
draw(ht1.2 + ht2.2, ht_gap = unit(0.5, "cm"))
dev.off()


ht1.2 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = "SPP1", title = names(object.list)[i], width = 10, height = 5, color.heatmap = "OrRd", font.size = 15)
ht2.2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = "SPP1", title = names(object.list)[i+1], width = 10, height = 5, color.heatmap = "OrRd", font.size = 15)

pdf("cellchat_pathway_heatmap_overall.pdf", width = 10, height = 5)
draw(ht1.2 + ht2.2, ht_gap = unit(0.5, "cm"))
dev.off()


gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Microglia", comparison = c(2,1))
gg1 + xlim(-0.25, 0.25) + ylim(-0.25, 0.25)


pathways.show <- c("SPP1") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


netVisual_bubble(cellchat, signaling = "SPP1", remove.isolate = T, comparison = c(2, 1), angle.x = 45)
