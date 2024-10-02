####21.09.13
####cellchat.NC comparison
library(cellchat.NC)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(mindr)
library(patchwork)

###Data process====
##create cellchat object from sccancer output
NC <- readRDS("scCANCER/res/SCRAN/NC.rds")
NC <- as.Seurat(NC, slot = "RNA")
data.input <- NC@assays$originalexp@data
#Idents(NC) <- NC$cellcycle
table(NC$old.ident)

labels <- Idents(NC)
meta <- data.frame(group = NC@meta.data$cellcycle, row.names = names(NC@meta.data$cellcycle))
unique(meta$group)
cellchat.NC <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat.NC

groupSize <- as.numeric(table(cellchat.NC@idents))
groupSize
cellchat <- cellchat.NC

##Load ligand receptor database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object

###PREPROCESSION====
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 16)

###find high expres_scransion ligand/receptor
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) #结果存在cellchat@LR$LRsig
cellchat <- projectData(cellchat, PPI.human)#将找到的受配体关系映射到PPI上 对@data.signaling中的表达进行校正

###infer interaction====
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)#如果不想用PPI的校正结果，raw.use=TURE
#filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net,"cellchat/res_scran/NC/net_lr.csv")

###infer 信号通路水平的细胞通讯网络（结果储存在@netP下面，有一个概率值和对应的p val）
cellchat <- computeCommunProbPathway(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "cellchat/res_scran/NC/netpathway.csv")

#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

###通讯网络系统分析====
##社会网络分析（通讯网络中的角色识别）
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
saveRDS(cellchat,"cellchat/res_scran/NC/cellchat_NC.rds")


###IR_2h===============
###Data process====
##create cellchat object from sccaIR_2her output
IR_2h <- readRDS("scCANCER/res/SCRAN/IR_2h.rds")
IR_2h <- as.Seurat(IR_2h, slot = "RNA")
data.input <- IR_2h@assays$originalexp@data
Idents(IR_2h) <- IR_2h$cellcycle
labels <- Idents(IR_2h)
meta <- data.frame(group = IR_2h@meta.data$cellcycle, row.names = names(IR_2h@meta.data$cellcycle))
unique(meta$group)
cellchat.IR_2h <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat.IR_2h
str(cellchat.IR_2h)
levels(cellchat.IR_2h)
groupSize <- as.numeric(table(cellchat.IR_2h@idents))
groupSize
cellchat <- cellchat.IR_2h

##Load ligand receptor database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object

###PREPROCESSION====
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 16)

###find high expres_scransion ligand/receptor
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) #结果存在cellchat@LR$LRsig
cellchat <- projectData(cellchat, PPI.human)#将找到的受配体关系映射到PPI上 对@data.signaling中的表达进行校正

###infer interaction====
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)#如果不想用PPI的校正结果，raw.use=TURE
#filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net,"cellchat/res_scran/IR_2h/net_lr.csv")

###infer 信号通路水平的细胞通讯网络（结果储存在@netP下面，有一个概率值和对应的p val）
cellchat <- computeCommunProbPathway(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "cellchat/res_scran/IR_2h/netpathway.csv")

#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

###通讯网络系统分析====
##社会网络分析（通讯网络中的角色识别）
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
saveRDS(cellchat,"cellchat/res_scran/IR_2h/cellchat_IR_2h.rds")

###IR_6h===============
###Data process====
##create cellchat object from sccaIR_6her output
IR_6h <- readRDS("scCANCER/res/SCRAN/IR_6h.rds")
IR_6h <- as.Seurat(IR_6h, slot = "RNA")
data.input <- IR_6h@assays$originalexp@data
Idents(IR_6h) <- IR_6h$cellcycle
labels <- Idents(IR_6h)
meta <- data.frame(group = IR_6h@meta.data$cellcycle, row.names = names(IR_6h@meta.data$cellcycle))
unique(meta$group)
cellchat.IR_6h <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat.IR_6h
str(cellchat.IR_6h)
levels(cellchat.IR_6h)
groupSize <- as.numeric(table(cellchat.IR_6h@idents))
groupSize
cellchat <- cellchat.IR_6h

##Load ligand receptor database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object

###PREPROCESSION====
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 16)

###find high expres_scransion ligand/receptor
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) #结果存在cellchat@LR$LRsig
cellchat <- projectData(cellchat, PPI.human)#将找到的受配体关系映射到PPI上 对@data.signaling中的表达进行校正

###infer interaction====
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)#如果不想用PPI的校正结果，raw.use=TURE
#filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net,"cellchat/res_scran/IR_6h/net_lr.csv")

###infer 信号通路水平的细胞通讯网络（结果储存在@netP下面，有一个概率值和对应的p val）
cellchat <- computeCommunProbPathway(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "cellchat/res_scran/IR_6h/netpathway.csv")

#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

###通讯网络系统分析====
##社会网络分析（通讯网络中的角色识别）
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
saveRDS(cellchat,"cellchat/res_scran/IR_6h/cellchat_IR_6h.rds")

###IR_6h===============
###Data process====
##create cellchat object from sccaIR_6her output
IR_6h <- readRDS("scCANCER/res/SCRAN/IR_6h.rds")
IR_6h <- as.Seurat(IR_6h, slot = "RNA")
data.input <- IR_6h@assays$originalexp@data
Idents(IR_6h) <- IR_6h$cellcycle
labels <- Idents(IR_6h)
meta <- data.frame(group = IR_6h@meta.data$cellcycle, row.names = names(IR_6h@meta.data$cellcycle))
unique(meta$group)
cellchat.IR_6h <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat.IR_6h
str(cellchat.IR_6h)
levels(cellchat.IR_6h)
groupSize <- as.numeric(table(cellchat.IR_6h@idents))
groupSize
cellchat <- cellchat.IR_6h

##Load ligand receptor database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object

###PREPROCESSION====
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 16)

###find high expres_scransion ligand/receptor
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) #结果存在cellchat@LR$LRsig
cellchat <- projectData(cellchat, PPI.human)#将找到的受配体关系映射到PPI上 对@data.signaling中的表达进行校正

###infer interaction====
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)#如果不想用PPI的校正结果，raw.use=TURE
#filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net,"cellchat/res_scran/IR_6h/net_lr.csv")

###infer 信号通路水平的细胞通讯网络（结果储存在@netP下面，有一个概率值和对应的p val）
cellchat <- computeCommunProbPathway(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "cellchat/res_scran/IR_6h/netpathway.csv")

#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

###通讯网络系统分析====
##社会网络分析（通讯网络中的角色识别）
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
saveRDS(cellchat,"cellchat/res_scran/IR_6h/cellchat_IR_6h.rds")






###Data load====
cellchat.NC <- readRDS("cellchat/res_scran/NC/cellchat_NC.rds")
cellchat.IR_6h <- readRDS("cellchat/res_scran/IR_6h/cellchat_IR_6h.rds")
object.list <- list(NC = cellchat.NC, IR_6h = cellchat.IR_6h)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
cellchat

###Compare ====
##the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
ggsave("cellchat/res_scran/compare/IR6h_NC/the total number of interactions and interaction strength.pdf", width = 6, height = 6, dpi = 600)
ggsave("cellchat/res_scran/compare/IR6h_NC/the total number of interactions and interaction strength.png", width = 6, height = 6, dpi = 600)

##Compare the number of interactions and interaction strength among different cell populations
#Differential number of interactions or interaction strength among different cell populations
pdf("cellchat/res_scran/compare/IR6h_NC/Differential number of interactions_1.pdf", width = 8, height = 6)
par(mfrow = c(1,2), xpd=TRUE)
gg1 <- netVisual_diffInteraction(cellchat, weight.scale = T)
gg2 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
gg1 + gg2
dev.off()

#heatmap
pdf("cellchat/res_scran/compare/IR6h_NC/Differential number of interactions_2.pdf", width = 8, height = 6)
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()
ggsave("cellchat/res_scran/compare/IR6h_NC/Differential number of interactions_2.pdf", width = 8, height = 6, dpi = 600)
ggsave("cellchat/res_scran/compare/IR6h_NC/Differential number of interactions_2.png", width = 8, height = 6, dpi = 600)

##Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
dev.new()
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
dev.off()
ggsave("cellchat/res_scran/compare/IR6h_NC/the major sources and targets.pdf", width = 8, height = 5, dpi = 600)
ggsave("cellchat/res_scran/compare/IR6h_NC/the major sources and targets.png", width = 8, height = 5, dpi = 600)

## Identify signaling changes associated with one cell group
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "G1")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "S")
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "G2M")
patchwork::wrap_plots(plots = list(gg1,gg2,gg3))
#save as cellchat/res_scran/compare/IR6h_NC/signaling chages.pdf

write.csv(gg1$data,"cellchat/res_scran/compare/IR6h_NC/signaling chages_G1.csv")
write.csv(gg2$data,"cellchat/res_scran/compare/IR6h_NC/signaling chages_S.csv")
write.csv(gg3$data,"cellchat/res_scran/compare/IR6h_NC/signaling chages_G2M.csv")

###Identify the conserved and context-specific signaling pathways
##Identify signaling networks with larger (or less) difference as well as signaling groups based on their functional/structure similarity
#Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
ggsave("cellchat/res_scran/compare/IR6h_NC/signaling groups based on their functional similarity.pdf")
ggsave("cellchat/res_scran/compare/IR6h_NC/signaling groups based on their functional similaritys.png")

#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
ggsave("cellchat/res_scran/compare/IR6h_NC/signaling groups based on their structure similarity.pdf")
ggsave("cellchat/res_scran/compare/IR6h_NC/signaling groups based on their structure similarity.png")
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
ggsave("cellchat/res_scran/compare/IR6h_NC/signaling groups based on their structure similarity_1.pdf")
ggsave("cellchat/res_scran/compare/IR6h_NC/signaling groups based on their structure similarity_1.png")

##Compute and visualize the pathway distance in the learned joint manifold
rankSimilarity(cellchat, type = "functional")
ggsave("cellchat/res_scran/compare/IR6h_NC/pathway distance.pdf", width = 3, height = 4, dpi = 600)
ggsave("cellchat/res_scran/compare/IR6h_NC/pathway distance.png", width = 3, height = 4, dpi = 600)

#Identify and visualize the conserved and context-specific signaling pathways
#Compare the overall information flow of each signaling pathway
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
ggsave("cellchat/res_scran/compare/IR6h_NC/conserved and context-specific signaling pathways.pdf")
ggsave("cellchat/res_scran/compare/IR6h_NC/conserved and context-specific signaling pathways.png")

#Compare outgoing (or incoming) signaling associated with each cell population
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pdf("cellchat/res_scran/compare/IR6h_NC/outgoing signaling pathways.pdf")
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

pdf("cellchat/res_scran/compare/IR6h_NC/incoming signaling pathways.pdf")
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

pdf("cellchat/res_scran/compare/IR6h_NC/overall signaling pathways.pdf")
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

###Identify the upgulated and down-regulated signaling ligand-receptor pairs====
##Identify dysfunctional signaling by comparing the communication probabities
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:3), comparison = c(1, 2), angle.x = 45)
ggsave("cellchat/res_scran/compare/IR6h_NC/ligand_receptor pairs_G1.pdf")
ggsave("cellchat/res_scran/compare/IR6h_NC/ligand_receptor pairs_G1.png")

netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:3), comparison = c(1, 2), angle.x = 45)
ggsave("cellchat/res_scran/compare/IR6h_NC/ligand_receptor pairs_G2M.pdf")
ggsave("cellchat/res_scran/compare/IR6h_NC/ligand_receptor pairs_G2M.png")

netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:3), comparison = c(1, 2), angle.x = 45)
ggsave("cellchat/res_scran/compare/IR6h_NC/ligand_receptor pairs_S.pdf")
ggsave("cellchat/res_scran/compare/IR6h_NC/ligand_receptor pairs_S.png")

# identify the upgulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset.
gg1 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:3),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in IR_6h", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:3),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in IR_6h", angle.x = 45, remove.isolate = T)
gg1 + gg2
ggsave("cellchat/res_scran/compare/IR6h_NC/ligand_receptor pairs_change_G1.pdf")
ggsave("cellchat/res_scran/compare/IR6h_NC/ligand_receptor pairs_change_G1.png")

gg1 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:3),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in IR_6h", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:3),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in IR_6h", angle.x = 45, remove.isolate = T)
gg1 + gg2
ggsave("cellchat/res_scran/compare/IR6h_NC/ligand_receptor pairs_change_G2M.pdf", height = 6, width = 8)
ggsave("cellchat/res_scran/compare/IR6h_NC/ligand_receptor pairs_change_G2M.png")

gg1 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:3),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in IR_6h", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:3),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in IR_6h", angle.x = 45, remove.isolate = T)
gg1 + gg2
ggsave("cellchat/res_scran/compare/IR6h_NC/ligand_receptor pairs_change_S.pdf", height = 5, width = 6)
ggsave("cellchat/res_scran/compare/IR6h_NC/ligand_receptor pairs_change_S.png")
#Identify dysfunctional signaling by using differential expres_scransion analysis
#(https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html#load-the-required-libraries)

#Visualize the upgulated and down-regulated signaling ligand-receptor pairs using Chord diagram
# Chord diagram
# par(mfrow = c(1,2), xpd=TRUE)
# netVisual_chord_gene(object.list[[2]], sources.use = 2, targets.use = c(1:3), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
# netVisual_chord_gene(object.list[[1]], sources.use = 2, targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))


###Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram====
#(https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html#load-the-required-libraries)


###Compare the signaling gene expres_scransion distribution between different datasets====
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NC", "IR_6h")) # set factor level
plotGeneExpression(cellchat, signaling = "GDF", split.by = "datasets", colors.ggplot = T)
ggsave("cellchat/res_scran/compare/IR6h_NC/signaling gene expres_scransion distribution between different datasets.pdf")
ggsave("cellchat/res_scran/compare/IR6h_NC/signaling gene expres_scransion distribution between different datasets.png")

cellchat@netP

saveRDS(cellchat, file = "cellchat/res_scran/compare/IR6h_NC/cellchat_comparisonAnalysis_IR_6h_vs_NC.rds")


####---------- by cluster ------------
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(mindr)
library(patchwork)

setwd('/media/liyaru/LYR10T/Radiation_A549/TH1/A549')

###Data process====
##create cellchat object from sccancer output
# NC <- readRDS("scCANCER/res/SCRAN/NC.rds")
NC <- readRDS("scCANCER/res/NC/NC_seurat.rds")
DimPlot(NC,group.by = 'RNA_snn_res.0.8',reduction = 'tsne')
#NC <- as.Seurat(NC, slot = "RNA")
data.input <- NC@assays$RNA@counts
meta <- data.frame(group_orig = NC@meta.data$RNA_snn_res.0.8,
                   group_num = as.numeric(NC@meta.data$RNA_snn_res.0.8), 
                   row.names = rownames(NC@meta.data))

cellchat.NC <- createCellChat(object = data.input, meta = meta, group.by = "group_num")
cellchat.NC

groupSize <- as.numeric(table(cellchat.NC@idents))
groupSize
cellchat <- cellchat.NC

##Load ligand receptor database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object

###PREPROCESSION====
cellchat <- subsetData(cellchat)
# future::plan("multiprocess", workers = 16)

###find high expres_scransion ligand/receptor
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) #结果存在cellchat@LR$LRsig
cellchat <- projectData(cellchat, PPI.human)#将找到的受配体关系映射到PPI上 对@data.signaling中的表达进行校正

###infer interaction====
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)#如果不想用PPI的校正结果，raw.use=TURE
#filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)

# write.csv(df.net,"cellchat/res_scran/NC/net_lr.csv")

###infer 信号通路水平的细胞通讯网络（结果储存在@netP下面，有一个概率值和对应的p val）
cellchat <- computeCommunProbPathway(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
# write.csv(df.netP, "cellchat/res_scran/NC/netpathway.csv")

#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

###通讯网络系统分析====
##社会网络分析（通讯网络中的角色识别）
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
# saveRDS(cellchat,"cellchat/res_scran/NC/cellchat_NC.rds")

groupSize <- as.numeric(table(cellchat@idents))

pdf('/media/liyaru/LYR10T/Radiation_A549/Revision/6_cellchat.pdf')
# par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

