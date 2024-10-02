###210902
###use scran to add cellcycle information in sccancer result
library(scran)
library(Seurat)
library(org.Hs.eg.db)
library(scatterplot3d)

setwd('/media/liyaru/LYR10T/Radiation_A549/TH1/A549')
###NC====
NC <- readRDS("scCANCER/res/NC/NC_seurat.rds")
NC <- DietSeurat(NC, graphs = "pca")
NC <- as.SingleCellExperiment(NC)

hm.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", 
                                package="scran"))
ensembl <- mapIds(org.Hs.eg.db, keys=rownames(NC), 
                  keytype="SYMBOL", column="ENSEMBL")
assigned <- cyclone(NC, pairs=hm.pairs, gene.names=ensembl)
head(assigned$scores)
table(assigned$phases)
saveRDS(assigned,'/media/liyaru/LYR10T/Radiation_A549/Revision/3_scran_list_NC.rds')

draw=cbind(assigned$score,assigned$phases) 
attach(draw)

cycle_color = c(rgb(247,229,163,maxColorValue = 255),
                rgb(190,186,218,maxColorValue = 255),
                rgb(141,211,199,maxColorValue = 255)
                )

# pdf("scCANCER/res/SCRAN/NC_cellcycle_2.pdf", height = 6, width = 8)
pdf("/media/liyaru/LYR10T/Radiation_A549/Revision/3_scran_3D_NC.pdf", height = 6, width = 8)
scatterplot3d(G1, S, G2M, angle=20,
              color = cycle_color[as.numeric(as.factor(assigned$phases))],
              grid=TRUE, box=FALSE)
dev.off()



save(assigned,file = 'cell_cycle_assigned.Rdata')
plot(assigned$score$G1, assigned$score$G2M, 
     xlab="G1 score", ylab="G2/M score", pch=16)
dev.off()
ggsave("scCANCER/res/SCRAN/NC_cellcycle_score.pdf")

NC$cellcycle <- assigned$phases
table(NC$cellcycle)
saveRDS(NC,"scCANCER/res/SCRAN/NC.rds")

NC <- readRDS("scCANCER/res/SCRAN/NC.rds")
NC$S.Score
NC$G2M.Score


###IR_2h====
IR_2h <- readRDS("scCANCER/res/IR_2h/IR_2h.rds")
IR_2h <- DietSeurat(IR_2h, graphs = "pca")
IR_2h <- as.SingleCellExperiment(IR_2h)

hm.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", 
                                package="scran"))
ensembl <- mapIds(org.Hs.eg.db, keys=rownames(IR_2h), 
                  keytype="SYMBOL", column="ENSEMBL")
assigned <- cyclone(IR_2h, pairs=hm.pairs, gene.names=ensembl)
head(assigned$scores)
table(assigned$phases)

saveRDS(assigned,'/media/liyaru/LYR10T/Radiation_A549/Revision/3_scran_list_2h.rds')

assigned <- readRDS('/media/liyaru/LYR10T/Radiation_A549/Revision/3_scran_list_2h.rds')

draw=cbind(assigned$score,assigned$phases) 
attach(draw) #attach的目的就是现在加载，之后直接引用即可

#pdf("scCANCER/res/SCRAN/IR_2h_cellcycle_2.pdf", height = 6, width = 8)
pdf("/media/liyaru/LYR10T/Radiation_A549/Revision/3_scran_3D_2h.pdf", height = 6, width = 8)
scatterplot3d(G1, S, G2M, angle=20,
              color = cycle_color[as.numeric(as.factor(assigned$phases))],
              grid=TRUE, box=FALSE)
dev.off()

save(assigned,file = 'cell_cycle_assigned.Rdata')
plot(assigned$score$G1, assigned$score$G2M, 
     xlab="G1 score", ylab="G2/M score", pch=16)
dev.off()
ggsave("scCANCER/res/SCRAN/IR_2h_cellcycle_score.pdf")

IR_2h$cellcycle <- assigned$phases
table(IR_2h$cellcycle)
saveRDS(IR_2h,"scCANCER/res/SCRAN/IR_2h.rds")

###IR_6h====
IR_6h <- readRDS("scCANCER/res/IR_6h/IR_6h.rds")
IR_6h <- DietSeurat(IR_6h, graphs = "pca")
IR_6h <- as.SingleCellExperiment(IR_6h)

hm.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", 
                                package="scran"))
ensembl <- mapIds(org.Hs.eg.db, keys=rownames(IR_6h), 
                  keytype="SYMBOL", column="ENSEMBL")
assigned <- cyclone(IR_6h, pairs=hm.pairs, gene.names=ensembl)
saveRDS(assigned,'/media/liyaru/LYR10T/Radiation_A549/Revision/3_scran_list_6h.rds')

assigned <- readRDS('/media/liyaru/LYR10T/Radiation_A549/Revision/3_scran_list_6h.rds')

head(assigned$scores)
table(assigned$phases)

draw=cbind(assigned$score,assigned$phases) 
attach(draw) #attach的目的就是现在加载，之后直接引用即可

pdf("/media/liyaru/LYR10T/Radiation_A549/Revision/3_scran_3D_6h.pdf", height = 6, width = 8)
#pdf("scCANCER/res/SCRAN/IR_6h_cellcycle_2.pdf", height = 6, width = 8)
scatterplot3d(G1, S, G2M, angle=20,
              color = cycle_color[as.numeric(as.factor(assigned$phases))],
              grid=TRUE, box=FALSE)
dev.off()

save(assigned,file = 'cell_cycle_assigned.Rdata')
plot(assigned$score$G1, assigned$score$G2M, 
     xlab="G1 score", ylab="G2/M score", pch=16)
dev.off()
ggsave("scCANCER/res/SCRAN/IR_6h_cellcycle_score.pdf")

IR_6h$cellcycle <- assigned$phases
table(IR_6h$cellcycle)
saveRDS(IR_6h,"scCANCER/res/SCRAN/IR_6h.rds")
table(IR_6h$cellcycle)
IR_2h$orig.ident
