library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(ggplot2)

###lOAD DATA====
# define a convenient function to load all the data and create a Seurat object
create_obj <- function(dir) {
  count.path <- list.files(path = dir, pattern = "*filtered_peak_bc_matrix.h5", full.names = TRUE)
  fragment.path <- list.files(path = dir, pattern = "*fragments.tsv.gz", full.names = TRUE)[1]
  counts <- Read10X_h5(count.path)
  md.path <- list.files(path = dir, pattern = "*singlecell.csv", full.names = TRUE)
  md <- read.table(file = md.path, stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)
  chrom_assay <- CreateChromatinAssay(counts = counts, 
                                      sep = c(":", "-"), 
                                      genome = 'GRCh38',
                                      fragments = fragment.path
                                      )
  obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = 'peaks',
    meta.data = md
  )
  return(obj)
}

N1 <- create_obj("scATAC_cellranger/A549_N1/outs/")
C2 <- create_obj("scATAC_cellranger/A549_C2/outs")
C3 <- create_obj("scATAC_cellranger/A549_C3/outs")

#Creating a common peak set
# 直接使用UnifyPeaks函数将多个不同对象的中peaks进行合并
combined.peaks <- UnifyPeaks(object.list = list(N1, C2, C3), mode = "reduce")

#Quantify peaks in each dataset
# 使用FeatureMatrix函数对每个数据集重新进行计数定量
N1.counts <- FeatureMatrix(
  fragments = CreateFragmentObject('scATAC_file/N1/N1_fragments.tsv.gz', cells = colnames(N1)),
  features = combined.peaks,
  sep = c(":", "-"),
)

C2.counts <- FeatureMatrix(
  fragments = CreateFragmentObject('scATAC_file/C2/C2_fragments.tsv.gz', cells = colnames(C2)),
  features = combined.peaks,
  sep = c(":", "-")
)

C3.counts <- FeatureMatrix(
  fragments = CreateFragmentObject('scATAC_file/C3/C3_fragments.tsv.gz', cells = colnames(C3)),
  features = combined.peaks,
  sep = c(":", "-")
)


# 使用CreateAssayObject函数新建一个assay对象存储计数的数据
N1[['peaks']] <- CreateAssayObject(counts = N1.counts)
C2[['peaks']] <- CreateAssayObject(counts = C2.counts)
C3[['peaks']] <- CreateAssayObject(counts = C3.counts)

###Merge objects====
# add information to identify dataset of origin
N1$sample <- 'NC'
C2$sample <- 'IR_immediate'
C3$sample <- 'IR_2h'

# merge all datasets, adding a cell ID to make sure cell names are unique
# 直接使用merge函数将多个数据集进行合并
combined <- merge(x = N1, y = list(C2, C3), add.cell.ids = c("NC", "IR_immediate", "IR_2h"))
table(combined$sample)
#IR_2h    IR_immediate      NC 
#4073         8596         8927 
saveRDS(combined,"Signac/integration/combined.rds")
# combined <- readRDS("Signac/integration/combined.rds")

# make sure to change to the assay containing common peaks
DefaultAssay(combined) <- "peaks"

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(
  combined,
  reduction.key = 'LSI_',
  reduction.name = 'lsi',
  irlba.work = 400
)

# RunChromVAR
combined <- RunChromVAR(
  object = combined,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(combined) <- 'chromvar'
combined <- computeDeviations(combined)
computeDeviations(combined,
                  annotations, background_peaks = getBackgroundPeaks(combined),
                  expectation = computeExpectations(combined))


saveRDS(combined, "Signac/integration/combined_chrom_VAR.rds")
# combined <- readRDS("Signac/integration/combined_chrom_VAR.rds")

Idents(combined) <- combined$orig.ident

#motif
differential.activity_N1_C3 <- FindMarkers(
  object = combined,
  ident.1 = 'NC',
  ident.2 = 'IR_2h',
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks')

MotifPlot(
  object = combined,
  motifs = rownames(differential.activity_N1_C3)[1:25],
  assay = 'peaks'
)
rownames(differential.activity_C3_N1)
ggsave("Signac/integration/TF/differential.motif_N1_C3_1.pdf")


write.csv(differential.activity_N1_C3,"Signac/integration/TF/differential.activity_N1_C3.csv")
differential.activity_N1_C3 <- read.csv("Signac/integration/TF/differential.activity_N1_C3.csv")
motif.name = as.data.frame(combined@assays$peaks@motifs@motif.names)
motif.name <- t(motif.name)
write.csv(motif.name, "Signac/integration/TF/motif_name.csv")
motif.name <- read.csv("Signac/integration/TF/motif_name.csv")
differential.activity_N1_C3 <- merge(differential.activity_N1_C3,motif.name, all = FALSE)
write.csv(differential.activity_N1_C3,"Signac/integration/TF/differential.activity_N1_C3.csv")

differential.activity_C3_N1 <- FindMarkers(
  object = combined,
  ident.1 = 'IR_2h',
  ident.2 = 'NC',
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
write.csv(differential.activity_C3_N1,"Signac/integration/TF/differential.activity_C3_N1.csv")
differential.activity_C3_N1 <- read.csv("Signac/integration/TF/differential.activity_C3_N1.csv")
differential.activity_C3_N1 <- merge(differential.activity_C3_N1,motif.name, by ='X')
write.csv(differential.activity_C3_N1,"Signac/integration/TF/differential.activity_C3_N1.csv")

MotifPlot(
  object = combined,
  motifs = head(rownames(differential.activity_C3_N1)),
  assay = 'peaks'
)
rownames(differential.activity_C3_N1)
ggsave("Signac/integration/TF/differential.motif_C3_N1.pdf")

###C2 VS N1====
differential.activity_N1_C2 <- FindMarkers(
  object = combined,
  ident.1 = 'NC',
  ident.2 = 'IR_immediate',
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)

MotifPlot(
  object = combined,
  motifs = rownames(differential.activity_N1_C2),
  assay = 'peaks'
)
ggsave("Signac/integration/TF/differential.motif_N1_C2.pdf")

write.csv(differential.activity_N1_C2,"Signac/integration/TF/differential.activity_N1_C2.csv")

differential.activity_N1_C2 <- read.csv("Signac/integration/TF/differential.activity_N1_C2.csv")
motif.name = as.data.frame(combined@assays$peaks@motifs@motif.names)
motif.name <- t(motif.name)
write.csv(motif.name, "Signac/integration/TF/motif_name.csv")
motif.name <- read.csv("Signac/integration/TF/motif_name.csv")
differential.activity_N1_C2 <- merge(differential.activity_N1_C2,motif.name, all = FALSE)
write.csv(differential.activity_N1_C2,"Signac/integration/TF/differential.activity_N1_C2.csv")

a <- read.csv("Signac/integration/TF/differential.activity_N1_C2.csv")

differential.activity_C2_N1 <- FindMarkers(
  object = combined,
  ident.1 = 'IR_immediate',
  ident.2 = 'NC',
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
write.csv(differential.activity_C2_N1,"Signac/integration/TF/differential.activity_C2_N1.csv")

MotifPlot(
  object = combined,
  motifs = rownames(differential.activity_C2_N1),
  assay = 'peaks'
)
ggsave("Signac/integration/TF/differential.motif_C2_N1.pdf")

differential.activity_C2_N1 <- read.csv("Signac/integration/TF/differential.activity_C2_N1.csv")
differential.activity_C2_N1 <- merge(differential.activity_C2_N1,motif.name, all = FALSE)
write.csv(differential.activity_C2_N1,"Signac/integration/TF/differential.activity_C2_N1.csv")


library(ggplot2)
library(ggrepel)
#load data
Dat1 <- read.csv("Signac/integration/TF/differential.activity_C3_N1.csv")
Dat2 <- read.csv("Signac/integration/TF/differential.activity_N1_C3.csv")
Dat2 <- Dat2[,-1]
Dat2$avg_log2FC <- -Dat2$avg_log2FC
Dat <- rbind(Dat1,Dat2)

Dat$threshold = factor(ifelse(Dat$p_val_adj < 0.05 & abs(Dat$avg_log2FC) >= 0.2, ifelse(Dat$avg_log2FC>= 0.2,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))

ggplot(Dat,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
  geom_text_repel(
    data = Dat[which(Dat$p_val_adj<0.05 & abs(Dat$avg_log2FC)>0.2),],
    aes(label = V1),
    size = 3,
    segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
  theme_bw()+
  theme(
    legend.title = element_blank()#不显示图例标题
  )+
  ylab('-log10 (p-adj)')+
  xlab('log2 (FoldChange)')+
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05
ggsave('Signac/integration/TF/IR2h_NC_colcano plot.pdf', width = 6, height = 4)

#Load data
Dat1 <- read.csv("Signac/integration/TF/differential.activity_C2_N1.csv")
Dat2 <- read.csv("Signac/integration/TF/differential.activity_N1_C2.csv")
Dat2$avg_log2FC <- -Dat2$avg_log2FC
Dat <- rbind(Dat1,Dat2)

Dat$threshold = factor(ifelse(Dat$p_val_adj < 0.05 & abs(Dat$avg_log2FC) >= 0.25, ifelse(Dat$avg_log2FC>= 0.25,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))

ggplot(Dat,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+
  geom_text_repel(
    data = Dat[which(Dat$p_val_adj<0.05 & abs(Dat$avg_log2FC)>0.25),],
    aes(label = V1),
    size = 3,
    segment.color = "black", show.legend = FALSE )+
  theme_bw()+
  theme(
    legend.title = element_blank()
  )+
  ylab('-log10 (p-adj)')+
  xlab('log2 (FoldChange)')+
  geom_vline(xintercept=c(-0.25,0.25),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)
ggsave('Signac/integration/TF/IRimmediate_NC_colcano plot_0.25.pdf', width = 6, height = 4)


####--------- cistopic-----------------------
library(cisTopic)
##Load data====
pathTo10X <- 'E:/A549/A549-ATAC-0208/outs/'
data_folder <- paste0(pathTo10X, 'filtered_peak_bc_matrix')
metrics <- paste0(pathTo10X, 'singlecell.csv')
cisTopicObject <- createcisTopicObjectFrom10Xmatrix(data_folder, metrics,  project.name='IR_2h')
cisTopicObject
##Building the models====
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(2, 5, 10:25, 30, 35, 40), seed=987, nCores=20, iterations = 500, addModels=FALSE)
saveRDS(cisTopicObject,"cisTOPIC/res/int/cisTopicObject_1.rds")
cisTopicObject<- readRDS("cisTOPIC/res/int/cisTopicObject_1.rds")
#add sample labels to cell data
A <- cisTopicObject@cell.data
A$sample <- gsub("^.*-", "", rownames(A))
A$sample <- gsub("1","NC", A$sample)
A$sample <- gsub("2","IR_immediate", A$sample)
A$sample <- gsub("3","IR_2h", A$sample)
table(A$sample)
cisTopicObject@cell.data <- A
cell.names <- cisTopicObject@cell.names


# # For WarpLDA
# cisTopicObject@calc.params[['runWarpLDAModels']]$seed <- seed
# cisTopicObject@calc.params[['runWarpLDAModels']]$iterations <- iterations
# cisTopicObject@calc.params[['runWarpLDAModels']]$alpha <- alpha
# cisTopicObject@calc.params[['runWarpLDAModels']]$alphaByTopic <- alphaByTopic
# cisTopicObject@calc.params[['runWarpLDAModels']]$beta <- beta
# 
# # For CGS
# cisTopicObject@calc.params[['runCGSModels']]$seed <- seed
# cisTopicObject@calc.params[['runCGSModels']]$burnin <- burnin
# cisTopicObject@calc.params[['runCGSModels']]$iterations <- iterations
# cisTopicObject@calc.params[['runCGSModels']]$alpha <- alpha
# cisTopicObject@calc.params[['runCGSModels']]$alphaByTopic <- alphaByTopic
# cisTopicObject@calc.params[['runCGSModels']]$beta <- beta

#Selection of the best model
pdf("cisTOPIC/res/int/model_topic_number.pdf", height = 9, width = 9)
par(mfrow=c(3,3))
cisTopicObject <- selectModel(cisTopicObject, type='maximum')
cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
cisTopicObject <- selectModel(cisTopicObject, type='derivative')
dev.off()

##Interpreting the models====
#Identification of cell states using the cell-cisTopic distributions
cisTopicObject <- runtSNE(cisTopicObject, target='cell', seed=123, pca=FALSE, method='Probability')
cisTopicObject <- runUmap(cisTopicObject, target='cell', seed=123, method='Probability')
cellassign <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
saveRDS(cisTopicObject,"cisTOPIC/res/int/cisTopicObject_2.rds")

set.seed(123)
library(Rtsne)
DR <- Rtsne(t(cellassign), pca=F)
DRdist <- dist(DR$Y)
library(densityClust)
dclust <- densityClust(DRdist,gaussian=T)
dclust <- findClusters(dclust, rho = 50, delta = 2.5)
#Distance cutoff calculated to 4.229945

# Check thresholds
pdf("cisTOPIC/res/int/theroshold_check.pdf", height = 6, width = 8)
options(repr.plot.width=6, repr.plot.height=6)
plot(dclust$rho,dclust$delta,pch=20,cex=0.6,xlab='rho', ylab='delta')
points(dclust$rho[dclust$peaks],dclust$delta[dclust$peaks],col="red",pch=20,cex=0.8)
text(dclust$rho[dclust$peaks]-2,dclust$delta[dclust$peaks]+1.5,labels=dclust$clusters[dclust$peaks])
abline(v=50)
abline(h=2.5)
dev.off()

# Add cluster information
densityClust <- dclust$clusters
densityClust <- as.data.frame(densityClust)
rownames(densityClust) <- cisTopicObject@cell.names
colnames(densityClust) <- 'densityClust'
densityClust[,1] <- as.factor(densityClust[,1])
cisTopicObject <- addCellMetadata(cisTopicObject, densityClust)

pdf("cisTOPIC/res/int/cluster.pdf",height = 4,width = 9)
par(mfrow=c(1,3))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c('nCounts', 'nAcc','densityClust'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)
dev.off()

plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c('densityClust'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)

pdf("cisTOPIC/res/int/cluster_sample.pdf",height = 5,width = 5)
par(mfrow=c(1,1))
plotFeatures(cisTopicObject, method='tSNE', target='cell', topic_contr=NULL, colorBy=c('sample'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)
dev.off()

#heatmap
pdf("cisTOPIC/res/int/heatmap.pdf",height = 6,width = 9)
cellTopicHeatmap(cisTopicObject, method='Z-score', colorBy=c('densityClust'), col.low = "dodgerblue", col.mid = "floralwhite", col.high = "brown1")
dev.off()

#heatmap
pdf("cisTOPIC/res/int/heatmap_sample.pdf",height = 6,width = 9)
cellTopicHeatmap(cisTopicObject, method='Z-score', colorBy=c('sample'), col.low = "dodgerblue", col.mid = "floralwhite", col.high = "brown1")
dev.off()

#heatmap
pdf("cisTOPIC/res/int/heatmap_sample_densityClust.pdf",height = 6,width = 9)
cellTopicHeatmap(cisTopicObject, method='Z-score', colorBy=c('sample',"densityClust"), col.low = "dodgerblue", col.mid = "floralwhite", col.high = "brown1")
dev.off()

pdf("cisTOPIC/res/int/heatmap_sample_densityClust_Probability.pdf",height = 6,width = 9)
cellTopicHeatmap(cisTopicObject, method='Probability', colorBy=c('sample',"densityClust"), col.low = "dodgerblue", col.mid = "floralwhite", col.high = "brown1")
dev.off()

pdf("cisTOPIC/res/int/heatmap_sample_Probability.pdf",height = 6,width = 9)
cellTopicHeatmap(cisTopicObject, method='Probability', colorBy=c('sample'), col.low = "dodgerblue", col.mid = "floralwhite", col.high = "brown1")
dev.off()

#To color the Umap by topic score:
pdf("cisTOPIC/res/int/color the Umap by topic score.pdf",height = 9,width = 12)
par(mfrow=c(3,4))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr='Z-score', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
dev.off()

# Plot
par(mfrow=c(1,1))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c("sample"), cex.legend = 0.4, factor.max=.75, dim=2, legend=TRUE, intervals=10)


# Compute cell rankings
library(AUCell)
aucellRankings <- AUCell_buildRankings(pred.matrix, plot=FALSE, verbose=FALSE)

# Check signature enrichment in cells (Reference time: 1 min)
cisTopicObject <- signatureCellEnrichment(cisTopicObject, aucellRankings, selected.signatures='all', aucMaxRank = 0.3*nrow(aucellRankings), plot=FALSE)

# Analysis of the regulatory topics
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

getBigwigFiles(cisTopicObject, path='cisTOPIC/res/int/cisTopics_asBW', seqlengths=seqlengths(txdb))

pdf("cisTOPIC/res/int/topic_region_score.pdf",height = 9,width = 9)
par(mfrow=c(2,6))
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=TRUE)
dev.off()

getBedFiles(cisTopicObject, path='cisTOPIC/res/int/cisTopics_asBed')

##Topic visualization====
cisTopicObject <- runtSNE(cisTopicObject, target='region', perplexity=200, check_duplicates=FALSE)
pdf("cisTOPIC/res/int/region-based-tSNEs-3.pdf",height = 5,width = 5)
par(mfrow=c(1,1))
plotFeatures(cisTopicObject, method='tSNE', target='region', topic_contr=NULL, colorBy=c('nCells'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown', intervals=10)
dev.off()

pdf("cisTOPIC/res/int/region-based-tSNEs-2.pdf",height = 9,width = 9)
par(mfrow=c(3,4))
plotFeatures(cisTopicObject, method='tSNE', target='region', topic_contr='Z-score', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown')
dev.off()

#Annotation to genes and GO terms
library(org.Hs.eg.db)
cisTopicObject <- annotateRegions(cisTopicObject, txdb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb='org.Hs.eg.db')

par(mfrow=c(1,1))
signaturesHeatmap(cisTopicObject, selected.signatures = 'annotation')
ggsave("cisTOPIC/res/int/annotation-1.pdf")5*8

pdf("cisTOPIC/res/int/annotation-2.pdf",height = 9,width = 9)
par(mfrow=c(1,1))
plotFeatures(cisTopicObject, method='tSNE', target='region', 
             topic_contr=NULL, colorBy=c('annotation'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)
dev.off()

  
  
  ###(Transcription factor) motif enrichment====
# cisTopicObject <- binarizedcisTopicsToCtx(cisTopicObject, liftOver=hg38_to_hg19_list, genome="hg19")
# cisTopicObject <- scoredRegionsToCtx(cisTopicObject, liftOver=hg38_to_hg19_list, genome="hg19")
# 
pathToFeather <- "cisTOPIC/hg19-regions-9species.all_regions.mc9nr.feather"
cisTopicObject <- topicsRcisTarget(cisTopicObject, liftOver=hg38_to_hg19_list, genome="hg19", pathToFeather, reduced_database=FALSE, nesThreshold=3, rocthr=0.005, maxRank=20000, nCores=5)

path <- paste0("Topic_motif_enr","_",1:12)

for (i in 1:12){
  A <- cisTopicObject@binarized.RcisTarget[[i]]
  write.csv(A,paste0("cisTOPIC/res/int/",path[i],".csv"))
}


cisTopicObject <- getCistromes(cisTopicObject, annotation = 'Both', nCores=5)
DT::datatable(Topic1_motif_enr[,-c("enrichedRegions", "TF_lowConf"), with=FALSE], escape = FALSE, filter="top", options=list(pageLength=5))

# Compute AUC rankings based on the predictive distribution
pred.matrix <- predictiveDistribution(cisTopicObject)

library(AUCell)
aucellRankings <- AUCell_buildRankings(pred.matrix, plot=FALSE, verbose=FALSE)

cisTopicObject <- getCistromeEnrichment(cisTopicObject, topic=1, TFname='E2F1', aucellRankings = aucellRankings, aucMaxRank = 0.05*nrow(aucellRankings), plot=FALSE)
cisTopicObject <- getCistromeEnrichment(cisTopicObject, topic=1, TFname='YY1', aucellRankings = aucellRankings, aucMaxRank = 0.05*nrow(aucellRankings), plot=FALSE)
cisTopicObject <- getCistromeEnrichment(cisTopicObject, topic=1, TFname='GATA1', aucellRankings = aucellRankings, aucMaxRank = 0.05*nrow(aucellRankings), plot=FALSE)
cisTopicObject <- getCistromeEnrichment(cisTopicObject, topic=1, TFname='DEAF1', aucellRankings = aucellRankings, aucMaxRank = 0.05*nrow(aucellRankings), plot=FALSE)
cisTopicObject <- getCistromeEnrichment(cisTopicObject, topic=1, TFname='POLR2A', aucellRankings = aucellRankings, aucMaxRank = 0.05*nrow(aucellRankings), plot=FALSE)
cisTopicObject <- getCistromeEnrichment(cisTopicObject, topic=1, TFname='NKX2-8', aucellRankings = aucellRankings, aucMaxRank = 0.05*nrow(aucellRankings), plot=FALSE)
cisTopicObject <- getCistromeEnrichment(cisTopicObject, topic=1, TFname='SETDB1', aucellRankings = aucellRankings, aucMaxRank = 0.05*nrow(aucellRankings), plot=FALSE)
cisTopicObject <- getCistromeEnrichment(cisTopicObject, topic=1, TFname='POU4F1', aucellRankings = aucellRankings, aucMaxRank = 0.05*nrow(aucellRankings), plot=FALSE)




cisTopicObject@cell.data$`Topic1_NKX2-8 (16p)`


pdf("cisTOPIC/res/int/Topic1_E2F1.pdf")
par(mfrow=c(1,1))
plotFeatures(cisTopicObject, method='tSNE', target='cell', topic_contr=NULL, colorBy=c('Topic1_E2F1 (17p)'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)
dev.off()
cisTopicObject@cell.data$
  
pdf("cisTOPIC/res/int/Topic1_YY1.pdf")
par(mfrow=c(1,1))
plotFeatures(cisTopicObject, method='tSNE', target='cell', topic_contr=NULL, colorBy=c('Topic1_YY1 (15p)'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)
dev.off()


pdf("cisTOPIC/res/int/Topic1_GATA1.pdf")
par(mfrow=c(1,1))
plotFeatures(cisTopicObject, method='tSNE', target='cell', topic_contr=NULL, colorBy=c('Topic1_GATA1 (53p)'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)
dev.off()

plotFeatures(cisTopicObject, method='tSNE', target='cell', topic_contr=NULL, colorBy=c('Topic1_DEAF1 (9p)'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)
plotFeatures(cisTopicObject, method='tSNE', target='cell', topic_contr=NULL, colorBy=c('Topic1_POLR2A (19p)'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)
plotFeatures(cisTopicObject, method='tSNE', target='cell', topic_contr=NULL, colorBy=c('Topic1_POU4F1 (29p)'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)
plotFeatures(cisTopicObject, method='tSNE', target='cell', topic_contr=NULL, colorBy=c('Topic1_SETDB1 (26p)'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)
plotFeatures(cisTopicObject, method='tSNE', target='cell', topic_contr=NULL, colorBy=c('Topic1_NKX2-8 (16p)'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)
# pdf("cisTOPIC/res/int/Topic1_DEAF1.pdf")
# pdf("cisTOPIC/res/int/Topic1_POLR2A.pdf")
# pdf("cisTOPIC/res/int/Topic1_POU4F1.pdf")
# pdf("cisTOPIC/res/int/Topic1_SETDB1.pdf")
# pdf("cisTOPIC/res/int/Topic1_NKX2-8.pdf")
par(mfrow=c(1,1))
plotFeatures(cisTopicObject, method='tSNE', target='cell', topic_contr=NULL, colorBy=c('Topic1_DEAF1 (9p)'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)
dev.off()

saveRDS(cisTopicObject,"cisTOPIC/res/int/cisTopicObject_5.rds")
cisTopicObject <- readRDS("cisTOPIC/res/int/cisTopicObject_5.rds")

A <-paste0(cisTopicObject@cell.data$sample,"_",cisTopicObject@cell.data$densityClust)
table(cisTopicObject@cell.data$densityClust)  
cellTopicHeatmap(cisTopicObject, method='Probability', colorBy=c('densityClust'))
ggsave("")