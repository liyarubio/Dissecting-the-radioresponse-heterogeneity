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
saveRDS(combined,"Signac/integration/combined_unfiltered.rds")
combined <- readRDS("Signac/integration/combined_unfiltered.rds")
# 使用NucleosomeSignal函数计算Nucleosome banding pattern
combined <- NucleosomeSignal(object = combined)
combined1 <- combined
###QC====
# 计算peaks中reads的比例
combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$passed_filters * 100
# 计算比对到“blacklist”区域的reads比率
combined$blacklist_ratio <- combined$blacklist_region_fragments / combined$peak_region_fragments
# 质控指标的可视化
p1 <- VlnPlot(combined, c('pct_reads_in_peaks', 'peak_region_fragments'), pt.size = 0.1)
p2 <- VlnPlot(combined, c('blacklist_ratio', 'nucleosome_signal'), pt.size = 0.1) & scale_y_log10()

p1 | p2
ggsave("Signac/integration/QC.pdf")
ggsave("Signac/integration/QC.png")

# # 根据核小体条带信号强度进行分组
# combined$nucleosome_group <- ifelse(combined$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
# # 使用FragmentHistogram函数可视化核小体条带分布模式
# FragmentHistogram(object = combined, group.by = 'nucleosome_group')
# ggsave("Signac/integration/QC_1.pdf")
# ggsave("Signac/integration/QC_1.png")
# 
# # create granges object with TSS positions
# gene.ranges <- genes(EnsDb.Hsapiens.v86)
# seqlevelsStyle(gene.ranges) <- 'UCSC'
# gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
# gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')
# 
# tss.ranges <- GRanges(
#   seqnames = seqnames(gene.ranges),
#   ranges = IRanges(start = start(gene.ranges), width = 2),
#   strand = strand(gene.ranges)
# )
# 
# seqlevelsStyle(tss.ranges) <- 'UCSC'
# tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')

# # to save time use the first 2000 TSSs
# # 使用TSSEnrichment函数计算TSS富集得分
# combined <- TSSEnrichment(object = combined, tss.positions = tss.ranges[1:2000])
# combined$high.tss <- ifelse(combined$TSS.enrichment > 2, 'High', 'Low')
# # 使用TSSPlot函数对TSS富集得分进行可视化
# TSSPlot(combined, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()

# make sure to change to the assay containing common peaks
DefaultAssay(combined) <- "peaks"
# 对合并后的数据进行归一化，降维与可视化
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(
  combined,
  reduction.key = 'LSI_',
  reduction.name = 'lsi',
  irlba.work = 400
)
combined <- RunUMAP(combined, dims = 2:30, reduction = 'lsi')
combined1 <- combined
combined <- combined1
###visulization
table(combined$sample)
combined$sample <- factor(combined$sample, levels = c('NC','IR_immediate','IR_2h'))
Idents(combined) <- table(combined$sample)
DimPlot(combined, pt.size = 0.5, group.by = "sample")
ggsave("Signac/integration/dimplot_sample.pdf", width = 8, height = 6)
ggsave("Signac/integration/dimplot_sample.png", width = 8, height = 6)

combined <- FindClusters(combined)
DepthCor(combined)
ggsave("Signac/integration/lsi_correlation.pdf")
# 对细胞执行基于图的聚类
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:30)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3)
combined$seurat_clusters
combined$cluster <- paste0("cluster",as.numeric(combined$seurat_clusters))
combined$cluster <- factor(combined$cluster, levels = paste0("cluster",1:11))
combined$sample_cluster <- paste0(combined$sample,"_cluster",as.numeric(combined$seurat_clusters))
a <- table(combined$sample_cluster)
write.csv(a,"Signac/integration/sample_cluster.csv")
saveRDS(combined,"Signac/integration/combined_cluster.rds")

#鉴定不同聚类群之间差异可及性peaks区域
# 使用FindMarkers函数进行差异分析
Idents(combined) <- combined$sample
IR2h_NC_da_down <- FindMarkers(
  object = combined,
  ident.1 = "NC",
  ident.2 = "IR_2h",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.csv(da_peaks_C3_N1,"Signac/integration/IR2h_NC_da_down.csv")

IR2h_NC_da_up <- FindMarkers(
  object = combined,
  ident.1 = "IR_2h",
  ident.2 = "NC",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.csv(da_peaks_C3_N1,"Signac/integration/IR2h_NC_da_up.csv")

# 差异分析结果可视化
# plot1 <- VlnPlot(
#   object = combined,
#   features = rownames(da_peaks_C3_N1)[2],
#   pt.size = 0.1,
#   idents = c("NC","IR_2h")
# )
# plot2 <- FeaturePlot(
#   object = combined,
#   features = rownames(da_peaks_C3_N1)[2],
#   pt.size = 0.1
# )
# 
# plot1 | plot2



IR2h_IR_immediate_da_down <- FindMarkers(
  object = combined,
  ident.1 = "IR_immediate",
  ident.2 = "IR_2h",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.csv(da_peaks_C3_N1,"Signac/integration/IR2h_IR_immediate_da_down.csv")

IR2h_IR_immediate_da_up <- FindMarkers(
  object = combined,
  ident.1 = "IR_2h",
  ident.2 = "IR_immediate",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.csv(da_peaks_C3_N1,"Signac/integration/IR2h_IR_immediate_da_up.csv")



# 提取显著的差异可及性peak区域信息
open_IR2h_NC <- rownames(IR2h_NC_da_down[IR2h_NC_da_down$avg_log2FC > 0.5, ])
open_IR2h_IR_immediate <- rownames(IR2h_IR_immediate_da_down[IR2h_IR_immediate_da_down$avg_log2FC >0.5, ])

# create granges object with TSS positions
gene.ranges <- genes(EnsDb.Hsapiens.v86)
seqlevelsStyle(gene.ranges) <- 'UCSC'
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

# 使用ClosestFeature函数对差异开放性peak区域进行注释
closest_genes_IR2h_NC <- ClosestFeature(combined,
  regions = open_IR2h_NC,
  annotation = gene.ranges,
  sep = c(':', '-')
)
write.csv(closest_genes_IR2h_NC,"Signac/integration/closest_genes_IR2h_NC_log2FC_0.5.csv")

closest_genes_IR2h_IR_immediate <- ClosestFeature(combined,
  regions = open_IR2h_IR_immediate,
  annotation = gene.ranges,
  sep = c(':', '-')
)
write.csv(closest_genes_IR2h_IR_immediate,"Signac/integration/closest_genes_IR2h_IR_immediate_log2FC_0.5.csv")

closest_genes_IR2h_NC_1 <- ClosestFeature(combined,
                                        regions = rownames(IR2h_NC_da_down),
                                        annotation = gene.ranges,
                                        sep = c(':', '-')
)
write.csv(closest_genes_IR2h_NC_1,"Signac/integration/closest_genes_IR2h_NC.csv")

closest_genes_IR2h_IR_immediate_1 <- ClosestFeature(combined,
                                                  regions = rownames(IR2h_IR_immediate_da_down),
                                                  annotation = gene.ranges,
                                                  sep = c(':', '-')
)
write.csv(closest_genes_IR2h_IR_immediate_1,"Signac/integration/closest_genes_IR2h_IR_immediate.csv")

# 差异分析结果可视化
# plot1 <- VlnPlot(
#   object = combined,
#   features = rownames(da_peaks_C3_N1)[2],
#   pt.size = 0.1,
#   idents = c("NC","IR_2h")
# )
# plot2 <- FeaturePlot(
#   object = combined,
#   features = rownames(da_peaks_C3_N1)[2],
#   pt.size = 0.1
# )
# 
# plot1 | plot2

###创建基因表达活性矩阵（gene activity matrix）====
# Extend coordinates upstream to include the promoter
genebodyandpromoter.coords <- Extend(x = gene.ranges, upstream = 2000, downstream = 0)
Fragments(combined)
# create a gene by cell matrix
gene.activities <- FeatureMatrix(
  fragments = Fragments(combined),
  features = genebodyandpromoter.coords,
  cells = colnames(combined)
)

# convert rownames from chromsomal coordinates into gene names
gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- gene.key[rownames(gene.activities)]

# add the gene activity matrix to the Seurat object as a new assay, and normalize it
combined[['RNA']] <- CreateAssayObject(counts = gene.activities)
combined <- NormalizeData(
  object = combined,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_RNA)
)
saveRDS(combined,"Signac/integration/combined_final.rds")

DefaultAssay(combined) <- 'peaks'

DimPlot(
  object = combined,
  pt.size = 0.1,
  split.by = 'sample',
  ncol = 3
)
# table(combined$sample)
ggsave("Signac/integration/dimplot_sample.pdf", height = 5, width = 8)


#使用CoveragePlot函数绘制peaks的覆盖图
# CoveragePlot(
#   object = combined,
#   region = rownames(da_peaks)[c(1,5)],
#   sep = c(":", "-"),
#   peaks = StringToGRanges(rownames(combined), sep = c(":", "-")),
#   annotation = gene.ranges,
#   extend.upstream = 20000,
#   extend.downstream = 20000,
#   ncol = 1
# )


combined$cluster <- factor(combined$cluster, levels = paste0('cluster',1:11))
####Analysis by cluster
combined <- readRDS("Signac/integration/combined_final.rds")
Idents(combined) <- combined$sample_cluster
da_peaks <- FindMarkers(combined, ident.1 = 'NC_cluster2', 'IR_2h_cluster4')
da_peaks <- da_peaks[which(da_peaks$p_val_adj<0.05),]
da_peaks_1 <- da_peaks[which(da_peaks$avg_log2FC > 0.5),]
closest_genes <- ClosestFeature(combined,
                                regions = rownames(da_peaks_1),
                                annotation = gene.ranges,
                                sep = c(':', '-'))
da_peaks_anno <- cbind(da_peaks_1,closest_genes)
write.csv(da_peaks_anno,'Signac/integration/NC_2_IR_2h_4.csv')

Idents(combined) <- combined$peak_region_fragments
VlnPlot(combined, features = 'peak_region_fragments', group.by = 'cluster', split.by = 'sample')
VlnPlot(combined, features = 'nfeature_peaks', group.by = 'cluster', split.by = 'sample')

combined$nCount_peaks

NC <- subset(combined, idents = 'NC')
Idents(NC) <- NC$cluster
da_peaks_4_1 <- FindMarkers(NC, ident.1 = 'cluster1', ident.2 = 'cluster4')
da_peaks <- da_peaks[which(da_peaks$p_val_adj<0.05),]
da_peaks_1 <- da_peaks[which(da_peaks$avg_log2FC > 0.25),]

closest_genes_1_4 <- ClosestFeature(combined,
                                  regions = rownames(da_peaks_1_4),
                                  annotation = gene.ranges,
                                  sep = c(':', '-'))
write.csv(closest_genes_NC,"Signac/integration/closest_genes_NC_log2FC_0.5.csv")


# 直接使用merge函数将多个数据集进行合并
combined <- readRDS("Signac/integration/combined_final.rds")
DefaultAssay(combined) <- 'peaks'
# Get a list of motif position frequency matrices from the JASPAR database
# 使用getMatrixSet函数从JASPAR数据库中提取Motif的PFM矩阵信息
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# Scan the DNA sequence of each peak for the presence of each motif
# 使用CreateMotifMatrix函数构建Motif矩阵对象
# Add the Motif combined to the assay
# 使用AddMotifcombined函数将Motif类添加到Seurat对象???
combined[['peaks']] <- AddMotifs(
  object = combined[['peaks']],
  pfm = pfm,
  genome = 'hg38',
)
# 使用RegionStats函数计算peaks区域的序列特???
combined <- RegionStats(
  object = combined,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c("-", "-")
)

# # 使用FindMarkers函数鉴定差异可及性peaks
# da_peaks <- FindMarkers(
#   object = combined,
#   ident.1 = '0',
#   ident.2 = '1',
#   only.pos = TRUE,
#   test.use = 'LR',
#   latent.vars = 'nCount_peaks'
# )
# 
# # Test the differentially accessible peaks for overrepresented motifs
# # 使用FindMotifs函数进行motif富集分析
# enriched.motifs <- FindMotifs(
#   object = combined,
#   features = head(rownames(da_peaks), 1000)
# )
# 
# # 查看motif富集分析的结果
# # sort by p-value and fold change
# enriched.motifs <- enriched.motifs[order(enriched.motifs[, 7], -enriched.motifs[, 6]), ]
# head(enriched.motifs)

# 使用RunChromVAR函数计算所有细胞中的motif activities
combined <- RunChromVAR(
  object = combined,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(combined) <- 'chromvar'
combined <- computeDeviations(combined)
computeDeviations(combined,
                  annotations, background_peaks = getBackgroundPeaks(combined),
                  expectation = computeExpectations(combined))

A@x
# # look at the activity of Mef2c, the top hit from the overrepresentation testing
# p2 <- FeaturePlot(
#   object = combined,
#   features = rownames(enriched.motifs)[[1]],
#   min.cutoff = 'q10',
#   max.cutoff = 'q90',
#   pt.size = 0.1
# )
# p1 + p2
# ggsave()
saveRDS(combined, "Signac/integration/combined_chrom_VAR.rds")
combined <- readRDS("Signac/integration/combined_chrom_VAR.rds")

Idents(combined) <- combined$orig.ident
differential.activity_all <- read.csv("Signac/integration/TF/differential.activity_all.csv")
motif_name <- read.csv("Signac/integration/TF/motif_name.csv")
TEMP_x <- substr(differential.activity_all$X, 1,8)
differential.activity_all$X <- TEMP_x
differential.activity_all <- merge(differential.activity_all, motif_name)
write.csv(differential.activity_all, "Signac/integration/TF/differential.activity_all.csv")

#测试不同细胞类型之间motif的差异活性得分，这和对不同细胞类型之间的差异可及性peaks进行富集测试的结果相类似???
differential.activity_N1_C3 <- FindMarkers(
  object = combined,
  ident.1 = 'NC',
  ident.2 = 'IR_2h',
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks')

# 使用MotifPlot函数对富集到的motif进行可视???
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

# 使用MotifPlot函数对富集到的motif进行可视???
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
# 使用MotifPlot函数对富集到的motif进行可视???
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

# 使用MotifPlot函数对富集到的motif进行可视???
MotifPlot(
  object = combined,
  motifs = rownames(differential.activity_C2_N1),
  assay = 'peaks'
)
ggsave("Signac/integration/TF/differential.motif_C2_N1.pdf")

differential.activity_C2_N1 <- read.csv("Signac/integration/TF/differential.activity_C2_N1.csv")
differential.activity_C2_N1 <- merge(differential.activity_C2_N1,motif.name, all = FALSE)
write.csv(differential.activity_C2_N1,"Signac/integration/TF/differential.activity_C2_N1.csv")

#加载包
library(ggplot2)
library(ggrepel)
#读取数据
Dat1 <- read.csv("Signac/integration/TF/differential.activity_C3_N1.csv")
Dat2 <- read.csv("Signac/integration/TF/differential.activity_N1_C3.csv")
Dat2 <- Dat2[,-1]
Dat2$avg_log2FC <- -Dat2$avg_log2FC
Dat <- rbind(Dat1,Dat2)
#确定是上调还是下调，用于给图中点上色）
Dat$threshold = factor(ifelse(Dat$p_val_adj < 0.05 & abs(Dat$avg_log2FC) >= 0.2, ifelse(Dat$avg_log2FC>= 0.2,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))

ggplot(Dat,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
  geom_text_repel(
    data = Dat[which(Dat$p_val_adj<0.05 & abs(Dat$avg_log2FC)>0.2),],
    aes(label = V1),
    size = 3,
    segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
  theme_bw()+#修改图片背景
  theme(
    legend.title = element_blank()#不显示图例标题
  )+
  ylab('-log10 (p-adj)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05
ggsave('Signac/integration/TF/IR2h_NC_colcano plot.pdf', width = 6, height = 4)

#读取数据
Dat1 <- read.csv("Signac/integration/TF/differential.activity_C2_N1.csv")
Dat2 <- read.csv("Signac/integration/TF/differential.activity_N1_C2.csv")
Dat2$avg_log2FC <- -Dat2$avg_log2FC
Dat <- rbind(Dat1,Dat2)
#确定是上调还是下调，用于给图中点上色）
Dat$threshold = factor(ifelse(Dat$p_val_adj < 0.05 & abs(Dat$avg_log2FC) >= 0.25, ifelse(Dat$avg_log2FC>= 0.25,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))

ggplot(Dat,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
  geom_text_repel(
    data = Dat[which(Dat$p_val_adj<0.05 & abs(Dat$avg_log2FC)>0.25),],
    aes(label = V1),
    size = 3,
    segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
  theme_bw()+#修改图片背景
  theme(
    legend.title = element_blank()#不显示图例标题
  )+
  ylab('-log10 (p-adj)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  geom_vline(xintercept=c(-0.25,0.25),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05
ggsave('Signac/integration/TF/IRimmediate_NC_colcano plot_0.25.pdf', width = 6, height = 4)


####--------- cistopic-----------------------
library(cisTopic)
##Load data====
pathTo10X <- 'E:/A549/scATAC_cellranger/A549_N1/outs/'
data_folder <- paste0(pathTo10X, 'filtered_peak_bc_matrix')
metrics <- paste0(pathTo10X, 'singlecell.csv')
cisTopicObject <- createcisTopicObjectFrom10Xmatrix(data_folder, metrics,  project.name='NC')

##Building the models====
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(2, 5, 10:25, 30, 35, 40), seed=987, nCores=20, iterations = 500, addModels=FALSE)
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
pdf("cisTOPIC/res/N1/model_topic_number.pdf", height = 9, width = 9)
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

set.seed(123)
library(Rtsne)
DR <- Rtsne(t(cellassign), pca=F)
DRdist <- dist(DR$Y)
library(densityClust)
dclust <- densityClust(DRdist,gaussian=T)
dclust <- findClusters(dclust, rho = 50, delta = 2.5)
#Distance cutoff calculated to 4.229945

# Check thresholds
pdf("cisTOPIC/res/N1/theroshold_check.pdf", height = 6, width = 8)
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

pdf("cisTOPIC/res/N1/cluster.pdf",height = 4,width = 9)
par(mfrow=c(1,3))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c('nCounts', 'nAcc','densityClust'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)
dev.off()

pdf("cisTOPIC/res/N1/cluster_1.pdf",height = 5,width = 5)
par(mfrow=c(1,1))
plotFeatures(cisTopicObject, method='tSNE', target='cell', topic_contr=NULL, colorBy=c('densityClust'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)
dev.off()

#heatmap
pdf("cisTOPIC/res/N1/heatmap.pdf",height = 6,width = 9)
cellTopicHeatmap(cisTopicObject, method='Z-score', colorBy=c('densityClust'), col.low = "dodgerblue", col.mid = "floralwhite", col.high = "brown1")
dev.off()

#To color the Umap by topic score:
pdf("cisTOPIC/res/N1/color the Umap by topic score.pdf",height = 9,width = 12)
par(mfrow=c(3,4))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr='Z-score', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
dev.off()

#Enrichment of epigenomic signatures in the cells
pred.matrix <- predictiveDistribution(cisTopicObject)

# Compute cell rankings
library(AUCell)
aucellRankings <- AUCell_buildRankings(pred.matrix, plot=FALSE, verbose=FALSE)

# Check signature enrichment in cells (Reference time: 1 min)
# cisTopicObject <- signatureCellEnrichment(cisTopicObject, aucellRankings, selected.signatures='all', aucMaxRank = 0.3*nrow(aucellRankings), plot=FALSE)

# Analysis of the regulatory topics
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

getBigwigFiles(cisTopicObject, path='cisTOPIC/res/N1/cisTopics_asBW', seqlengths=seqlengths(txdb))

pdf("cisTOPIC/res/N1/topic_region_score.pdf",height = 9,width = 9)
par(mfrow=c(2,6))
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=TRUE)
dev.off()

getBedFiles(cisTopicObject, path='cisTOPIC/res/N1/cisTopics_asBed')

##Topic visualization====
cisTopicObject <- runtSNE(cisTopicObject, target='region', perplexity=200, check_duplicates=FALSE)
pdf("cisTOPIC/res/N1/region-based-tSNEs-3.pdf",height = 5,width = 5)
par(mfrow=c(1,1))
plotFeatures(cisTopicObject, method='tSNE', target='region', topic_contr=NULL, colorBy=c('nCells'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown', intervals=10)
dev.off()
pdf("cisTOPIC/res/N1/region-based-tSNEs-2.pdf",height = 9,width = 9)
par(mfrow=c(3,4))
plotFeatures(cisTopicObject, method='tSNE', target='region', topic_contr='Z-score', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown')
dev.off()

#Annotation to genes and GO terms
library(org.Hs.eg.db)
cisTopicObject <- annotateRegions(cisTopicObject, txdb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb='org.Hs.eg.db')

par(mfrow=c(1,1))
signaturesHeatmap(cisTopicObject, selected.signatures = 'annotation')
#ggsave("cisTOPIC/res/N1/annotation-1.pdf")5*8

pdf("cisTOPIC/res/N1/annotation-2.pdf",height = 9,width = 9)
par(mfrow=c(1,1))
plotFeatures(cisTopicObject, method='tSNE', target='region', 
             topic_contr=NULL, colorBy=c('annotation'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)
dev.off()

#identifying enriched GO terms
#Running rGREAT and RcisTarget with hg38
# url and file name for a chain file
library(R.utils)
url <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
hg38ToHg19.chain <- "data/hg38ToHg19.over.chain"
download.file(url, destfile = paste0(hg38ToHg19.chain, ".gz"))
gunzip(paste0(hg38ToHg19.chain, ".gz"))

# Import chain file
library(rtracklayer)
hg38ToHg19.chain <- import.chain(hg38ToHg19.chain)

# Obtain liftOver dictionary (as list)
hg38_coord <- cisTopicObject@region.ranges
hg38_to_hg19_list <- liftOver(hg38_coord, hg38ToHg19.chain)

# Run GREAT based on liftover to hg19 coordinates
cisTopicObject <- GREAT(cisTopicObject, genome='hg19', liftOver=hg38_to_hg19_list, fold_enrichment=2, geneHits=1, sign=0.05, request_interval=10)
cisTopicObject <- binarizedcisTopicsToCtx(cisTopicObject, liftOver=hg38_to_hg19_list, genome="hg19")
cisTopicObject <- scoredRegionsToCtx(cisTopicObject, liftOver=hg38_to_hg19_list, genome="hg19")

# pathToFeather <- "hg19-regions-1M-9species.all_regions.mc9nr.feather"
# cisTopicObject <- topicsRcisTarget(cisTopicObject, genome='hg19', pathToFeather, reduced_database=FALSE, nesThreshold=3, rocthr=0.005, maxRank=20000, nCores=1)
# cisTopicObject<- getCistromes(cisTopicObject, annotation = 'Both', nCores=5)


ontologyDotPlot(cisTopicObject, top=5, topics=c(1:12), var.y='name', order.by='Binom_Adjp_BH')
ggsave("cisTOPIC/res/N1/GREAT_GO.pdf", height = 10, width = 12)

saveRDS(cisTopicObject,"cisTOPIC/res/N1/cisTopicObject.rds")

