####Time210826
setwd("E://A549")
library(Signac)
library(Seurat)
library(dplyr)
library(igraph)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
set.seed(1234)

##Pre-processing workflow
# 读取peak/cell matrix矩阵
counts <- Read10X_h5(filename = "scATAC_file/N1/N1_filtered_peak_bc_matrix.h5")
# 查看counts信息
head(counts)
# 查看peaks和细胞的数目
dim(counts)
# 读取细胞注释信息
metadata <- read.csv(file = "scATAC_file/N1/N1_singlecell.csv",
                     header = TRUE,
                     row.names = 1
)
# 查看metadata信息
names(metadata)
head(metadata)[,1:5]

#加入fragments信息
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'GRCh38',
  fragments = 'scATAC_file/N1/N1_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

# 构建Seurat对象
N1 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = 'peaks',
  meta.data = metadata
)
# 查看Seurat对象
N1
N1[['peaks']]

#see the genomic ranges associated with each feature in the object
granges(N1)

#annotation
#extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
class(annotations)
# change to UCSC style since the data was mapped to GRCh38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "GRCh38"

write.table(annotation, "annotation/annotation_GRCh38_extracted.txt")

# add the gene information to the object
Annotation(N1) <- annotations

##QC
# compute nucleosome signal score per cell
N1 <- NucleosomeSignal(object = N1)
# compute TSS enrichment score per cell
N1 <- TSSEnrichment(object = N1, fast = FALSE)
# add blacklist ratio and fraction of reads in peaks
N1$pct_reads_in_peaks <- N1$peak_region_fragments / N1$passed_filters * 100
N1$blacklist_ratio <- N1$blacklist_region_fragments / N1$peak_region_fragments

N1$high.tss <- ifelse(N1$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(N1, group.by = 'high.tss') + NoLegend()
ggsave("Signac/N1/N1_TSS_enrichment.png",width = 8, height = 6, dpi = 600)
# 根据核小体条带信号强度进行分组
N1$nucleosome_group <- ifelse(N1$nucleosome_signal > 1, 'NS > 1', 'NS < 1')
FragmentHistogram(object = N1, group.by = 'nucleosome_group')
ggsave("Signac/N1/N1_核小体条带分布模式.png",width = 8, height = 6, dpi = 600)

#quality control visualization
p1 <- VlnPlot(object = N1, c('pct_reads_in_peaks', 'peak_region_fragments','TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'), pt.size = 0.1, ncol = 5)
ggsave("Signac/N1/N1_QC.png",width = 16, height = 8, dpi = 600)
N1

# 根据不同QC指标进行数据过滤
#N1 <- subset(
#  x = N1,
#  subset = peak_region_fragments > 500 &
#    peak_region_fragments < 15000 &
#    pct_reads_in_peaks > 15 &
#   blacklist_ratio < 0.005 &
#    nucleosome_signal < 4 &
#    TSS.enrichment > 2
#)
N1
# 使用RunTFIDF函数进行数据归一化
N1 <- RunTFIDF(N1)
# 使用FindTopFeatures函数选择可变的features
N1 <- FindTopFeatures(N1, min.cutoff = 'q0')
# 使用RunSVD函数进行线性降维
N1 <- RunSVD(N1)
DepthCor(N1)
ggsave("Signac/N1/N1_DepthCor.png",width = 8, height = 6, dpi = 600)

# 使用RunUMAP函数进行UMAP非线性降维
N1 <- RunUMAP(object = N1, reduction = 'lsi', dims = 2:30)
# 对细胞执行基于图的聚类
N1 <- FindNeighbors(object = N1, reduction = 'lsi', dims = 2:30)
N1 <- FindClusters(object = N1, verbose = FALSE, algorithm = 3)
# 使用DimPlot函数进行数据可视化
DimPlot(object = N1, label = TRUE, pt.size = 1, label.size = 6, )
ggsave("Signac/N1/N1_cluster.pdf",width = 8, height = 6, dpi = 600)
ggsave("Signac/N1/N1_cluster.pdf",width = 8, height = 6, dpi = 600)

##创建基因活性矩阵。
#之前的聚类区域所用的features是peaks，为了展示不同分群基因活性的差异，首先要创建一个类似RNA表达的矩阵。
#用基因加上游2000bp区域的比对片段数代表该基因的活性。
gene.activities <- GeneActivity(N1)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
N1[['RNA']] <- CreateAssayObject(counts = gene.activities)
N1[['RNA']] <- CreateAssayObject(counts = gene.activities)
N1 <- NormalizeData(
  object = N1,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(N1$nCount_RNA)
)

saveRDS(N1, file = "Signac/N1/N1.rds")

N1 <- readRDS("E:/A549/Signac/N1/N1_nonfiltered.rds")
N1$seurat_clusters

###scRNA_NB analysis with seurat
# Load data
RNA.NC <- readRDS("scCANCER/res/NC/NC_seurat.rds")
DefaultAssay(N1) <- "RNA"
# 使用FindTransferAnchors函数识别整合的anchors
transfer.anchors <- FindTransferAnchors(
  reference = RNA.NC,
  query = N1,
  reduction = 'cca'
)

# 使用TransferData函数基于识别出的anchors进行数据映射整合
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = RNA.NC$group,
  weight.reduction = N1[['lsi']],
  dims = 2:30
)

#N1 <- AddMetaData(object = N1, metadata = predicted.labels)
N1$group <- paste0("cluster",as.numeric(N1$seurat_clusters))
N1$group <- factor(N1$group, levels = paste0("cluster",1:10))
table(N1$group)
N1[[]]

Idents(N1) <- "group"

# 数据可视化
plot1 <- DimPlot(
  object = RNA.NC,
  group.by = 'group',
  pt.size =1,
  label.size = 5,
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

N1$seurat_clusters
DimPlot(
  object = N1,
  group.by = 'group',
  pt.size =1.5,
  label.size = 5,
  label = TRUE,
  repel = TRUE)
ggsave("Signac/N1/umap_220704.pdf")
ggsave("Signac/N1/umap_220704.png")

table(N1$predicted.id)
table(RNA.NC$group)

plot1 + plot2
ggsave("Signac/N1/N1_NC_group.pdf")
ggsave("Signac/N1/N1_NC_group.png")
a <- table(N1$predicted.id)
write.csv(a, "Signac/N1/N1_NC_cellcycle_summary.csv")

# 数据过滤 
#N1.filtered <- subset(N1, subset = prediction.score.max > 0.5) 
hist(N1$prediction.score.max) 


# switch back to working with peaks instead of gene activities
DefaultAssay(N1) <- 'RNA'

### 使用FindMarkers函数进行差异分析====
Idents(N1) <- N1$group
DefaultAssay(N1) <- "peaks"
da_peaks <- FindAllMarkers(
  object = N1,
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.csv(da_peaks,"Signac/N1/peak_markers_220704.csv")

Idents(N1) <- N1$group
N1@assays$
DefaultAssay(N1) <- "RNA"
da_genes <- FindAllMarkers(
  object = N1,
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.csv(da_genes,"Signac/N1/gene_markers_220704.csv")


CoveragePlot(object = N1, features = "DDX10", region = "chr11-108775380-108775896", extend.upstream = 10000, extend.downstream = 10000)
ggsave("Signac/N1/DDX10.pdf",height = 4, width = 6)
ggsave("Signac/N1/DDX10.png")
