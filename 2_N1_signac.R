#### Analysis of scATAC-seq data of NC
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
# read peak/cell matrix
counts <- Read10X_h5(filename = "scATAC_file/N1/N1_filtered_peak_bc_matrix.h5")
# check counts
head(counts)
# check peaks and cell numbers
dim(counts)
# read metadata
metadata <- read.csv(file = "scATAC_file/N1/N1_singlecell.csv",
                     header = TRUE,
                     row.names = 1
)
# check metadata
names(metadata)
head(metadata)[,1:5]

#add fragments information
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'GRCh38',
  fragments = 'scATAC_file/N1/N1_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

# construct Seurat object
N1 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = 'peaks',
  meta.data = metadata
)
# check Seurat object
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

N1 <- RunTFIDF(N1)
# FindTopFeatures
N1 <- FindTopFeatures(N1, min.cutoff = 'q0')
# RunSVD
N1 <- RunSVD(N1)
DepthCor(N1)
ggsave("Signac/N1/N1_DepthCor.png",width = 8, height = 6, dpi = 600)

# UMAP
N1 <- RunUMAP(object = N1, reduction = 'lsi', dims = 2:30)

# cluster
N1 <- FindNeighbors(object = N1, reduction = 'lsi', dims = 2:30)
N1 <- FindClusters(object = N1, verbose = FALSE, algorithm = 3)
# DimPlot visulization
DimPlot(object = N1, label = TRUE, pt.size = 1, label.size = 6, )
ggsave("Signac/N1/N1_cluster.pdf",width = 8, height = 6, dpi = 600)

##gene activity matix
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

### FindMarkers====
Idents(N1) <- N1$group
DefaultAssay(N1) <- "peaks"
da_peaks <- FindAllMarkers(
  object = N1,
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
write.csv(da_peaks,"Signac/N1/peak_markers.csv")


CoveragePlot(object = N1, features = "RAD51B", region = "chr14-61832750-67833500", extend.upstream = 10000, extend.downstream = 10000)
ggsave("Signac/N1/RAD51B.pdf",height = 4, width = 6)
