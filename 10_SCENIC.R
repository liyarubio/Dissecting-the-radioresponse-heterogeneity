###SCENIC downstream analysis
###20211023 
# Required packages:
library(SCopeLoomR)
library(AUCell)
library(SCENIC)

# For some of the plots:
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(Seurat)
library(htmltools)
packageVersion("SCENIC")

#1.2.4

setwd("E:/A549")
scenicLoomPath <- file.path("SCENIC/res/NC_SCENIC.loom")
file.exists(scenicLoomPath)

#motifEnrichmentFile <- file.path("pbmc10k_data__reg_mtf.csv.gz")
list.files()


loom <- open_loom(scenicLoomPath)
# Read information from loom file:
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)

# embeddings <- get_embeddings(loom)
# cellClusters <- get_clusterings(loom)

close_loom(loom)

## To check whether it was read properly:
length(regulons)
head(names(regulons))
regulonAUC
#length(regulonsAucThresholds)
#plot(embeddings$`SCENIC AUC UMAP`)


####------Scoring the network activity----------------
regulonAUC

####------Regulators for known cell types or clusters--
## To start from clusters/cell types from Scanpy:

####----------RNA-----------------------------
NC <- readRDS("scCANCER/res/SCRAN/NC.rds")
NC <- as.Seurat(NC)
NC$clusters <- paste0("cluster",as.numeric(NC$seurat_clusters))
Idents(NC) <- NC$clusters
DimPlot(NC)

cellClusters <- as.data.frame(NC@meta.data$clusters)
rownames(cellClusters) <- rownames(NC@meta.data)
colnames(cellClusters) <- "0"
head(cellClusters)


selectedResolution <- "0" # select resolution
# Split the cells by cluster:
cellsPerCluster <- split(rownames(cellClusters), cellClusters[,selectedResolution]) 

regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

####--------heatmap----------------------------
# plot:
options(repr.plot.width=8, repr.plot.height=10) # To set the figure size in Jupyter
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6))) # row font size
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6))) # row font size


class(regulonActivity_byCellType)
regulonActivity_byCellType <- as.data.frame(regulonActivity_byCellType)
regulonActivity_byCellType[1:10,]
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[1:10,], name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6))) # row font size


regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later


topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators)
viewTable(topRegulators, options = list(pageLength = 10))



####-------Cell-type specific regulators-----------------
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellClusters[colnames(regulonAUC), selectedResolution])
head(rss)
dim(rss)

## Showing regulons and cell types with any RSS > 0.01 
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter


for (i in 1:10){
  plotRSS_oneSet(rss, setName = paste0("cluster",i)) # cluster ID
}

plotRSS_oneSet(rss, setName = "cluster10") # cluster ID


#SCENIC::plotRSS_oneSet

rss
#setName = "FALSE"
setName = "cluster1"
n=10

library(ggplot2)
library(ggrepel)
rssThisType <- sort(rss[, setName], decreasing = TRUE)
thisRss <- data.frame(regulon = names(rssThisType), rank = seq_along(rssThisType), 
                      rss = rssThisType)
thisRss$regulon[(n + 1):nrow(thisRss)] <- NA
ggplot(thisRss, aes(x = rank, y = rss)) + 
  geom_point(color = "blue", size = 1) + 
  ggtitle(setName) + 
  geom_label_repel(aes(label = regulon), 
                   max.overlaps=20,
                   box.padding = 0.35, point.padding = 0.5, 
                   segment.color = "grey50", na.rm = TRUE) + 
  theme_classic()                               


#Cell states based on the GRN activity
#cat ( names ( embeddings ),  sep = "\n" )


####-------RNA embedding------------------------
library(Seurat)
NC_1 <- readRDS("scCANCER/res/NC/NC_seurat.rds")
NC_1@meta.data <- NC@meta.data
NC<- NC_1
NC.embed <- Embeddings(NC, reduction = "umap")
#embedding 细胞带有样本信息

m <- NC@meta.data
regulonAUC

# Overview of these embeddings (see below for details)

options(repr.plot.width=10, repr.plot.height=8) # To set the figure size in Jupyter
# par(mfrow=c(2, ceiling(length(names(embeddings))/2)))

rownames(regulonAUC)
TF <- rownames(regulonAUC)[1:10]

regulonsToPlot <- TF
regulonAUC[regulonsToPlot,]

pdf("SCENIC/res/pic/NC/TF.pdf")
t <- AUCell::AUCell_plotTSNE(NC.embed, exprMat_log, regulonAUC[regulonsToPlot,], 
                             plots=c("AUC"), cex = .5)
dev.off()

options(repr.plot.width=10, repr.plot.height=8) # To set the figure size in Jupyter
par(mfrow=c(2,3))
tfsToPlot <- rownames(regulonAUC)[1:10]
AUCell::AUCell_plotTSNE(NC.embed, exprMat_log[tfsToPlot,], plots=c("Expression"), cex = .5)
AUCell::AUCell_plotTSNE(NC.embed, exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5)



AUCell::AUCell_plotTSNE(NC.embed, exprMat_log, regulonAUC["ALX1(+)",], plots=c("AUC"), cex = .5)

# Density plot to detect most likely stable states
options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter
library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(NC.embed, .5)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=6, drawlabels=FALSE)


###--------Binarized regulon activity-----------
head(as.data.frame(regulonAucThresholds))
t <- names(regulonAucThresholds)
t <- as.numeric(t)
names(t) <- regulonAucThresholds
head(as.data.frame(t))
regulonAucThresholds <- t

pdf("SCENIC/res/pic/NC/TF_ON_OFF.pdf")
options(repr.plot.width=10, repr.plot.height=4) # To set the figure size in Jupyter
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(NC.embed, exprMat_log, 
                        regulonAUC[regulonsToPlot,], thresholds = regulonAucThresholds[regulonsToPlot],
                        plots=c("AUC", "histogram", "binary"), cex = .5)
dev.off()


# This function will be included in the next version of AUCell
binarizeAUC <- function(auc, thresholds)
{
  thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
  regulonsCells <- setNames(lapply(names(thresholds), 
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(getAUC(auc)[x,]>trh))
                                   }),names(thresholds))
  
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"  
  
  return(binaryRegulonActivity)
}
binaryRegulonActivity <- binarizeAUC(regulonAUC, regulonAucThresholds)
dim(binaryRegulonActivity)
binaryRegulonActivity[1:5,1:3]

# Subset some cells to make the plots faster:
nCells <- 1000
set.seed(123)
cellsSelected <- sample(colnames(regulonAUC), nCells) 
binAct_subset <- binaryRegulonActivity[, which(colnames(binaryRegulonActivity) %in% cellsSelected)]
dim(binAct_subset)

options(repr.plot.width=12, repr.plot.height=10) # To set the figure size in Jupyter
# binAct_subset <- binAct_subset[regulonOrder,]
ComplexHeatmap::Heatmap(binAct_subset, name="Binarized activity", 
                        col = c("white", "black"),
                        cluster_rows = TRUE, cluster_columns=TRUE,
                        show_column_names = FALSE,
                        row_names_gp=grid::gpar(fontsize=6)) # row font size



###------Binarized heatmap -----------
# percentage of cells of that cell type/cluster with the regulon active)
minPerc <- .3
binaryRegulonActivity
cellInfo_binarizedCells <- cellClusters[which(rownames(cellClusters)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$`0`), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
# TF开关为开的细胞 且比例>minPerc 
#t <-as.data.frame(rowSums(regulonActivity_byCellType_Binarized>minPerc))

binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
ComplexHeatmap::Heatmap(binaryActPerc_subset, name="Regulon activity (%)", col = c("white","pink","red"))


####--------The networks in detail: TFs, targets and motifs--------------------------
length(regulons)
sum(lengths(regulons)>=10)
viewTable(cbind(nGenes=lengths(regulons)), options=list(pageLength=10))

#Genes in a regulon
grep("CTCF", names(regulons), value=T) # paste0("^","EBF1","_")
# "CTCF(+)" 
regulons[["P21(+)"]]
gene <- "P53"

#otential TFs upstream of a given gene can be found with
names(regulons)[which(sapply(regulons, function(x) "P53" %in% x))]

dim(regulons_incidMat)

genes <- c("ITGAX","TNFRSF1B","CIB1") 
incidMat_subset <- regulons_incidMat[,genes]
incidMat_subset <- incidMat_subset[rowSums(incidMat_subset)>0,]
incidMat_subset

#Motifs supporting the regulons
motifEnrichmentFile <- file.path("reg.csv")
file.exists(motifEnrichmentFile)

motifEnrichment <- data.table::fread(motifEnrichmentFile, header=T, skip=1)[-3,]
colnames(motifEnrichment)[1:2] <- c("TF", "MotifID")


tableSubset <- motifEnrichment[TF=="BATF"]
viewMotifs(tableSubset, colsToShow = c("logo", "NES", "TF" ,"Annotation"), options=list(pageLength=5))

head(tableSubset)



####-------------IR2h------------------
library(SCopeLoomR)
library(AUCell)
library(SCENIC)

# For some of the plots:
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(Seurat)
library(htmltools)
packageVersion("SCENIC")

#1.2.4
### Initialize settings
# setwd("E:/A549")
setwd('/media/liyaru/LYR10T/Radiation_A549/TH1/A549')
scenicLoomPath <- file.path("SCENIC/res/IR2h_SCENIC.loom")
file.exists(scenicLoomPath)

loom <- open_loom(scenicLoomPath)
# Read information from loom file:
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)

# embeddings <- get_embeddings(loom)
# cellClusters <- get_clusterings(loom)

close_loom(loom)

## To check whether it was read properly:
length(regulons)
head(names(regulons))
regulonAUC
#length(regulonsAucThresholds)
#plot(embeddings$`SCENIC AUC UMAP`)


####------Scoring the network activity----------------
regulonAUC

####------Regulators for known cell types or clusters--
## To start from clusters/cell types from Scanpy:

####----------RNA-----------------------------
NC <- readRDS("scCANCER/res/int_scran_orig_cluster/2_IR_2h.rds")

# LYR edit : new meld group
m = NC@meta.data
meld = fread('/media/liyaru/LYR10T/Radiation_A549/Revision/1.meld_merge.csv')
table(meld$response_new)

hist(meld$chd_likelihood, breaks = 50, col = "grey", 
     freq = F,)
# 添加密度曲线
lines(density(meld$chd_likelihood),col= "red",lwd=2,)

table(meld$response_new,meld$class)
prop.table(table(meld$response_new,meld$class),margin = 2)

quantile(meld$chd_likelihood, probs = c(0.3333, 0.6667))

meld = meld[meld$class == 'IR_2h']
meld$X = paste0(meld$X,'_3')

# res_new = meld$response_new #new
res_new = meld$response #old
names(res_new) = meld$X
table(res_new)

NC = AddMetaData(NC,metadata = res_new,col.name = 'response_new')
table(NC$response_new)

cellClusters <- as.data.frame(NC@meta.data$response_new)
rownames(cellClusters) <- gsub("_3","",rownames(NC@meta.data))
colnames(cellClusters) <- "0"
head(cellClusters)

selectedResolution <- "0" # select resolution
# Split the cells by cluster:
cellsPerCluster <- split(rownames(cellClusters), cellClusters[,selectedResolution]) 

regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

####--------heatmap----------------------------
# plot:
options(repr.plot.width=8, repr.plot.height=10) # To set the figure size in Jupyter
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6))) # row font size
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6))) # row font size


class(regulonActivity_byCellType)
regulonActivity_byCellType <- as.data.frame(regulonActivity_byCellType)
regulonActivity_byCellType[1:60,]
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[1:60,], name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6))) # row font size


regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later


topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- as.character(topRegulators$CellType)
# topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]

top_regu <- topRegulators %>%
  group_by(CellType) %>%
  top_n(n=20,wt=RelativeActivity)

regulon = top_regu$Regulon %>% as.character()
regulon10 = as.data.frame(regulonActivity_byCellType_Scaled)[regulon,]
pheatmap(regulon10,cluster_rows = T)


####-------Cell-type specific regulators-----------------
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellClusters[colnames(regulonAUC), selectedResolution])
head(rss)
dim(rss)

topRegulators <- reshape2::melt(rss)
colnames(topRegulators) <- c("Regulon", "CellType", "RSS")

top_regu <- topRegulators %>%
  group_by(CellType) %>%
  top_n(n=20,wt=RSS)

regulon = top_regu$Regulon %>% as.character()

regulon_rss10 = as.data.frame(regulonActivity_byCellType_Scaled)[regulon,]
pheatmap(regulon_rss10,cluster_rows = F)

rss10 = as.data.frame(rss)[regulon,]
pheatmap(rss10,cluster_rows = F)


rss
#setName = "FALSE"
setName = "Hi"
n=20

fwrite(as.data.frame(rss),
       '/media/liyaru/LYR10T/Radiation_A549/Revision/2.rss_2h.csv',row.names = T)

library(ggplot2)
library(ggrepel)
rssThisType <- sort(rss[, setName], decreasing = TRUE)
thisRss <- data.frame(regulon = names(rssThisType), rank = seq_along(rssThisType), 
                      rss = rssThisType)
thisRss$regulon[(n + 1):nrow(thisRss)] <- NA
ggplot(thisRss, aes(x = rank, y = rss)) + 
  geom_point(color = "blue", size = 1) + 
  ggtitle(setName) + 
  geom_label_repel(aes(label = regulon), 
                   max.overlaps=99,
                   box.padding = 0.35, point.padding = 0.5, 
                   segment.color = "grey50", na.rm = TRUE) + 
  theme_classic()                      


####--------IR6h------------------------
library(SCopeLoomR)
library(AUCell)
library(SCENIC)

# For some of the plots:
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(Seurat)
library(htmltools)
packageVersion("SCENIC")

#1.2.4
### Initialize settings
# setwd("E:/A549")
setwd('/media/liyaru/LYR10T/Radiation_A549/TH1/A549')
scenicLoomPath <- file.path("SCENIC/res/IR6h_SCENIC.loom")
file.exists(scenicLoomPath)

loom <- open_loom(scenicLoomPath)
# Read information from loom file:
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)

# embeddings <- get_embeddings(loom)
# cellClusters <- get_clusterings(loom)

close_loom(loom)

## To check whether it was read properly:
length(regulons)
head(names(regulons))
regulonAUC
#length(regulonsAucThresholds)
#plot(embeddings$`SCENIC AUC UMAP`)


####------Scoring the network activity----------------
regulonAUC

####------Regulators for known cell types or clusters--
## To start from clusters/cell types from Scanpy:

####----------RNA-----------------------------
NC <- readRDS("scCANCER/res/int_scran_orig_cluster/2_IR_6h.rds")

cellClusters <- as.data.frame(NC@meta.data$group)
rownames(cellClusters) <- gsub("_2","",rownames(NC@meta.data))
colnames(cellClusters) <- "0"
head(cellClusters)


selectedResolution <- "0" # select resolution
# Split the cells by cluster:
cellsPerCluster <- split(rownames(cellClusters), cellClusters[,selectedResolution]) 

regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

####--------heatmap----------------------------
# plot:
options(repr.plot.width=8, repr.plot.height=10) # To set the figure size in Jupyter
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6))) # row font size
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6))) # row font size


class(regulonActivity_byCellType)
regulonActivity_byCellType <- as.data.frame(regulonActivity_byCellType)
regulonActivity_byCellType[1:60,]
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[1:60,], name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6))) # row font size


regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later


topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators)
viewTable(topRegulators, options = list(pageLength = 10))



####-------Cell-type specific regulators-----------------
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellClusters[colnames(regulonAUC), selectedResolution])
head(rss)
dim(rss)

## Showing regulons and cell types with any RSS > 0.01 
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter


# for (i in 1:10){
#   plotRSS_oneSet(rss, setName = paste0("cluster",i)) # cluster ID
# }

plotRSS_oneSet(rss, setName = "cluster12") # cluster ID


#SCENIC::plotRSS_oneSet

rss
#setName = "FALSE"
setName = "cluster12"
n=20

library(ggplot2)
library(ggrepel)
rssThisType <- sort(rss[, setName], decreasing = TRUE)
thisRss <- data.frame(regulon = names(rssThisType), rank = seq_along(rssThisType), 
                      rss = rssThisType)
thisRss$regulon[(n + 1):nrow(thisRss)] <- NA
ggplot(thisRss, aes(x = rank, y = rss)) + 
  geom_point(color = "blue", size = 1) + 
  ggtitle(setName) + 
  geom_label_repel(aes(label = regulon), 
                   max.overlaps=20,
                   box.padding = 0.35, point.padding = 0.5, 
                   segment.color = "grey50", na.rm = TRUE) + 
  theme_classic()                               


regulons$`KDM5A(+)`
regulons$`KDM5B(+)`

GRNBoost_linkList <- importArboreto(file.path("SCENIC/res/IR6h_adj.tsv"))
head(GRNBoost_linkList)

KDM5A <- GRNBoost_linkList[which(GRNBoost_linkList$TF=="KDM5A"),]
gene_KDM5A <- regulons$`KDM5A(+)`
KDM5A_net <- KDM5A[which(KDM5A$Target %in% gene_KDM5A),]
fwrite(KDM5A_net,"SCENIC/res/IR_6h/KDM5A.txt")

KDM5B <- GRNBoost_linkList[which(GRNBoost_linkList$TF=="KDM5B"),]
gene_KDM5B <- regulons$`KDM5B(+)`
KDM5B_net <- KDM5B[which(KDM5B$Target %in% gene_KDM5B),]
fwrite(KDM5B_net,"SCENIC/res/IR_6h/KDM5B.txt")

ZNF133 <- GRNBoost_linkList[which(GRNBoost_linkList$TF=="ZNF133"),]
gene_ZNF133 <- regulons$`ZNF133(+)`
ZNF133_net <- ZNF133[which(ZNF133$Target %in% gene_ZNF133),]
fwrite(ZNF133_net,"SCENIC/res/IR_6h/ZNF133.txt")

GATA6 <- GRNBoost_linkList[which(GRNBoost_linkList$TF=="GATA6"),]
gene_GATA6 <- regulons$`GATA6(+)`
GATA6_net <- GATA6[which(GATA6$Target %in% gene_GATA6),]
fwrite(GATA6_net,"SCENIC/res/IR_6h/GATA6.txt")

ZBTB7B <- GRNBoost_linkList[which(GRNBoost_linkList$TF=="ZBTB7B"),]
gene_ZBTB7B <- regulons$`ZBTB7B(+)`
ZBTB7B_net <- ZBTB7B[which(ZBTB7B$Target %in% gene_ZBTB7B),]
fwrite(ZBTB7B_net,"SCENIC/res/IR_6h/ZBTB7B.txt")







