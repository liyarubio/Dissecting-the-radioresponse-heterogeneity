###SCENIC downstream analysis
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

# new meld group
m = NC@meta.data
meld = fread('/media/liyaru/LYR10T/Radiation_A549/Revision/1.meld_merge.csv')
table(meld$response_new)

hist(meld$chd_likelihood, breaks = 50, col = "grey", 
     freq = F,)

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


# ####-------Cell-type specific regulators-----------------
# rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellClusters[colnames(regulonAUC), selectedResolution])
# head(rss)
# dim(rss)
# 
# topRegulators <- reshape2::melt(rss)
# colnames(topRegulators) <- c("Regulon", "CellType", "RSS")
# 
# top_regu <- topRegulators %>%
#   group_by(CellType) %>%
#   top_n(n=20,wt=RSS)
# 
# regulon = top_regu$Regulon %>% as.character()
# 
# regulon_rss10 = as.data.frame(regulonActivity_byCellType_Scaled)[regulon,]
# pheatmap(regulon_rss10,cluster_rows = F)
# 
# rss10 = as.data.frame(rss)[regulon,]
# pheatmap(rss10,cluster_rows = F)
# 
# 
# rss
# #setName = "FALSE"
# setName = "Hi"
# n=20
# 
# fwrite(as.data.frame(rss),
#        '/media/liyaru/LYR10T/Radiation_A549/Revision/2.rss_2h.csv',row.names = T)
# 
# library(ggplot2)
# library(ggrepel)
# rssThisType <- sort(rss[, setName], decreasing = TRUE)
# thisRss <- data.frame(regulon = names(rssThisType), rank = seq_along(rssThisType), 
#                       rss = rssThisType)
# thisRss$regulon[(n + 1):nrow(thisRss)] <- NA
# ggplot(thisRss, aes(x = rank, y = rss)) + 
#   geom_point(color = "blue", size = 1) + 
#   ggtitle(setName) + 
#   geom_label_repel(aes(label = regulon), 
#                    max.overlaps=99,
#                    box.padding = 0.35, point.padding = 0.5, 
#                    segment.color = "grey50", na.rm = TRUE) + 
#   theme_classic()                      