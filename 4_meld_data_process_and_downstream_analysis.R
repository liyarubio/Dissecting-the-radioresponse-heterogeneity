####output matrix for meld
library(Seurat)
###NC_IR6h====
IR6h_NC <- readRDS("scCANCER/res/NC_IR6h/SEURAT/IR6h_NC.rds")
a <- GetAssayData(IR6h_NC,assay = "RNA",slot = "count")
a <- as.data.frame(a)
a[1:5,1:5]
write.csv(a,"meld/res/expr_count.csv")
b <- IR6h_NC@meta.data
write.csv(b,"meld/code/res/meta_data.csv")


###NC_IR2h====
IR2h_NC <- readRDS("scCANCER/res/NC_IR2h/SEURAT/IR2h_NC.rds")
a <- GetAssayData(IR2h_NC,assay = "RNA",slot = "count")
a <- as.data.frame(a)
a[1:5,1:5]
write.csv(a,"meld/code/res/expr_count_IR2h_NC.csv")

b <- IR2h_NC@meta.data
write.csv(b,"meld/code/res/meta_data_IR2h_NC.csv")

####-----------downstream 2h----------------
library(Seurat)
library(magrittr)
library(dplyr)
library(stringr) 


IR_2h <- readRDS("scCANCER/res/IR_2h/IR_2h.rds")
meld_2h <- read.csv("meld/code/res/2h_metadata.csv")

A <- IR_2h@meta.data
table(colnames(IR_2h))
meld_2h <- meld_2h[which(meld_2h$orig.ident == "IR_2h"),]
meld_2h$X <- gsub("_2","",meld_2h$X)
colnames(meld_2h)
meld_2h_1 <- meld_2h[,c("X", 'chd_likelihood')]

meld_2h_1[meld_2h_1$chd_likelihood >= 0.7,'response'] <- 'Hi'
meld_2h_1[meld_2h_1$chd_likelihood <= 0.5,'response'] <- 'Low'
meld_2h_1[meld_2h_1$chd_likelihood >0.5&meld_2h_1$chd_likelihood <0.7,'response'] <- 'Mid'

plot(density(meld_2h_1$chd_likelihood), main = "Density Plot of Data", xlab = "Data", ylab = "Density")


A$X <- rownames(A)
A1 <- merge(A,meld_2h_1,by = 'X', all = TRUE)

rownames(A1) <- A1$X
colnames(A1)
A2 <- A1[,-1]
rownames(A2)
colnames(A2)

table(rownames(A1) %in% colnames(IR_2h))

IR_2h@meta.data <- A2
saveRDS(IR_2h,"meld/res/IR_2h/IR_2h.rds")

Idents(IR_2h) <- IR_2h$response
response_markers <- FindAllMarkers(IR_2h,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(response_markers,"meld/res/IR_2h/DE.csv")
colnames(response_markers)
mypal <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))
gene <- unique(response_markers$gene)

# expression matrix
exp <- GetAssayData(IR_2h, slot = "counts")
exp <- log10(exp + 1)
head(IR_2h$response)
new_celltype <- sort(IR_2h$response)
head(new_celltype)
cs_data <- as.matrix(exp[gene, names(new_celltype)])
ac=data.frame(cluster=new_celltype)
rownames(ac)=colnames(cs_data)
library(pheatmap)
pheatmap(cs_data,show_colnames =F,show_rownames = T,
         cluster_rows = F,
         cluster_cols = F,
         annotation_col=ac,
         border_color = "black")


library(ComplexHeatmap)
library(paletteer) 
color = paletteer_d("ggsci::nrc_npg")[c(1,3,4)]
names(color) <- levels(new_celltype)
top_anno <- HeatmapAnnotation(cluster=anno_block(gp=gpar(fill=color),
                                                 labels = levels(new_celltype),
                                                 labels_gp = gpar(cex=0.5,color='white',fontsize = 18)))

Heatmap(cs_data,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = T,
        column_split = new_celltype,
        top_annotation = top_anno, 
        column_title = NULL,
        heatmap_legend_param = list(
          title='Expression',
          title_position='leftcenter-rot'),
        col = colorRampPalette(c("white","#66CCFF","#333366"))(100),
        border = "black",
        row_names_gp = gpar(fontsize = 8))

DimPlot(IR_2h)

####-----------downstream 6h----------------
library(Seurat)
library(magrittr)
library(dplyr)
library(stringr) 


IR_6h <- readRDS("scCANCER/res/IR_6h/IR_6h.rds")
meld_6h <- read.csv("meld/code/res/6h_metadata.csv")

A <- IR_6h@meta.data
table(colnames(IR_6h))
meld_6h <- meld_6h[which(meld_6h$orig.ident == "IR_6h"),]
meld_6h$X <- gsub("_3","",meld_6h$X)
colnames(meld_6h)
meld_6h_1 <- meld_6h[,c("X", 'chd_likelihood')]

meld_6h_1[meld_6h_1$chd_likelihood >= 0.7,'response'] <- 'Hi'
meld_6h_1[meld_6h_1$chd_likelihood <= 0.5,'response'] <- 'Low'
meld_6h_1[meld_6h_1$chd_likelihood >0.5&meld_6h_1$chd_likelihood <0.7,'response'] <- 'Mid'

plot(density(meld_6h_1$chd_likelihood), main = "Density Plot of Data", xlab = "Data", ylab = "Density")

# LYR edit ： meld distribution
meld_6h_1$class = 'IR_6h'
meld_2h_1$class = 'IR_2h'

meld_merge = rbind(meld_2h_1,meld_6h_1)

q = quantile(meld_merge$chd_likelihood)

meld_merge[meld_merge$chd_likelihood >= q[4],'response_new'] <- 'High'
meld_merge[meld_merge$chd_likelihood <= q[2],'response_new'] <- 'Low'
meld_merge[meld_merge$chd_likelihood >q[2] & meld_merge$chd_likelihood < q[4],'response_new'] <- 'Mid'
table(meld_merge$response_new)
fwrite(meld_merge,'/media/liyaru/LYR10T/Radiation_A549/Revision/1.meld_merge.csv')

pdf('/media/liyaru/LYR10T/Radiation_A549/Revision/1.meld_distribution.pdf',width = 5,height = 4)
ggplot(meld_merge, aes(x = chd_likelihood))+ 
  geom_density(aes(fill = class,alpha=0.4))+
  scale_fill_manual(values = c(rgb(141,180,226,maxColorValue = 255),rgb(253,214,181,maxColorValue =255)))+
  theme_bw()
dev.off()

ggplot(meld_merge, aes(x = chd_likelihood))+ 
  geom_density()

summary(meld_merge$chd_likelihood)

A$X <- rownames(A)
A1 <- merge(A,meld_6h_1,by = 'X', all = TRUE)

rownames(A1) <- A1$X
colnames(A1)
A2 <- A1[,-1]
rownames(A2)
colnames(A2)

table(rownames(A1) %in% colnames(IR_6h))

IR_6h@meta.data <- A2
saveRDS(IR_6h,"meld/res/IR_6h/IR_6h.rds")
Idents(IR_6h) <- IR_6h$response
response_markers <- FindAllMarkers(IR_6h,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(response_markers,"meld/res/IR_6h/DE.csv")
colnames(response_markers)
mypal <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))
# gene <- unique(response_markers$gene)
gene1 <- read.csv("meld/res/IR_6h/IR6h_hi_lo_up.csv")
gene <- gene1[which(gene1$avg_log2FC>0.5),]$X

# expression matrix
IR_6h1 <- subset(IR_6h, idents = "Mid", invert = TRUE)
exp <- GetAssayData(IR_6h1, slot = "counts")
exp <- log10(exp + 1)
head(IR_6h$response)
new_celltype <- sort(IR_6h1$response)
head(new_celltype)
cs_data <- as.matrix(exp[gene, names(new_celltype)])
ac=data.frame(cluster=new_celltype)
rownames(ac)=colnames(cs_data)

library(ComplexHeatmap)
library(paletteer) 
color = paletteer_d("ggsci::nrc_npg")[c(1,3)]
names(color) <- levels(new_celltype)
# top_anno <- HeatmapAnnotation(cluster=anno_block(gp=gpar(fill=color),
#                                                  labels = levels(new_celltype),
#                                                  labels_gp = gpar(cex=0.5,color='white',fontsize = 18)))

Heatmap(cs_data,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = T,
        column_split = new_celltype,
        top_annotation = top_anno, 
        column_title = NULL,
        heatmap_legend_param = list(
          title='Expression',
          title_position='leftcenter-rot'),
        col = colorRampPalette(c("white","#66CCFF","#333366"))(100),
        border = "black",
        row_names_gp = gpar(fontsize = 8))

DimPlot(IR_6h)