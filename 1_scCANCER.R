####scCANCER for scRNA-seq analysis
library(scCancer)
###NC====
dataPath <- "scRNA_file/A549NB/"
savePath <- "scCANCER/res/NC/"
statPath <- "scCANCER/res/NC/"
sampleName <- "NC"
authorName <- "Chen-Lab@TH"

# Run scStatistics
stat.results <- runScStatistics(
  dataPath = dataPath,
  savePath = savePath,
  sampleName = sampleName,
  authorName = authorName
)

# Run scAnnotation
anno.results <- runScAnnotation(
  dataPath = dataPath,
  statPath = statPath,
  savePath = savePath,
  authorName = authorName,
  sampleName = sampleName,
  geneSet.method = "average"
)

###IR_2h====
dataPath <- "scRNA_file/A549CC/"
savePath <- "scCANCER/res/IR_2h/"
statPath <- "scCANCER/res/IR_2h/"
sampleName <- "IR_2h"
authorName <- "Chen-Lab@TH"

# Run scStatistics
stat.results <- runScStatistics(
  dataPath = dataPath,
  savePath = savePath,
  sampleName = sampleName,
  authorName = authorName
)

# Run scAnnotation
anno.results <- runScAnnotation(
  dataPath = dataPath,
  statPath = statPath,
  savePath = savePath,
  authorName = authorName,
  sampleName = sampleName,
  geneSet.method = "average"
)

###IR_6h====
dataPath <- "scRNA_file/A549CD/"
savePath <- "scCANCER/res/IR_6h/"
statPath <- "scCANCER/res/IR_6h/"
sampleName <- "IR_6h"
authorName <- "Chen-Lab@TH"

# Run scStatistics
stat.results <- runScStatistics(
  dataPath = dataPath,
  savePath = savePath,
  sampleName = sampleName,
  authorName = authorName
)

# Run scAnnotation
anno.results <- runScAnnotation(
  dataPath = dataPath,
  statPath = statPath,
  savePath = savePath,
  authorName = authorName,
  sampleName = sampleName,
  geneSet.method = "average"
)

####combination by different integration methods====
###IR2h_NC====
# The paths of all sample's "runScAnnotation" results
single.savePaths <- c("scCANCER/res/NC/", "scCANCER/res/IR_2h/")
sampleNames <- c("NC", "IR_2h")    # The labels for all samples
savePath <- paste0("scCANCER/res/NC_IR2h/")       # A path to save the results
combName <- "NC_IR2h"                 # A label of the combined samples
authorName <- "Chen-Lab@TH"                # The author name to mark the report
comb.method <- "SeuratMNN"               # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")

# Run scCombination
comb.results <- runScCombination(
  single.savePaths = single.savePaths, 
  sampleNames = sampleNames, 
  savePath = savePath, 
  combName = combName,
  authorName = authorName,
  comb.method = comb.method)

###IR6h_NC====
# The paths of all sample's "runScAnnotation" results
single.savePaths <- c("scCANCER/res/NC/", "scCANCER/res/IR_6h/")
sampleNames <- c("NC", "IR_6h")    # The labels for all samples
savePath <- paste0("scCANCER/res/NC_IR6h/")       # A path to save the results
combName <- "NC_IR6h"                 # A label of the combined samples
authorName <- "Chen-Lab@TH"                # The author name to mark the report
comb.method <- "SeuratMNN"               # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")

# Run scCombination
comb.results <- runScCombination(
  single.savePaths = single.savePaths, 
  sampleNames = sampleNames, 
  savePath = savePath, 
  combName = combName,
  authorName = authorName,
  comb.method = comb.method)

###DEG====
seu <- readRDS('A549/scCANCER/res/NC_IR6h/IR6h_NC.rds')
DimPlot(seu,reduction = "tsne",split.by = "orig.ident",group.by ="integrated_snn_res.0.8" )

Idents(seu) <- seu$orig.ident
seu <- subset(seu,ident='IR_6h')

DefaultAssay(seu) <- 'RNA'
seu@meta.data$group = paste0('cluster',
                             (as.integer(seu@meta.data$integrated_snn_res.0.8)))

m = FindMarkers(seu,group.by = "group",ident.1 = "cluster5",ident.2 = "cluster12",
                logfc.threshold=0)
m$gene <- rownames(m)
fwrite(m,'DEG/C5_C12_DEG_IR6h.csv',row.names = T)

m <- fread('DEG/C5_C12_DEG_IR6h.csv')
# volcano plot
{
  Dat<- as.data.frame(m)
  Dat$gene <- Dat$V1
  fc = 0.5

  Dat$threshold <- 'Other'
  Dat[Dat$p_val_adj < 0.05 & Dat$avg_log2FC > fc,'threshold'] <- 'Up'
  Dat[Dat$p_val_adj < 0.05 & Dat$avg_log2FC < -fc,'threshold'] <- 'Down'
  table(Dat$threshold)
  Dat$threshold <- factor(Dat$threshold,levels = c('Up','Down','Other'))
  
  p <- ggplot(Dat,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
    geom_point()+
    scale_color_manual(values=c("#CD3333","#000080", "#808080"))+
    geom_text_repel(
      data = Dat[(Dat$p_val_adj<0.05 & abs(Dat$avg_log2FC) > 1),],
      aes(label = gene),
      size = 3,
      segment.color = "black", show.legend = FALSE,
      max.overlaps=40)+
    theme_bw()+
    theme(
      legend.title = element_blank(),
    )+
    ylab('-log10 (p-adj)')+
    xlab('log2 (FoldChange)')+
    geom_vline(xintercept=c(-0.5,0.5),lty=3,col="black",lwd=0.5) +
    geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)
  
}

c5 <- m[m$p_val_adj < 0.05 & m$avg_log2FC > 0.5,]
fwrite(c5,'DEG/C5_up_genes_vs_C12_IR6h.csv')
