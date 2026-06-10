library(dyno)
library(tidyverse)
library(Matrix)
library(Seurat)
library(dynwrap)
library(sampling)

####-----------NC--------------------
###prepare data====
##NC
NC <- readRDS("scCANCER/res/NC/NC.rds")
NC <- as.Seurat(NC)
colnames(NC@meta.data)
NC$group <- paste0("cluster",as.numeric(NC$seurat_clusters))
NC_meta <- NC@meta.data[,c("orig.ident", "nCount_RNA","nFeature_RNA","group","cellcycle")]
NC@meta.data <- NC_meta
NC$group_orig <- NC$group

NC_2 <- RenameCells(NC, new.names=paste0(colnames(NC),"_1"))
rownames(NC_2@meta.data)
IR6h_NC <- readRDS("scCANCER/res/NC_IR6h/IR6h_NC.rds")
Idents(IR6h_NC) <- IR6h_NC$orig.ident
NC_1 <- subset(IR6h_NC, idents = "NC")
rownames(NC_1)
NC_1@meta.data <- NC_2@meta.data
rownames(NC_1)
NC <- NC_1
saveRDS(NC,"scCANCER/res/int_scran_orig_cluster/NC.rds")
# NC <- readRDS("scCANCER/res/int_scran_orig_cluster/NC.rds")
DefaultAssay(NC) <- 'RNA'
all.genes <- rownames(NC)
NC <- NC %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(features = all.genes) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  FindNeighbors(dims = 1:30) %>%
  #FindClusters(resolution=c(0.3,0.5,0.8,1.0,1.2)) %>%
  RunTSNE(dims = 1:30) %>% RunUMAP(dims = 1:30)
saveRDS(NC, "scCANCER/res/int_scran_orig_cluster/2_NC.rds")
NC <- readRDS("scCANCER/res/int_scran_orig_cluster/2_NC.rds")

###downsample====
meta <- NC@meta.data
table(meta$cellcycle)
t <- sampling::strata(meta, stratanames='cellcycle', size= c(640,120,240), method="srswor", pik,description=FALSE)
data <- getdata(meta,t)
table(data$cellcycle)
NC_downsample <- subset(NC,cells = rownames(data)) 
saveRDS(NC_downsample, "Velocity/res/NC/NC_downsample.rds")
NC <- NC_downsample
NC <- readRDS("Velocity/res/NC/NC_downsample.rds")

NC_meta <- NC@meta.data
NC_meta$cellcycle <- factor(NC@meta.data$cellcycle, levels = c("G1","S","G2M"))
# NC_meta[order(NC@meta.data$cellcycle),]
NC_meta
# unique(NC_meta[c(641:760),]$cellcycle)
NC_meta <- NC_meta[c(1:640,761:1000,641:760),]
NC@meta.data <- NC_meta
###prepare data====  
dataset_NC <- wrap_expression(
  counts = t(NC@assays$RNA@counts),
  expression = t(NC@assays$RNA@data))

t <- rownames(NC)
# counts = t(NC@assays$RNA@counts)
# dim(counts)


#add prior information
start_id <- rownames(NC_meta[which(NC_meta$cellcycle == "G1"),])
end_id <- rownames(NC_meta[which(NC_meta$cellcycle == "G2M"),])
dataset_NC <- add_prior_information(dataset_NC, start_id = start_id)

#adding group information
table(NC$cellcycle)
# 
# a <- c(dataset_NC$grouping[271:910],dataset_NC$grouping[1:270],dataset_NC$grouping[911:1000])
# 


dataset_NC <- add_grouping(
  dataset_NC,
  NC$cellcycle
)

table(dataset_NC$grouping)

dataset_NC <- dynwrap::add_dimred(dataset_NC,
  NC@reductions$umap@cell.embeddings)


# dataset_NC1 <- dataset_NC
# ###Selecting the best methods for a dataset_NC====
# guidelines_shiny(dataset_NC = dataset_NC)
# 
# answers <- dynguidelines::answer_questions(
#   multiple_disconnected = FALSE,
#   expect_topology = NULL,  
#   expected_topology = NULL,
#   n_cells = 15000,
#   n_features = 2500,
#   time = "20h",
#   memory = "100GB",
#   prior_information = NULL,
#   docker = TRUE
# )
# guidelines <- guidelines_shiny(answers = answers)
# 
# methods_selected <- guidelines$methods_selected
# methods_selected <- "slingshot"
# method <- get_ti_methods(as_tibble = FALSE)[[1]]$fun

###Inferring trajectories====
set.seed(1)#To make the execution of a method reproducible, fix the seed
# model_NC_slingshot <- infer_trajectory(dataset_NC, method = "slingshot")
# model_NC_angle <- infer_trajectory(dataset_NC, method = "angle")
# model_NC_paga <- infer_trajectory(dataset_NC, method = "paga")
model_NC_paga_tree <- infer_trajectory(dataset_NC, method = "paga_tree")

# # model_NC_slingshot <- infer_trajectory(dataset_NC, method = "slingshot")
# plot_dimred(
#   model_NC_slingshot,
#   color_density = "grouping",
#   expression_source = dataset_NC$expression, 
#   grouping = dataset_NC$grouping,
#   plot_milestone_network = TRUE,
# )
# ggsave("Velocity/res/NC/slingshot_combined_sample_cellcycle.pdf")
# 
# plot_dimred(
#   model_NC_paga,
#   color_density = "grouping",
#   expression_source = dataset_NC$expression, 
#   grouping = dataset_NC$grouping,
#   plot_milestone_network = TRUE,
# )
# ggsave("Velocity/res/NC/paga_combined_sample_cellcycle.pdf")
# 
# model_NC_paga_tree <- infer_trajectory(dataset_NC, method = "paga_tree")

plot_dimred(
  model_NC_paga_tree,
  color_density = "grouping",
  expression_source = dataset_NC$expression, 
  grouping = dataset_NC$grouping,
  plot_milestone_network = TRUE,
)

ggsave("Velocity/res/NC/paga_tree_cellcycle.pdf")

# model_rooted <- model_NC_paga_tree %>% add_root(root_milestone_id = c("7","4"))

plot_dimred(
  model_NC_paga_tree,
  color_density = "grouping",
  expression_source = dataset_NC$expression, 
  grouping = group_onto_nearest_milestones(model_NC_paga_tree),
  plot_milestone_network = TRUE
)

model_NC_paga_tree <- add_root(
  model_NC_paga_tree,
  root_milestone_id = "6",
  flip_edges = TRUE
)

plot_dimred(
  model_NC_paga_tree,
  color_density = "grouping",
  expression_source = dataset_NC$expression, 
  grouping = dataset_NC$grouping,
  plot_milestone_network = TRUE
)
ggsave("Velocity/res/NC/paga_tree_cellcycle_adapted.pdf")

###heatmap====
plot_heatmap(model_NC_angle, expression_source = dataset_NC)
ggsave("Velocity/res/NC/angle_heatmap.pdf")
plot_heatmap(model_NC_slingshot, expression_source = dataset_NC)
ggsave("Velocity/res/NC/slingshot_heatmap.pdf")
plot_heatmap(model_NC_paga_tree, expression_source = dataset_NC)
ggsave("Velocity/res/NC/paga_tree_heatmap.pdf")
plot_heatmap(model_NC_paga, expression_source = dataset_NC)
ggsave("Velocity/res/NC/paga_heatmap.pdf")




# plot_dimred(
#   model_NC_paga_tree,
#   color_density = "grouping",
#   expression_source = dataset_NC$expression, 
#   grouping = dataset_NC$grouping,
#   plot_milestone_network = TRUE,
# )
# 
# plot_dimred(
#   model_NC_paga_tree,
#   color_density = "grouping",
#   expression_source = dataset_NC$expression, 
#   grouping = group_onto_nearest_milestones(model_NC_paga_tree),
#   plot_milestone_network = TRUE,
# )

# plot_dendro(model_NC_paga_tree,grouping = dataset_NC$grouping)

# disconnected_dataset_NC <- dyntoy::generate_trajectory(model = "dataset_NC", num_cells = 1000)
# plot_graph(disconnected_dataset_NC)

###save model====
saveRDS(model_NC_angle,"Velocity/res/NC/model_NC_angle.rds")
saveRDS(model_NC_slingshot,"Velocity/res/NC/model_NC_slingshot.rds")
saveRDS(model_NC_paga,"Velocity/res/NC/model_NC_paga.rds")
saveRDS(model_NC_paga_tree,"Velocity/res/NC/model_NC_paga_tree.rds")


####------ IR 2h-------------------
library(dyno)
library(tidyverse)
library(Matrix)
library(Seurat)
library(dynwrap)
library(sampling)

###prepare data====
##IR_2h
IR_2h <- readRDS("scCANCER/res/IR_2h/IR_2h.rds")
IR_2h <- as.Seurat(IR_2h)
colnames(IR_2h@meta.data)
IR_2h$group <- paste0("cluster",as.numeric(IR_2h$seurat_clusters))
IR_2h_meta <- IR_2h@meta.data[,c("orig.ident", "IR_2hount_RNA","nFeature_RNA","group","cellcycle")]
IR_2h@meta.data <- IR_2h_meta
IR_2h$group_orig <- IR_2h$group

IR_2h_2 <- RenameCells(IR_2h, new.names=paste0(colnames(IR_2h),"_1"))
rownames(IR_2h_2@meta.data)
IR6h_IR_2h <- readRDS("scCANCER/res/IR_2h_IR6h/IR6h_IR_2h.rds")
Idents(IR6h_IR_2h) <- IR6h_IR_2h$orig.ident
IR_2h_1 <- subset(IR6h_IR_2h, idents = "IR_2h")
rownames(IR_2h_1)
IR_2h_1@meta.data <- IR_2h_2@meta.data
rownames(IR_2h_1)
IR_2h <- IR_2h_1
saveRDS(IR_2h,"scCANCER/res/int_scran_orig_cluster/IR_2h.rds")

DefaultAssay(IR_2h) <- 'RNA'
all.genes <- rownames(IR_2h)
IR_2h <- IR_2h %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(features = all.genes) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  FindNeighbors(dims = 1:30) %>%
  #FindClusters(resolution=c(0.3,0.5,0.8,1.0,1.2)) %>%
  RunTSNE(dims = 1:30) %>% RunUMAP(dims = 1:30)
saveRDS(IR_2h, "scCANCER/res/int_scran_orig_cluster/2_IR_2h.rds")
IR_2h <- readRDS("scCANCER/res/int_scran_orig_cluster/2_IR_2h.rds")
###downsaple====
meta <- IR_2h@meta.data
table(meta$cellcycle)
t <- sampling::strata(meta, stratanames='cellcycle', size= c(270,640,90), method="srswor", pik,description=FALSE)
data <- getdata(meta,t)
table(data$cellcycle)
IR_2h_downsample <- subset(IR_2h,cells = rownames(data)) 
saveRDS(IR_2h_downsample, "Velocity/res/IR_2h/IR_2h_downsample.rds")
IR_2h <- IR_2h_downsample
table(IR_2h$cellcycle)

IR_2h <- readRDS("Velocity/res/IR_2h/IR_2h_downsample.rds")
IR_2h$cellcycle <- factor(IR_2h$cellcycle, levels = c("G1","S","G2M"))
IR_2h_meta <- IR_2h@meta.data

# IR_2h_meta[order(IR_2h@meta.data$cellcycle),]
IR_2h_meta
unique(IR_2h_meta[c(641:760),]$cellcycle)

IR_2h_meta <- IR_2h_meta[c(271:910,1:270,911:1000),]
IR_2h@meta.data <- IR_2h_meta
###prepare data====  
dataset_IR_2h <- wrap_expression(
  counts = t(IR_2h@assays$RNA@counts),
  expression = t(IR_2h@assays$RNA@data))

t <- rownames(IR_2h)
# counts = t(IR_2h@assays$RNA@counts)
# dim(counts)


#add prior information
start_id <- rownames(IR_2h@meta.data[which(IR_2h@meta.data$cellcycle == "G1"),])
end_id <- rownames(IR_2h@meta.data[which(IR_2h@meta.data$cellcycle == "G2M"),])
dataset_IR_2h <- add_prior_information(dataset_IR_2h, start_id = start_id)

#adding group information
IR_2h$cellcycle

dataset_IR_2h <- add_grouping(
  dataset_IR_2h,
  IR_2h$cellcycle
)

dataset_IR_2h <- dynwrap::add_dimred(dataset_IR_2h,
                                     IR_2h@reductions$umap@cell.embeddings)


# dataset_IR_2h$group_ids <- c("G1","S","G2M")

# dataset_IR_2h1 <- dataset_IR_2h
# ###Selecting the best methods for a dataset_IR_2h====
# guidelines_shiny(dataset_IR_2h = dataset_IR_2h)
# 
# answers <- dynguidelines::answer_questions(
#   multiple_disconnected = FALSE,
#   expect_topology = NULL,
#   expected_topology = NULL,
#   n_cells = 15000,
#   n_features = 2500,
#   time = "20h",
#   memory = "100GB",
#   prior_information = NULL,
#   docker = TRUE
# )
# guidelines <- guidelines_shiny(answers = answers)
# 
# methods_selected <- guidelines$methods_selected
# methods_selected <- "slingshot"
# method <- get_ti_methods(as_tibble = FALSE)[[1]]$fun

###Inferring trajectories====
set.seed(1)#To make the execution of a method reproducible, fix the seed
# model_IR_2h_slingshot <- infer_trajectory(dataset_IR_2h, method = "slingshot")
# model_IR_2h_angle <- infer_trajectory(dataset_IR_2h, method = "angle")
# model_IR_2h_paga <- infer_trajectory(dataset_IR_2h, method = "paga")
model_IR_2h_paga_tree <- infer_trajectory(dataset_IR_2h, method = "paga_tree")

# model_IR_2h_slingshot <- infer_trajectory(dataset_IR_2h, method = "slingshot")
table(dataset_IR_2h$grouping)

IR_2h$cellcycle <- a

class(a)
A <- dynplot:::add_cell_coloring(dataset_IR_2h$cell_info, color_cells = "grouping", grouping = dataset_IR_2h$grouping, trajectory = model_NC_paga_tree)

plot_dimred(
  model_IR_2h_slingshot,
  color_density = "grouping",
  color_cells = ,
  expression_source = dataset_IR_2h$expression, 
  plot_milestone_network = TRUE,
)

ggsave("Velocity/res/IR_2h/slingshot_combined_cellcycle.pdf")

plot_dimred(
  model_IR_2h_paga,
  color_density = "grouping",
  expression_source = dataset_IR_2h$expression, 
  grouping = dataset_IR_2h$grouping,
  plot_milestone_network = TRUE
)
ggsave("Velocity/res/IR_2h/paga_combined_cellcycle.pdf")


# model_IR_2h_paga_tree <- infer_trajectory(dataset_IR_2h, method = "paga_tree")

plot_dimred(
  model_IR_2h_paga_tree,
  color_density = "grouping",
  # color_cells = "grouping",
  expression_source = dataset_IR_2h$expression, 
  grouping = dataset_IR_2h$grouping,
  plot_milestone_network = TRUE
)


ggsave("Velocity/res/IR_2h/paga_tree_cellcycle.pdf")

###heatmap====
plot_heatmap(model_IR_2h_angle, expression_source = dataset_IR_2h)
ggsave("Velocity/res/IR_2h/angle_heatmap.pdf")
plot_heatmap(model_IR_2h_slingshot, expression_source = dataset_IR_2h)
ggsave("Velocity/res/IR_2h/slingshot_heatmap.pdf")
plot_heatmap(model_IR_2h_paga_tree, expression_source = dataset_IR_2h)
ggsave("Velocity/res/IR_2h/paga_tree_heatmap.pdf")
plot_heatmap(model_IR_2h_paga, expression_source = dataset_IR_2h)
ggsave("Velocity/res/IR_2h/paga_heatmap.pdf")




# plot_dimred(
#   model_IR_2h_paga_tree,
#   color_density = "grouping",
#   expression_source = dataset_IR_2h$expression, 
#   grouping = dataset_IR_2h$grouping,
#   plot_milestone_network = TRUE,
# )
# 
# plot_dimred(
#   model_IR_2h_paga_tree,
#   color_density = "grouping",
#   expression_source = dataset_IR_2h$expression, 
#   grouping = group_onto_nearest_milestones(model_IR_2h_paga_tree),
#   plot_milestone_network = TRUE,
# )

# plot_dendro(model_IR_2h_paga_tree,grouping = dataset_IR_2h$grouping)

# disconnected_dataset_IR_2h <- dyntoy::generate_trajectory(model = "dataset_IR_2h", num_cells = 1000)
# plot_graph(disconnected_dataset_IR_2h)

###save model====
saveRDS(model_IR_2h_angle,"Velocity/res/IR_2h/model_IR_2h_angle.rds")
saveRDS(model_IR_2h_slingshot,"Velocity/res/IR_2h/model_IR_2h_slingshot.rds")
saveRDS(model_IR_2h_paga,"Velocity/res/IR_2h/model_IR_2h_paga.rds")
saveRDS(model_IR_2h_paga_tree,"Velocity/res/IR_2h/model_IR_2h_paga_tree.rds")

dynplot::plot_edge_flips(model_NC_paga_tree,model_IR_2h_paga_tree)



####------- IR 6h--------------------
library(dyno)
library(tidyverse)
library(Matrix)
library(Seurat)
library(dynwrap)
library(sampling)

###prepare data====
##IR_6h
IR_6h <- readRDS("scCANCER/res/IR_6h/IR_6h.rds")
IR_6h <- as.Seurat(IR_6h)
colnames(IR_6h@meta.data)
IR_6h$group <- paste0("cluster",as.numeric(IR_6h$seurat_clusters))
IR_6h_meta <- IR_6h@meta.data[,c("orig.ident", "IR_6hount_RNA","nFeature_RNA","group","cellcycle")]
IR_6h@meta.data <- IR_6h_meta
IR_6h$group_orig <- IR_6h$group

IR_6h_2 <- RenameCells(IR_6h, new.names=paste0(colnames(IR_6h),"_1"))
rownames(IR_6h_2@meta.data)
IR6h_IR_6h <- readRDS("scCANCER/res/IR_6h_IR6h/IR6h_IR_6h.rds")
Idents(IR6h_IR_6h) <- IR6h_IR_6h$orig.ident
IR_6h_1 <- subset(IR6h_IR_6h, idents = "IR_6h")
rownames(IR_6h_1)
IR_6h_1@meta.data <- IR_6h_2@meta.data
rownames(IR_6h_1)
IR_6h <- IR_6h_1
saveRDS(IR_6h,"scCANCER/res/int_scran_orig_cluster/IR_6h.rds")
IR_6h <- readRDS("scCANCER/res/int_scran_orig_cluster/IR_6h.rds")
DefaultAssay(IR_6h) <- 'RNA'
all.genes <- rownames(IR_6h)
IR_6h <- IR_6h %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(features = all.genes) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  FindNeighbors(dims = 1:30) %>%
  #FindClusters(resolution=c(0.3,0.5,0.8,1.0,1.2)) %>%
  RunTSNE(dims = 1:30) %>% RunUMAP(dims = 1:30)
saveRDS(IR_6h, "scCANCER/res/int_scran_orig_cluster/2_IR_6h.rds")
IR_6h <- readRDS("scCANCER/res/int_scran_orig_cluster/2_IR_6h.rds")
###downsaple====
meta <- IR_6h@meta.data
table(meta$cellcycle)
t <- sampling::strata(meta, stratanames='cellcycle', size= c(270,640,90), method="srswor", pik,description=FALSE)
data <- getdata(meta,t)
table(data$cellcycle)
IR_6h_downsample <- subset(IR_6h,cells = rownames(data)) 
saveRDS(IR_6h_downsample, "Velocity/res/IR_6h/IR_6h_downsample.rds")
IR_6h <- IR_6h_downsample
table(IR_6h$cellcycle)

IR_6h <- readRDS("Velocity/res/IR_6h/IR_6h_downsample.rds")
IR_6h$cellcycle <- factor(IR_6h$cellcycle, levels = c("G1","S","G2M"))
IR_6h_meta <- IR_6h@meta.data

# IR_6h_meta[order(IR_6h@meta.data$cellcycle),]
IR_6h_meta
unique(IR_6h_meta[c(641:760),]$cellcycle)

IR_6h_meta <- IR_6h_meta[c(271:910,1:270,911:1000),]
IR_6h@meta.data <- IR_6h_meta
###prepare data====  
dataset_IR_6h <- wrap_expression(
  counts = t(IR_6h@assays$RNA@counts),
  expression = t(IR_6h@assays$RNA@data))

t <- rownames(IR_6h)
# counts = t(IR_6h@assays$RNA@counts)
# dim(counts)


#add prior information
start_id <- rownames(IR_6h@meta.data[which(IR_6h@meta.data$cellcycle == "G1"),])
end_id <- rownames(IR_6h@meta.data[which(IR_6h@meta.data$cellcycle == "G2M"),])
dataset_IR_6h <- add_prior_information(dataset_IR_6h, start_id = start_id)

#adding group information
IR_6h$cellcycle

dataset_IR_6h <- add_grouping(
  dataset_IR_6h,
  IR_6h$cellcycle
)

dataset_IR_6h <- dynwrap::add_dimred(dataset_IR_6h,
                                     IR_6h@reductions$umap@cell.embeddings)


# dataset_IR_6h$group_ids <- c("G1","S","G2M")

# dataset_IR_6h1 <- dataset_IR_6h
# ###Selecting the best methods for a dataset_IR_6h====
# guidelines_shiny(dataset_IR_6h = dataset_IR_6h)
# 
# answers <- dynguidelines::answer_questions(
#   multiple_disconnected = FALSE,
#   expect_topology = NULL,
#   expected_topology = NULL,
#   n_cells = 15000,
#   n_features = 2500,
#   time = "20h",
#   memory = "100GB",
#   prior_information = NULL,
#   docker = TRUE
# )
# guidelines <- guidelines_shiny(answers = answers)
# 
# methods_selected <- guidelines$methods_selected
# methods_selected <- "slingshot"
# method <- get_ti_methods(as_tibble = FALSE)[[1]]$fun

###Inferring trajectories====
set.seed(1)#To make the execution of a method reproducible, fix the seed
# model_IR_6h_slingshot <- infer_trajectory(dataset_IR_6h, method = "slingshot")
# model_IR_6h_angle <- infer_trajectory(dataset_IR_6h, method = "angle")
# model_IR_6h_paga <- infer_trajectory(dataset_IR_6h, method = "paga")
model_IR_6h_paga_tree <- infer_trajectory(dataset_IR_6h, method = "paga_tree")

# model_IR_6h_slingshot <- infer_trajectory(dataset_IR_6h, method = "slingshot")
table(dataset_IR_6h$grouping)

IR_6h$cellcycle <- a

class(a)
A <- dynplot:::add_cell_coloring(dataset_IR_6h$cell_info, color_cells = "grouping", grouping = dataset_IR_6h$grouping, trajectory = model_NC_paga_tree)

plot_dimred(
  model_IR_6h_slingshot,
  color_density = "grouping",
  color_cells = ,
  expression_source = dataset_IR_6h$expression, 
  plot_milestone_network = TRUE,
)

ggsave("Velocity/res/IR_6h/slingshot_combined_cellcycle.pdf")

plot_dimred(
  model_IR_6h_paga,
  color_density = "grouping",
  expression_source = dataset_IR_6h$expression, 
  grouping = dataset_IR_6h$grouping,
  plot_milestone_network = TRUE
)
ggsave("Velocity/res/IR_6h/paga_combined_cellcycle.pdf")


# model_IR_6h_paga_tree <- infer_trajectory(dataset_IR_6h, method = "paga_tree")

plot_dimred(
  model_IR_6h_paga_tree,
  color_density = "grouping",
  # color_cells = "grouping",
  expression_source = dataset_IR_6h$expression, 
  grouping = dataset_IR_6h$grouping,
  plot_milestone_network = TRUE
)


ggsave("Velocity/res/IR_6h/paga_tree_cellcycle.pdf")

###heatmap====
plot_heatmap(model_IR_6h_angle, expression_source = dataset_IR_6h)
ggsave("Velocity/res/IR_6h/angle_heatmap.pdf")
plot_heatmap(model_IR_6h_slingshot, expression_source = dataset_IR_6h)
ggsave("Velocity/res/IR_6h/slingshot_heatmap.pdf")
plot_heatmap(model_IR_6h_paga_tree, expression_source = dataset_IR_6h)
ggsave("Velocity/res/IR_6h/paga_tree_heatmap.pdf")
plot_heatmap(model_IR_6h_paga, expression_source = dataset_IR_6h)
ggsave("Velocity/res/IR_6h/paga_heatmap.pdf")




# plot_dimred(
#   model_IR_6h_paga_tree,
#   color_density = "grouping",
#   expression_source = dataset_IR_6h$expression, 
#   grouping = dataset_IR_6h$grouping,
#   plot_milestone_network = TRUE,
# )
# 
# plot_dimred(
#   model_IR_6h_paga_tree,
#   color_density = "grouping",
#   expression_source = dataset_IR_6h$expression, 
#   grouping = group_onto_nearest_milestones(model_IR_6h_paga_tree),
#   plot_milestone_network = TRUE,
# )

# plot_dendro(model_IR_6h_paga_tree,grouping = dataset_IR_6h$grouping)

# disconnected_dataset_IR_6h <- dyntoy::generate_trajectory(model = "dataset_IR_6h", num_cells = 1000)
# plot_graph(disconnected_dataset_IR_6h)

###save model====
saveRDS(model_IR_6h_angle,"Velocity/res/IR_6h/model_IR_6h_angle.rds")
saveRDS(model_IR_6h_slingshot,"Velocity/res/IR_6h/model_IR_6h_slingshot.rds")
saveRDS(model_IR_6h_paga,"Velocity/res/IR_6h/model_IR_6h_paga.rds")
saveRDS(model_IR_6h_paga_tree,"Velocity/res/IR_6h/model_IR_6h_paga_tree.rds")

dynplot::plot_edge_flips(model_NC_paga_tree,model_IR_6h_paga_tree)










