####06.18
####scCANCER for scRNA-seq analysis of A549NB
library(scCancer)
####NC====
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

####IR_2h====
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

####IR_6h====
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
Integration_methods <-  c("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")
for (i in 1:6){ 
  # The paths of all sample's "runScAnnotation" results
  single.savePaths <- c("scCANCER/res/NC/", "scCANCER/res/IR_2h/","scCANCER/res/IR_6h/")
  sampleNames <- c("NC", "IR_2h", "IR_6h")    # The labels for all samples
  savePath <- paste0("scCANCER/res/comnination/",Integration_methods[i],"/")       # A path to save the results
  combName <- "A549_com"                 # A label of the combined samples
  authorName <- "Chen-Lab@TH"                # The author name to mark the report
  comb.method <- Integration_methods[i]               # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")
  
  # Run scCombination
  comb.results <- runScCombination(
    single.savePaths = single.savePaths, 
    sampleNames = sampleNames, 
    savePath = savePath, 
    combName = combName,
    authorName = authorName,
    comb.method = comb.method
  )
}

####IR2h_NC====
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

####IR6h_NC====
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

