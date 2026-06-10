#####download TCGA clinical information
library(TCGAbiolinks)
tmp<-getGDCprojects()
tmp$project_id
clinical <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")  
write.csv(clinical,"RSI/data/clinical_information/clinical_LUAD.csv",row.names = TRUE)

table(clinical$treatments_radiation_treatment_or_therapy)
A <- clinical[which(clinical$treatments_radiation_treatment_or_therapy == "yes"),]
sample <- A$bcr_patient_barcode

### Load data
exp <- read.table("RSI/data/convert_exp_FPKM.txt", header = TRUE)
rownames(exp) <- exp$Tag
a <- gsub("-",".",sample)
a %in% substr(colnames(exp),1,12)
a_50 <- a[50]
a_61 <- a[61]
b <- a[-c(50,61)]
b %in% substr(colnames(exp),1,12)
exp1 <- exp[,substr(colnames(exp),1,12) %in% b]  
c=as.numeric(substr(colnames(exp1),14,15))
table(c)
colnames(exp1)
exp2 <- exp1[,-(108:121)]
dim(exp2)
b %in% substr(colnames(exp2),1,12)
###107 radiation patients

####gene expression analysis
exp4 <- as.data.frame((exp2))
exp4 <- log(exp4+1)
dim(exp4)
head(colnames(exp4))
dim(meta)
head(rownames(meta))
rownames(meta) <- gsub("-",".",rownames(meta))

exp6 <- as.data.frame(t(exp4))
exp6$RSI <- -0.0098009*exp6$AR+0.0128283*exp6$JUN+0.0254552*exp6$STAT1-0.0017589*exp6$PRKCB-0.0038171*exp6$RELA+0.1070213*exp6$ABL1-0.0002509*exp6$SUMO1-0.0092431*exp6$CDK1-0.0204469*exp6$HDAC1-0.0441683*exp6$IRF1
median(exp6[["RSI"]])
exp6$Radisensitivity <- ifelse(exp6$RSI<median(exp6[["RSI"]]),'Radiosensitive','Radioresistant') 
table(exp6$Radisensitivity)
my_comparisons <- list(c("Radiosensitive","Radioresistant"))
genes <- c("KDM4B","CDKN1A")
for (i in genes){
  ggdotplot(exp6, x = "Radisensitivity",y = i, add = "boxplot", add.params = list(color = "Radisensitivity"), title = i, binwidth = 0.1, color = "Radisensitivity", fill = "Radisensitivity",
            palette = c("#00AFBB", "#E7B800"))+ 
    stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")+
    theme(text=element_text(size=15))
  ggsave(paste0("RSI/res3/expression2/",i,"_1.pdf"))
  ggsave(paste0("RSI/res3/expression2/",i,"_1.png"))
}


