# #在当前工作目录下读取单个xml文件用于test：
# library(XML)#用于读取xml文件的R包
# options(stringsAsFactors = F)
# xmltest <- xmlParse("RSI/data/clinical_information/TCGA-LUAD/00e01b05-d939-49e6-b808-e1bab0d6773b.BCR XML")
# 
# #形成根目录列表数据
# xmltop <- xmlRoot(xmltest)
# class(xmltop) #查看类
# xmlName(xmltop) #查看根目录名
# xmlSize(xmltop) #查看根目录存在节点数
# 
# xmltop[1]
# xmltop[2] #查看节点的详细信息
# 
# CHOL_cl1_df <- xmlToDataFrame(xmltop [2])#把xml文件中Root[2]包含内容转化为数据框
# CHOL_cl1_df[,1:6]
# head(t(CHOL_cl1_df))#t()转置；head(t())查看转置后的前六行
# 
# #单个样本test无误后，开始导入全部xml
# CHOL_cl_xmls <- dir("RSI/data/clinical_information/TCGA-LUAD/",
#                     pattern = "*XML",
#                     recursive = T)
# CHOL_cl_xmls[1:4] #查看list的前四列，522个样本（522列）已全部读取
# 
# #自定义一个把每一个xml文件转换为dataframe的函数
# cldf <- function(x){
#   xmlresult <- xmlParse(file.path("RSI/data/clinical_information/TCGA-LUAD/",x))
#   xmltop2 <- xmlRoot(xmlresult)
#   CHOLcldf <- xmlToDataFrame(xmltop2[2])
#   return(t(CHOLcldf))
# }
# 
# #test一个xml文件看自定义函数是否成功
# head(cldf("00e01b05-d939-49e6-b808-e1bab0d6773b.BCR XML"))
# 
# #函数test无误开始临床数据合并
# cl <- lapply(CHOL_cl_xmls,cldf) #把所有的xml读取进来
# CHOL_cl <- t(do.call(cbind,cl)) #数据合并
# CHOL_cl_df <- data.frame(CHOL_cl)#转置后是矩阵，需要恢复成数据框
# #输出整理后的clinical文件(csv)
# write.csv(CHOL_cl_df,"RSI/data/clinical_information/clinical.csv",row.names = TRUE)
# 
# #查看包含哪些临床信息
# colnames(CHOL_cl_df)
# #分离做过辐射处理的患者
# A <- CHOL_cl_df[which(CHOL_cl_df$radiations != ""),]
# write.csv(A,"RSI/data/clinical_information/clinical_radiation.csv",row.names = TRUE)
# table(A$bcr_patient_barcode)



#####利用R工具下载TCGA临床信息
library(TCGAbiolinks)
tmp<-getGDCprojects()
tmp$project_id
clinical <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")  
write.csv(clinical,"RSI/data/clinical_information/clinical_LUAD.csv",row.names = TRUE)
# library(DT)
# datatable(clinical, filter = 'top', 
#           options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
#           rownames = FALSE)
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
###以上共获得107个radiation patients样本

rownames(exp2)

exp3 <- as.data.frame(t(exp2))
exp3 <- log(exp3+1)

# gene_sig <- c('JUN','HDAC1','RELA','PRKCB','SUMO1','ABL1','STAT1','AR','PAK2','IRF1')    
exp3$RSI <- -0.0098009*exp3$AR+0.0128283*exp3$JUN+0.0254552*exp3$STAT1-0.0017589*exp3$PRKCB-0.0038171*exp3$RELA+0.1070213*exp3$ABL1-0.0002509*exp3$SUMO1-0.0092431*exp3$PAK2-0.0204469*exp3$HDAC1-0.0441683*exp3$IRF1
median(exp3[["RSI"]])

test <- c("CDKN1A","MDM2","PTTG1","FDXR","DKK1","TP53I3","PHLDA3","ALDH3A1","TOB1","NINJ1","ELF3","IGFL2-AS1","KDM4B","NEAT1",
          "BBC3","CYFIP2","CDKN3","ARL6IP1","PURPL","GAS6-AS1","PFKFB3","CCND1","NUPR1","AVPI1","IFI16","BAX","ZMAT3","NRAV")
genes <- NULL  
# test <- c("KDM4B","CCND1")
for (i in test){
  targetgene <- i
  median(exp3[[targetgene]])
  high <- exp3[which(exp3[[targetgene]] >= median(exp3[[targetgene]])),]
  low <- exp3[which(exp3[[targetgene]] < median(exp3[[targetgene]])),]
  high$classify <- "High"
  low$classify <- "Low"
  
  exp_plot <- rbind(high,low)
  exp_plot$classify 
  
  my_comparisons <- list(c("High","Low"))
  
  ggdotplot(exp_plot, x = "classify",y = "RSI",add = "mean_sd", add.params = list(color = "black"), fill = "classify", binwidth = 0.005, color = "classify",
            palette = c("#00AFBB", "#E7B800"))+ 
    stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")+
    theme(text=element_text(size=16))
  ggsave(paste0("RSI/res3/cluster12/",targetgene,".pdf"), height = 4,width = 5)
  ggsave(paste0("RSI/res3/cluster12/",targetgene,".png"),height = 4,width = 5)
}


###生存分析
# 加载R包
library(survival)
library(survminer)

data(lung)
View(lung)
# 整理数据
targetgene <- "CDKN1A"
median(exp3[[targetgene]])
high <- exp3[which(exp3[[targetgene]] >= median(exp3[[targetgene]])),]
low <- exp3[which(exp3[[targetgene]] < median(exp3[[targetgene]])),]
high$classify <- "High"
low$classify <- "Low"
rownames(high)
B <- A[-c(50,61),]
b %in% gsub("-",".",B$bcr_patient_barcode)
B$CDKN1A <- "Low"
B[gsub("-",".",B$bcr_patient_barcode) %in% substr(rownames(high),1,12),"CDKN1A"]="High"
table(B$CDKN1A)
rownames(B) <- B$bcr_patient_barcode
B$ajcc_pathologic_stage
meta=B[,c("vital_status","days_to_last_follow_up",
                                 "days_to_death",
                                 "race",
                                 "gender",
                                 "age_at_index",
                                 "ajcc_pathologic_stage")]

exp4 <- as.data.frame((exp2))
exp4 <- log(exp4+1)
dim(exp4)
head(colnames(exp4))
dim(meta)
head(rownames(meta))
rownames(meta) <- gsub("-",".",rownames(meta))

#1、计算生存时间
meta$days_to_death[is.na(meta$days_to_death)] <- 0   #缺失值标记为0
meta$days_to_last_follow_up[is.na(meta$days_to_last_follow_up)] <- 0
meta$days=as.numeric(meta[,2])+as.numeric(meta[,3])

#时间以月份记，保留两位小数
meta$time=round(meta$days/30,2)

#2、根据生死定义活着是0，死的是1
meta$event=ifelse(meta$vital_status=='Alive',0,1)
table(meta$event)

#3 年龄分组(部分样本缺失，考虑可能的影响应该不大)
meta$age_at_index[is.na(meta$age_at_index)] <- 0
meta$age_at_index=as.numeric(meta$age_at_index)
meta$age_group=ifelse(meta$age_at_index>median(meta$age_at_index),'older','younger')
table(meta$age_group)

#4 癌症阶段
table(meta$ajcc_pathologic_stage)

#5 race 人种
table(meta$race)

#6 性别 gender
table(meta$gender)

save(exp4,meta,file="RSI/data/tosur.RData")

#利用这两个文件接下来就可以开始生存分析了

setwd('/media/liyaru/LYR10T/Radiation_A549/TH1/A549')

rm(list=ls())
load("RSI/data/tosur.RData")
library(survival)
library(survminer)

sfit <- survfit(Surv(time, event)~gender, data=meta)
print(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)

ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)

#对性别和年龄生存分析拼图
sfit1=survfit(Surv(time, event)~gender, data=meta)
sfit2=survfit(Surv(time, event)~age_group, data=meta)
splots <- list()
splots[[1]] <- ggsurvplot(sfit1,pval =TRUE, data = meta, risk.table = TRUE)
splots[[2]] <- ggsurvplot(sfit2,pval =TRUE, data = meta, risk.table = TRUE)
arrange_ggsurvplots(splots, print = TRUE,  ncol = 2, nrow = 1, risk.table.height = 0.4)

#基因表达的生存分析
# (1) 直接指定感兴趣基因，可以是1个或多个
exprSet=exp4  #套流程
test <- c("CDKN1A","MDM2","PTTG1","FDXR","DKK1","TP53I3","PHLDA3","ALDH3A1","TOB1","NINJ1","ELF3","IGFL2-AS1","KDM4B","NEAT1",
          "BBC3","CYFIP2","CDKN3","ARL6IP1","PURPL","GAS6-AS1","PFKFB3","CCND1","NUPR1","AVPI1","IFI16","BAX","ZMAT3","NRAV")
colnames(exprSet) <- substr(colnames(exprSet),1,12)

  g = test[20] # 随便选一个
  meta$gene = ifelse(as.numeric(exprSet[g,])>median(as.numeric(exprSet[g,])),'high','low')
  sfit1=survfit(Surv(time, event)~gene, data=meta)
  ggsurvplot(sfit1,palette = c("#E7B800", "#2E9FDF"),
             risk.table =TRUE,pval =TRUE,
             conf.int =TRUE,xlab ="Time in months", 
             ggtheme =theme_light(), 
             ncensor.plot = TRUE)

#文章图
  g=c("NINJ1")
  meta$NINJ1 = ifelse(as.numeric(exprSet[g,])>median(as.numeric(exprSet[g,])),'high','low')
  sfit <- survfit(Surv(time, event)~NINJ1, data=meta)
  print(sfit)
  ggsurvplot(sfit, conf.int=F, pval=TRUE)
  
  ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
             risk.table =TRUE,pval =TRUE,
             conf.int =TRUE,xlab ="Time in months", 
             ggtheme =theme_light(), 
             ncensor.plot = TRUE)

  # g=c("IGFL2-AS1")
  # meta$gene = ifelse(as.numeric(exprSet[g,])>median(as.numeric(exprSet[g,])),'high','low')
  # sfit <- survfit(Surv(time, event)~gene, data=meta)
  # print(sfit)
  # ggsurvplot(sfit, conf.int=F, pval=TRUE)
  # 
  # ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
  #            risk.table =TRUE,pval =TRUE,
  #            conf.int =TRUE,xlab ="Time in months", 
  #            ggtheme =theme_light(), 
  #            ncensor.plot = TRUE)




####基因差异表达分析
###辐射与未辐射
exp1 <- exp[,-(517:575)]
colnames(exp1)
rownames(exp1) <- exp1$Tag
exp2 <- exp1[,-1]
colnames(exp1)
##转置表格
exp3 <- as.data.frame(t(exp2))
exp3 <- log(exp3+1)
exp5 <- exp3
rownames(exp5)
exp5$Radiation <- ifelse(substr(rownames(exp5),1,12) %in% b, 'YES','NO') 
my_comparisons <- list(c("YES","NO"))

for (i in test){
  ggboxplot(exp5, x = "Radiation",y = i, add = "mean_sd", add.params = list(color = "black"), fill = "Radiation", binwidth = 0.1, color = "Radiation",
          palette = c("#00AFBB", "#E7B800"))+ 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")+
  theme(text=element_text(size=15))
 ggsave(paste0("RSI/res3/expression/",i,".pdf"))
 ggsave(paste0("RSI/res3/expression/",i,".png"))
}

###高敏感低敏感
#总数据
exp5$RSI <- -0.0098009*exp5$AR+0.0128283*exp5$JUN+0.0254552*exp5$STAT1-0.0017589*exp5$PRKCB-0.0038171*exp5$RELA+0.1070213*exp5$ABL1-0.0002509*exp5$SUMO1-0.0092431*exp5$CDK1-0.0204469*exp5$HDAC1-0.0441683*exp5$IRF1
median(exp5[["RSI"]])
exp5$Radisensitivity <- ifelse(exp5$RSI<median(exp5[["RSI"]]),'high','low') 
table(exp5$Radisensitivity)
my_comparisons <- list(c("high","low"))
for (i in test){
  ggboxplot(exp5, x = "Radisensitivity",y = i, add = "mean_sd", add.params = list(color = "black"), binwidth = 0.1, color = "Radisensitivity",
            palette = c("#00AFBB", "#E7B800"))+ 
    stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")+
    theme(text=element_text(size=15))
  ggsave(paste0("RSI/res3/expression1/",i,".pdf"))
  ggsave(paste0("RSI/res3/expression1/",i,".png"))
}

#辐射数据
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


group_list <- ifelse(exp5$RSI>median(exp5[["RSI"]]),'high','low') 
group_list <- factor(group_list, levels = c('high','low'))
#DESeq2方法做差异分析
library(edgeR)
colData <- data.frame(row.names = colnames(exp4),condition=group_list)
dge <- DGEList(counts=exp2,group=group_list)
dge$samples$lib.size <- colSums(dge$counts)
dge <- calcNormFactors(dge)
design <- model.matrix(~0+group_list)
row.names(design)<- colnames(dge)
colnames(design)<-levels(group_list)
dge <- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge,design)
dge <- estimateGLMTagwiseDisp(dge,design)
fit <- glmFit(dge,design)
fit2 <- glmLRT(fit, contrast = c(-1,1))

DEG2 <- topTags(fit2, n=row(exp2))
DEG2 <- as.data.frame(DEG2)
logFC_cutoff2 <- with(DEG2,mean(abs(logFC)) + 2*sd(abs(logFC)))
DEG2$change = as.factor(
  ifelse(DEG2$PValue < 0.05 & abs(DEG2$logFC) > logFC_cutoff2,
         ifelse(DEG2$logFC > logFC_cutoff2, "UP", "DOWN"),"NOT")
)
head(DEG2)


dds <-DESeqDataSetFromMatrix(
  countData = exp4,
  colData = colData,
  design = ~ condition)

