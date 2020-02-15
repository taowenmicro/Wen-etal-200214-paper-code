#清空内存#######
rm(list=ls()) 
library("phyloseq")
library("ggplot2")
library("dada2")
library("tidyverse")
# ###制作phyloseq对象并保存
# aa = "D:/Shared_Folder/pro2_3_bac/data_processing_gg135/a9_usearch_otu_table/otu_table_tax_phylosep.txt"
# aab = import_qiime(aa)
# aab
# #导入mapping文件
# metadata <- import_qiime_sample_data("./mapping.txt")
# #colnames(tax_table(China_phyloseq))<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# tree = import_qiime(treefilename = "./data_processing_gg135/a9_tree/rep_set.tree")
# 
# ps <- merge_phyloseq(metadata, aab,tree)
# ps
# saveRDS(ps, "./ps.rds")
ps6 = readRDS("~/Desktop/DATA_get_wilt/a2_result_gg135/ps6.rds")
ps6
##########做网络筛选OTU  0.0005#########

net_table  = otu_table(ps6)
dim(net_table)
head(net_table)
psD_rhi <- ps6 %>%
  subset_taxa(
    sample_data(ps6)$SampleType %in% c("D_rhi")
  )
psD_rhi

psH_rhi <- ps6 %>%
  subset_taxa(
    sample_data(ps6)$SampleType %in% c("H_rhi")
  )
psH_rhi
D_rhi <- as.data.frame(otu_table(psD_rhi))
head(D_rhi)
psD_bulk <- ps6 %>%
  subset_taxa(
    sample_data(ps6)$SampleType %in% c("D_bulk")
  )
psD_bulk

psH_bulk <- ps6 %>%
  subset_taxa(
    sample_data(ps6)$SampleType %in% c("H_bulk")
  )
psH_bulk



net_tax  = tax_table(ps6)
head(net_tax )
net_tax = as.data.frame(net_tax)


dir.create("./result_and_script//a9_network/data")

write.table(as.data.frame(otu_table(psD_rhi)),"./result_and_script//a9_network/data/net_D_rhi.txt",sep = "\t",col.names =NA)
write.table(as.data.frame(otu_table(psH_rhi)),"./result_and_script//a9_network/data/net_H_rhi.txt",sep = "\t",col.names =NA)
write.table(as.data.frame(otu_table(psD_bulk)),"./result_and_script//a9_network/data/net_D_bulk.txt",sep = "\t",col.names =NA)
write.table(as.data.frame(otu_table(psH_bulk)),"./result_and_script//a9_network/data/net_H_bulk.txt",sep = "\t",col.names =NA)
##########做网络筛选OTU  0.0005#########





####使用稀有OTU做分析网络#################
ps1 = filter_taxa(ps, function(x) sum(x > 3) > (0.20*length(x)), TRUE);ps1
ps2  = transform_sample_counts(ps1, function(x) x / sum(x) );ps2
ps2 <- ps2 %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "mitochondria" &
      Class   != "Chloroplast"
  )

ps2
ps3 = filter_taxa(ps2, function(x) mean(x) < 0.00005, TRUE);ps3###0.001
#ps3 = filter_taxa(ps3, function(x) mean(x) > 0.00001, TRUE);ps3###0.001
net_table  = otu_table(ps1)[row.names(otu_table(ps3))]
dim(net_table)
head(net_table)

net_tax  = tax_table(ps)[row.names(otu_table(ps3))]
head(net_tax )
net_tax = as.data.frame(net_tax)

net_table <- as.data.frame(net_table)
net_B80<- select(net_table, matches("B80[0-9]"))
net_B80Z<- select(net_table, matches("B80Z"));head(net_B80Z)
net_LF <- select(net_table, matches("LF[0-9]"));head(net_LF)
net_LFZ<- select(net_table, matches("LFZ"))
dir.create("./a9_network/data")

write.table(net_B80,"./a9_network/data//net_B80_0005.txt",sep = "\t",col.names =NA)
write.table(net_B80Z ,"./a9_network/data//net_B80Z_0005.txt",sep = "\t",col.names =NA)
write.table(net_LF,"./a9_network/data//net_LF_0005.txt",sep = "\t",col.names =NA)
write.table(net_LFZ,"./a9_network/data//net_LFZ_0005.txt",sep = "\t",col.names =NA)

#############

# sh ~/Desktop/Shared_Folder/sparcc/run_sparcc.sh  -f ./data/net_LFZ_0005.txt -o ./spracc_LFZ -n 10
# sh ~/Desktop/Shared_Folder/sparcc/run_sparcc.sh  -f ./data/net_LF_0005.txt -o ./spracc_LF -n 10
# sh ~/Desktop/Shared_Folder/sparcc/run_sparcc.sh  -f ./data/net_B80Z_0005.txt -o ./spracc_B80Z -n 10
# sh ~/Desktop/Shared_Folder/sparcc/run_sparcc.sh  -f ./data/net_B80_0005.txt -o ./spracc_B80 -n 10
ps4 = filter_taxa(ps, function(x) sum(x > 3) > (0.50*length(x)), TRUE);ps4
ps8  = transform_sample_counts(ps4, function(x) log(1 + x) );ps8
###去除真核生物，古细菌，叶绿体和线粒体#######
ps5 <- ps8 %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "mitochondria" &
      Class   != "Chloroplast"
  )

ps5 = filter_taxa(ps5, function(x) mean(x) > 2 , TRUE);ps5#
pspre = ps5
pspre = ps2
library(caret)
#构建一个年龄区间列
pspre = ps3
dataMatrix <- data.frame(SampleType = sample_data(pspre)$SampleType, t(otu_table(pspre)))#合并样品分组
# 随机抽取八个分组
trainingMice <- c("B80","LF")
#按照分组提取样品
inTrain <- which(sample_data(pspre)$SampleType %in% trainingMice)
length(inTrain)#
inTrain
#mapping文件中的样品数量多于otu表格中，未进行清洗，此处加工
#inTrain=inTrain[1:221]
training <- dataMatrix[inTrain,]
dim(training)
##将因子矫正
training$SampleType = as.character(training$SampleType)
training$SampleType = as.factor(training$SampleType)

testing <- dataMatrix[-inTrain,]
testing$SampleType = as.character(testing$SampleType)
testing$SampleType = as.factor(testing$SampleType)

# 作预测建模
iris.rf = randomForest(SampleType ~ ., training , ntree=1000, 
                       nPerm=10,  proximity=TRUE, importance=TRUE)  
print(iris.rf)  
# 验证预测模型
iris.pred = predict(iris.rf,  testing )  
# 输出预测与观测对应表
table(observed=testing[,"SampleType"], predicted=iris.pred) 



