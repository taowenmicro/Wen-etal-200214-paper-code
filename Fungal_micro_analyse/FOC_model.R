

ps7= ps_add_out 
ps7

ps7 <- ps7 %>%
  subset_taxa(
    #Kingdom == "Bacteria" &
    # Genus  == "Fusarium"
    #Species == "Fusarium_oxysporum"
  )
ps7

ps7 <- ps7 %>%
  subset_taxa(
    #Kingdom == "Bacteria" &
    # Genus  == "Fusarium"
    Species == "Fusarium_oxysporum"
  )
ps7

ps7 <- ps7 %>%
  subset_taxa(
    #Kingdom == "Bacteria" &
    # Genus  == "Fusarium"
    Species %in%c("Fusarium_oxysporum","Fusarium_keratoplasticum") 
  )

ps7 <- ps7 %>%
  subset_taxa(
    #Kingdom == "Bacteria" &
    # Genus  == "Fusarium"
    Species %in%c("Fusarium_keratoplasticum") 
  )



####进行phyloseq下游分析
# ps7 <- prune_samples(sample_sums(ps7) >=500,ps7);ps7
# ps7
ps7 = filter_taxa(ps7, function(x) sum(x ) > 0 , TRUE);ps7
ps7  = transform_sample_counts(ps7, function(x) x / sum(x) );ps7



mapping = as.data.frame(sample_data(ps7))
mapping$SampleType

table(mapping$SampleType)
table(mapping$Description)
vegan_otu <-  function(physeq){
  OTU <-  otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <-  t(OTU)
  }
  return(as(OTU,"matrix"))
}
otutab = as.data.frame(t(vegan_otu(ps7)))
dim(otutab)
head(otutab)
otutab_t[is.na(otutab_t)] = 0
mapping7 = as.data.frame(sample_data(ps7))
mapping7$SampleType
library(randomForest)

# otutab need transposition for randomForest function
otutab_t = as.data.frame(t(otutab))
otutab_t[is.na(otutab_t)] = 0

# Set classification info.
otutab_t$group = factor(mapping7$SampleType,levels= c("D_bulk","H_bulk"))
colnames(otutab_t) = paste("OTU",colnames(otutab_t),sep = "")
# set random seed for reproducible
set.seed(315)
set.seed(23)
set.seed(21)
set.seed(2)
set.seed(100)
# RandomForest Classification
model_add= randomForest(OTUgroup ~ ., data=otutab_t, importance=TRUE, proximity=TRUE)

print(model_add)

#预测项目owe
model = model_add
# model = model_orig
ps6 = readRDS("./add_sample_for_predict/owe//a9_usearch_otu_table/ps_owe.rds")
ps6
map = as.data.frame(sample_data(ps6))
map
map$SampleType = paste(map$SampleTypeDH,map$zone,sep = "_")
sample_data(ps6) = map

result = predict_rand(ps7,ps6,model = model)
result[[2]]
result[[1]]

##No217

ps6 = readRDS("./add_sample_for_predict/No217//a9_usearch_otu_table/ps_NCBI2.rds")
ps6
map = as.data.frame(sample_data(ps6))
map
map$SampleType = paste(map$SampleTypeDH,map$zone,sep = "_")
sample_data(ps6) = map

result = predict_rand(ps7,ps6,model = model)
result[[2]]
result[[1]]


###No218：这个研究采样自香蕉20厘米的根际，意外的是将发病和健康预测相反
ps6 = readRDS("./add_sample_for_predict/No218//a9_usearch_otu_table/ps_NCBI2.rds")
ps6
map = as.data.frame(sample_data(ps6))
map
map$SampleType = paste(map$SampleTypeDH,map$zone,sep = "_")
sample_data(ps6) = map

result = predict_rand(ps7,ps6,model = model)
result[[2]]
result[[1]]


###No25：刘红军师兄样本,土体大部分预测为发病

ps6 = readRDS("./add_sample_for_predict/No_25//a9_usearch_otu_table/ps_NCBI2.rds")
ps6
map = as.data.frame(sample_data(ps6))
map
map$SampleType = paste(map$SampleTypeDH,map$zone,sep = "_")
sample_data(ps6) = map

result = predict_rand(ps7,ps6,model = model)
result[[2]]
result[[1]]


###No224
ps6 = readRDS("./add_sample_for_predict/No224//a9_usearch_otu_table/ps_NCBI2.rds")
ps6
map = as.data.frame(sample_data(ps6))
map
map$SampleType = paste(map$SampleTypeDH,map$zone,sep = "_")
sample_data(ps6) = map

result = predict_rand(ps7,ps6,model = model)
result[[2]]
result[[1]]

###liuhongjun
ps6 = readRDS("./add_sample_for_predict/liuhongjun//a9_usearch_otu_table/ps_NCBI2.rds")
ps6
map = as.data.frame(sample_data(ps6))
map
map$SampleType = paste(map$SampleTypeDH,map$zone,sep = "_")
sample_data(ps6) = map

result = predict_rand(ps7,ps6,model = model)
result[[2]]
result[[1]]