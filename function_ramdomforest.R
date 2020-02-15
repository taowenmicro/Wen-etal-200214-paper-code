

### 构造一个函数
#合并ps对象
#有两种模式，model =1 是取共有部分进行合并
#model =2 是保留x中的全部OTU，y中如果有其他OTU，全部舍弃，如果没有x中的这些otu使用0填充。
#注意进化树可不能合并，必须重新运算，这里没有必要花费很长的时间来运算这个。
#model = 3,xy全部OTU进行合并

merge_ps = function(ps1,ps2,model = 1){
  library(randomForest)
  library(stats)
  vegan_otu <-  function(physeq){
    OTU <-  otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU <-  t(OTU)
    }
    return(as(OTU,"matrix"))
  }

  vegan_tax <-  function(physeq){
    tax <-  tax_table(physeq)

    return(as(tax,"matrix"))
  }


  #需要分三次合并
  #otu_table
  if(model == 1){
    otu_table1 = as.data.frame(t(vegan_otu(ps1)))
    dim(otu_table1)
    otu_table2 = as.data.frame(t(vegan_otu(ps2)))
    dim(otu_table2)
    otu_table3 = merge(otu_table1,otu_table2,by = "row.names",all = TRUE)
    dim(otu_table3)
    row.names(otu_table3) = otu_table3$Row.names
    otu_table3$Row.names = NULL
    otu_table3 = as.matrix(otu_table3)
    otu_table3[is.na(otu_table3)] <- 0
    #tax_table

    tax_table1 = as.data.frame((vegan_tax(ps1)))
    dim(tax_table1)
    tax_table2 = as.data.frame((vegan_tax(ps2)))
    dim(tax_table2)
    head(tax_table1)
    head(tax_table2)
    tax_table3 = rbind(tax_table1,tax_table2)
    dim(tax_table3)

  }
  if(model == 2){
    otu_table1 = as.data.frame(t(vegan_otu(ps1)))
    dim(otu_table1)
    otu_table2 = as.data.frame(t(vegan_otu(ps2)))
    dim(otu_table2)
    otu_table3 = merge(otu_table1,otu_table2,by = "row.names",all.x = TRUE)
    dim(otu_table3)
    row.names(otu_table3) = otu_table3$Row.names
    otu_table3$Row.names = NULL
    otu_table3 = as.matrix(otu_table3)
    otu_table3[is.na(otu_table3)] <- 0
    #tax_table

    tax_table1 = as.data.frame((vegan_tax(ps1)))
    dim(tax_table1)
    tax_table3 = tax_table1
  }
  if(model == 3){
    otu_table1 = as.data.frame(t(vegan_otu(ps1)))
    dim(otu_table1)
    otu_table2 = as.data.frame(t(vegan_otu(ps2)))
    dim(otu_table2)
    otu_table3 = merge(otu_table1,otu_table2,by = "row.names",all = TRUE)
    dim(otu_table3)
    row.names(otu_table3) = otu_table3$Row.names
    otu_table3$Row.names = NULL
    otu_table3 = as.matrix(otu_table3)
    otu_table3[is.na(otu_table3)] <- 0
    #tax_table
    tax_table1 = as.data.frame((vegan_tax(ps1)))
    dim(tax_table1)
    tax_table2 = as.data.frame((vegan_tax(ps2)))
    dim(tax_table2)
    head(tax_table1)
    head(tax_table2)
    tax_table3 = rbind(tax_table1,tax_table2)
    dim(tax_table3)
  }


  mapping1 = as.data.frame(sample_data(ps1))
  head(mapping1)
  mapping1 = mapping1[,c("SampleType","Description")]
  mapping2 = as.data.frame(sample_data(ps2))
  head(mapping2)
  # mapping2$BarcodeSequence = NULL
  # mapping2$LinkerPrimerSequence = NULL
  mapping2$SampleType = paste(mapping2$fianl_SampleType,mapping2$zone,sep = "_")
  mapping2 = mapping2[,c("SampleType","Description")]
  mapping3 = rbind(mapping1,mapping2)
  head(mapping3)


  ps_add_out =phyloseq(otu_table(as.matrix(otu_table3),taxa_are_rows = TRUE),
                       sample_data(mapping3),
                       tax_table(as.matrix(tax_table3)))
  return(ps_add_out)
}

merge_ps_tax = function(ps1,ps2,model = 1){
  library(randomForest)
  library(stats)
  vegan_otu <-  function(physeq){
    OTU <-  otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU <-  t(OTU)
    }
    return(as(OTU,"matrix"))
  }
  
  vegan_tax <-  function(physeq){
    tax <-  tax_table(physeq)
    
    return(as(tax,"matrix"))
  }
  
  
  #需要分三次合并
  #otu_table
  if(model == 1){
    otu_table1 = as.data.frame(t(vegan_otu(ps1)))
    dim(otu_table1)
    otu_table2 = as.data.frame(t(vegan_otu(ps2)))
    dim(otu_table2)
    otu_table3 = merge(otu_table1,otu_table2,by = "row.names",all = TRUE)
    dim(otu_table3)
    row.names(otu_table3) = otu_table3$Row.names
    otu_table3$Row.names = NULL
    otu_table3 = as.matrix(otu_table3)
    otu_table3[is.na(otu_table3)] <- 0
    #tax_table
    
    tax_table1 = as.data.frame((vegan_tax(ps1)))
    dim(tax_table1)
    tax_table2 = as.data.frame((vegan_tax(ps2)))
    dim(tax_table2)
    head(tax_table1)
    head(tax_table2)
    tax_table3 = rbind(tax_table1,tax_table2)
    dim(tax_table3)
    
  }
  if(model == 2){
    otu_table1 = as.data.frame(t(vegan_otu(ps1)))
    dim(otu_table1)
    otu_table2 = as.data.frame(t(vegan_otu(ps2)))
    dim(otu_table2)
    otu_table3 = merge(otu_table1,otu_table2,by = "row.names",all.x = TRUE)
    dim(otu_table3)
    row.names(otu_table3) = otu_table3$Row.names
    otu_table3$Row.names = NULL
    otu_table3 = as.matrix(otu_table3)
    otu_table3[is.na(otu_table3)] <- 0
    #tax_table
    
    tax_table1 = as.data.frame((vegan_tax(ps1)))
    dim(tax_table1)
    tax_table3 = tax_table1
  }
  if(model == 3){
    otu_table1 = as.data.frame(t(vegan_otu(ps1)))
    dim(otu_table1)
    otu_table2 = as.data.frame(t(vegan_otu(ps2)))
    dim(otu_table2)
    otu_table3 = merge(otu_table1,otu_table2,by = "row.names",all = TRUE)
    dim(otu_table3)
    row.names(otu_table3) = otu_table3$Row.names
    otu_table3$Row.names = NULL
    otu_table3 = as.matrix(otu_table3)
    otu_table3[is.na(otu_table3)] <- 0
    #tax_table
    tax_table1 = as.data.frame((vegan_tax(ps1)))
    dim(tax_table1)
    tax_table2 = as.data.frame((vegan_tax(ps2)))
    dim(tax_table2)
    head(tax_table1)
    head(tax_table2)
    tax_table3 = rbind(tax_table1,tax_table2)
    dim(tax_table3)
  }
  
  
  mapping1 = as.data.frame(sample_data(ps1))
  head(mapping1)
  mapping1 = mapping1[,c("SampleType","Description")]
  mapping2 = as.data.frame(sample_data(ps2))
  head(mapping2)
  # mapping2$BarcodeSequence = NULL
  # mapping2$LinkerPrimerSequence = NULL
  mapping2$SampleType = paste(mapping2$fianl_SampleType,mapping2$zone,sep = "_")
  mapping2 = mapping2[,c("SampleType","Description")]
  mapping3 = rbind(mapping1,mapping2)
  head(mapping3)
  
  
  ps_add_out =phyloseq(otu_table(as.matrix(otu_table3),taxa_are_rows = TRUE),
                       sample_data(mapping3),
                       tax_table(as.matrix(tax_table3)))
  ps_add_out = transform_sample_counts(ps_add_out, function(x) x / sum(x) )
  return(ps_add_out)
}


# ps_cs = merge_ps(ps1 = ps7,ps2 = ps6,model = 2)

predict_rand <- function(ps7,ps6,model) {


  ps6  = transform_sample_counts(ps6, function(x) x / sum(x) );ps6
  vegan_otu <-  function(physeq){
    OTU <-  otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU <-  t(OTU)
    }
    return(as(OTU,"matrix"))
  }

  vegan_tax <-  function(physeq){
    tax <-  tax_table(physeq)

    return(as(tax,"matrix"))
  }
  library(tidyverse)
  psX <- ps6 %>%
    subset_taxa(
      row.names(tax_table(ps6)) %in% row.names(tax_table(ps7))
    )
  psX

  otu  = as.data.frame(t(vegan_otu(psX)))

  head(otu)

  otu_add = matrix(0,nrow = length(setdiff(row.names(tax_table(ps7)), row.names(tax_table(ps6)))),
                   ncol = dim(otu)[2])
  dim(otu_add)
  otu_add = as.data.frame(otu_add)
  row.names(otu_add) = setdiff(row.names(tax_table(ps7)), row.names(tax_table(ps6)))
  colnames(otu_add) = colnames(otu)
  dim(otu)
  otu_fil = rbind(otu,otu_add)
  dim(otu_fil)
  ps8 <- psX


  otutab = otu_fil
  dim(otutab)
  mapping8 = as.data.frame(sample_data(ps8))
  # mapping8 $SampleType = mapping8$fianl_SampleType
  # sample_data(ps6)=  mapping8
  # saveRDS(ps6,"./add_sample_for_predict/No_166//ps_NCBI10.rds")
  # install.packages("randomForest")
  library(randomForest)

  # otutab need transposition for randomForest function
  otutab_t = as.data.frame(t(otutab))

  # as.factor(mapping8$SampleType)
  # Set classification info.
  otutab_t$group = factor(mapping8$SampleType)
  # otutab_t$group = factor(mapping8$SampleType)
  colnames(otutab_t) = paste("OTU",colnames(otutab_t),sep = "")

  otutab.pred = stats::predict(model, otutab_t )
  pre_tab = table(observed=otutab_t[,"OTUgroup"], predicted=otutab.pred)
  pre_tab
  # save prediction result
  predict = data.frame(group = otutab_t[,"OTUgroup"], predicted=otutab.pred)



  summary(predict)

  return(list(otutab.pred,pre_tab,summary(predict),otutab_t))
}


predict_rand_tax <- function(ps7,ps6,model,taxGlomRank = "Genus") {



  vegan_otu <-  function(physeq){
    OTU <-  otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU <-  t(OTU)
    }
    return(as(OTU,"matrix"))
  }

  vegan_tax <-  function(physeq){
    tax <-  tax_table(physeq)

    return(as(tax,"matrix"))
  }
  library(tidyverse)
  psX <- ps6 %>%
    subset_taxa(
      row.names(tax_table(ps6)) %in% row.names(tax_table(ps7))
    )
  psX

  otu  = as.data.frame(t(vegan_otu(psX)))

  head(otu)

  otu_add = matrix(0,nrow = length(setdiff(row.names(tax_table(ps7)), row.names(tax_table(ps6)))),
                   ncol = dim(otu)[2])
  dim(otu_add)
  otu_add = as.data.frame(otu_add)
  row.names(otu_add) = setdiff(row.names(tax_table(ps7)), row.names(tax_table(ps6)))
  colnames(otu_add) = colnames(otu)
  dim(otu)
  otu_fil = rbind(otu,otu_add)
  dim(otu_fil)
  ps8 <- psX
  
  ps8 = tax_glom(ps8, taxrank = taxGlomRank)


  otutab = otu_fil
  dim(otutab)
  mapping8 = as.data.frame(sample_data(ps8))
  # mapping8 $SampleType = mapping8$fianl_SampleType
  # sample_data(ps6)=  mapping8
  # saveRDS(ps6,"./add_sample_for_predict/No_166//ps_NCBI10.rds")
  # install.packages("randomForest")
  library(randomForest)

  # otutab need transposition for randomForest function
  otutab_t = as.data.frame(t(otutab))

  # as.factor(mapping8$SampleType)
  # Set classification info.
  otutab_t$group = factor(mapping8$SampleType)
  # otutab_t$group = factor(mapping8$SampleType)
  colnames(otutab_t) = paste("OTU",colnames(otutab_t),sep = "")

  otutab.pred = stats::predict(model, otutab_t )
  pre_tab = table(observed=otutab_t[,"OTUgroup"], predicted=otutab.pred)
  pre_tab
  # save prediction result
  predict = data.frame(group = otutab_t[,"OTUgroup"], predicted=otutab.pred)



  summary(predict)

  return(list(pre_tab,summary(predict)))
}
