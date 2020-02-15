library(phyloseq)
library(tidyverse)
library(pheatmap)


ps_add_out = readRDS("./ps_all_OTU.rds")
ps_add_out



###多样性分析
library("phyloseq")
library(tidyverse)
head(a3)
mapping = as.data.frame(sample_data(ps_add_out))
head(mapping)
mapping$ID = row.names(mapping) 
sample_data(ps_add_out) = mapping
ps1 <- ps_add_out %>%
  subset_taxa(
    row.names(tax_table(ps_add_out)) %in% a3$id
  )
ps1

ps = ps1
ps

vegan_otu <-  function(physeq){
  OTU <-  otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <-  t(OTU)
  }
  return(as(OTU,"matrix"))
}
otu_table = as.data.frame(t(vegan_otu(ps)))
head(otu_table)


# otu_table = as.data.frame(t(vegan_otu(ps_Genus)))
head(otu_table)
design = as.data.frame(sample_data(ps))
## 计算相对丰度，计算每个物种丰度均值，按照均值排序
OTU = as.matrix(otu_table)
norm = t(t(OTU)/colSums(OTU,na=TRUE)) #* 100 # normalization to total 100
norma = norm %>% 
  t() %>% as.data.frame()
#数据分组计算平均值
iris.split <- split(norma,as.factor(design$SampleType))

iris.apply <- lapply(iris.split,function(x)colMeans(x,na.rm = TRUE))
# 组合结果
norm2 <- do.call(rbind,iris.apply)%>% 
  t() 
norm2 = as.data.frame(norm2)
norm2$mean=apply(norm2,1,mean)
norm2$ID = row.names(norm2)
colnames(norm2)
##按照mean、列进行排序desc设置从大到小排序
norm3<- arrange(norm2, desc(mean))
row.names(norm3) = norm3$ID
norm3$ID = NULL
### 提取前30个属
wt = norm3
wt$mean = NULL

### 添加OTU注释信息
vegan_tax <-  function(physeq){
  tax <-  tax_table(physeq)
  
  return(as(tax,"matrix"))
}
tax_table = as.data.frame(vegan_tax(ps))
head(tax_table)
wt
wt_tax = merge(wt,tax_table,by = "row.names",all = F)
head(wt_tax)
row.names(wt_tax) = wt_tax$Row.names

res = wt_tax
ID = rep("A",length(res$Row.names))
i = 1
for (i in 1:length(res$Row.names)) {
  
  if( is.na(as.character(res$Species[i])) |res$Species[i] == "unidentified"){
    if( is.na(as.character(res$Genus[i]))|res$Genus[i] == "unidentified"){
      if( is.na(as.character(res$Family[i]))|res$Family[i] == "unidentified"){
        if( is.na(as.character(res$Order[i]))|res$Order[i] == "unidentified"){
          if( is.na(as.character(res$Class[i]))|res$Class[i] == "unidentified"){
            if( is.na(as.character(res$Phylum[i]))|res$Phylum[i] == "unidentified"){
              if( is.na(as.character(res$Kingdom[i]))|res$Kingdom[i] == "unidentified"){
                
                
                
              }else{ID[i] = as.character(res$Kingdom[i])}
              
              
            }else{ID[i] = as.character(res$Phylum[i])}
            
          }else{ID[i] = as.character(res$Class[i])}
          
        }else{ID[i] = as.character(res$Order[i])}
        
      }else{ID[i] = as.character(res$Family[i])}
      
    }else{ID[i] = as.character(res$Genus[i])}
    
  }else{ID[i] = as.character(res$Species[i])}
  
}

ID
res$ID = ID



wt = res[,c(colnames(wt))]
head(wt)
color = colorRampPalette(c("navy", "white", "firebrick3"))(60)

wt2<-sqrt(wt)
wt2[wt2>0.5]<-0.5
wt2<-sqrt(wt2)




#开始绘图
#number_format设置保留小数点位数
#display_numbers是否显示数字
#注意可以自己设定显示内容：pheatmap(test, display_numbers = matrix(ifelse(test > 5, "*", ""), nrow(test)))#还可以自己设定要显示的内容；
# pheatmap(wt2,fontsize=6,cellwidth = 10, cellheight = 10,cluster_rows = FALSE,
#          color = colorRampPalette(c("navy", "white", "firebrick3"))(60),
#          display_numbers = TRUE, filename = paste(path,"heatmap.pdf"))


p = pheatmap(wt2,fontsize=6,cellwidth = 12, cellheight =6,cluster_rows = FALSE,cluster_cols = FALSE,
             color = colorRampPalette(c("navy", "white", "firebrick3"))(60),
             display_numbers = FALSE,labels_row  = res$ID,labels_col  = c("Disease","Health"))


p
path = getwd()
FileName2 <- paste(path,"/","mean_laading",".pdf", sep = "")
ggsave(FileName2, p, width = 6, height =12, device = cairo_pdf, family = "Times New Roman" )


#-----------------------------------------------全部-----------------样本--------------------



# otu_table = as.data.frame(t(vegan_otu(ps_Genus)))
head(otu_table)
design = as.data.frame(sample_data(ps))
## 计算相对丰度，计算每个物种丰度均值，按照均值排序
OTU = as.matrix(otu_table)
norm = t(t(OTU)/colSums(OTU,na=TRUE)) #* 100 # normalization to total 100
norma = norm %>% 
  t() %>% as.data.frame()
#数据分组计算平均值
iris.split <- split(norma,as.factor(design$SampleType))

iris.apply <- lapply(iris.split,function(x)colMeans(x,na.rm = TRUE))
# 组合结果
norm2 <- do.call(rbind,iris.apply)%>% 
  t() 
norm2 = as.data.frame(norm2)
norm2$mean=apply(norm2,1,mean)
norm2$ID = row.names(norm2)
colnames(norm2)
##按照mean、列进行排序desc设置从大到小排序
norm3<- arrange(norm2, desc(mean))
row.names(norm3) = norm3$ID
norm3$ID = NULL
### 提取前30个属
head(norm)
norm = as.data.frame(norm)
dim(norm)

dd = merge(norm,a3,by = "row.names",all = FALSE)
dim(dd)
ddd<- arrange(dd, desc(MeanDecreaseAccuracy))
head(ddd)
norm = ddd[1:280]
head(norm)
row.names(norm) =norm$Row.names
norm$Row.names = NULL
wt = norm
head(wt)
# wt$mean = NULL

### 添加OTU注释信息
vegan_tax <-  function(physeq){
  tax <-  tax_table(physeq)
  
  return(as(tax,"matrix"))
}
tax_table = as.data.frame(vegan_tax(ps))
head(tax_table)
head(wt)
wt_tax = merge(wt,tax_table,by = "row.names",all = TRUE)
head(wt_tax)
row.names(wt_tax) = wt_tax$Row.names
wt_tax$Row.names = NULL
#------paixu--------------------------------

dd = merge(wt_tax,a3,by = "row.names",all = FALSE)
head(dd)
ddd<- arrange(dd, desc(MeanDecreaseAccuracy))
head(ddd)
wt_tax1= ddd[1:287]
head(wt_tax1)
row.names(wt_tax1) =wt_tax1$Row.names
wt_tax1$Row.names = NULL





mapping =as.data.frame( sample_data(ps))
head(mapping)
res = wt_tax1
res$Row.names = row.names(res)
head(res)
# res$Species[i] == "unidentified"
# res$Genus[i] == "unidentified"
# 


ID = rep("A",length(res$Row.names))
i = 1
for (i in 1:length(res$Row.names)) {
  
  if( is.na(as.character(res$Species[i])) |res$Species[i] == "unidentified"){
    if( is.na(as.character(res$Genus[i]))|res$Genus[i] == "unidentified"){
      if( is.na(as.character(res$Family[i]))|res$Family[i] == "unidentified"){
        if( is.na(as.character(res$Order[i]))|res$Order[i] == "unidentified"){
          if( is.na(as.character(res$Class[i]))|res$Class[i] == "unidentified"){
            if( is.na(as.character(res$Phylum[i]))|res$Phylum[i] == "unidentified"){
              if( is.na(as.character(res$Kingdom[i]))|res$Kingdom[i] == "unidentified"){
                
                
                
              }else{ID[i] = as.character(res$Kingdom[i])}
              
              
            }else{ID[i] = as.character(res$Phylum[i])}
            
          }else{ID[i] = as.character(res$Class[i])}
          
        }else{ID[i] = as.character(res$Order[i])}
        
      }else{ID[i] = as.character(res$Family[i])}
      
    }else{ID[i] = as.character(res$Genus[i])}
    
  }else{ID[i] = as.character(res$Species[i])}
  
}

ID

ID
res$ID = ID



wt = res[,c(colnames(wt))]
head(wt)
color = colorRampPalette(c("navy", "white", "firebrick3"))(60)

wt2<-sqrt(wt)
wt2[wt2>3]<-3
wt2<-sqrt(wt2)
wt2<-sqrt(wt2)
wt2<-sqrt(wt2)


#开始绘图
#number_format设置保留小数点位数
#display_numbers是否显示数字
#注意可以自己设定显示内容：pheatmap(test, display_numbers = matrix(ifelse(test > 5, "*", ""), nrow(test)))#还可以自己设定要显示的内容；
# pheatmap(wt2,fontsize=6,cellwidth = 10, cellheight = 10,cluster_rows = FALSE,
#          color = colorRampPalette(c("navy", "white", "firebrick3"))(60),
#          display_numbers = TRUE, filename = paste(path,"heatmap.pdf"))

# annotation_row = data.frame(tax = mapping$SampleType)
# rownames(annotation_row) = rownames(wt2)

head(wt2)

aa = t(wt2)
ps
head(aa)
head(mapping)


index = merge(aa,mapping,by = "row.names",all = TRUE)
dim(index)
head(index)

asa<- arrange(index, SampleType)
head(asa)

#
annotation_col = data.frame(asa$SampleType) 
rownames(annotation_col) =asa$Row.names



asa$SampleType
dim(asa)
fia = asa[,1:41]
head(fia)
row.names(fia) = fia$Row.names
fia$Row.names=  NULL





fia =t(fia)
p = pheatmap(fia,fontsize=6,cellwidth = 0.5, cellheight =6,cluster_rows = FALSE,cluster_cols = FALSE,
             color = colorRampPalette(c( "white", "firebrick3"))(60),labels_row  = res$ID,annotation_col = annotation_col,
             display_numbers = FALSE,labels_col  = c("Disease","Health"))

# 
p


path = getwd()
FileName2 <- paste(path,"/","all_laading1",".pdf", sep = "")
ggsave(FileName2, p, width = 10, height =12, device = cairo_pdf, family = "Times New Roman" )








