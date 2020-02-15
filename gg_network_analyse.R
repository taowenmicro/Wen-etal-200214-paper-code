r.threshold=0.4
p.threshold=0.05
d = "all"
ps = ps1 = readRDS("../ps_model_40.rds")
layouts = c("H_bulk","D_bulk")
y = matrix(1:2409,nrow = 14,ncol = length(layouts))
layout = layouts[1]
layout = layouts[2]
# 存储otu-sample矩阵的文件名
otu_sample_file <- paste("./data/",layout,"_",d,"_net",".txt",sep = "")
# sparcc cor 文件
r_sparcc_file<-paste("sparcc_net_result_matrix/",layout,"_",d,"/basis_corr/cor_sparcc.out",sep = "")
# sparcc p-value 文件
p_sparcc_file<-paste("sparcc_net_result_matrix/",layout,"_",d,"/pvals/pvals_two_sided.txt",sep = "")
# 文件读取
otu <- read.delim(otu_sample_file,row.names=1)
head(otu)
# 转置otu表
otu <-t(otu)
dim(otu)

# 读取r值矩阵
r_sparcc <- read.table(r_sparcc_file,row.names=1)
dim(r_sparcc)

# str(r_sparcc)
r_sparcc = as.matrix(r_sparcc)
# 读取p值矩阵
p_sparcc <- read.table(p_sparcc_file,row.names=1)
# str(p_sparcc)
p_sparcc  = as.matrix(p_sparcc)


occor.p<-as.matrix(p_sparcc)
##R值
occor.r<-as.matrix(r_sparcc)
# 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
occor.r[occor.p>p.threshold|abs(occor.r)<r.threshold] = 0 

head(occor.r)
dim(occor.r)




library(network)
library(ggplot2)
library(sna)
library(ergm)
library(igraph)
g <- network(occor.r, directed=FALSE,vertex.attrnames=TRUE)
summary(g)
net  = g
m <- as.matrix.network.adjacency(net)  # get sociomatrix
plotcord 

source("node_pro.R")
nodepro_result<-node_pro(igraph)
head(netplot)

# netplot$elements






# plotcord <- data.frame(gplot.layout.fruchtermanreingold(m, NULL))
plotcord<- data.frame(gplot.layout.circle(m, NULL))
# plotcord[[ 5]] <- data.frame(gplot.layout.circrand(m, NULL));names(plotcord[[5]]) = "circrand"
# for (ii in 1:18) {
plotcor = plotcord
dim(plotcor)
colnames(plotcor) = c("X1", "X2")
head(plotcor)
plotcor$elements <- colnames(occor.r)
edglist <- as.matrix.network.edgelist(net)
edglist = as.data.frame(edglist)
# aaaa = as.matrix.network(net)
# 构建igraph对象构建邻接矩阵
dim(occor.r)
igraph <- graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=TRUE)
E(igraph)$weight
edglist$weight = E(igraph)$weight

edges <- data.frame(plotcor[edglist[, 1], ], plotcor[edglist[, 2], ])
dim(edges)
edges$weight = E(igraph)$weight
##这里将边权重根据正负分为两类
aaa = rep("a",length(edges$weight))

for (i in 1:length(edges$weight)) {
  if (edges$weight[i]> 0) {
    aaa[i] = "+"
  }
  if (edges$weight[i]< 0) {
    aaa[i] = "-"
  }
  
}
#添加到edges中
edges$wei_label = aaa
colnames(edges) <- c("X1", "Y1","OTU_1", "X2", "Y2","OTU_2","weight","wei_label")
edges$midX <- (edges$X1 + edges$X2)/2
edges$midY <- (edges$Y1 + edges$Y2)/2
head(edges)
edges$wei_label
dim(edges)
# library(ggrepel)

#设置标签、
row.names(plotcor) = gsub("X","",plotcor$elements)
dim(plotcor)

### 合并细菌真菌ps对象
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
otu_table_16S = as.data.frame(t(vegan_otu(ps)))
head(otu_table_16S)

tax_table_16S = as.data.frame(vegan_tax(ps))

dim(tax_table_16S)

head(plotcor)
row.names(plotcor) = row.names(tax_table_16S)
#合并注释文件
dim(plotcor)
res = merge(plotcor,tax_table_16S,by = "row.names",all = F)


match(plotcor$elements,row.names(tax_table_16S))

dim(res)
head(res)
row.names(res) = res$Row.names
res$Row.names = NULL
head(res,n = 12)
res$Species[i] == "unidentified"
res$Genus[i] == "unidentified"



ID = rep("A",length(res$X1))
i = 1
for (i in 1:length(res$X1)) {
  
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




## 计算丰度均值
head(otu_table_16S)
sub_design = as.data.frame(sample_data(ps))

count  = vegan_otu(ps)
count2 = as.data.frame(count)
head(count2)
# count2$ID <- sub_design$SampleType
iris.split <- split(count2,as.factor(sub_design$SampleType))
#数据分组计算平均值
iris.apply <- lapply(iris.split,function(x)colMeans(x))
# 组合结果
iris.combine <- do.call(rbind,iris.apply)
ven2 = t(iris.combine)
head(ven2)

index = merge(res,ven2,by = "row.names",all = F)
dim(res)
dim(index)
row.names(index) = index$Row.names
index$Row.names = NULL


plotcord = index
head(plotcord)
library("dplyr")
#求取平均丰丰度
plotcord = mutate(plotcord, mean = (D_bulk+H_bulk)/2)

plotcord$mean

dim(plotcord)
row.names(plotcord) = plotcord$elements
plotcord$id = c(1:length(plotcord$X1))
plotcord$angle <- 90 - 360 * plotcord$id / length(plotcord$X1)
plotcord$angle[21:40] = plotcord$angle[21:40] +180
# plotcord$angle[23:44] = plotcord$angle[23:44] +180



plotcord$hjust <-  -0.15
plotcord$hjust[21:40] = 1.15
# plotcord$hjust[23:44] = 1.15


source("node_pro.R")
nodepro_result<-node_pro(igraph)
head(nodepro_result)

ind = merge(plotcord,nodepro_result,by = "row.names",all = FALSE)

dim(ind)
head(ind)
library(tidyverse)
netplot<- arrange(ind, desc(igraph.cen.degree))
head(netplot)

xa = netplot$elements[1:10]

xa

ss = rep("more",length(netplot$elements))
cc = rep(1,length(netplot$elements))
i  = 1

for (i in 1:length(netplot$elements)) {
  if( netplot$elements[i] %in% xa){
    ss[i] = "more"
  }else{
    
    ss[i] = "less"
  }
}


netplot$fil =  ss
netplot$fil = as.factor(netplot$fil)



head(edges)
library(ggrepel)
library(ggtree)

x=  c(-1.25,0,0,1.25)
x= x*1.5
y=  c(0,-1.25,1.25,0)
y = y*1.5
mer = data.frame(x = x,y =y)

geom_label2()
pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(wei_label)),
                                data = edges, size = 0.5) +
  geom_point(aes(X1, X2,size = mean,fill = fil),  pch = 21, data = netplot) + scale_colour_brewer(palette = "Set1") +
  geom_point(aes(x = x,y = y),data = mer,color = "white")+
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  labs( title = paste(layout,"network",sep = "_"))+
  
  # geom_text_repel(aes(X1, X2,label=Family),size=4, data = plotcord)+
  geom_text(aes(X1, X2,label=ID, angle = angle,  hjust=hjust),size=4, data = netplot)+
  scale_color_manual(values = c("#E41A1C","green","yellow" ))+
  # theme(plot.margin=unit(rep(3,4),'lines'))+
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())


pnet

# ?geom_segment


path = getwd()
FileName2 <- paste(path,"/",layout,".pdf", sep = "")



ggsave(FileName2, pnet, width = 12, height =10, device = cairo_pdf, family = "Times New Roman" )



occor.r = as.matrix(occor.r)
# # 构建igraph对象构建邻接矩阵
igraph <- graph_from_adjacency_matrix(occor.r,weighted=TRUE,diag=FALSE,mode = "directed")
igraph

?graph_from_adjacency_matrix
###网络边的赋值及其设置
igraph.weight <- E(igraph)$weight# 将igraph weight属性赋值到igraph.weight,用于后边做图
E(igraph)$weight <- NA
igraph<-remove.edge.attribute(igraph,"weight")#把边值删除
netpro_result<-net_pro(igraph)
colnames(netpro_result)<-layout
y = as.data.frame(y)
colnames(y) = layouts
# head(y)
y[layout] = netpro_result[,1]
row.names(y) = row.names(netpro_result)
y

FileName2 <- paste(path,"/","network_total_index",".csv", sep = "")
write.csv(y,tablename)

### 下面开始做网络模块分析




occor.r = as.matrix(occor.r)
# # 构建igraph对象构建邻接矩阵
igraph <- graph_from_adjacency_matrix(occor.r,weighted=TRUE,diag=FALSE,mode = "directed")
igraph

?graph_from_adjacency_matrix
###网络边的赋值及其设置
igraph.weight <- E(igraph)$weight# 将igraph weight属性赋值到igraph.weight,用于后边做图
E(igraph)$weight <- NA
igraph<-remove.edge.attribute(igraph,"weight")#把边值删除
netpro_result<-net_pro(igraph)
colnames(netpro_result)<-layout
y = as.data.frame(y)
colnames(y) = layouts
# head(y)
y[layout] = netpro_result[,1]
row.names(y) = row.names(netpro_result)
y

FileName2 <- paste(path,"/","network_total_index",".csv", sep = "")
write.csv(y,tablename)

### 下面开始做网络模块分析

