igraph
##方式二：使用模块性给节点上色
# 模块性 modularity
fc <- cluster_fast_greedy(igraph,weights =NULL)# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
fc
modularity <- modularity(igraph,membership(fc))



# 按照模块为节点配色
comps <- membership(fc)
colbar <- rainbow(max(comps))
V(igraph)$color <- colbar[comps] 

layout

###改变网络节点边框宽度
mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}

add.vertex.shape("fcircle", clip=igraph.shape.noclip,
                 plot=mycircle, parameters=list(vertex.frame.color=1,
                                                vertex.frame.width=1))


# 
# plotname <- paste(e,"/co-occurrence_net",d,c,".pdf",sep = "")
# pdf(file = plotname,width = 30,height = 20)

l <- layout_with_fr(igraph)
# l <- layout_on_sphere(igraph)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

set.seed(123)
#节点颜色代表门，节点大小代表平均丰度，连线颜色代表正负相关
plot(igraph,main=paste(layout, "network",sep = ""), 
     layout=l,
     vertex.shape="fcircle", 
     # vertex.frame.color="white",
     vertex.label=NA ,
     vertex.frame.color="#984EA3",
     vertex.color="#984EA3",
     vertex.size =2,
     vertex.frame.alpha =0.5,
     edge.width=0.5,
     edge.lty=2,
     edge.curved=F,
     margin=c(0,0,0,0)
     #vertex.frame.width=5
)

legend(1,1.5, legend=c("pos.cor","neg.cor"),col=c("red","blue"),lty=1,lwd=2, bty="n",cex=1)
# dev.off()

path = "./cluster_analyse"
dir.create("./cluster_analyse")

plotname1 <- paste(path,"/co-occurrence_net_",d,layout,"_cluster",".pdf",sep = "")
pdf(file = plotname1,width = 8,height = 6)
ceb <- cluster_edge_betweenness(igraph,weights =NULL) 
#?cluster_edge_betweenness
dendPlot(ceb, mode="hclust")
# l <- layout_with_fr(igraph)
#l <- layout_on_sphere(igraph)
set.seed(123)
# plot(net.bg, rescale=F, layout=l*1.0)
#节点颜色代表门，节点大小代表平均丰度，连线颜色代表正负相关

plot(ceb,igraph,main=paste(layout,d ,"Co-occurrence network",sep= ""), 
     layout=l*1,#vertex.shape="fcircle",
     vertex.frame.color="black",
     vertex.label=NA,
     vertex.size =5,
     edge.width=1,edge.lty=1*2,
     edge.curved=F,margin=c(0,0,0,0),
     vertex.frame.width=5)
# 
# plot(ceb,igraph,main=paste(layout,d ,"Co-occurrence network cluster",sep= ""), 
#      layout=l,
#      vertex.shape="fcircle", 
#      vertex.frame.color="black",
#      vertex.label=NA,
#      vertex.label.cex=2,
#      vertex.label.color="white",
#      # vertex.size =igraph.size,
#      edge.width=1,
#      edge.lty=1,
#      edge.curved=F,
#      margin=c(0,0,0,0),
#      vertex.frame.width=5
# )
legend(1,1.5, legend=c("pos.cor","neg.cor"),col=c("red","blue"),lty=1,lwd=2, bty="n",cex=1)


dev.off()



# ceb <- cluster_edge_betweenness(igraph,weights =NULL) 
#?cluster_edge_betweenness
# dendPlot(ceb, mode="hclust")
# l <- layout_with_fr(igraph)
#l <- layout_on_sphere(igraph)

# 模块性 modularity
fc <- cluster_fast_greedy(igraph,weights =NULL)# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
modularity <- modularity(igraph,membership(fc))


fc

# 按照模块为节点配色
comps <- membership(fc)
#查看关注OTU存在那个模块中
comps["SH489839.07FU_JX376499_reps_singleton"]
set.seed(123)
## 按照某个模块构建子图
modules.g<-induced_subgraph(igraph,comps==4)
E(modules.g)$weight<-NA
# plotname5 <- paste(e,"/co-occurrence_net_sub_cluster",d,c,".pdf",sep = "")
# pdf(file = plotname5,width = 20,height = 10)
plot(modules.g,main=paste(layout,d,"modules_network",sep = ""), 
     #layout=l,
     vertex.shape="fcircle", 
     vertex.frame.color="black",
     vertex.label.cex=1,
     vertex.label.color="white",
     
     # vertex.label=node.label,
     # vertex.size =igraph.size,
     edge.width=1,
     edge.lty=1,
     edge.curved=F,
     margin=c(0,0,0,0),
     vertex.frame.width=5
)
# legend(1,1, legend=otu.tax.levels,col=levels(node.col), pch=19,cex=1,bty="n")
legend(1,1.5, legend=c("pos.cor","neg.cor"),col=c("red","blue"),lty=1,lwd=2, bty="n",cex=1)
# legend(1,.4, legend=node_lg1,bty="n",pch = 19,cex=1)
# legend(1.1,.4, legend=node_lg2, pch = 19,bty="n")

# dev.off()


##健康网络主要分为三个三个模块，现在分别做可视化

plotcord <- data.frame(gplot.layout.fruchtermanreingold(m, NULL))

# for (ii in 1:18) {
plotcor = plotcord
colnames(plotcor) = c("X1", "X2")
head(plotcor)
plotcor$elements <- colnames(occor.r)
edglist <- as.matrix.network.edgelist(net)
edglist = as.data.frame(edglist)
# aaaa = as.matrix.network(net)
# 构建igraph对象构建邻接矩阵
igraph <- graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
E(igraph)$weight
edglist$weight = E(igraph)$weight

edges <- data.frame(plotcor[edglist[, 1], ], plotcor[edglist[, 2], ])
head(edges)
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

head(tax_table_16S)

head(plotcor)

res = merge(plotcor,tax_table_16S,by = "row.names",all = F)
dim(res)
head(res)
row.names(res) = res$Row.names
res$Row.names = NULL
head(res)

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
#求取平均丰丰度
plotcord = mutate(plotcord, mean = (D_bulk+H_bulk)/2)

plotcord$mean

head(plotcord)

ab = fc[1]$`1`
ab = fc[2]$`2`
ab = fc[3]$`3`
ab = fc[4]$`4`
ab = fc[5]$`5`
ab = fc[6]$`6`
ab = fc[7]$`7`
plotcord_sub<- filter(plotcord, elements%in%ab)
dim(plotcord_sub)

head(edges)
edges_sub<- filter(edges, OTU_2%in%ab&OTU_1%in%ab)
dim(edges_sub)


library(ggrepel)
pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(wei_label)),
                                data = edges_sub, size = 0.5) +
  geom_point(aes(X1, X2,fill = Kingdom,size = mean),  pch = 21, data = plotcord_sub) + scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  labs( title = paste(layout,"network",sep = "_"))+
  geom_text_repel(aes(X1, X2,label=Family),size=4, data = plotcord_sub)+
  scale_color_manual(values = c("#377EB8","#E41A1C" ))+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())


pnet

FileName2 <- paste(path,"/",layout,"mod3",".pdf", sep = "")



ggsave(FileName2, pnet, width = 6, height =4, device = cairo_pdf, family = "Times New Roman" )
y
igraph

nodepro_result<-node_pro(igraph)
head(nodepro_result)
##合并丰度和注释文件
cs = plotcord
row.names(cs) = cs$elements

node = merge(nodepro_result,cs,by = "row.names",all = F)
head(node)
filename2 <- paste(path,"/igraph.node.pro",layout,d,".csv",sep = "")
write.csv(node,filename2)

### 提取hub节点-----------------------------------------------------------------------------------------------
##基于以下四种方法计算节点的网络属性
net.cn <- closeness(igraph)
net.bn <- betweenness(igraph)
net.pr <- page_rank(igraph)$vector
net.hs <- hub_score(igraph)$vector
#基于hub值排序
net.hs.sort <- sort(net.hs,decreasing = TRUE)
# ?sort
#选择前五个keystone 种类
net.hs.top5 <- head(net.hs.sort,n=5)

net.hs.top5

?hub_score

