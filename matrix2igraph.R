matrix2igraph<-function(matr,r.threshold,p.threshold){
  occor<-corAndPvalue(matr,method = c( "spearman"))
  # multiple test the p values
  mtadj<-mt.rawp2adjp(unlist(occor$p),proc="BH")
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor,dim(matr)[2])
  ##R值
  occor.r<-occor$cor
  # 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
  occor.r[occor.p>p.threshold|abs(occor.r)<r.threshold] = 0 
  
  # 构建igraph对象
  igraph <- graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
  igraph
  # NOTE:可以设置weighted=NULL,但是此时要注意此函数只能识别相互作用矩阵内正整数，所以应用前请确保矩阵正确。
  # 可以按下面命令转换数据
  # occor.r[occor.r!=0] <- 1
  # igraph <- graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=NULL,diag=FALSE)
  
  # 是否去掉孤立顶点，根据自己实验而定
  # remove isolated nodes，即去掉和所有otu均无相关性的otu 可省略，前期矩阵已处理过
  bad.vs <- V(igraph)[degree(igraph) == 0]
  igraph <- delete.vertices(igraph, bad.vs)
  igraph
}