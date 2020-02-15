net_pro<-function(igraph,outdir){
  # network property
  # 边数量 The size of the graph (number of edges)
  num.edges <- length(E(igraph)) # length(curve_multiple(igraph))
  num.edges
  # 顶点数量 Order (number of vertices) of a graph
  num.vertices <- length(V(igraph))# length(diversity(igraph, weights = NULL, vids = 	V(igraph)))
  num.vertices
  # 连接性(connectance) 网络中物种之间实际发生的相互作用数之和（连接数之和）占总的潜在相互作用数（连接数）的比例，可以反映网络的复杂程度
  connectance <- edge_density(igraph,loops=FALSE)# 同 graph.density;loops如果为TRUE,允许自身环（self loops即A--A或B--B）的存在
  connectance
  # 平均度(Average degree)
  average.degree <- mean(igraph::degree(igraph))# 或者为2M/N,其中M 和N 分别表示网络的边数和节点数。
  average.degree
  # 平均路径长度(Average path length)
  average.path.length <- average.path.length(igraph) # 同mean_distance(igraph) # mean_distance calculates the average path length in a graph
  average.path.length
  # 直径(Diameter)
  diameter <- diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
  diameter
  # 群连通度 edge connectivity / group adhesion
  edge.connectivity <- edge_connectivity(igraph)
  edge.connectivity
  # 聚集系数(Clustering coefficient)：分局域聚类系数和全局聚集系数，是反映网络中节点的紧密关系的参数，也称为传递性。整个网络的全局聚集系数C表征了整个网络的平均的“成簇性质”。
  clustering.coefficient <- transitivity(igraph) 
  clustering.coefficient
  no.clusters <- no.clusters(igraph)
  no.clusters
  # 度中心性(Degree centralization)
  centralization.degree <- centralization.degree(igraph)$centralization
  centralization.degree
  # 介数中心性(Betweenness centralization)
  centralization.betweenness <- centralization.betweenness(igraph)$centralization 
  centralization.betweenness
  # 紧密中心性(Closeness centralization)
  centralization.closeness <- centralization.closeness(igraph)$centralization
  centralization.closeness
  
  num.pos.edges<-sum(igraph.weight>0)# number of postive correlation
  num.neg.edges<-sum(igraph.weight<0)# number of negative correlation
  
  igraph.network.pro <- rbind(num.edges,num.pos.edges,num.neg.edges,num.vertices,connectance,average.degree,average.path.length,diameter,edge.connectivity,clustering.coefficient,no.clusters,centralization.degree,centralization.betweenness,centralization.closeness)
  rownames(igraph.network.pro)<-c("num.edges","num.pos.edges","num.neg.edges","num.vertices","connectance","average.degree","average.path.length","diameter","edge.connectivity","clustering.coefficient","no.clusters","centralization.degree","centralization.betweenness","centralization.closeness")
  colnames(igraph.network.pro)<- "value"
  igraph.network.pro
}