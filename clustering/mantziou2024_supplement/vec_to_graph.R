# transform vectorised adjacency to graph
graph_reform<-function(vec,n_nodes){ 
  t<-matrix(rep(0,n_nodes*n_nodes),nrow = n_nodes,ncol = n_nodes)
  t[upper.tri(t)]<-vec
  t<-t(t)+t
  g<-graph_from_adjacency_matrix(t,mode = "undirected")
  return(g)
}

# transform vector to adjacency
adjac_reform_undir<-function(vec,num_nodes){ 
  t<-matrix(rep(0,num_nodes*num_nodes),nrow = num_nodes,ncol = num_nodes)
  t[upper.tri(t)]<-vec
  t<-t(t)+t
  return(t)
}

#function that vectorises adjacency matrices, excluding diagonal for undirected networks
mat_vect<-function(mat){
  return(mat[which(upper.tri(mat),arr.ind = TRUE)])
}
