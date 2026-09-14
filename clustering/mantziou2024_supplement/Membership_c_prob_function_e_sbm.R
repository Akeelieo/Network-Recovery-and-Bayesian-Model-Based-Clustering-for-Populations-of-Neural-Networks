mem_c_prob<-function(w_c,th_c,adj_c){
  mem_c<-adj_c%*%log(th_c/(1-th_c))+sum(log(1-th_c))
  return(log(w_c)+mem_c)
}
