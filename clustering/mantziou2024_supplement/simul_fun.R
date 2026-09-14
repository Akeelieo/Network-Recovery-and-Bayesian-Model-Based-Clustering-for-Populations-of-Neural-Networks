
simulate_data_e<-function(G_m,p,q,n,c,seed_count){ #G_m list of representatives for each cluster, p,q vector of false pos/false neg for each cluster, n number of data/networks to be simulated,c number of clusters
  G_simul_data<-lapply(1:n, function(x) array(NA,dim = length(G_m[[1]])))
  count<-1
  for (k in 1:c){
    for (i in 1:(n%/%c)){ 
      for (j in 1:length(G_m[[1]])){
        if (G_m[[k]][j]==1){
          set.seed(seed_count) 
          G_simul_data[[count]][j]<-rbinom(1,1,1-q[k])
        }else{
          set.seed(seed_count) 
          G_simul_data[[count]][j]<-rbinom(1,1,p[k])
        }
        seed_count<-seed_count+1
      }
      count<-count+1
    }
  }
  if (n%%c !=0){
    for(i in 1:(n%%c)){
      for (j in 1:length(G_m[[1]])){
        if (G_m[[c]][j]==1){
          set.seed(seed_count) 
          G_simul_data[[count]][j]<-rbinom(1,1,1-q[c])
        }else{
          set.seed(seed_count) 
          G_simul_data[[count]][j]<-rbinom(1,1,p[c])
        }
        seed_count<-seed_count+1
      }
      count<-count+1
    }
  }
  return(G_simul_data)
}

# difference to above function: set.seed() is placed outside the for loop and is increased by 1 after whole vector/network is generated with current seed
simulate_data_e_alter_seed<-function(G_m,p,q,n,c,seed_count){ #G_m list of representatives for each cluster, p,q vector of false pos/false neg for each cluster, n number of data/networks to be simulated,c number of clusters
  G_simul_data<-lapply(1:n, function(x) array(NA,dim = length(G_m[[1]])))
  count<-1
  for (k in 1:c){
    for (i in 1:(n%/%c)){
      set.seed(seed_count) 
      for (j in 1:length(G_m[[1]])){
        if (G_m[[k]][j]==1){
          G_simul_data[[count]][j]<-rbinom(1,1,1-q[k])
        }else{
          G_simul_data[[count]][j]<-rbinom(1,1,p[k])
        }
      }
      seed_count<-seed_count+1
      count<-count+1
    }
  }
  if (n%%c !=0){
    for(i in 1:(n%%c)){
      set.seed(seed_count)
      for (j in 1:length(G_m[[1]])){
        if (G_m[[c]][j]==1){
          G_simul_data[[count]][j]<-rbinom(1,1,1-q[c])
        }else{
          G_simul_data[[count]][j]<-rbinom(1,1,p[c])
        }
      }
      seed_count<-seed_count+1
      count<-count+1
    }
  }
  return(G_simul_data)
}


#Generate the prob vector (for proposals of Gm) using the data within the clusters (for case of C=3 only)
prob_vec_fun<-function(z_true,sim_data){
  sum_1<-sum_2<-sum_3<-0
  for (i in which(z_true==1)){
    sum_1<-sum_1+sim_data[[i]]
  }
  for (i in which(z_true==2)){
    sum_2<-sum_2+sim_data[[i]]
  }
  for (i in which(z_true==3)){
    sum_3<-sum_3+sim_data[[i]]
  }
  return(list(sum_1/length(which(z_true==1)),sum_2/length(which(z_true==2)),sum_3/length(which(z_true==3))))
}

# prob vec for general case (any size of C)
prob_vec_fun_gen<-function(z_init,sim_data,C_max){
  prob_vec_list<-list()
  for (c in 1:C_max){
    prob_vec_list[[c]]<-Reduce('+',sim_data[which(z_init==c)])/length(which(z_init==c))
  }
  return(prob_vec_list)
}

#generate edges of graphs with prob that depends on the node membership
gener_edges_sbm<-function(named_adj_vec,theta,seed_val){ 
  for(i in 1:length(named_adj_vec)){
    pr<-theta[match(names(named_adj_vec[i]),names(theta))]
    set.seed(seed_val)
    named_adj_vec[i]<-rbinom(1,1,pr)
    seed_val<-seed_val+1
  }
  return(named_adj_vec)
}
