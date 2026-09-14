
#log-likelihood of data
LogLikelih_mix_e<-function(centr,p,q,data_list,index){ #centr should be a vector of upper/lower trianle adjac
  sum_vec <- Reduce('+',data_list[index]) 
  Likel<-centr*(log(((1-q)*(1-p))/(p*q))*sum_vec+length(index)*log(q/(1-p)))+log(p/(1-p))*sum_vec+length(index)*log(1-p)
  return(sum(Likel))
}

#log-priors for Gm and p,q
Logpriors_mix_e<-function(centr,p,q,theta_c,memb_c,n_nodes,posit,a_0,b_0,c_0,d_0){ #centr should be vectors
  #prior of centroid Gm is an SBM
  centr<-lab_adj(centr,memb_c,posit) 
  theta_vec<-theta_c[match(names(centr),names(theta_c))]
  G_m_prior<-centr%*%log(theta_vec/(1-theta_vec))+sum(log(1-theta_vec))

  #Beta priors for p,q
  p_prior<-(a_0-1)*log(p)+(b_0-1)*log(1-p)
  q_prior<-(c_0-1)*log(q)+(d_0-1)*log(1-q)
  return(sum(G_m_prior,p_prior,q_prior))
}

#Log-posterior for Gm and p,q
Logpost_mix_e<-function(centr,p,q,data_list,index,theta_c,memb_c,n_nodes,posit,a_0,b_0,c_0,d_0){ #centr should be a vector
  return(LogLikelih_mix_e(centr,p,q,data_list,index)+Logpriors_mix_e(centr,p,q,theta_c,memb_c,n_nodes,posit,a_0,b_0,c_0,d_0))
}

#log-likelihood of data when UPDATING CENTROID for  OUTLIER cluster model
LogLikelih_outlier_clust<-function(centr,p_vec,q_vec,data_list){ #centr should be a vector of upper/lower trianle adjac
  log_vec_1<-log(((1-p_vec)*(1-q_vec))/(p_vec*q_vec))# returns a vector of length N_data
  log_vec_2<-log(q_vec/(1-p_vec))# returns a vector of length N_data
  log_vec_3<-log(p_vec/(1-p_vec))# returns a vector of length N_data
  log_vec_4<-log(1-p_vec)# returns a vector of length N_data
  sum_vec_1<-0
  sum_vec_2<-0
  for (i in 1:length(data_list)){
    sum_vec_1<-sum_vec_1+(data_list[[i]]*log_vec_1[i])+log_vec_2[i] #result is a vec
    sum_vec_2<-sum_vec_2+(data_list[[i]]*log_vec_3[i])+log_vec_4[i] #result is a vec
  }
  Likel<-centr%*%sum_vec_1+sum(sum_vec_2)
  return(Likel)
}
#log-priors for Gm 
Logpriors_outlier_clust<-function(centr,theta_c,memb_c,n_nodes,posit){ 
  #prior of centroid Gm is an SBM
  centr<-lab_adj(centr,memb_c,posit) 
  theta_vec<-theta_c[match(names(centr),names(theta_c))]
  G_m_prior<-centr%*%log(theta_vec/(1-theta_vec))+sum(log(1-theta_vec))
  return(G_m_prior)
}


Logpost_outlier_clust<-function(centr,p_vec,q_vec,data_list,theta_c,memb_c,n_nodes,posit){
  return(LogLikelih_outlier_clust(centr,p_vec,q_vec,data_list)+Logpriors_outlier_clust(centr,theta_c,memb_c,n_nodes,posit))
}