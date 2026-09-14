MH_update_rep <- function(ind_cl,G_m_current,G_data,perturb,c_current,lab,p_current,q_current,theta_current,n_nodes,a_0,b_0,c_0,d_0){
  G_m_new<-abs(G_m_current-rbinom(length(G_data[[1]]),1,perturb))
  G_m_new<-lab_adj(G_m_new,c_current,lab) #give membership labels to proposed centroid
  mhr=exp(Logpost_mix_e(G_m_new,p_current,q_current,G_data,ind_cl,theta_current,c_current,n_nodes,lab,a_0,b_0,c_0,d_0)-
            Logpost_mix_e(G_m_current,p_current,q_current,G_data,ind_cl,theta_current,c_current,n_nodes,lab,a_0,b_0,c_0,d_0))
  prob<-min(1, mhr)
  gen<-rbinom(1, 1, prob)
  if(gen==1){
    G_m_current<-G_m_new 
  }
  return(G_m_current)
}

MH_update_rep_out <- function(G_m_current,G_data,perturb,c_current,lab,p_current,q_current,theta_current,n_nodes){
  G_m_new<-abs(G_m_current-rbinom(length(G_data[[1]]),1,perturb))
  G_m_new<-lab_adj(G_m_new,c_current,lab) #give membership labels to proposed centroid
  mhr=exp(Logpost_outlier_clust(G_m_new,p_current,q_current,G_data,theta_current,c_current,n_nodes,lab)-
            Logpost_outlier_clust(G_m_current,p_current,q_current,G_data,theta_current,c_current,n_nodes,lab))
  prob<-min(1, mhr)
  gen<-rbinom(1, 1, prob)
  if(gen==1){
    G_m_current<-G_m_new 
  }
  return(G_m_current)
}

MH_update_p <- function(p_current,rw_steps,G_m_current,q_current,G_data,ind_cl,theta_current,c_current,n_nodes,lab,a_0,b_0,c_0,d_0){
  y<-p_current+runif(1,min = -rw_steps,max = rw_steps) # candidate proposal y, check the below constraints
  if (y>0.5){
    p_new<-(1-y)
  }else if (y<0){
    p_new<-(-y)
  }else if (y>0 & y<0.5){
    p_new<-y
  }
  mhr=exp(Logpost_mix_e(G_m_current,p_new,q_current,G_data,ind_cl,theta_current,c_current,n_nodes,lab,a_0,b_0,c_0,d_0)-
            Logpost_mix_e(G_m_current,p_current,q_current,G_data,ind_cl,theta_current,c_current,n_nodes,lab,a_0,b_0,c_0,d_0))
  
  prob<-min(1, mhr)
  gen<-rbinom(1, 1, prob)
  if (gen==1){
    p_current<-p_new
  }
  return(p_current)
}

MH_update_q <- function(p_current,rw_steps,G_m_current,q_current,G_data,ind_cl,theta_current,c_current,n_nodes,lab,a_0,b_0,c_0,d_0){
  y<-q_current+runif(1,min = -rw_steps,max = rw_steps) # candidate proposal y, check the below constraints
  if (y>0.5){
    q_new<-(1-y)
  }else if (y<0){
    q_new<-(-y)
  }else if (y>0 & y<0.5){
    q_new<-y
  }
  mhr=exp(Logpost_mix_e(G_m_current,p_current,q_new,G_data,ind_cl,theta_current,c_current,n_nodes,lab,a_0,b_0,c_0,d_0)-
            Logpost_mix_e(G_m_current,p_current,q_current,G_data,ind_cl,theta_current,c_current,n_nodes,lab,a_0,b_0,c_0,d_0))
  
  prob<-min(1, mhr)
  gen<-rbinom(1, 1, prob)
  if (gen==1){
    q_current<-q_new
  }
  return(q_current)
}

