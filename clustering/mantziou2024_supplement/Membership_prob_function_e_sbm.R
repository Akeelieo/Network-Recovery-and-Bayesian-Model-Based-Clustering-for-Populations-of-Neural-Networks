
Membership_prob_e<-function(taus_c,p_c,q_c,G_m_c,G_dat){
  mem<-sum(G_m_c*(log(1-q_c)*G_dat+log(q_c)*(1-G_dat))+(1-G_m_c)*(log(p_c)*G_dat+log(1-p_c)*(1-G_dat)))
  return(log(taus_c)+mem)
}



