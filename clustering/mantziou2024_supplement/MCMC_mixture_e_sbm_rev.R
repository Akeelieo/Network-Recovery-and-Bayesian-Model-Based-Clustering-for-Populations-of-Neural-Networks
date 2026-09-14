#############################################################################################################################
########### Main code for MCMC algorithm for mixture of measurement error models with SBM network representatives ###########
############################################ as presented in Section 4.3  ############################################
#############################################################################################################################

# load libraries
library(stats)
library(boot)
library(igraph)
library(coda)
library(nnet)
library(Matrix)
library(MCMCprecision)
library(Hmisc)
library(rlist)
library(plyr)
library(stringr)

# source R scripts containing functions called in MCMC 
#setwd("/YOUR_DIRECTORY")
source("Log_post_e_sbm_rev.R")
source("Membership_prob_function_e_sbm.R")
source("normalizing_log.R")
source("label_fun.R")
source("Membership_c_prob_function_e_sbm.R")
source("MH_updates.R")

MCMC_mix_err_sbm<-function(C,B,iter,burn_in,n_nodes,N_data,G_data,prob_vec,a_0,b_0,c_0,d_0,e_0,f_0,z_init,c_init,p_init,q_init,theta_init,perturb,rw_steps){ 
  ### START ARGUMENTS DESCRIPTION ###
  # C: number of clusters, B: number of blocks, iter: number of iterations of MCMC, burn_in: number of iterations to burn,
  # N_data: number of network observations, n_nodes: number of nodes of network data,
  # G_data: list of vectors/upper triangle matrices for symmetric networks, 
  # prob_vec: list of vector probabilities to produce proposal representative network in each cluster
  # a_0, b_0,c_0,d_0,e_0,f_0 each a vector of C length, denoting hyperparameters for Beta prior, for different clusters for p,q respectively,
  # z_init: vector of initial cluster membership of networks, c_init: list of vectors, one for each cluster, denoting initial block membership of the nodes for representative of each cluster
  # p_init, q_init: scalar values of p and q for each cluster, theta_init: initial val for theta,
  # perturb: vector with elements the perturbations of the representatives'edges,
  # rw_steps: vector with elements step sizes of RW proposals for p and q
  ### END OF ARGUMENTS DESCRIPTION ###
  
  #chains for posterior samples                                          
  chain_G_m <- lapply(1:C, function(x) matrix(NA, nrow = iter-burn_in, ncol = length(G_data[[1]]))) #list of chains for Gm for each cluster c
  chain_p <- matrix(NA, nrow = iter-burn_in, ncol = C)
  chain_q <- matrix(NA, nrow = iter-burn_in, ncol = C)
  chain_taus <- matrix(NA, nrow = iter-burn_in, ncol = C)
  chain_z <- matrix(NA, nrow = iter-burn_in, ncol = N_data)
  chain_c <- lapply(1:C, function(x) matrix(NA, nrow = iter-burn_in, ncol = n_nodes)) # list of chains for c membership of representatives' nodes, for each cluster c
  chain_w <- lapply(1:C, function(x) matrix(NA, nrow = iter-burn_in, ncol = B)) # list of chains for w prob of membership of representatives' nodes for each cluster c
  chain_theta <- lapply(1:C, function(x) matrix(NA, nrow = iter-burn_in, ncol = (B*(B+1))/2)) # list of chains for theta block probabilities for each cluster c
  
  #initialization of current values
  G_m_current <- lapply(1:C,function(x) rbinom(length(G_data[[1]]),1,prob_vec[[x]])) # initialize the representatives of each cluster, producing them from Bernoulli with prob=sum(data in cluster)/length(of data in cluster)
  p_current <- rep(p_init,C)  
  q_current <- rep(q_init,C) 
  z_current <- z_init
  c_current <- c_init # list with C vectors of the block membership of the representatives' nodes in each c cluster
  theta_current <- lapply(1:C, function(x) rep(theta_init,B*(B+1)/2)) # B(B+1)/2 corresponds to the number of elements in upper triangl and the diagonal
  sam <- 1 #counter for samples after burn_in
  w_current <- list()
  for (c in 1:C){ # Gibbs step Initialization of w_c for each cluster c using initialized c_current (weights for block membership of nodes)
    h2 <- tabulate(c_current[[c]])
    if(length(h2)<B){h2 <- c(h2,rep(0,B-length(h2)))} # if a block has remained with no allocation, fill with 0s missing due to tabulate fun
    w_current[[c]] <- rdirichlet(1,1+h2) #hyperarameter 1
  }
  
  # create names for the representative entries and for theta entries, corresponding to the block membership of the nodes, for each cluster
  toy_mat <- matrix(nrow = n_nodes,ncol = n_nodes)
  # <- which(upper.tri(toy_mat),arr.ind = TRUE) # gives the positions of the upper triangle elements in matrix of size n_nodes, i.e. the pairs of nodes(possible edges) for undirected networks
  #lab <- cbind(rep(1:n_nodes,n_nodes),rep(1:n_nodes,each=n_nodes)) # lab for directed case
  #colnames(lab) <- c("row","col")
  k <- 1
  lab_theta <- array()
  for (i in 1:B){
    for (j in i:B){
      lab_theta[k] <- paste(i,j,sep = "") # labels for theta that are combinations of the blocks
      k <- k+1
    }
  }
  for (k in 1:C){ # give names to elements of vectors in list (theta vector for each cluster) where these names remain the same throughout
    names(theta_current[[k]]) <- lab_theta # give names to each vector of theta probabilities for each cluster c
  }
  
  G_m_current <- lab_adj_list(G_m_current,c_current,lab) # call function lab_adj_list that returns a list with vectorized ajdacencies with names for each element, given the current membership for each cluster

  for (i in 1:iter){
    print(c("iter: ",i))
    # Gibbs step: Update τ1,...,τc weights for clusters c
    h <- tabulate(z_current)  
    if(length(h)<C){h <- c(h,rep(0,C-length(h)))} 
    taus_current <- rdirichlet(1,1+h) #update taus using z_current
    
    # MH Steps and Gibbs steps
    
    for (c in 1:C){
      v <- rmultinom(1,1,c(0.5,0.3,0.4,0.5,0.3,0.3,rep(c(0.2,0.2,0.3,0.4,0.4,0.2),2))) #to see if update Gm or q or q, for c cluster, and with which proposal to update

      # MH Step: Update Gm,c
      if (1 %in% v[1:5]){ # Update representative with one of the 6 possible proposals
        ind_cl <- which(z_current==c) # find which graphs belong to c class, gives their index to identify them in list of G_data
        G_m_current[[c]] <- MH_update_rep(ind_cl,G_m_current[[c]],G_data,perturb[which(v==1)],c_current[[c]],lab,
                                          p_current[c],q_current[c],theta_current[[c]],n_nodes,a_0[c],b_0[c],c_0[c],d_0[c])
        
      }
        
      if (v[6]==1){
        ind_cl <- which(z_current==c) #find which graphs belong to c class, gives their index to identify them in list of G_data
        G_m_new <- rbinom(length(G_data[[1]]),1,prob_vec[[c]]) # generate edges of proposal with diff probs according to prob_vec
        G_m_new <- lab_adj(G_m_new,c_current[[c]],lab) #give membership labels to proposed centroid
        d1 <- sum(dbinom(G_m_current[[c]],1,prob_vec[[c]],log = TRUE)) #proposal distr density for current state
        d2 <- sum(dbinom(G_m_new,1,prob_vec[[c]],log = TRUE)) #proposal distr density for proposal vector
        mhr = exp(Logpost_mix_e(G_m_new,p_current[c],q_current[c],G_data,ind_cl,theta_current[[c]],c_current[[c]],n_nodes,lab,a_0[c],b_0[c],c_0[c],d_0[c])+d1-
                  Logpost_mix_e(G_m_current[[c]],p_current[c],q_current[c],G_data,ind_cl,theta_current[[c]],c_current[[c]],n_nodes,lab,a_0[c],b_0[c],c_0[c],d_0[c])-d2)
        prob <- min(1, mhr)
        gen <- rbinom(1, 1, prob)
        if (gen==1){
          G_m_current[[c]] <- G_m_new
        }
      }
      
      # MH Step: Update p_c
      if(1 %in% v[7:12]){ 
        ind_cl <- which(z_current==c) #find which graphs belong to c class, gives their index to identify them in list of G_data
        p_current[c] <- MH_update_p(p_current[c],rw_steps[which(v==1)-6],G_m_current[[c]],q_current[c],G_data,ind_cl,theta_current[[c]],
                                    c_current[[c]],n_nodes,lab,a_0[[c]],b_0[[c]],c_0[[c]],d_0[[c]])
      }
      
      # MH Step: Update q_c
      if(1 %in% v[13:18]){ 
        ind_cl <- which(z_current==c) #find which graphs belong to c class, gives their index to identify them in list of G_data
        q_current[c] <- MH_update_q(p_current[c],rw_steps[which(v==1)-12],G_m_current[[c]],q_current[c],G_data,ind_cl,theta_current[[c]],
                                    c_current[[c]],n_nodes,lab,a_0[[c]],b_0[[c]],c_0[[c]],d_0[[c]])
      }

      
      # Gibbs step: Update w_c (weights for block membership of nodes) for each cluster c 
      h2 <- tabulate(c_current[[c]])
      if(length(h2)<B){h2 <- c(h2,rep(0,B-length(h2)))} 
      w_current[[c]] <- rdirichlet(1,0.5+h2) #hyperarameter 1
      
      # Gibbs step: Update theta_c (probabil of edge for given block membership of nodes) for each cluster c 
      for (b in 1:length(theta_current[[1]])){ #update each theta value for possible pairs of blocks
        ind <- names(theta_current[[c]][b])
        adj_ind <- G_m_current[[c]][which(names(G_m_current[[c]]) %in% ind)]
        n_ind <- length(adj_ind)
        adj_ind <- sum(adj_ind)
        theta_current[[c]][b] <- rbeta(1,adj_ind+e_0[c],f_0[c]+n_ind-adj_ind) 
      }
      
      # Gibbs step: Update c_c (membership of nodes of representative of each cluster c)
      pr<-array() #array of probs for node i to belong to group b
      for(n in 1:n_nodes){
        adj_local <- G_m_current[[c]][c(which(match(lab[,"col"],as.character(n))==1),which(match(lab[,"row"],as.character(n))==1))] # identify entries of adjacency that involve node i
        adj_local <- G_m_current[[c]][union(which(match(lab[,"col"],as.character(n))==1),which(match(lab[,"row"],as.character(n))==1))] # for directed case
        j_nodes <- setdiff(seq(1,n_nodes,1),n) # define the j nodes with i!=j
        j_nodes <- c(setdiff(1:n_nodes, n), setdiff(1:n_nodes, n))        
        c_local <- c_current[[c]][j_nodes] # find the membership of the j nodes
        for(b in 1:B){
          nam <- ifelse(c_local<b,paste(c_local,b,sep = ""),paste(b,c_local,sep = "")) # paste the membership of the j nodes with b of current iteration
          theta_local <- theta_current[[c]][match(nam,names(theta_current[[c]]))] # create a theta_local that corresponds to the entries of the adjacency given the memberships of the j nodes and b
          pr[b] <- mem_c_prob(w_current[[c]][b],theta_local,adj_local)
        }
        l_2 <- logxpy_2(pr[1],pr[2]) # normalize to sum to 1 for 2 Blocks (B=2)
        pr <- c(exp(pr[1]-l_2),exp(pr[2]-l_2))
        c_current[[c]][n] <- which(rmultinom(1,1,pr)!=0,arr.ind=TRUE)[,"row"]
      }
      G_m_current[[c]] <- lab_adj(G_m_current[[c]],c_current[[c]],lab) # update labels for nodes membership for adjacency
    } # end of iterations over c
    
    
    #Gibbs step: Update z1,...,zc
    
    p <- matrix(NA,nrow=N_data,ncol=C) # matrix with each row prob
    for (l in 1:N_data){
      for (m in 1:C){
        p[l,m] <- Membership_prob_e(taus_current[m],p_current[m],q_current[m],G_m_current[[m]],G_data[[l]])
      }
      p_temp <- exp(p[l,]-logxpy_C(p[l,]))    
      p[l, ] <- p_temp
    }
    p <- t(apply(p,1,function(x) {x/sum(x)})) # normalize the probability matrix in order for the rows to sum to 1
    z_current <- t(rMultinom(p,1))
    
    
    # Keep updates in chain if iteration>burn_in
    if(i>burn_in){
      chain_taus[sam, ] <- taus_current
      for (k in 1:C){
        chain_G_m[[k]][sam, ] <- G_m_current[[k]] 
        chain_c[[k]][sam, ] <- c_current[[k]]
        chain_w[[k]][sam, ] <- w_current[[k]]
        chain_theta[[k]][sam, ] <- theta_current[[k]]
      }
      chain_p[sam, ] <- p_current
      chain_q[sam, ] <- q_current
      chain_z[sam, ] <- z_current
      sam <- sam+1
    }
    
  } # end of MCMC iterations i
  return(list(chain_taus,chain_z,chain_p,chain_q,chain_G_m,chain_c,chain_w,chain_theta)) 
}

