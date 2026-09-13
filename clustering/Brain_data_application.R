# This script contains code for brain data application in Section 6.2

# source R script with main code for MCMC with SFM extension
setwd("/YOUR_DIRECTORY")
source("vec_to_graph.R")
source("MCMC_outlier_cluster_sbm.R")

library(fst)
library(stringr)
library(igraph)

# read data in folder neurodata_edge_lists (downloaded from : https://neurodata.io/mri/)

# Create graphs, adjacencies and vectorised adjacencies from edgelist brain data
ind<-1
indiv<-c(27:56) # individuals in study
nodes<-seq(1:200)
all_brain_vec_list<-list()
all_brain_adj_list<-list()
all_brain_graph_list<-list()
for (i in indiv){
  for (j in 1:10){ # corresponding to measurements for each individual
    el<-read.csv(paste("/YOUR_DIRECTORY/neurodata_edge_lists/sub-00254",i,"_ses-",j,"_dwi_CPAC200.ssv",sep = ""))
    el_mat<-str_split_fixed(el[,1], " ",3) #split to columns
    el_mat<-apply(el_mat,2,as.integer) #make the entries integers
    g<-simplify(graph_from_data_frame(as.data.frame(el_mat)[ ,1:2],directed = FALSE,vertices= nodes))
    all_brain_graph_list[[ind]]<-g
    all_brain_adj_list[[ind]]<-get.adjacency(g,sparse = TRUE)
    all_brain_vec_list[[ind]]<-mat_vect(as.matrix(all_brain_adj_list[[ind]]))
    ind<-ind+1 
  }
}


# Visualize network data
par(mfrow=c(1,2))
g_1<-graph_reform(all_brain_vec_list[[1]],200)
coord<-layout.auto(g_1)
g_2<-graph_reform(all_brain_vec_list[[2]],200)
g_120<-graph_reform(all_brain_vec_list[[120]],200)
plot(g_120,vertex.size=2,edge.color="slategrey", vertex.color="black",vertex.label=NA,layout=coord)


#### DISTANCE MATRICES FOR JACCARD, WAVELETS, L_2 DISTANCE ####

#Function that Creates a matrix of Jaccard distances between graphs
jaccard_mat<-function(adj_lst){
  jaccard_mt<-matrix(rep(0,length(adj_lst)*length(adj_lst)) ,nrow = length(adj_lst),ncol = length(adj_lst))
  for (i in 1:length(adj_lst)){
    for (j in i:length(adj_lst)){
      print(c(i,j))
      jaccard_mt[i,j]<-sum(abs(adj_lst[[i]]-adj_lst[[j]]))/(sum(pmax(adj_lst[[i]],adj_lst[[j]])))
      jaccard_mt[j,i]<-jaccard_mt[i,j]
    }
  }
  return(jaccard_mt)
}

jacc_mat<-jaccard_mat(all_brain_adj_list)


#Function that Creates matrix for l2 distance measure
l2_mat<-function(adj_lst){
  l2_mt<-matrix( rep(0,length(adj_lst)*length(adj_lst)) ,nrow = length(adj_lst),ncol = length(adj_lst))
  for (i in 1:length(adj_lst)){
    for (j in i:length(adj_lst)){
      print(c(i,j))
      tr<-eigen(laplacian_matrix(graph_from_adjacency_matrix(adj_lst[[i]],mode = "undirected"),normalized = TRUE),only.values = TRUE)$values
      tr2<-eigen(laplacian_matrix(graph_from_adjacency_matrix(adj_lst[[j]],mode = "undirected"),normalized = TRUE),only.values = TRUE)$values
      l2_mt[i,j]<-sum(abs(exp(-0.1*tr2)-exp(-0.1*tr))^2)
      l2_mt[j,i]<-l2_mt[i,j]
    }
  }
  return(l2_mt)
}

l2_mat<-l2_mat(all_brain_adj_list)

#Function that Creates matrix for Wavelets distance measure
wave_mat<-function(adj_lst,n_nodes,tau){
  wave_mt<-matrix( rep(0,length(adj_lst)*length(adj_lst)),nrow = length(adj_lst),ncol = length(adj_lst))
  for (i in 1:length(adj_lst)){
    for (j in i:length(adj_lst)){
      print(c(i,j))
      tr<-eigen(laplacian_matrix(graph_from_adjacency_matrix(adj_lst[[i]],mode = "undirected"),normalized = TRUE))
      xt1<-diag(exp(-tau*(tr$values))) 
      delta1<-(tr$vectors)%*%xt1%*%t(tr$vectors)
      
      tr2<-eigen(laplacian_matrix(graph_from_adjacency_matrix(adj_lst[[j]],mode = "undirected"),normalized = TRUE))
      xt2<-diag(exp(-tau*(tr2$values))) 
      delta2<-(tr2$vectors)%*%xt2%*%t(tr2$vectors)
      
      delta<-delta1-delta2
      wave_mt[i,j]<-sum(diag(t(delta)%*%delta))/n_nodes
      wave_mt[j,i]<-wave_mt[i,j]
    }
  }
  return(wave_mt)
}

wave_mat<-wave_mat(all_brain_adj_list,200,1.2)


####################################################################################
############## INITIALIZATION FOR OUTLIER CLUST ALGO AND SBM CENTROID ##############
####################################################################################

#### INITIALIZATION BASED ON COMBINATION OF METRICS AND K-MEMOID CLUST ###

library(kmed)
jacc_memoid<-fastkmed(jacc_mat,ncluster = 2,iterate = 50)$cluster
l2_memoid<-fastkmed(l2_mat,ncluster = 2,iterate = 50)$cluster
wave_memoid<-fastkmed(wave_mat,ncluster = 2,iterate = 50)$cluster
#majority vote of communities for the three metrics, to determine a single initialization of cluster membership
combin_comm_memoid<-c()
for(i in 1:300){
  t<-c(jacc_memoid[i],l2_memoid[i],wave_memoid[i])
  if(length(which(t==1))>=2){
    combin_comm_memoid[i]<-1
  }else{
    combin_comm_memoid[i]<-2
  }
}

z_init<-combin_comm_memoid


prob_vec<-Reduce("+",all_brain_vec_list)/length(all_brain_vec_list)

a_0<-b_0<-c_0<-d_0<-rep(0.5,2)
e_0<-f_0<-0.5
set.seed(12)
centr_est<-rbinom(19900,1,prob_vec)

# transform representative from vector to graph
centr_est_g<-graph_reform(centr_est,length(nodes))

##SBM estimation of representative 
library("blockmodels")
library("ramify")
centr_est_adj<-adjac_reform_undir(centr_est,length(nodes))#get adjacency from vector
sbm_1<-BM_bernoulli("SBM",centr_est_adj)
sbm_1$estimate()
K<-2
sbm_1_memb<-argmax(sbm_1$memberships[[K]]$Z)

c_init<-sbm_1_memb

perturb<-c(0.0001,0.0006,0.0003,0.001,0.0002)
rw_steps<-c(0.01,0.005,0.001,0.0001,0.0005,0.05)


# Run on cluster in parallel the MCMC for different initialisation of p,q,theta

# Index from 1 to 3
job_ind<-as.numeric(Sys.getenv("SOME_ARRAYID"))

p_init_mul<-q_init_mul<-c(0.01,0.018,0.35)
theta_init_mul<-c(0.2,0.5,0.7)

MCMC_mix_err_sbm(2,2,1000000,50000,0,200,300,all_brain_vec_list,prob_vec,a_0,
                 b_0,c_0,d_0,e_0,f_0,z_init,c_init,p_init_mul[job_ind],q_init_mul[job_ind],theta_init_mul[job_ind],
                 perturb,rw_steps)

# Visualisations

# function to bind results from MCMC sub iterations saved from MCMC algorithm
bind_fun<-function(name_list,param,direct){
  new_list<-list()
  for(i in 1:length(name_list)){
    new_list[[i]]<-readRDS(paste(direct,param,name_list[i],".rds",sep = ""))
  }
  return(do.call(rbind,new_list))
}

# Visualise Figure 24, Section 6.2
#### MEMBERSHIP OF NETWORKS Z ####
#### LAST 100,000, DIFFERENT INTIALIZATIONS ####

dir<-"/YOUR_DIRECTORY"
par<-"z_mem_"
iter_val_1<-c(950000,1e+06)
iter_val_2<-c("950000_init_0.018","1e+06_init_0.018")
iter_val_3<-c("950000_init_0.35","1e+06_init_0.35")
z_val_1<-bind_fun(iter_val_1,par,dir)
z_val_2<-bind_fun(iter_val_2,par,dir)
z_val_3<-bind_fun(iter_val_3,par,dir)

library(ggplot2)
library(reshape2)
z_memb_df<-function(z,iter_start,iter_end){
  memb_prop<-matrix(NA,nrow = iter_end-iter_start+1,ncol = 2)
  ind<-1
  for(i in iter_start:iter_end){
    test<-tabulate(z[,i]) #prop of time individ i in "1" and in "2"
    if(length(test)<2){test<-c(test,rep(0,2-length(test)))}
    memb_prop[ind,]<-test
    ind<-ind+1
  }
  z_df<-as.data.frame(memb_prop)
  names(z_df)<-c("1","2")
  melted_df<-melt(z_df)
  melted_df<-data.frame(melted_df,c(seq(1,iter_end-iter_start+1,1),seq(1,iter_end-iter_start+1,1)))
  names(melted_df)<-c("Cluster","value","Observation")
  return(melted_df)
}

start_list<-seq(1,300,10)
for(j in 1:30){
  ind_z<-z_memb_df(z_val_1,start_list[j],start_list[j]+9)
  assign(paste("individ_",j,"_z",sep = ""),ggplot(ind_z, aes(fill=Cluster, y=value/100000, x=Observation)) + 
           geom_bar(position="stack", stat="identity")+ylab("proportion")+
           theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16))+
           scale_x_continuous(expand = c(0,0.05),breaks=seq(1, 10, by = 1))+scale_y_continuous(expand = c(0,0.01))+scale_fill_manual(values = c("1" = "grey", "2" = "blue"))
  )
}

library(patchwork)
combined <- individ_1_z+individ_2_z+individ_3_z+individ_4_z+individ_5_z+individ_6_z+individ_7_z+
  individ_8_z+individ_9_z+individ_10_z+individ_11_z+individ_12_z+individ_13_z+individ_14_z+individ_15_z+individ_16_z+individ_17_z+
  individ_18_z+individ_19_z+individ_20_z+individ_21_z+individ_22_z+individ_23_z+individ_24_z+individ_25_z+individ_26_z+individ_27_z+individ_28_z+
  individ_29_z+individ_30_z& theme(legend.position = "bottom",legend.text=element_text(size=16),legend.title=element_text(size=16))

combined[[16]] <-combined[[16]] + theme(axis.title.y = element_text(size = 16))
combined[[28]] <-combined[[28]] + theme(axis.title.x = element_text(size = 16))

combined + plot_layout(guides = "collect",ncol=5)

# Visualise Figure 25, Section 6.2

# membership of individuals (from first initialisation)
memb_prop<-matrix(NA,nrow = 300,ncol = 2)
for(i in 1:300){
  test<-tabulate(z_val_1[,i]) #prop of time individ i in "1" and in "2"
  if(length(test)<2){test<-c(test,rep(0,2-length(test)))}
  memb_prop[i,]<-test
}
z_df_1<-as.data.frame(memb_prop)
# z_df_2<-as.data.frame(memb_prop)
# z_df_3<-as.data.frame(memb_prop)

ind_memb_1<-apply(z_df_1, 1,which.max)
# ind_memb_2<-apply(z_df_2, 1,which.max)
# ind_memb_3<-apply(z_df_3, 1,which.max)

ind_memb_1_cl_1<-which(ind_memb_1==1)
ind_memb_1_cl_2<-which(ind_memb_1==2)

# features of graphs of individuals according to cluster membership

features<-matrix(NA,nrow = 3,ncol = length(ind_memb_1_cl_1))
ind<-1
for (i in ind_memb_1_cl_1){
  gr_temp<-all_brain_graph_list[[i]]
  features[1,ind]<-transitivity(gr_temp)
  features[2,ind]<-average.path.length(gr_temp, directed = FALSE)
  features[3,ind]<-mean_distance(gr_temp,directed = FALSE)
  ind<-ind+1
}
features1<-features

features<-matrix(NA,nrow = 3,ncol = length(ind_memb_1_cl_2))
ind<-1
for (i in ind_memb_1_cl_2){
  gr_temp<-all_brain_graph_list[[i]]
  features[1,ind]<-transitivity(gr_temp)
  features[2,ind]<-average.path.length(gr_temp, directed = FALSE)
  features[3,ind]<-mean_distance(gr_temp,directed = FALSE)
  ind<-ind+1
}
features2<-features
features_all<-cbind(features1,features2)

col_class<-c(rep(1,length(ind_memb_1_cl_1)),rep(2,length(ind_memb_1_cl_2)))

# from features keep the transitivity (clustering coeff) and average path
df_features_all<-data.frame(x=features_all[1,],y=features_all[2,],cluster=as.factor(col_class))


#gray scale ggplot
ggplot(df_features_all, aes(x = x, y = y, color = cluster)) +geom_point() +scale_color_manual(breaks = c("1", "2"),values=c("black", "gray60"))+theme_classic()+
  xlab("clustering coefficient")+ylab("average shortest path length")+theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),
                                                                            axis.title.x=element_text(size=13,hjust = 0.5),axis.title.y=element_text(size=13,hjust = 0.5),
                                                                            legend.text=element_text(size=13))

#### MEMBERSHIP OF NODES C for network representative ####
#### LAST 100,000 ####
par<-"c_mem_"
iter_val_1<-c(950000,1e+06)
iter_val_2<-c("950000_init_0.018","1e+06_init_0.018")
iter_val_3<-c("950000_init_0.35","1e+06_init_0.35")
c_val_1<-bind_fun(iter_val_1,par,dir)
c_val_2<-bind_fun(iter_val_2,par,dir)
c_val_3<-bind_fun(iter_val_3,par,dir)

memb_prop<-matrix(NA,nrow = 200,ncol = 2)
for(i in 1:200){
  test<-tabulate(c_val_1[,i]) #prop of time individ i in "1" and in "2"
  if(length(test)<2){test<-c(test,rep(0,2-length(test)))}
  memb_prop[i,]<-test
}
memb_prop_df<-as.data.frame(memb_prop)

memb_prop<-matrix(NA,nrow = 200,ncol = 2)
for(i in 1:200){
  test<-tabulate(c_val_2[,i]) #prop of time individ i in "1" and in "2"
  if(length(test)<2){test<-c(test,rep(0,2-length(test)))}
  memb_prop[i,]<-test
}
memb_prop_df_2<-as.data.frame(memb_prop)

memb_prop<-matrix(NA,nrow = 200,ncol = 2)
for(i in 1:200){
  test<-tabulate(c_val_3[,i]) #prop of time individ i in "1" and in "2"
  if(length(test)<2){test<-c(test,rep(0,2-length(test)))}
  memb_prop[i,]<-test
}
memb_prop_df_3<-as.data.frame(memb_prop)

# Visualise results for representative networks, Figure 23, Section 6.2

######### Posterior mode for network representative (3 initialisations)

# get posterior mode 
library(mgcv)
t<-uniquecombs(as.matrix(rep_1_mem)) # rep_1_mem: saved results from MCMC last 50,000 iterations

t_ind<-which.max(table(attr(t,"index")))

if(length(table(attr(t,"index")))!=1){
  post_mode_1<-t[t_ind,]
}else{
  post_mode_1<-t
}

# repeat above for posterior mode from the other two initialisations

g_1<-graph_reform(as.vector(unlist(post_mode_1[[1]])),200) # from init 1
g_2<-graph_reform(as.vector(unlist(post_mode_2[[1]])),200) # from init 2
g_3<-graph_reform(as.vector(unlist(post_mode_3[[1]])),200) # from init 3

#get edgelists 
el_1<-get.edgelist(g_1)
el_2<-get.edgelist(g_2)
el_3<-get.edgelist(g_3)

nodes<-seq(1:200)
g_1_rev<-simplify(graph_from_data_frame(as.data.frame(el_1),directed = FALSE,vertices= nodes))
g_2_rev<-simplify(graph_from_data_frame(as.data.frame(el_2),directed = FALSE,vertices= nodes))
g_3_rev<-simplify(graph_from_data_frame(as.data.frame(el_3),directed = FALSE,vertices= nodes))

dif_1<-difference(g_2_rev,g_1_rev)
dif_2<-difference(g_1_rev,g_2_rev)
E(dif_1)$color<-"magenta"
E(dif_2)$color<-"black"
new_g<-graph.union(dif_1, dif_2, byname=TRUE)

dif_3<-difference(g_3_rev,g_1_rev)
dif_4<-difference(g_1_rev,g_3_rev)
E(dif_3)$color<-"magenta"
E(dif_4)$color<-"darkslategrey"
new_g_2<-graph.union(dif_3, dif_4, byname=TRUE)

#### Plot centroid networks with coloured nodes according to membership #### 
node_memb<-apply(memb_prop_df, 1,which.max)
node_memb_2<-apply(memb_prop_df_2, 1,which.max)
node_memb_3<-apply(memb_prop_df_3, 1,which.max)


################# GRAY SCALES ################# 
 
#keep layout fixed
set.seed(1)
lt<-layout.fruchterman.reingold(g_1)

E(new_g)$color<-ifelse(is.na(E(new_g)$color_2),"gray48","black")
E(new_g)$width<-ifelse(is.na(E(new_g)$color_2),1,1.2) # or 1.5 for black

E(new_g_2)$color<-ifelse(is.na(E(new_g_2)$color_2),"gray48","black")
E(new_g_2)$width<-ifelse(is.na(E(new_g_2)$color_2),1,1.2)

#representative of 1st initial value
V(g_1_rev)$color<-ifelse(node_memb==1,"white","grey36")
plot(g_1_rev,vertex.size=5,vertex.label.color="black",
     vertex.label.font=2, edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,vertex.label=NA,layout=lt)
#representative of 2nd initial value (NOT IN COMMON EDGES WITH REF NETWORK FROM 1st INIT VALUE)
V(new_g)$color<-ifelse(node_memb_2==1,"white","grey36")
plot(new_g,vertex.size=5,vertex.label.color="black",
     vertex.label.font=2, edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,vertex.label=NA,layout=lt)
#representative of 3rd initial value (NOT IN COMMON EDGES WITH REF NETWORK FROM 1st INIT VALUE)
V(new_g_2)$color<-ifelse(node_memb_3==1,"white","grey36")
plot(new_g_2,vertex.size=5,vertex.label.color="black",
     vertex.label.font=2, edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,vertex.label=NA,layout=lt)
