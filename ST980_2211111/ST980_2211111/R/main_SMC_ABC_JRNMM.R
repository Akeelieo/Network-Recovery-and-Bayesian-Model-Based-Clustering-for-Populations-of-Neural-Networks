#-----------------------------------------------------------------------------------
# Author: Akeel Shah
# Date:   2025-03-28
#
# Description: nSMC-ABC method for network inference and parameter estimation
#              in the stochastic multi-population JRNMM, proposed in the paper:            
#
#              Network inference via approximate Bayesian computation. 
#              Illustration on a stochastic multi-population neural mass model, 
#              by S. Ditlevsen, M. Tamborrino and I. Tubikanec
#-------------------------------------------------------------------------------

#--------------------------------------------
#load files
#--------------------------------------------

source(file="required_packages.R")
source(file="matrices_SplittingJRNMM.R")
source(file="functions_SMC_ABC_JRNMM.R")
source(file="prepare_EEG_functions.R")        # <-- new

#--------------------------------------------
#set file names
#--------------------------------------------

#--------------------------------------------
# set file names (from command line)
#--------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) args <- c("chb01_03", "during")

record        <- args[1]                                    # e.g. "chb01_15"
period        <- args[2]                                    # "before" or "during"
seizure_index <- if (length(args) >= 3) as.integer(args[3]) else 1

edf_dir  <- "."      # folder holding the .edf files

# summary file follows from the record: "chb01_15" -> "chb01-summary.txt"
patient      <- sub("_.*$", "", record)
summary_path <- file.path(edf_dir, paste0(patient, "-summary.txt"))

# fail fast, before the expensive pilot run
stopifnot(length(args) >= 2,
          period %in% c("before", "during"),
          file.exists(file.path(edf_dir, paste0(record, ".edf"))),
          file.exists(summary_path))


# pull this record's seizure interval out of the summary file
sm  <- readLines(summary_path)
hit <- grep(paste0("File Name: ", record, ".edf"), sm, fixed = TRUE)
stopifnot(length(hit) == 1)

blk <- sm[hit:min(hit + 12, length(sm))]
st  <- as.numeric(sub("\\D*(\\d+)\\s*seconds.*", "\\1",
                      grep("Seizure.*Start Time", blk, value = TRUE)[seizure_index]))
en  <- as.numeric(sub("\\D*(\\d+)\\s*seconds.*", "\\1",
                      grep("Seizure.*End Time",   blk, value = TRUE)[seizure_index]))

T_window <- en - st
stopifnot(is.finite(T_window), T_window > 0, st >= T_window)

cat("record", record, period, "| seizure", st, "-", en,
    "| window", T_window, "s\n")

#derived: both paths carry the record and period, so runs cannot overwrite
#each other and every results folder is self-describing
filename_data   <-file.path("RefData",record,
                            ifelse(period=="during","duringSeizure","beforeSeizure"),
                            "Data")
#include the seizure index so multiple seizures within the same record/period
#write to distinct folders instead of overwriting one another
filename_results<-paste0("ABC_Results_",record,"_",period,"_s",seizure_index)

dir.create(filename_results,recursive=TRUE,showWarnings=FALSE)

#-------------------------------------------------------------------------------
# PREPARE SMC-ABC
#-------------------------------------------------------------------------------

#number of cores used for parallel computation.
#On the cluster, use the SLURM allocation (SLURM_CPUS_PER_TASK) rather than
#detectCores(), which reports the whole physical node and would oversubscribe
#a shared node. Falls back to detectCores()-1 when run outside SLURM (e.g. locally).
ncl<-as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset=NA))
if(is.na(ncl)) ncl<-max(1, detectCores()-1)
N<-4 #number of populations
M<-500 #number of kept samples per iteration

#Prior distribution for continuous model parameters
Pr_cont<-matrix(0,nrow=N+6,ncol=2)
for(j in 1:N){
  Pr_cont[j,]<-c(1,15) #A's
}
Pr_cont[N+1,]<-c(100,3000) #L
Pr_cont[N+2,]<-c(0.5,1) #c
Pr_cont[N+3,]<-c(100,15000) #sigL
Pr_cont[N+4,]<-c(100,15000) #sigR
Pr_cont[N+5,]<-c(1,200) #muL
Pr_cont[N+6,]<-c(1,200) #muR

#stay-probability of the discrete perturbation kernel
stay_prob<-0.9

#--------------------------------------------
# Read reference data and compute their summaries 
#--------------------------------------------

#Time grid
T<-T_window
h<-2*10^-3
grid<-seq(from=0,to=T,by=h)

prep<-prepare_EEG_data(record=record,
                       period=period,
                       summary_path=summary_path,
                       edf_dir=edf_dir,
                       out_root="RefData",
                       T_window=T_window,
                       hsim=h,
                       seizure_index=seizure_index,
                       force=FALSE)      # set TRUE to regenerate from the .edf

X<-prep$X

#sanity checks before a long run
stopifnot(nrow(X)==N, ncol(X)==length(grid))
cat("Data:",record,period,"| channel sd:",round(apply(X,1,sd),3),"\n")

#Determine summary parameters, summaries and summary weights
summaries_parameters<-ABC_summaries_parameters(X,T,h)
summaries<-ABC_summaries(X,T,h,summaries_parameters)
summaries_weights<-ABC_summaries_weights(summaries,summaries_parameters,N)

#--------------------------------------------
#Determine (known) model & network parameters 
#--------------------------------------------

#A's (unknown)
A<-rep(Inf,N)
#B
B<-rep(22.0,N)
#a
a<-rep(100.0,N)
#b
b<-rep(20.0,N)
#C
C<-rep(70.0,N)
#mu
mu<-rep(Inf,N)
#v0
v0<-rep(6.0,N)
#r
r<-rep(0.56,N)
#vmax
vmax<-rep(5.0,N)
#sig
sig<-rep(Inf,N)

Theta<-matrix(0,nrow=N,ncol=9)
Theta[,1]<-A
Theta[,2]<-B
Theta[,3]<-a
Theta[,4]<-b
Theta[,5]<-C
Theta[,6]<-mu
Theta[,7]<-v0
Theta[,8]<-r
Theta[,9]<-vmax

ab_vec<-c()
for(j in 1:N){
  ab_vec<-c(ab_vec,c(a[j],a[j],b[j]))
}
Gamma<-diag(ab_vec,3*N,3*N)
#Sigma is built per proposal, since sigma is unknown

Rho<-matrix(Inf,nrow=N,ncol=N) #all rho's are unknown

L<-Inf #L unkown
c<-Inf #c unknown
K<-matrix(Inf,nrow=N,ncol=N) #all K's are unknown

#starting value X0
startv<-rep(0,6*N) 
#zero-vector for simulation of multivariate normal
meanVec<-rep(0,6*N)

#--------------------------------------------
# exponential matrix and covariance matrix (after Cholesky decomposition) for splitting-simulation
#--------------------------------------------

dm_sim<-exponential_matrix_JRNMM(N,Gamma,h)
#cm is built per proposal inside the ABC functions, since sigma is unknown

#-------------------------------------------------------------------------------
# START SMC-ABC
#-------------------------------------------------------------------------------

#prepare parallel computation
cl<-makeCluster(ncl)
registerDoSNOW(cl)

#--------------------------------------------
# Start pilot run for SMC-ABC
#--------------------------------------------

n_pilot<-10^4
merge_d <- foreach(i=1:n_pilot, .combine='rbind',
                   .packages=c('SplittingJRNMM','mvnfast','expm')) %dopar% {
                     ABC_pilot(N,T,h,grid,startv,Theta,Rho,K,dm_sim,meanVec,Gamma,
                               summaries,summaries_parameters,summaries_weights,Pr_cont)
                   }

#determine delta_1 as the median of the n_pilot calculated distances
delta_1<-median(merge_d) 

#--------------------------------------------
# Start round 1 for SMC-ABC
#--------------------------------------------

start_time<-Sys.time()

r<-1

#carry out for-loop in parallel  
merge_d<-foreach(i=1:M,.combine='rbind',.packages = c('SplittingJRNMM','mvnfast')) %dopar% {
  ABC_SMC_r1(N,T,h,grid,startv,Theta,Rho,K,dm_sim,meanVec,Gamma,summaries,summaries_parameters,summaries_weights,Pr_cont,delta_1)
}

#access the kept particles
dim_merge<-N*(N-1)+N+7
merge_d<-merge_d[,1:dim_merge]

#sort them with respect to the distances
sort_d<-merge_d[order(merge_d[,1]),]

#keep all the sorted values
Dvec<-sort_d[,1]

Amat<-matrix(0,nrow=M,ncol=N)
for(j in 1:N){
  Amat[,j]<-sort_d[,j+1]
}

Lvec<-sort_d[,N+2]
cvec<-sort_d[,N+3]
sig1vec<-sort_d[,N+4]
sig2vec<-sort_d[,N+5]
mu1vec<-sort_d[,N+6]
mu2vec<-sort_d[,N+7]

pmat<-matrix(0,nrow=M,ncol=(N*(N-1)))
for(j in 1:(N*(N-1))){
  pmat[,j]<-sort_d[,j+(N+7)]
}

#Initialize the weights
weights_c<-rep(1,M)

#Normalise the initialized weights
norm_weights_c<-weights_c/sum(weights_c)

#prepare hat_p_vec for the Bernoulli network parameters
hat_p_vec<-rep(0,N*(N-1))
for(k in 1:(N*(N-1))){
  hat_p_vec[k]<-sum(pmat[,k])/M 
} 

#determine covariance matrix hat_Sigma_r from theta_c,kept,r-1
theta_c_kept_r_1<-matrix(0,nrow=M,ncol=(N+6))
for(j in 1:N){
  theta_c_kept_r_1[,j]<-Amat[,j]
}
theta_c_kept_r_1[,N+1]<-Lvec
theta_c_kept_r_1[,N+2]<-cvec
theta_c_kept_r_1[,N+3]<-sig1vec
theta_c_kept_r_1[,N+4]<-sig2vec
theta_c_kept_r_1[,N+5]<-mu1vec
theta_c_kept_r_1[,N+6]<-mu2vec
sigma_kernel<-2*cov.wt(theta_c_kept_r_1,wt=norm_weights_c)$cov

#initialize the acceptance rate
ar<-1

#--------------------------------------------
# Start rounds r>1 for SMC-ABC 
#--------------------------------------------

while(ar>0.001){ #repeat until stopping criterion reached 
  r<-r+1
  
  #determine delta_r
  if(ar<0.01){
    delta_r<-Dvec[[M*0.75]]
  } else {
    delta_r<-median(Dvec) 
  }
  
  #carry out for-loop in parallel
  merge_d<-foreach(i=1:M,.combine='rbind',.packages = c('SplittingJRNMM','mvnfast')) %dopar% {
    ABC_SMC(N,T,h,grid,startv,Theta,Rho,K,dm_sim,meanVec,Gamma,summaries,summaries_parameters,summaries_weights,Amat,Lvec,cvec,sig1vec,sig2vec,mu1vec,mu2vec,norm_weights_c,sigma_kernel,hat_p_vec,stay_prob,Pr_cont,delta_r)
  }
  
  #access the kept particles
  dim_merge<-N*(N-1)+N+8
  merge_d<-merge_d[,1:dim_merge]
  
  #sort them with respect to the distances
  sort_d<-merge_d[order(merge_d[,1]),]
  
  #keep all the sorted values
  Dvec<-sort_d[,1]
  
  Amat_new<-matrix(0,nrow=M,ncol=N)
  for(j in 1:N){
    Amat_new[,j]<-sort_d[,j+1]
  }
  
  Lvec_new<-sort_d[,N+2]
  cvec_new<-sort_d[,N+3]
  sig1vec_new<-sort_d[,N+4]
  sig2vec_new<-sort_d[,N+5]
  mu1vec_new<-sort_d[,N+6]
  mu2vec_new<-sort_d[,N+7]
  
  pmat_new<-matrix(0,nrow=M,ncol=(N*(N-1)))
  for(j in 1:(N*(N-1))){
    pmat_new[,j]<-sort_d[,j+(N+7)]
  }
  
  #determine acceptance rate
  num_sim_dist<-sort_d[,N*(N-1)+N+8]
  ar<-M/sum(num_sim_dist)
  
  #update the weights for the kept samples of that iteration
  for(j in 1:M){
    
    #Multivariate normal perturbation kernel
    sum_c<-0
    x_kernel<-c(Amat_new[j,],Lvec_new[j],cvec_new[j],sig1vec_new[j],sig2vec_new[j],mu1vec_new[j],mu2vec_new[j])
    
    for(l in 1:M){
      mu_kernel<-c(Amat[l,],Lvec[l],cvec[l],sig1vec[l],sig2vec[l],mu1vec[l],mu2vec[l])
      sum_c<-sum_c+norm_weights_c[l]*continuous_Kernel(x_kernel,mu_kernel,sigma_kernel)
    }
    
    #Prior
    prior_val_As<-1
    for(jj in 1:N){
      prior_val_As<-prior_val_As*dunif(Amat_new[j,jj],min=Pr_cont[jj,1],max=Pr_cont[jj,2])
    }
    
    #update weights
    weights_c[j]<-(prior_val_As*dunif(Lvec_new[j],min=Pr_cont[N+1,1],max=Pr_cont[N+1,2])*dunif(cvec_new[j],min=Pr_cont[N+2,1],max=Pr_cont[N+2,2])*dunif(sig1vec_new[j],min=Pr_cont[N+3,1],max=Pr_cont[N+3,2])*dunif(sig2vec_new[j],min=Pr_cont[N+4,1],max=Pr_cont[N+4,2])*dunif(mu1vec_new[j],min=Pr_cont[N+5,1],max=Pr_cont[N+5,2])*dunif(mu2vec_new[j],min=Pr_cont[N+6,1],max=Pr_cont[N+6,2]))/sum_c
  }
  
  #normalize the newly calculated weights
  norm_weights_c<-weights_c/sum(weights_c)
  
  #update Amat, Lvec, cvec and pmat
  Amat<-Amat_new
  Lvec<-Lvec_new
  cvec<-cvec_new
  sig1vec<-sig1vec_new
  sig2vec<-sig2vec_new
  mu1vec<-mu1vec_new
  mu2vec<-mu2vec_new
  pmat<-pmat_new
  
  #update hat_p_vec 
  #estimate hat_p_r^k, k=1,...,dn, from theta_d,kept,r-1
  hat_p_vec<-rep(0,N*(N-1))
  for(k in 1:(N*(N-1))){
    hat_p_vec[k]<-sum(pmat[,k])/M 
  } 
  
  #determine covariance matrix hat_Sigma_r from theta_c,kept,r-1
  theta_c_kept_r_1<-matrix(0,nrow=M,ncol=(N+6))
  for(j in 1:N){
    theta_c_kept_r_1[,j]<-Amat[,j]
  }
  theta_c_kept_r_1[,N+1]<-Lvec
  theta_c_kept_r_1[,N+2]<-cvec
  theta_c_kept_r_1[,N+3]<-sig1vec
  theta_c_kept_r_1[,N+4]<-sig2vec
  theta_c_kept_r_1[,N+5]<-mu1vec
  theta_c_kept_r_1[,N+6]<-mu2vec
  sigma_kernel<-2*cov.wt(theta_c_kept_r_1,wt=norm_weights_c)$cov
  
  #modify some input, according to what is unknown
  Theta[,1]<-rep(Inf,N) #A's are unknown
  Theta[,6]<-rep(Inf,N) #mu's are unknown
  Rho<-matrix(Inf,nrow=N,ncol=N) #all rho's are unknown
  K<-matrix(Inf,nrow=N,ncol=N) #all K's are unknown
  
  #print current iteration and corresponding acceptance rate
  print(r)
  print(ar)
}

#-------------------------------------------------------------------------------
# End SMC-ABC 
#-------------------------------------------------------------------------------

stopCluster(cl)

#--------------------------------------------
# Store kept samples and weights of last iteration 
#--------------------------------------------

#L,c
write(t(Lvec_new), file = paste(filename_results,"/Lvec.txt",sep=""),ncolumns = M,sep = " ")
write(t(cvec_new), file = paste(filename_results,"/cvec.txt",sep=""),ncolumns = M,sep = " ")

#sigmas and mus
write(t(sig1vec_new), file = paste(filename_results,"/sig1vec.txt",sep=""),ncolumns = M,sep = " ")
write(t(sig2vec_new), file = paste(filename_results,"/sig2vec.txt",sep=""),ncolumns = M,sep = " ")
write(t(mu1vec_new), file = paste(filename_results,"/mu1vec.txt",sep=""),ncolumns = M,sep = " ")
write(t(mu2vec_new), file = paste(filename_results,"/mu2vec.txt",sep=""),ncolumns = M,sep = " ")

#A's
for(j in 1:N){
  write(t(Amat_new[,j]), file = paste(filename_results,"/A",j,"vec.txt",sep=""),ncolumns = M,sep = " ")
}

#Rho's
counter<-0
for(j in 1:N){
  for(k in 1:N){
    if(j!=k){
      counter<-counter+1
      write(t(pmat_new[,counter]), file = paste(filename_results,"/p",j,k,"vec.txt",sep = ""),ncolumns = M,sep = " ")
    }
  }
}

#normalized weights  
write(t(norm_weights_c),file = paste(filename_results,"/norm_weights_c.txt",sep = ""),ncolumns = M,sep = " ")