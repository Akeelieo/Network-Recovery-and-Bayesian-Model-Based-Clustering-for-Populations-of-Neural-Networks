#function that gives label/names to elements of vectorized adjacency given the SBM membership of nodes
#list format for arguments
lab_adj_list<-function(vect_adj_list,memb_list,posit){#posit: positions of elements in upper triangle of adjacency
  for (j in 1:length(vect_adj_list)){
    lab1<-posit
    map1<-setNames(memb_list[[j]],seq(1,length(memb_list[[1]]),1)) #label the membership by the label of the node
    lab1[]<-map1[lab1] #substitute membership of node within set of labeled pairs of nodes (lab)
    ind<-which(lab1[,"row"]>lab1[,"col"]) #identify the pairs where i>j and reverse the sequence
    lab1[ind, ]<- c(lab1[ind,"col"],lab1[ind,"row"])
    names(vect_adj_list[[j]])<-paste(lab1[,"row"],lab1[,"col"],sep = "")
  }
  return(vect_adj_list)
}

#vector format for arguments
lab_adj<-function(vect_adj,memb,posit){#posit: positions of elements in upper triangle of adjacency
  map1<-setNames(memb,seq(1,length(memb),1))
  posit[]<-map1[posit]
  ind<-which(posit[,"row"]>posit[,"col"]) #identify the pairs where i>j and reverse the sequence
  posit[ind, ]<- c(posit[ind,"col"],posit[ind,"row"])
  names(vect_adj)<-paste(posit[,"row"],posit[,"col"],sep = "")
  return(vect_adj)
}
