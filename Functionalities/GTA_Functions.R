#### Description ####

# R script for calculation of metrics of graph theoretical analysis


#### Libraries ####
library(tidyverse)


#### Functions ####

WeightedCharPathLength<-function(input){
  edges<-data.frame(get.edgelist(input),as.numeric(E(input)$weight))
  colnames(edges)<-c("from","to","weight")
  edges$length<-1/edges$weight
  nodes<-data.frame(ID_long=V(input)$name,label_proper=V(input)$label_proper)
  n_nodes<-nrow(nodes)
  length<-matrix(Inf,nrow=n_nodes,ncol=n_nodes)
  for (i in 1:nrow(edges)){
    length[which(as.character(nodes$ID_long)==edges[i,"from"]),
           which(as.character(nodes$ID_long)==edges[i,"to"])]<-edges[i,"length"]
    length[which(as.character(nodes$ID_long)==edges[i,"to"]),
           which(as.character(nodes$ID_long)==edges[i,"from"])]<-edges[i,"length"]
  }
  distance<-matrix(Inf,ncol=n_nodes,nrow=n_nodes)
  for (i in 1:n_nodes){
    distance[i,i]<-0
  }
  path<-matrix(0,ncol=n_nodes,nrow=n_nodes)
  colnames(length)<-rownames(length)<-nodes$ID_long
  colnames(distance)<-rownames(distance)<-nodes$ID_long
  colnames(path)<-rownames(path)<-nodes$ID_long
  for (i in 1:n_nodes){
    unreached_id<-rep(T,n_nodes)
    length_unreached<-length
    closest_remaining<-i
    while (T){
      unreached_id[closest_remaining]<-F
      length_unreached[,closest_remaining]<-0
      for (j in closest_remaining){
        unreached_existing_id<-which(length_unreached[i,]>0)
        for (k in unreached_existing_id){
          old_distance<-distance[j,k]
          new_distance<-distance[i,j]+length_unreached[j,k]
          if (new_distance<old_distance){
            distance[j,k]<-new_distance
            path[i,k]<-path[i,j]+1
          }
        }
      }
      if (sum(unreached_id)==0){
        break
      }
      min_distance<-min(distance[i,which(unreached_id)])
      if (is.infinite(min_distance)){
        break
      }
      closest_remaining<-which(distance[i,]==min_distance)
    }
  }
  return(list(distance,path))
}