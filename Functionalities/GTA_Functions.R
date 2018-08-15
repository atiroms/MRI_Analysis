#### Description ####

# R script for calculation of metrics of graph theoretical analysis


#### Libraries ####
library(igraph)


#### Weighted Graph Functions ####

iGraph2Nodes<-function(input){
  nodes<-data.frame(ID_long=V(input)$name,label_proper=V(input)$label_proper)
  return(nodes)
}


iGraph2Edges<-function(input){
  edges<-data.frame(get.edgelist(input),as.numeric(E(input)$weight))
  colnames(edges)<-c("from","to","weight")
  return(edges)
}


iGraph2WeightMat<-function(input){
  nodes<-iGraph2Nodes(input)
  edges<-iGraph2Edges(input)
  n_nodes<-nrow(nodes)
  weight<-matrix(0,nrow=n_nodes,ncol=n_nodes)
  for (i in 1:nrow(edges)){
    weight[which(as.character(nodes$ID_long)==edges[i,"from"]),
           which(as.character(nodes$ID_long)==edges[i,"to"])]<-edges[i,"weight"]
    weight[which(as.character(nodes$ID_long)==edges[i,"to"]),
           which(as.character(nodes$ID_long)==edges[i,"from"])]<-edges[i,"weight"]
  }
  colnames(weight)<-rownames(weight)<-nodes$ID_long
  return(weight)
}


iGraph2LengthMat<-function(input){
  nodes<-iGraph2Nodes(input)
  edges<-iGraph2Edges(input)
  edges$length<-1/edges$weight
  n_nodes<-nrow(nodes)
  length<-matrix(Inf,nrow=n_nodes,ncol=n_nodes)
  for (i in 1:nrow(edges)){
    length[which(as.character(nodes$ID_long)==edges[i,"from"]),
           which(as.character(nodes$ID_long)==edges[i,"to"])]<-edges[i,"length"]
    length[which(as.character(nodes$ID_long)==edges[i,"to"]),
           which(as.character(nodes$ID_long)==edges[i,"from"])]<-edges[i,"length"]
  }
  colnames(length)<-rownames(length)<-nodes$ID_long
  return(length)
}


WeightedDistance<-function(input_igraph=NULL,input_length=NULL){
  if (is.null(input_length)){
    length<-iGraph2LengthMat(input_igraph)
  }else{
    length<-input_length
  }
  n_nodes<-ncol(length)  
  distance<-matrix(Inf,ncol=n_nodes,nrow=n_nodes)
  diag(distance)<-0
  path<-matrix(0,ncol=n_nodes,nrow=n_nodes)
  colnames(distance)<-rownames(distance)<-colnames(length)
  colnames(path)<-rownames(path)<-colnames(length)
  for (i in 1:n_nodes){
    unreached_id<-rep(T,n_nodes)
    length_unreached<-length
    closest_remaining<-i
    while (T){
      unreached_id[closest_remaining]<-F
      length_unreached[,closest_remaining]<-0
      for (j in closest_remaining){
        unreached_existing_id<-which(length_unreached[j,]>0)
        for (k in unreached_existing_id){
          old_distance<-distance[i,k]
          new_distance<-distance[i,j]+length_unreached[j,k]
          if (new_distance<old_distance){
            distance[i,k]<-new_distance
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


WeightedCharPath<-function(input_igraph=NULL,input_distancemat=NULL){
  if (is.null(input_distancemat)){
    distance<-WeightedDistance(input_igraph)[[1]]
  }else{
    distance<-input_distancemat
  }
  diag(distance)<-NA
  output<-mean(distance,na.rm=T)
  return(output)
}


WeightedGlobalEfficiency<-function(input_igraph=NULL,input_distancemat=NULL){
  if (is.null(input_distancemat)){
    distance<-WeightedDistance(input_igraph)[[1]]
  }else{
    distance<-input_distancemat
  }
  diag(distance)<-NA
  output<-mean(1/distance,na.rm=T)
  return(output)
}


WeightedEccentricity<-function(input_igraph=NULL,input_distancemat=NULL){
  if (is.null(input_distancemat)){
    distance<-WeightedDistance(input_igraph)[[1]]
  }else{
    distance<-input_distancemat
  }
  diag(distance)<-NA
  output<-data.frame(matrix(ncol=2,nrow=0))
  for (i in 1:nrow(distance)){
    output<-rbind(output,cbind(ID_long=colnames(distance)[i],value=max(distance[,i],na.rm=T)))
  }
  return(output)
}


WeightedRadius<-function(input_igraph=NULL,input_distancemat=NULL){
  if (is.null(input_distancemat)){
    distance<-WeightedDistance(input_igraph)[[1]]
  }else{
    distance<-input_distancemat
  }
  eccentricity<-WeightedEccentricity(input_distancemat=distance)
#  eccentricity$value<-as.numeric(levels(eccentricity$value))[eccentricity$value]
  output<-min(as.numeric.factor(eccentricity$value))
  return(output)
}


WeightedDiameter<-function(input_igraph=NULL,input_distancemat=NULL){
  if (is.null(input_distancemat)){
    distance<-WeightedDistance(input_igraph)[[1]]
  }else{
    distance<-input_distancemat
  }
  eccentricity<-WeightedEccentricity(input_distancemat=distance)
#  eccentricity$value<-as.numeric(levels(eccentricity$value))[eccentricity$value]
  output<-max(as.numeric.factor(eccentricity$value))
  return(output)
}


WeightedClustCoef<-function(input){
  weight<-iGraph2WeightMat(input)
  k<-rowSums(weight!=0)
  cyc<-weight^(1/3)
  cyc3<-diag(cyc %*% cyc %*% cyc)
  k[cyc3==0]<-Inf
  c<-cyc3/(k*(k-1))
  return(c)
}


WeightedTransitivity<-function(input){
  weight<-iGraph2WeightMat(input)
  k<-rowSums(weight!=0)
  cyc<-weight^(1/3)
  cyc3<-diag(cyc %*% cyc %*% cyc)
  k[cyc3==0]<-Inf
  t<-sum(cyc3)/sum((k*(k-1)))
  return(t)
}


WeightedLocalEfficiency<-function(input){
  #algorithm of Wang 2016
  weight<-iGraph2WeightMat(input)
  length<-iGraph2LengthMat(input)
  length[length==Inf]<-0
  adjacency<-weight>0
  ot<-1/3
  n_nodes<-ncol(weight)
  e<-rep(0,n_nodes)
  cbrt_weight<-weight^ot
  cbrt_length<-length^ot
  for (u in 1:n_nodes){
    V<-intersect(which(adjacency[u,]),which(adjacency[,u]))
    sw<-cbrt_weight[u,V]+t(cbrt_weight[V,u])
    di<-WeightedDistance(input_length=cbrt_length[V,V])[[1]]
    di<-1/di
    diag(di)<-0
    se<-di+t(di)
    numer<-(sum(sum((t(sw)%*%sw)*se)))/2
    if (numer!=0){
      sa<-adjacency[u,V]+t(adjacency[V,u])
      denom<-sum(sa)^2-sum(sa^2)
      e[u]<-numer/denom
    }
  }
  return(e)
}


WeightedModularity<-function(input){
  
}


WeightedDegreeDist<-function(input){
  
}


WeightedAssortabilityCoef<-function(input){
  
}


WeightedSmallWorldness<-function(input){
  
}


WeightedStrength<-function(input){
  
}


WeightedClosenessCentrality<-function(input){
  
}

#### Binary Graph Functions ####