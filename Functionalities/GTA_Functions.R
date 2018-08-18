#### Description ####

# R script for calculation of metrics of graph theoretical analysis


#### Libraries ####
library(igraph)


#### Weighted Graph Functions ####


#### Basic functions ####
iGraph2Nodes<-function(input){
  nodes<-data.frame(node=V(input)$name,label_proper=V(input)$label_proper)
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
    weight[which(as.character(nodes$node)==edges[i,"from"]),
           which(as.character(nodes$node)==edges[i,"to"])]<-edges[i,"weight"]
    weight[which(as.character(nodes$node)==edges[i,"to"]),
           which(as.character(nodes$node)==edges[i,"from"])]<-edges[i,"weight"]
  }
  colnames(weight)<-rownames(weight)<-nodes$node
  return(weight)
}

iGraph2LengthMat<-function(input){
  nodes<-iGraph2Nodes(input)
  edges<-iGraph2Edges(input)
  edges$length<-1/edges$weight
  n_nodes<-nrow(nodes)
  length<-matrix(Inf,nrow=n_nodes,ncol=n_nodes)
  for (i in 1:nrow(edges)){
    length[which(as.character(nodes$node)==edges[i,"from"]),
           which(as.character(nodes$node)==edges[i,"to"])]<-edges[i,"length"]
    length[which(as.character(nodes$node)==edges[i,"to"]),
           which(as.character(nodes$node)==edges[i,"from"])]<-edges[i,"length"]
  }
  colnames(length)<-rownames(length)<-nodes$node
  return(length)
}


#### Basic Measures ####

# Strength / Weighted Degree
WeightedStrength<-function(input){
  weight<-iGraph2WeightMat(input)
  strength<-rowSums(weight)
  output<-list(NULL,strength)
  names(output)<-c("graph","node")
  return(output)
}

# Distance / Shortest Path Length
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
  output<-list(distance,path)
  names(output)<-c("distance","binary_shortest_path_length")
  return(output)
}

# Weighted Geometric Means of Triangles


#### Measures of Integration ####

# Characteristic Path Length / Average Path Length
WeightedCharPath<-function(input_igraph=NULL,input_distancemat=NULL){
  if (is.null(input_distancemat)){
    distance<-WeightedDistance(input_igraph)[[1]]
  }else{
    distance<-input_distancemat
  }
  diag(distance)<-NA
  charpath<-mean(distance,na.rm=T)
  output<-list(charpath,NULL)
  names(output)<-c("graph","node")
  return(output)
}

# (not in Rubinov NeuroImage 2010)
WeightedEccentricity<-function(input_igraph=NULL,input_distancemat=NULL){
  if (is.null(input_distancemat)){
    distance<-WeightedDistance(input_igraph)[[1]]
  }else{
    distance<-input_distancemat
  }
  diag(distance)<-NA
  eccentricity<-NULL
  for (i in 1:nrow(distance)){
    eccentricity<-c(eccentricity,max(distance[i,],na.rm=T))
  }
  names(eccentricity)<-rownames(distance)
  output<-list(NULL,eccentricity)
  names(output)<-c("graph","node")
  return(output)
}

# (not in Rubinov NeuroImage 2010)
WeightedRadius<-function(input_igraph=NULL,input_distancemat=NULL){
  if (is.null(input_distancemat)){
    distance<-WeightedDistance(input_igraph)[[1]]
  }else{
    distance<-input_distancemat
  }
  eccentricity<-WeightedEccentricity(input_distancemat=distance)[[2]]
  #  eccentricity$value<-as.numeric(levels(eccentricity$value))[eccentricity$value]
  radius<-min(eccentricity)
  output<-list(radius,NULL)
  names(output)<-c("graph","node")
  return(output)
}

# (not in Rubinov NeuroImage 2010)
WeightedDiameter<-function(input_igraph=NULL,input_distancemat=NULL){
  if (is.null(input_distancemat)){
    distance<-WeightedDistance(input_igraph)[[1]]
  }else{
    distance<-input_distancemat
  }
  eccentricity<-WeightedEccentricity(input_distancemat=distance)[[2]]
  #  eccentricity$value<-as.numeric(levels(eccentricity$value))[eccentricity$value]
  dianeter<-max(eccentricity)
  output<-list(diameter,NULL)
  names(output)<-c("graph","node")
  return(output)
}

# Global Efficiency
WeightedGlobalEfficiency<-function(input_igraph=NULL,input_distancemat=NULL){
  if (is.null(input_distancemat)){
    distance<-WeightedDistance(input_igraph)[[1]]
  }else{
    distance<-input_distancemat
  }
  diag(distance)<-NA
  efficiency<-rowMeans(1/distance,na.rm=T)
  ave<-mean(efficiency)
  output<-list(ave,efficiency)
  names(output)<-c("graph","node")
  return(output)
}


#### Measures of Segregation ####

# Weighted clustering coefficient
WeightedClustCoef<-function(input){
  weight<-iGraph2WeightMat(input)
  k<-rowSums(weight!=0)
  cyc<-weight^(1/3)
  cyc3<-diag(cyc %*% cyc %*% cyc)
  k[cyc3==0]<-Inf
  c<-cyc3/(k*(k-1))
  ave<-mean(c)
  output<-list(ave,c)
  names(output)<-c("graph","node")
  return(output)
}


# Weighted Transitivity
WeightedTransitivity<-function(input){
  weight<-iGraph2WeightMat(input)
  k<-rowSums(weight!=0)
  cyc<-weight^(1/3)
  cyc3<-diag(cyc %*% cyc %*% cyc)
  k[cyc3==0]<-Inf
  t<-sum(cyc3)/sum((k*(k-1)))
  output<-list(t,NULL)
  names(output)<-c("graph","node")
  return(output)
}

# Weighted local efficiency
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
  names(e)<-rownames(weight)
  ave<-mean(e)
  output<-list(ave,e)
  names(output)<-c("graph","node")
  return(output)
}


# Weighted modularity
WeightedModularity<-function(input,gamma_v=1){
  # Newman 2006 algorithm
  weight<-iGraph2WeightMat(input)
  n_nodes<-ncol(weight)
  K<-colSums(weight)
  m<-sum(K)
  B<-as.matrix(weight-gamma_v*(K%*%t(K))/m)
  Ci<-rep(1,n_nodes)
  cn<-1
  U<-c(1,0)
  ind<-1:n_nodes
  Bg<-B
  Ng<-n_nodes
  while (U[1]){
    eig<-eigen(Bg)
    D<-eig$values
    V<-eig$vectors
    i1<-which.max(D)
    v1<-V[,i1]
    #    v1<--V[,i1]
    S<-rep(1,Ng)
    S[v1<0]<--1
    q<-t(S) %*% Bg %*% S
    if (q>1e-10){
      qmax<-q
      diag(Bg)<-0
      indg<-rep(T,Ng)
      Sit<-S
      while (sum(indg,na.rm=T)!=0){
        Qit<-rep(qmax,length(Sit))-4*Sit*(Bg %*% Sit)
        imax<-which.max(Qit*indg)
        qmax<-max(Qit*indg,na.rm=T)
        Sit[imax]<--Sit[imax]
        indg[imax]<-NaN
        if (qmax>q){
          q<-qmax
          S<-Sit
        }
      }
      
      if (abs(sum(S))==Ng){
        U<-U[-1]
      }else{
        cn<-cn+1
        Ci[ind[S==1]]<-U[1]
        Ci[ind[S==-1]]<-cn
        U<-c(cn,U)
      }
    }else{
      U<-U[-1]
    }
    ind<-which(Ci==U[1])
    bg<-B[ind,ind]
    Bg<-bg-diag(rowSums(bg))
    Ng<-length(ind)
  }
  names(Ci)<-rownames(weight)
  s<-replicate(n_nodes,Ci)
  Q<-(!(s-t(s)))*B/m
  Q<-sum(Q)
  output<-list(Q,Ci)
  names(output)<-c("graph","community_of_node")
  return(output)
}

# Community Structure using Louvain algorithm (not in Rubinov NeuroImage 2010)

# Community structure using link-based algorithm (not in Rubinov NeuroImage 2010)

# Community structure using clique-percolation (not in Rubinov NeuroImage 2010)


#### Means of centrality ####

# Weighted Closeness Centrality
WeightedClosenessCentrality<-function(input){
  
}

# Weighted Betweenness Centrality
WeightedBetweennessCentrality<-function(input){
  
}

# Weighted Eigenvector Centrality (not in Rubinov NeuroImage 2010)
WeightedEigenvectorCentrality<-function(input){
  
}

# Weighted Within-module degree z-score

# Weighted Participation coefficient


#### Network Motifs ####
# Anatomical and functional motifs

# Motif z-score

# motif fingerprint

#### Measures of resilience ####

# Cumulative weighted degree distribution
CumWeightedDegreeDist<-function(input){
  
}

# Average Weighted Neighbor Degree
AverageWeightedNeighborDegree<-function(input){
  
}

# Weighted Assortativity Coefficient
WeightedAssortativityCoef<-function(input){
  weight<-iGraph2WeightMat(input)
  strength<-rowSums(weight)
  id<-which(upper.tri(weight)*weight>0,arr.ind = T)
  K<-nrow(id)
  stri<-strength[id[,1]]
  strj<-strength[id[,2]]
  r<-(sum(stri*strj)/K-(sum(0.5*(stri+strj))/K)^2)/(sum(0.5*(stri^2+strj^2))/K-(sum(0.5*(stri+strj))/K)^2)
  output<-list(r,NULL)
  names(output)<-c("graph","node")
  return(output)
}

#### Other Concepts ####

# Weighted network Small-worldness
WeightedNetworkSmallWorldness<-function(input){
  
}