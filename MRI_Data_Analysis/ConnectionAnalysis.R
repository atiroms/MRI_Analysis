#### Description ####

# R script to analyze relationship between structural/functional connection and clinical data.
# Inputs can be functional correlation from rsfMRI data, or Jackknife esimate of structural covariance from T1 data.
# 


#### Parameters ####

#parent_dir <- "D:/atiroms"
parent_dir <- "C:/Users/atiro"

script_dir <- file.path(parent_dir,"GitHub/MRI_Analysis")
input_dir <- file.path(parent_dir,"DropBox/MRI/Statistics/Connection")
output_dir <- file.path(input_dir,"Connection_data")

#connection_file <- "W1_HO_FC.csv"
#connection_file <- "W1_Power_FC.csv"
#connection_file <- "W1_DK_FC.csv"
connection_file <- "W1_DK_Male_TS1_FC.csv"

roi_subset<- ""
#roi_subset<- "cortex"
#roi_subset<- "subcortex"
#roi_subset<- "cerebellum"
#roi_subset<- "global"
#roi_subset<- "misc"

#subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1)
#subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1, Sex=1)
#subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1, Sex=2)
subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1, Sex=1,W1_Tanner_Stage=1)


p_uncorrected<-0.001

n_components<-10
#n_components<-30
#n_components<-5
tsne_dims<-2
tsne_perplexity<-30
tsne_max_itr<-1000

cost<-seq(0.04,0.25,0.01)


#### Libraries ####

library(multcomp)
library(FactoMineR)
library(ica)
library(Rtsne)
library(tidyverse)
#library(ggpubr)
library(igraph)
library(qgraph)


#### Functionalities ####

source(file.path(script_dir,"Functionalities/Functions.R"))
source(file.path(script_dir,"Functionalities/GTA_Functions.R"))
source(file.path(script_dir,"Functionalities/Figures.R"))


#### Data Loading ####
source(file.path(script_dir,"Functionalities/LoadClinicalData.R"))
HeatmapPlot(clinical_data,"Clinical Data","Clinical Measure",
            colnames(clinical_data)[-1],scale_data = T)

connection_data <-read.csv(file.path(input_dir,connection_file))
connection_data$flag<-F
for (i in subject_id){
  connection_data[which(connection_data$ID_pnTTC==i),"flag"]<-T
}
connection_data<-connection_data[which(connection_data$flag),-ncol(connection_data)]
if (roi_subset!=""){
  connection_data<-connection_data[which(ConvertID(connection_data$from,roi_data,"ID_long","group")==roi_subset),]
  connection_data<-connection_data[which(ConvertID(connection_data$to,roi_data,"ID_long","group")==roi_subset),]
}

connections<-connection_data[which(connection_data$ID_pnTTC==connection_data[1,"ID_pnTTC"]),2:5]
connections$edge<-paste(connections$from,connections$to,sep="_")
colnames(connections)<-c("from","from_label","to","to_label","edge")
n_connections<-nrow(connections)
rois<-data.frame(label=c(as.character(connections$from),as.character(connections$to)))
rois<-unique(rois$label)
rois<-as.character(rois[order(rois)])
n_rois<-length(rois)


#### GLM Analysis  ####

#not yet checked after update
GLMroutine<-function(input_mri_data,input_covar,id_covar,n_expvar){
  output<-data.frame(matrix(ncol=2+5*n_expvar,nrow=n_connections))
  collabel<-colnames(input_covar)[id_covar+1]
  input_covar<-data.frame(input_covar[,id_covar+1])
  colnames(input_covar)<-collabel
  for (i in 1:n_connections){
    edge_data<-input_mri_data[which(input_mri_data$from==connections[i,"from"]),]
    edge_data<-input_mri_data[which(input_mri_data$to==connections[i,"to"]),"r"]
    if (length(id_covar)==1){
      glmfit<-lm(edge_data~input_covar[,1])
    }else if (length(id_covar)==2){
      glmfit<-lm(edge_data~input_covar[,1]+input_covar[,2])
    }else if (length(id_covar)==3){
      glmfit<-lm(edge_data~input_covar[,1]+input_covar[,2]+input_covar[,3])
    }else if (length(id_covar)==4){
      glmfit<-lm(edge_data~input_covar[,1]+input_covar[,2]+input_covar[,3]+input_covar[,4])
    }
    if (length(id_covar)>=2){
      vifactor<-vif(glmfit)
    }else{
      vifactor<-NA
    }
    stats<-c(AIC(glmfit),BIC(glmfit))
    for (j in 1:n_expvar){
      contrast<-matrix(0L,nrow=1, ncol=length(id_covar)+1)
      contrast[1,j+1]<-1
      ttest <- summary(glht(glmfit, linfct = contrast))$test
      stats <-c(stats, ttest$coefficients[1],ttest$sigma[1],ttest$tstat[1],ttest$pvalues[1],vifactor[j])
    }
    output[i,]<-stats
  }
  collabel<-NULL
  for (j in 1:n_expvar){
    collabel<-c(collabel,paste(colnames(input_covar)[j],c("beta","sigma","t","p","VIF"),sep="_"))
  }
  collabel<-c("AIC","BIC",collabel)
  model_name<-paste(colnames(input_covar),collapse="_")
  collabel<-paste(collabel,"of",model_name,"model",sep="_")
  colnames(output)<-collabel
  return(output)
}

DoGLM_FC<-function(covariates_label=c("W1_Tanner_Stage","W1_Age_at_MRI")){
  dirname<-ExpDir("GLM_FC")
  n_covariates<-length(covariates_label)
  output<-connections
  clinical_data_subset<-clinical_data
  for (i in 1:n_covariates){
    clinical_data_subset<-clinical_data_subset[which(!is.na(clinical_data_subset[,covariates_label[i]])),]
  }
  subject_id_subset<-clinical_data_subset$ID_pnTTC
  
  covariates_data_subset<-data.frame(ID_pnTTC=clinical_data_subset$ID_pnTTC)
  for (i in 1:n_covariates){
    covariates_data_subset<-cbind(covariates_data_subset,clinical_data_subset[,covariates_label[i]])
  }
  for (i in 2:ncol(covariates_data_subset)){
    ave<-mean(covariates_data_subset[,i])
    covariates_data_subset[,i]<-covariates_data_subset[,i]-ave
  }
  colnames(covariates_data_subset)[-1]<-covariates_label
  
  connection_data_subset<-data.frame(matrix(ncol=ncol(connection_data),nrow=0))
  for (i in subject_id_subset){
    connection_data_subset<-rbind(connection_data_subset,connection_data[which(connection_data$ID_pnTTC==i),])
  }
  colnames(connection_data_subset)<-colnames(connection_data)
  
  for (i in n_covariates:1){
    n_expvar<-i
    for (j in 1:dim(combn(n_covariates,i))[2]){
      id_covar<-combn(n_covariates,i)[,j]
      output<-cbind(output,GLMroutine(connection_data_subset, covariates_data_subset,id_covar,n_expvar))
    }
  }
  
  best_model<-data.frame(matrix(ncol=1, nrow=(ncol(connection_data_subset)-1)))
  for (i in c('AIC', 'BIC')){
    xic<-output[, grep(i, names(output))]
    for (j in 1:(ncol(connection_data_subset)-1)){
      best_model[j,1]<-which.min(xic[j,])
    }
    colnames(best_model)<-paste(i,"best_model",sep="_")
    output<-cbind(output,best_model)
  }
  write.csv(output, file.path(dirname,"GLM.csv"),row.names=F)
  return(output)
}


#### Principal Component Analysis ####

# not yet updated
DoPCA_FC<-function(){
  dirname<-ExpDir("PCA_FC")
  data<-connection_data[-1]
  indexcolumn<-connection_data[1]
  data<-data.matrix(data)
  pca <-PCA(data,scale.unit = TRUE, ncp = n_components, graph = FALSE)
  varfactor<-data.frame(pca$var$coord)
  varfactor<-cbind(connections,ConvertID(connections_matrix[,1],roiid_data,input_roiid_type,"label_proper"),ConvertID(connections_matrix[,2],roiid_data,input_roiid_type,"label_proper"),varfactor)
  colnames(varfactor)<-c("ROI_ID","ROI_proper",sprintf("Dim_%02d",1:n_components))
  rownames(varfactor)<-NULL
  indfactor<-data.frame(pca$ind$coord)
  colnames(indfactor)<-sprintf("Dim_%02d",1:n_components)
#  pairs(indfactor)
  indfactor<-cbind(indexcolumn,indfactor)
  varianceaccounted<-pca$eig
  
  write.csv(varfactor, file.path(dirname,"VariableFactor.csv"),row.names=F)
  write.csv(indfactor, file.path(dirname,"IndividualFactor.csv"),row.names=F)
  write.csv(varianceaccounted, file.path(dirname,"VarianceAccounted.csv"))
  ComponentClinicalCorr(indfactor[-1],dirname)
  graph<-ComponentPlot(indfactor[-1],"PCA")
  return(list(varfactor,indfactor,varianceaccounted,graph))
}


#### Independent Component Analysis ####

# not yet updated
DoICA_FC<-function(){
  dirname<-ExpDir("ICA_FC")
  data<-connection_data[-1]
  indexcolumn<-connection_data[1]
  data<-data.matrix(data)
  ica <-icafast(data, nc=n_components,center=TRUE,maxit=100,tol=1e-6,alg="par",fun="logcosh",alpha=1)
  varfactor<-data.frame(ica$M)
  varfactor<-cbind(connections,ConvertID(connections_matrix[,1],roiid_data,input_roiid_type,"label_proper"),ConvertID(connections_matrix[,2],roiid_data,input_roiid_type,"label_proper"),varfactor)
  colnames(varfactor)<-c("ROI_ID","ROI_proper",sprintf("Dim_%02d",1:n_components))
  rownames(varfactor)<-NULL
  indfactor<-data.frame(ica$S)
  colnames(indfactor)<-sprintf("Dim_%02d",1:n_components)
#  pairs(indfactor)
  indfactor<-cbind(indexcolumn,indfactor)
  varianceaccounted<-ica$vafs
  
  write.csv(varfactor, file.path(dirname,"VariableFactor.csv"),row.names=F)
  write.csv(indfactor, file.path(dirname,"IndividualFactor.csv"),row.names=F)
  write.csv(varianceaccounted, file.path(dirname,"VarianceAccounted.csv"))
  ComponentClinicalCorr(indfactor[-1],dirname)
  graph<-ComponentPlot(indfactor[-1],"ICA")
  
  return(list(varfactor,indfactor,varianceaccounted,graph))
}


#### t-SNE Analysis ####

# not yet updated
DoTSNE_FC<-function(){
  dirname<-ExpDir("tSNE_FC")
  indexcolumn<-connection_data[1]
  data<-data.matrix(connection_data[-1])
  tsne_fc<-Rtsne(data, dims = tsne_dims, perplexity=tsne_perplexity, verbose=TRUE, max_iter = tsne_max_itr)
  
  indfactor<-data.frame(tsne_fc$Y)
  colnames(indfactor)<-sprintf("Dim_%02d",1:tsne_dims)
  indfactor<-cbind(indexcolumn,indfactor)
  

  write.csv(indfactor, file.path(dirname,"IndividualFactor.csv"),row.names=F)
  graph<-TwoComponentPlot(indfactor[-1], "t-SNE_FC")
  return(list(tsne_fc,graph))
}


#### Graph object construction ####

Edges2Graph<-function(input){
  edges<-data.frame(matrix(ncol=3,nrow=n_connections))
  edges[,1:2]<-connections[,c("from","to")]
  for (i in 1:n_connections){
    edges[i,3]<-as.numeric(input[intersect(which(input$from==connections[i,"from"]),
                               which(input$to==connections[i,"to"])),"r"])
  }
  colnames(edges)<-c("from","to","weight")
  nodes<-data.frame(id=rois)
  nodes$label_proper<-ConvertID(nodes$id,roi_data,"ID_long","label_proper")
  output <- graph.data.frame(d = edges, vertices = nodes, directed = F)
  return(output)
}


#### Binary Graph Calculation ####

# Subset edges according to desired cost
SubsetEdges<-function(input_graph, input_cost){
  n_edges4cost<-as.integer(n_rois*(n_rois-1)/2*input_cost)
  edges2delete<-head(order(E(input_graph)$weight),(n_connections-n_edges4cost))
  output<-delete.edges(input_graph,edges2delete)
  return(output)
}


# Calculate binary graph metrics
BinaryMetrics<-function(input_graph){
  metrics<-data.frame(matrix(nrow=0,ncol=3))
  colnames(metrics)<-c("node","metric","value")
  ## graph-level metrics
  # characteristic path length
  metrics<-rbind(metrics,cbind(node="graph",metric="characteristic path length",
                               value=average.path.length(input_graph)))
  # global efficiency
  eff<-1/(shortest.paths(input_graph))
  eff[!is.finite(eff)]<-0
  metrics<-rbind(metrics,cbind(node="graph",metric="global efficiency",
                               value=mean(eff,na.rm=TRUE)))
  # global clustering coefficient
  metrics<-rbind(metrics,cbind(node="graph",metric="global clustering coefficient",
                               value=transitivity(input_graph)))
  # average clustering coefficient
  metrics<-rbind(metrics,cbind(node="graph",metric="average clustering coefficient",
                               value=transitivity(input_graph,type="average")))
  # local efficiency
  # modularity
  # small-worldness
  suppressWarnings(metrics<-rbind(metrics,cbind(node="graph",metric="small-world index",
                                                value=smallworldIndex(input_graph)$index)))
  
  ## node-level metrics
  # degree centrality
  metrics<-rbind(metrics,cbind(node=rois,metric="degree centrality",
                               value=centr_degree(input_graph)$res))
  # betweenness centrality
  metrics<-rbind(metrics,cbind(node=rois,metric="betweenness centrality",
                               value=centr_betw(input_graph)$res))
  # eigenvector centrality
  metrics<-rbind(metrics,cbind(node=rois,metric="eigenvector centrality",
                               value=eigen_centrality(input_graph)$vector))
  
  rownames(metrics)<-NULL
  return(metrics)
}

# Iterate over costs
ItrCost<-function(input_graph){
  output<-data.frame()
  for (i in cost){
    subgraph<-SubsetEdges(input_graph,i)
    E(subgraph)$weight<-1
    metrics<-BinaryMetrics(subgraph)
    metrics<-cbind(cost=i,metrics)
    output<-rbind(output,metrics)
  }
  output$value<-as.numeric(output$value)
  metric_list<-output[which(output$cost==cost[1]),c("node","metric")]
  average<-data.frame()
  for (i in 1:nrow(metric_list)){
    average<-rbind(average,
                   cbind(cost="average",node=as.character(metric_list[i,"node"]),metric=as.character(metric_list[i,"metric"]),
                         value=mean(output[intersect(which(output$node==metric_list[i,"node"]),
                                                     which(output$metric==metric_list[i,"metric"])),"value"])))
  }
  output<-rbind(output,average)
  rownames(output)<-NULL
  return(output)
}


#### Weighted Graph Calculations ####

AddMetric<-function(input){
  output<-data.frame(matrix(nrow=0,ncol=4))
  if (!is.null(input$graph)){
    output_add<-cbind(node="graph",node_label=NA,metric=input$name[[1]],value=input$graph)
    output<-rbind(output,output_add)
  }
  if (!is.null(input$node)){
    output_add<-cbind(node=names(input$node),node_label=ConvertID(names(input$node),roi_data,"ID_long","label_proper"),metric=input$name[[1]],value=input$node)
    output<-rbind(output,output_add)
  }
  colnames(output)<-c("node","node_label","metric","value")
  return(output)
}

WeightedMetric<-function(input_igraph){
  metrics<-data.frame(matrix(nrow=0,ncol=4))
  distance<-WeightedDistance(input_igraph)$distance
  
  metrics<-rbind(metrics,AddMetric(WeightedCharPath(input_distance=distance)))
  metrics<-rbind(metrics,AddMetric(WeightedEccentricity(input_distance = distance)))
  metrics<-rbind(metrics,AddMetric(WeightedRadius(input_distance = distance)))
  metrics<-rbind(metrics,AddMetric(WeightedDiameter(input_distance = distance)))
  metrics<-rbind(metrics,AddMetric(WeightedGlobalEfficiency(input_distance = distance)))
  metrics<-rbind(metrics,AddMetric(WeightedClustCoef(input = input_igraph)))
  metrics<-rbind(metrics,AddMetric(WeightedTransitivity(input = input_igraph)))
  metrics<-rbind(metrics,AddMetric(WeightedLocalEfficiency(input = input_igraph)))
  metrics<-rbind(metrics,AddMetric(WeightedModularity(input = input_igraph)))
  metrics<-rbind(metrics,AddMetric(WeightedStrength(input = input_igraph)))
  metrics<-rbind(metrics,AddMetric(WeightedClosenessCentrality(input_distance = distance)))
  metrics<-rbind(metrics,AddMetric(WeightedBetweennessCentrality(input = input_igraph)))
  metrics<-rbind(metrics,AddMetric(WeightedEigenvectorCentrality(input = input_igraph)))
  metrics<-rbind(metrics,AddMetric(WeightedNeighborDegree(input = input_igraph)))
  metrics<-rbind(metrics,AddMetric(WeightedAssortativityCoef(input = input_igraph)))
  
  colnames(metrics)<-c("node","node_label","metric","value")
  rownames(metrics)<-NULL
  return(metrics)
}


#### Graph Theoretical Analysis ####

DoGTA<-function(){
  dirname<-ExpDir("GTA")
#  output_binary<-data.frame()
  output_weighted<-data.frame()
  for (i in 1:n_subject){
    subject_graph<-Edges2Graph(connection_data[which(connection_data$ID_pnTTC==subject_id[i]),])
#    subject_metric_binary<-ItrCost(subject_graph)
#    output_binary<-rbind(output_binary,
#                         cbind(ID_pnTTC=rep(subject_id[i],nrow(subject_metric_binary)),
#                               subject_metric_binary))
    E(subject_graph)$weight<-abs(E(subject_graph)$weight)
    subject_metric_weighted<-WeightedMetric(subject_graph)
    output_weighted<-rbind(output_weighted,
                           cbind(ID_pnTTC=rep(subject_id[i],nrow(subject_metric_weighted)),
                                 subject_metric_weighted))
  }
#  write.csv(output_binary, file.path(dirname,"GTA_binary.csv"),row.names=F)
  write.csv(output_weighted, file.path(dirname,"GTA_weighted.csv"),row.names=F)
#  output<-list(output_binary,output_weighted)
#  names(output)<-c("Binary","Weighted")
  output<-output_weighted
  return(output)
}
