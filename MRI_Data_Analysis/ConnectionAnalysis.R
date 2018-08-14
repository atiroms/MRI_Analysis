#### Description ####

# R script to analyze relationship between structural/functional connection and clinical data.
# Inputs can be functional correlation from rsfMRI data, or Jackknife esimate of structural covariance from T1 data.
# 


#### Parameters ####

parent_dir <- "D:/atiroms"
#parent_dir <- "C:/Users/atiro"

script_dir <- file.path(parent_dir,"GitHub/MRI_Analysis")
input_dir <- file.path(parent_dir,"DropBox/MRI/Statistics/Connection")
output_dir <- file.path(input_dir,"Connection_data")

connection_file <- "W1_HO_FC.csv"
#connection_file <- "W1_HO_Power.csv"
#connection_file <- "W1_HO_DK.csv"

roi_subset<- ""
#roi_subset<- "cortex"
#roi_subset<- "subcortex"
#roi_subset<- "cerebellum"
#roi_subset<- "global"
#roi_subset<- "misc"

subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1)
#subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1, Sex=1)
#subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1, Sex=2)
#subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1, Sex=1,W1_Tanner_Stage=1)


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

connections<-connection_data[which(connection_data$ID_pnTTC==connection_data[1,ID_pnTTC]),2:5]
connections$edge<-paste(connections$from,connections$to,sep="_")
colnames(connections)<-c("from","from_label","to","to_label","edge")
n_connections<-nrow(connections)
roi<-data.frame(label=c(connections$from,connections$to))
roi<-as.character(unique(roi$label))
roi<-roi[order(roi)]
n_rois<-length(roi)


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


#### Graph object creation ####

Edges2Graph<-function(input){
  edges<-data.frame(matrix(ncol=3,nrow=n_connections))
  edges[,1:2]<-connections[,c("from","to")]
  for (i in 1:n_connections){
    edges[,3]<-as.numeric(input[intersect(which(input$from==connections[i,"from"]),
                               which(input$to==connections[i,"to"])),"r"])
  }
  colnames(edges)<-c("from","to","weight")
  nodes<-data.frame(id=roi)
  nodes$label_proper<-ConvertID(nodes$id,roi_data,"ID_long","label_proper")
  output <- graph.data.frame(d = edges, vertices = nodes, directed = F)
  return(output)
}


#### Subset edges according to desired cost ####

SubsetEdges<-function(input_graph, input_cost){
  n_nodes<-length(V(input_graph))
  n_edges<-n_nodes*(n_nodes-1)/2
  n_edges4cost<-as.integer(n_nodes*(n_nodes-1)/2*input_cost)
  edges2delete<-head(order(E(input_graph)$weight),(n_edges-n_edges4cost))
  output<-delete.edges(input_graph,edges2delete)
  return(output)
}


#### Iterate over costs ####

ItrCost<-function(input_graph, input_cost=cost){
  n_nodes<-length(V(input_graph))
  output<-data.frame()
  for (i in input_cost){
    subgraph<-SubsetEdges(input_graph,i)
    E(subgraph)$weight<-1
    ## graph-level metrics
    metric<-data.frame(matrix(nrow=1,ncol=0))
    # characteristic path length
    metric<-cbind(metric,averagepathlength=average.path.length(subgraph))
    # global efficiency
    eff<-1/(shortest.paths(subgraph))
    eff[!is.finite(eff)]<-0
    metric<-cbind(metric,globalefficiency=mean(eff,na.rm=TRUE))
    # global clustering coefficient
    metric<-cbind(metric,globalclusteringcoefficient=transitivity(subgraph))
    # average clustering coefficient
    metric<-cbind(metric,averageclusteringcoefficient=transitivity(subgraph,type="average"))
    # local efficiency
    # modularity
    # small-worldness
    suppressWarnings(metric<-cbind(metric,smallworldindex=smallworldIndex(subgraph)$index))
    
    ## node-level metrics
    # degree centrality
    nodemetric<-data.frame(matrix(ncol=n_nodes,nrow=1))
    nodemetric[1,]<-centr_degree(subgraph)$res
    colnames(nodemetric)<-paste("degreecentrality",regions,sep="_")
    metric<-cbind(metric,nodemetric)
    # betweenness centrality
    nodemetric<-data.frame(matrix(ncol=n_nodes,nrow=1))
    nodemetric[1,]<-centr_betw(subgraph)$res
    colnames(nodemetric)<-paste("betweennesscentrality",regions,sep="_")
    metric<-cbind(metric,nodemetric)
    # eigenvector centrality
    nodemetric<-data.frame(matrix(ncol=n_nodes,nrow=1))
    nodemetric[1,]<-eigen_centrality(subgraph)$vector
    colnames(nodemetric)<-paste("eigenvectorcentrality",regions,sep="_")
    metric<-cbind(metric,nodemetric)
    
    output<-rbind(output,metric)
  }
  output<-data.frame(cost=input_cost,output)
  
  average<-data.frame(matrix(ncol=ncol(output),nrow=1))
  colnames(average)<-colnames(output)
  for (i in 2:ncol(output)){
    average[1,i]<-mean(output[,i])
  }
  output<-rbind(average,output)
  return(output)
}


#### Weighted Graph Theory Metric Calculation ####
WeightedMetric<-function(){
  
}


#### Graph Theoretical Analysis, subject-wise ####

DoGTA<-function(){
  dirname<-ExpDir("GTA")
  output_unweighted<-data.frame()
  output_weighted<-data.frame()
  for (i in 1:n_subject){
    subject_graph<-Edges2Graph(connection_data[which(connection_data$ID_pnTTC==subject_id[i]),])
    subject_metric_unweighted<-ItrCost(subject_graph)
    output_unweighted<-rbind(output_unweighted,subject_metric_unweighted[1,-1])
    subject_metric_weighted<-WeightedMetric(subject_graph)
    output_weighted<-rbind(output_weighted,wubject_metric_weighted[1,-1])
  }
  output_unweighted<-data.frame(ID_pnTTC=subject_id,output_unweighted)
  output_weighted<-data.frame(ID_pnTTC=subject_id,output_weighted)
  write.csv(output_unweighted, file.path(dirname,"GTA_unweighted.csv"),row.names=F)
  write.csv(output_weighted, file.path(dirname,"GTA_weighted.csv"),row.names=F)
  return(output)
}
