#**************************************************
# Description =====================================
#**************************************************
# R script to analyze relationship between structural/functional connection and clinical data.
# Inputs can be functional correlation from rsfMRI data, or Jackknife esimate of structural covariance from T1 data.


#**************************************************
# Parameters ======================================
#**************************************************
# parameters for fc_corr()
path_exp <- "DropBox/MRI/pnTTC/Puberty/Stats/func_XCP"
#dir_in <- c("13_fc_temp","14_fc_t1w","15_fc_temponly","16_fc_36p_1mm","17_fc_36p_2mm",
dir_in <- c("13_fc_temp","14_fc_t1w","15_fc_temponly","17_fc_36p_2mm",
            "18_fc_36p_native","19_fc_aroma_2mm","20_fc_acompcor_2mm")
dir_out <- "22_fc_corr"
subset_subj <- list(list("column"="W1_5sub","value"=1))


#**************************************************
# Libraries =======================================
#**************************************************
library(ggplot2)
library(GGally)


#**************************************************
# Create path list ================================
#**************************************************
func_path<-function(list_path_root = c("D:/atiroms","C:/Users/atiro","/home/atiroms"),
                    path_exp_=path_exp,
                    dir_in_=dir_in,
                    dir_out_=dir_out){
  path_root<-NA
  for(p in list_path_root){
    if(file.exists(p)){
      path_root<-p
    }
  }
  if(is.na(path_root)){
    print("Error: root path could not be found.")
  }
  path_script <- file.path(path_root,"GitHub/MRI_Analysis")
  path_common <- file.path(path_root,"DropBox/MRI/pnTTC/Puberty/Stats/CommonData")
  path_in     <- file.path(path_root,path_exp_,dir_in_)
  path_out    <- file.path(path_root,path_exp_,dir_out_)
  output <- list("script"=path_script,"input"=path_in,"output"=path_out,
                 "common"=path_common,"dir_in"=dir_in_,"dir_out"=dir_out_)
  return(output)
}

paths<-func_path()


#**************************************************
# Original library ================================
#**************************************************
source(file.path(paths$script,"Functionalities/Functions.R"))
source(file.path(paths$script,"Functionalities/Graphs.R"))


#**************************************************
# FC-FC correlation ===============================
#**************************************************
# for comparison of preprocessing methods

fc_corr<-function(paths_=paths,subset_subj_=subset_subj){
  print("Starting to calculate FC-FC correlation.")
  data_clinical<-func_clinical_data(paths_,subset_subj_)
  nullobj<-func_createdirs(paths_,copy_log=F)
  figs<-list()
  for (id_subj in data_clinical$list_id_subj){
    print(paste("Starting to calculate",as.character(id_subj),sep=" "))
    list_path_file_input<-NULL
    id_study_first<-0
    for (id_study in seq(length(paths_$dir_in))){
      file_input<-paste("fc_",sprintf("%05d", id_subj),"_rp.csv",sep="")
      path_file_input<-file.path(paths_$input[id_study],"output",file_input)
      if (file.exists(path_file_input)){
        list_path_file_input<-c(list_path_file_input,path_file_input)
        if(id_study_first==0){
          id_study_first<-id_study
        }
      }else{
        list_path_file_input<-c(list_path_file_input,NA)
      }
    }
    if(id_study_first==0){
      print("No input available.")
    }
    
    for (id_study in seq(length(paths_$dir_in))){
      path_file_input<-list_path_file_input[id_study]
      if(!is.na(path_file_input)){
        df_fc<-read.csv(path_file_input)
        if(id_study==id_study_first){
          df_fc_allstudy<-data.frame(matrix(ncol=length(paths_$dir_in)+2,nrow=nrow(df_fc)))
          colnames(df_fc_allstudy)<-c("from","to",paths_$dir_in)
        }
        df_fc_allstudy[,c("from","to",paths_$dir_in[id_study])]<-df_fc[,c("from","to","r")]
      }
    }
    
    fig<-ggpairs(df_fc_allstudy[,c(-1,-2)],
                 lower=list(continuous=wrap("points",alpha=0.1, size=0.001)),
                 #lower=list(continuous=wrap("points", size=0.001)),
                 title=paste(sprintf("%05d",id_subj),"fc_corr",sep="_"))
    ggsave(paste(sprintf("%05d",id_subj),"fc_corr.eps",sep="_"),plot=fig,device=cairo_ps,
           path=file.path(paths$output,"output"),dpi=300,height=10,width=10,limitsize=F)
    figs<-c(figs,list(fig))
    print(paste("Finished calculating subject",as.character(id_subj),sep=" "))
  }
  names(figs)<-as.character(data_clinical$list_id_subj)
  print("Finished calculating FC-FC correlation.")
}






#**************************************************
# OBSOLETE ========================================
#**************************************************






#### Parameters ####

#parent_dir <- "D:/atiroms"
#parent_dir <- "C:/Users/atiro"
list_parent_dir <- c("D:/atiroms","C:/Users/atiro")
for(dir in list_parent_dir){
  if(file.exists(dir)){
    parent_dir<-dir
  }
}


script_dir <- file.path(parent_dir,"GitHub/MRI_Analysis")
#input_dir <- file.path(parent_dir,"DropBox/MRI/Statistics/Connection")
input_dir <- file.path(parent_dir,"DropBox/MRI/pnTTC/Puberty/Statistics")
output_dir <- file.path(input_dir,"Connection_data")

connection_file <- "W1_HO_FC.csv"
#connection_file <- "W1_Power_FC.csv"
#connection_file <- "W1_DK_FC.csv"
#connection_file <- "W1_DK_Male_TS1_FC.csv"
#connection_file <- "W1_DK_Male_Subcortex_FC.csv"

#roi_subset<- NULL
#roi_subset<- "cortex"
#roi_subset<- "subcortex"
#roi_subset<- "cerebellum"
#roi_subset<- "global"
#roi_subset<- "misc"
#roi_subset <- c("cortex","subcortex")
roi_subset <- c("cortex","subcortex","cerebellum")

#for Power Atlas
#roi_subset<-c("Uncertain","Default mode","Sensory/somatomotor Hand","Sensory/somatomotor Mouth",
#              "Fronto-parietal Task Control","Cingulo-opercular Task Control","Subcortical",
#              "Salience","Auditory","Visual","Dorsal attention","Ventral attention",
#              "Memory retrieval?","Cerebellar")
#roi_subset<-c("Default mode","Sensory/somatomotor Hand","Sensory/somatomotor Mouth",
#              "Fronto-parietal Task Control","Cingulo-opercular Task Control","Subcortical",
#              "Salience","Auditory","Visual","Dorsal attention","Ventral attention",
#              "Memory retrieval?","Cerebellar")

#subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1)
#subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1, Sex=1)
#subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1, Sex=2)
#subject_subset <- data.frame(W1_T1QC_rsfMRIexist=1, Sex=1,W1_Tanner_Stage=1)
subject_subset <- data.frame(W1_T1QC_rsfMRIexist_CONNvoxelQC20=1, Sex=1)
#subject_subset <- data.frame(W1_T1QC_rsfMRIexist_CONNvoxelQC20=1, Sex=2)
#subject_subset <- data.frame(W1_T1QC_rsfMRIexist_CONNvoxelQC20=1,Sex=1,W1_Tanner_Stage=1)


covariate_label<-c("W1_Tanner_Stage","W1_Age_at_MRI")

p_uncorrected<-0.001
p_corrected<-0.05

n_components<-10
#n_components<-30
#n_components<-5
tsne_dims<-2
tsne_perplexity<-30
tsne_max_itr<-1000

cost<-seq(0.04,0.25,0.01)


#### Libraries ####

#library(multcomp)
library(FactoMineR)
library(ica)
library(Rtsne)
library(tidyverse)
#library(ggpubr)
library(igraph)
library(qgraph)


#### Functionalities ####

source(file.path(script_dir,"Functionalities/Functions.R"))
source(file.path(script_dir,"Functionalities/GLM_Functions.R"))
source(file.path(script_dir,"Functionalities/GTA_Functions.R"))
source(file.path(script_dir,"Functionalities/LI_Functions.R"))
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
connection_data<-connection_data[which(connection_data$flag),
                                 -which(colnames(connection_data)=="flag")]

nodes<-c(as.character(unique(connection_data$from)),
         as.character(unique(connection_data$to)))
nodes<-unique(nodes)
nodes<-nodes[order(nodes)]
nodes<-data.frame(node=nodes,
                  node_label=ConvertID(nodes,roi_data,"ID_long","label_proper"),
                  node_group=ConvertID(nodes,roi_data,"ID_long","group"))

if (!is.null(roi_subset)){
  nodes$include<-F
  for (i in roi_subset){
    nodes[which(nodes$node_group==i),"include"]<-T
  }
  included_nodes<-as.character(nodes[which(nodes$include),"node"])
  connection_data$from_flag<-F
  connection_data$to_flag<-F
  for (i in included_nodes){
    connection_data[which(connection_data$from==i),"from_flag"]<-T
    connection_data[which(connection_data$to==i),"to_flag"]<-T
  }
  connection_data<-connection_data[intersect(which(connection_data$from_flag),
                                             which(connection_data$to_flag)),
                                   -c(which(colnames(connection_data)=="from_flag"),
                                      which(colnames(connection_data)=="to_flag"))]
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

DoGLM_FC<-function(){
  dirname<-ExpDir("GLM_FC")
  connection_data_tidy<-connection_data[,-which(colnames(connection_data)=="p")]
  connection_data_tidy<-rename(connection_data_tidy,value=r)
  glm<-CommonGLM(connection_data_tidy,covariate_label,F,dirname,"GLM_FC.csv")
  
  models_expvars<-glm[intersect(which(glm[,"from"]==glm[1,"from"]),
                                which(glm[,"to"]==glm[1,"to"])),
                      c("model","exp_var")]
  fig<-NULL
  glm_ordered<-NULL
  for (i in 1:nrow(models_expvars)){
    id_obs<-which(glm[,"model"]==models_expvars[i,"model"])
    id_obs<-intersect(id_obs,which(glm[,"exp_var"]==models_expvars[i,"exp_var"]))
    glm_subset<-glm[id_obs,]
    glm_subset<-cbind(glm_subset,MultCompCorr(glm_subset))
    for (j in rois){
      id_obs<-union(which(glm_subset$from==j),which(glm_subset$to==j))
      id_obs<-id_obs[order(id_obs)]
      glm_subsubset<-glm_subset[id_obs,]
      pvalues<-MultCompCorr(glm_subsubset)
      for (k in 1:length(id_obs)){
        for (l in colnames(pvalues)){
          if (is.null(glm_subset[id_obs[k],paste("seed",l,sep="_")])){
            glm_subset[id_obs[k],paste("seed",l,sep="_")]<-pvalues[k,l]
          }else if (is.na(glm_subset[id_obs[k],paste("seed",l,sep="_")])){
            glm_subset[id_obs[k],paste("seed",l,sep="_")]<-pvalues[k,l]
          }else{
            glm_subset[id_obs[k],paste("seed",l,sep="_")]<-min(glm_subset[id_obs[k],paste("seed",l,sep="_")],
                                                               pvalues[k,l])
          }
        }
      }
    }
    nodes_edges<-GLM_FC2Graph(glm_subset,rois)
    fig_title<-paste("GLM Beta of Model:",models_expvars[i,"model"],
                     ", Explanatory Variable:",models_expvars[i,"exp_var"],sep=" ")
    fig<-c(fig,list(CircularPlot(nodes_edges,
                                 pvalue_type="seed_p_Benjamini_Hochberg",
                                 input_title=fig_title)))
#    fig<-c(fig,list(CircularPlot(nodes_edges,pvalue_type="p","GLM Beta Values")))
    glm_ordered<-rbind(glm_ordered,glm_subset)
  }
  write.csv(glm_ordered, file.path(dirname,"GLM_FC_ordered.csv"),row.names=F)
  output<-list(glm_ordered,fig)
  names(output)<-c("GLM","Figures")
  return(output)
}


GLM_FC_Replot<-function(file_path,pvalue_type="seed_p_Benjamini_Hochberg"){
  glm<-read.csv(file_path)
  models_expvars<-glm[intersect(which(glm[,"from"]==glm[1,"from"]),
                                which(glm[,"to"]==glm[1,"to"])),
                      c("model","exp_var")]
  rois<-c(as.character(unique(glm$from)),as.character(unique(glm$to)))
  rois<-unique(rois)
  rois<-rois[order(rois)]
  fig<-NULL
  for (i in 1:nrow(models_expvars)){
    id_obs<-which(glm[,"model"]==models_expvars[i,"model"])
    id_obs<-intersect(id_obs,which(glm[,"exp_var"]==models_expvars[i,"exp_var"]))
    glm_subset<-glm[id_obs,]
    nodes_edges<-GLM_FC2Graph(glm_subset,rois)
    fig_title<-paste("GLM Beta of Model:",models_expvars[i,"model"],
                     ", Explanatory Variable:",models_expvars[i,"exp_var"],sep=" ")
    fig<-c(fig,list(CircularPlot(nodes_edges,
                                 pvalue_type=pvalue_type,
                                 input_title=fig_title)))
  }
  return(fig)
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

Edges2iGraph<-function(input){
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
                   cbind(cost="average",node=as.character(metric_list[i,"node"]),
                         metric=as.character(metric_list[i,"metric"]),
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
    output_add<-cbind(node=names(input$node),
                      node_label=ConvertID(names(input$node),roi_data,"ID_long","label_proper"),
                      metric=input$name[[1]],value=input$node)
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

DoGTA<-function(absolute=T,threshold=NA){
  dirname<-ExpDir("GTA")
#  output_binary<-data.frame()
  for (i in 1:n_subject){
    print(paste("Calculating subject No.",i,", ID_pnTTC:",subject_id[i]))
    Sys.sleep(0.01)
    flush.console()
    subject_graph<-Edges2iGraph(connection_data[which(connection_data$ID_pnTTC==subject_id[i]),])
#    subject_metric_binary<-ItrCost(subject_graph)
#    output_binary<-rbind(output_binary,
#                         cbind(ID_pnTTC=rep(subject_id[i],nrow(subject_metric_binary)),
#                               subject_metric_binary))
    if (absolute){
      E(subject_graph)$weight<-abs(E(subject_graph)$weight)
    }
    subject_metric_weighted<-WeightedMetric(subject_graph)
    output_weighted_add<-cbind(ID_pnTTC=rep(subject_id[i],nrow(subject_metric_weighted)),
                               subject_metric_weighted)
    write.csv(output_weighted_add, file.path(dirname,sprintf("GTA_weighted_%05d.csv",i)),row.names=F)
  }
  output_weighted<-data.frame()
  for (i in 1:n_subject){
    output_weighted_add<-read.csv(file.path(dirname,sprintf("GTA_weighted_%05d.csv",i)))
    output_weighted<-rbind(output_weighted,output_weighted_add)
  }
  
#  write.csv(output_binary, file.path(dirname,"GTA_binary.csv"),row.names=F)
  write.csv(output_weighted, file.path(dirname,"GTA_weighted.csv"),row.names=F)
#  output<-list(output_binary,output_weighted)
#  names(output)<-c("Binary","Weighted")
  glm<-CommonGLM(output_weighted,covariate_label,global_covariate=F,dirname=dirname,"GLM_GTA.csv")
  output_weighted_tidy<-output_weighted[,c("ID_pnTTC","node","metric","value")]
  li<-CommonLI(output_weighted_tidy,"node",dirname,"LI_GTA.csv")
  li_tidy<-li[,c("ID_pnTTC","node","metric","L_ROI_ID","R_ROI_ID","Laterality_Index")]
  colnames(li_tidy)[6]<-"value"
  glm_li<-CommonGLM(li_tidy,covariate_label,global_covariate=F,dirname=dirname,"GLM_LI_GTA.csv")
  output<-list(output_weighted,glm,li,glm_li)
  names(output)<-c("Weighted_GTA","GLM_of_GTA","LI_of_GTA","GLM_of_LI_of_GTA")
  return(output)
}
