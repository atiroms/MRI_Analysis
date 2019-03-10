#**************************************************
# Description =====================================
#**************************************************
# R script to analyze relationship between structural/functional connection and clinical data.
# Inputs can be functional correlation from rsfMRI data, or Jackknife esimate of structural covariance from T1 data.


#**************************************************
# Parameters ======================================
#**************************************************
# parameters for glm_fc()
list_covar<-c("W1_Tanner_Max","W1_Age_at_MRI")

# parameters for fc_corr()
#path_exp <- "DropBox/MRI/pnTTC/Puberty/Stats/func_XCP"
path_exp <- "DropBox/MRI/pnTTC/Puberty/Stats/func_XCP/test_5sub"
#dir_in <- c("13_fc_temp","14_fc_t1w","15_fc_temponly","17_fc_36p_2mm",
#            "18_fc_36p_native","19_fc_aroma_2mm","20_fc_acompcor_2mm")
#dir_in <- c("17_fc_36p_2mm","19_fc_aroma_2mm","20_fc_acompcor_2mm")
#dir_out <- "22_fc_corr"
#dir_out <- "23_fc_corr_heatmap"
#dir_in <- "25_fc_acompcor_2mm"
#dir_out <- "26_glm_fc_acompcor"
#dir_out <- "27_glm_fc_acompcor2"
dir_in<-"28_fc_acompcor"
dir_out<-"29_glm_fc_acompcor"

#dir_in <- "07_fc_acompcor"
#dir_out <- "11_glm_fc_acompcor"
subset_subj <- list(list("column"="W1_5sub","value"=1))
#subset_subj <- list(list("column"="W1_5sub","value"=1),list("column"="Sex","value"=1))

list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")
#list_atlas<-c("aal116")

thr_pvalue <- 0.05


#**************************************************
# Libraries =======================================
#**************************************************
library(ggplot2)
library(GGally)
library(igraph)
library(qgraph)


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
source(file.path(paths$script,"functionality/function.R"))
source(file.path(paths$script,"functionality/glm_function.R"))
source(file.path(paths$script,"functionality/graph.R"))
source(file.path(paths$script,"functionality/gta_function.R"))



#**************************************************
# GTA functionalities =============================
#**************************************************

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


# Subset edges according to desired cost
SubsetEdges<-function(input_graph, input_cost){
  n_edges4cost<-as.integer(n_rois*(n_rois-1)/2*input_cost)
  edges2delete<-head(order(E(input_graph)$weight),(n_connections-n_edges4cost))
  output<-delete.edges(input_graph,edges2delete)
  return(output)
}


# add metric to output list in weighted GTA
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


#**************************************************
# Binary GTA ======================================
#**************************************************

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


#**************************************************
# Weighted GTA ======================================
#**************************************************

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

gta<-function(absolute=T,threshold=NA){
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



#**************************************************
# General linear model of FCs =====================
#**************************************************
glm_fc<-function(paths_=paths,subset_subj_=subset_subj,list_covar_=list_covar,
                 thr_pvalue_=thr_pvalue,list_atlas_=list_atlas){
  print("Starting glm_fc()")
  data_clinical<-func_clinical_data(paths_,subset_subj_)
  nullobj<-func_createdirs(paths_,copy_log=T)
  
  for (atlas in list_atlas_){
    print(paste("  Starting to calculate for atlas: ",atlas,sep=""))
    df_fc<-read.csv(file.path(paths_$input,"output",paste("fc_",atlas,".csv",sep="")))
    # Convert NaN's to zero, delete p column and change column name r to value
    df_fc$r[which(is.nan(df_fc$r))]<-0
    df_fc<-df_fc[,-which(colnames(df_fc)=="p")]
    colnames(df_fc)[colnames(df_fc)=="r"]<-"value"
    # Calculate GLM
    print("    Starting to calculate GLM of FCs.")
    data_glm<-func_glm(df_mri=df_fc,data_clinical,list_covar=list_covar_)
    print("    Finished calculating GLM of FCs.")
    print("    Starting to save GLM of FCs.")
    write.csv(data_glm$glm,file.path(paths_$output,"output",paste("glm_",atlas,".csv",sep="")),row.names = F)
    write.csv(data_glm$ic,file.path(paths_$output,"output",paste("ic_",atlas,".csv",sep="")),row.names = F)
    write.csv(data_glm$min_ic,file.path(paths_$output,"output",paste("min_ic_",atlas,".csv",sep="")),row.names = F)
    write.csv(data_glm$vif,file.path(paths_$output,"output",paste("vif_",atlas,".csv",sep="")),row.names = F)
    print("    Finished saving GLM of FCs.")
    
    print("    Starting to calculate seed-level multiple comparison correction and graphs.")
    
    list_roi <- sort(unique(c(as.character(unique(df_fc$from)),as.character(unique(df_fc$to)))))
    dict_roi <- func_dict_roi(paths_)
    dict_roi <- dict_roi[is.element(dict_roi$id,list_roi),]
    
    for (i in seq(length(data_glm$list_model))){
      model<-data_glm$list_model[[i]]
      for (var_exp in model){
        name_model<-paste(list_covar[model],collapse="_")
        name_var_exp<-list_covar[var_exp]
        print(paste("      Starting to calculate for model: ",name_model,", expvar: ",name_var_exp,sep=""))
        df_glm_subset<-data_glm$glm[intersect(which(data_glm$glm[,"model"]==name_model),
                                              which(data_glm$glm[,"var_exp"]==name_var_exp)),]
        # Add columns for global-level multiple comparison-corrected p values
        df_glm_subset<-cbind(df_glm_subset,mltcomp_corr(df_glm_subset))
        
        for (j in list_roi){
          id_obs<-sort(union(which(df_glm_subset$from==j),which(df_glm_subset$to==j)))
          df_glm_subsubset<-df_glm_subset[id_obs,]  # subset of df_glm_subset, starts or ends ad ROI j
          pvalues<-mltcomp_corr(df_glm_subsubset)  # multiple comparison-corrected p values
          for (k in seq(length(id_obs))){  # iterate over connections which starts / ends at ROI j
            for (l in colnames(pvalues)){
              if (is.null(df_glm_subset[id_obs[k],paste("seed",l,sep="_")])){
                df_glm_subset[id_obs[k],paste("seed",l,sep="_")]<-pvalues[k,l]
              }else if (is.na(df_glm_subset[id_obs[k],paste("seed",l,sep="_")])){
                df_glm_subset[id_obs[k],paste("seed",l,sep="_")]<-pvalues[k,l]
              }else{
                df_glm_subset[id_obs[k],paste("seed",l,sep="_")]<-min(df_glm_subset[id_obs[k],paste("seed",l,sep="_")],
                                                                      pvalues[k,l])
              }
            }
          }
        }
        
        
        # For each model / expvar pair, convert df into nodes/edges data
        graph<-glm_fc2graph(df_glm_subset,list_roi)
        for (j in 1:nrow(graph$node)){
          graph$node[j,"label"]<-as.character(dict_roi[which(dict_roi$id==graph$node[j,"label"]),"label"])
        }
        
        fig_circular<-graph_circular(input=graph,type_pvalue="seed_p_benjamini_hochberg",thr_pvalue=thr_pvalue_)
        fig_circular<-fig_circular +
          ggtitle(paste("GLM Beta\nModel: ",name_model,"\nExplanatory Variable: ",name_var_exp,sep=" ")) +
          theme(plot.title = element_text(hjust = 0.5))
        
        ggsave(paste("glm_graph_",atlas,"_model-",name_model,"_expvar-",name_var_exp,".eps",sep=""),
               plot=fig_circular,device=cairo_ps,path=file.path(paths_$output,"output"),
               dpi=300,height=10,width=10,limitsize=F)
        
        print(paste("      Finished calculating for model: ",name_model,", expvar: ",name_var_exp,sep=""))
      }
    }
    print("    Finished calculating seed-level multiple comparison correction and graphs.")
    print(paste("  Finished calculating for atlas: ",atlas,sep=""))
  }
  
  print("Finished all calculations of glm_fc()")
}



#**************************************************
# FC-FC correlation ===============================
#**************************************************
# for comparison of preprocessing methods

fc_corr<-function(paths_=paths,subset_subj_=subset_subj){
  print("Starting to calculate FC-FC correlation.")
  df_clinical<-func_clinical_data(paths_,subset_subj_)
  nullobj<-func_createdirs(paths_,copy_log=F)
  figs<-list()
  for (id_subj in df_clinical$list_id_subj){
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
                 upper=list(continuous=custom_corr_heatmap),
                 #lower=list(continuous=wrap("points",alpha=0.01,size=0.001,stroke = 0, shape = ".")),
                 lower=list(continuous=custom_smooth),
                 diag=list(continuous=custom_densityDiag),
                 title=paste(sprintf("%05d",id_subj),"fc_corr",sep="_"))
    ggsave(paste(sprintf("%05d",id_subj),"fc_corr.eps",sep="_"),plot=fig,device=cairo_ps,
           path=file.path(paths$output,"output"),dpi=300,height=10,width=10,limitsize=F)
    figs<-c(figs,list(fig))
    print(paste("Finished calculating subject",as.character(id_subj),sep=" "))
  }
  names(figs)<-as.character(df_clinical$list_id_subj)
  print("Finished calculating FC-FC correlation.")
  return(figs)
}

