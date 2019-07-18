#**************************************************
# Description =====================================
#**************************************************
# R script to analyze relationship between structural/functional connection and clinical data.
# Inputs can be functional correlation from rsfMRI data, or Jackknife esimate of structural covariance from T1 data.


#**************************************************
# Parameters ======================================
#**************************************************
# parameters for gta_bin() and gta_weight()
path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"
#path_exp <- "Dropbox/MRI/pnTTC/Puberty/Stats/func_XCP/test_5sub"

dir_in<-"54_fc_acompcor"
#dir_out<-"55_gta_bin"
dir_out<-"55_fingerprint"
#dir_out<-"59_pca_fc"

list_wave <- c(1,2)

subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
                             list("key"="W1_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)),
                    "2"=list(list("key"="W2_T1QC","value"=1),
                             list("key"="W2_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)))

list_covar<-list("tanner"=list("1"="W1_Tanner_Max",
                               "2"="W2_Tanner_Max",
                               "label"="Tanner stage"),
                 "age"=list("1"="W1_Age_at_MRI",
                            "2"="W2_Age_at_MRI",
                            "label"="Age"),
                 "sex"=list("1"="Sex",
                            "2"="Sex",
                            "label"="Sex"))

#list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")
list_atlas<-"aal116"
#list_atlas<-"schaefer400"
#list_atlas<-c("glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")
#thr_pvalue <- 0.05

list_cost<-seq(0.15,0.40,0.01)
absolute<-T
threshold<-NA


#**************************************************
# Libraries =======================================
#**************************************************
library(ggplot2)
library(GGally)
library(igraph)
library(qgraph)
library(FactoMineR)
library(missMDA)
library(ggrepel)
library(colorRamps)
library(tidyverse)


#**************************************************
# Create path list ================================
#**************************************************
func_path<-function(list_path_root = c("D:/atiroms","C:/Users/atiro","/home/atiroms","C:/Users/NICT_WS"),
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
  path_common <- file.path(path_root,"DropBox/MRI_img/pnTTC/puberty/common")
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
# Principal component analysis of FC ==============
#**************************************************
pca_fc<-function(paths_=paths,
                 list_atlas_=list_atlas,
                 list_wave_=list_wave,
                 list_covar_=list_covar,
                 subset_subj_=subset_subj){
  print("Starting pca_fc().")
  nullobj<-func_createdirs(paths_)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar_,rem_na_clin=F)
  df_clin<-data_clin$df_clin
  colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
  
  for (atlas in list_atlas_){
    # Load and examine FC data
    print(paste("Loding FC of atlas: ",atlas,sep=""))
    df_conn<-read.csv(file.path(paths_$input,"output",paste("atl-",atlas,"_fc.csv",sep="")))
    df_edge<-df_conn[which(df_conn$ID_pnTTC==df_conn[1,"ID_pnTTC"]),]
    df_edge<-df_edge[which(df_edge$ses==df_edge[1,"ses"]),c("from","to")]
    n_edge<-dim(df_edge)[1]
    list_node<-sort(unique(c(as.character(unique(df_edge$from)),as.character(unique(df_edge$to)))))
    n_node<-length(list_node)
    
    # Create list of subjects who meet subsetting condition and whose MRI data exist
    list_ses_exist <- sort(unique(df_conn$ses))
    list_id_subj_exist<-list()
    for (ses in list_ses_exist){
      df_conn_ses<-df_conn[df_conn$ses==ses,]
      id_subj_exist<-unique(df_conn_ses$ID_pnTTC)
      id_subj_subset<-df_clin[df_clin$ses==ses,"ID_pnTTC"]
      id_subj_exist<-intersect(id_subj_exist,id_subj_subset)
      list_id_subj_exist[[as.character(ses)]]<-sort(id_subj_exist)
    }
    
    # Cbind FC data (Fisher-z transform of FC) as input for PCA function
    df_conn_cbind<-data.frame(matrix(nrow=n_edge,ncol=0))
    df_clin_exist<-data.frame(matrix(nrow=0,ncol=ncol(df_clin)))
    colnames(df_clin_exist)<-colnames(df_clin)
    for (ses in list_ses_exist){
      for (id_subj in list_id_subj_exist[[ses]]){
        df_conn_subj<-df_conn[which(df_conn$ID_pnTTC==id_subj),]
        df_conn_subj<-df_conn_subj[which(df_conn_subj$ses==ses),]
        df_conn_cbind<-cbind(df_conn_cbind,df_conn_subj[["z_r"]])
        df_clin_exist<-rbind(df_clin_exist,df_clin[df_clin$ses==ses & df_clin$ID_pnTTC==id_subj,])
      }
    }
    colnames(df_conn_cbind)<-as.character(seq(ncol(df_conn_cbind)))
    rownames(df_conn_cbind)<-NULL
    
    # Calculate PCA of FC
    print("Starting to calculate PCA of FC.")
    # Transpose connection dataframe (rows >> data for each subject/session, columns >> data for each edge)
    df_conn<-t(df_conn_cbind)
    data_pca<-func_pca(df_src=df_conn,df_var=df_edge,df_indiv=df_clin_exist)
    write.csv(data_pca$df_fac_var,file.path(paths_$output,"output",paste("atl-",atlas,"_pca_variable_factor.csv",sep="")),row.names=F)
    write.csv(data_pca$df_fac_indiv,file.path(paths_$output,"output",paste("atl-",atlas,"_pca_individual_factor.csv",sep="")),row.names=F)
    write.csv(data_pca$df_var_accounted,file.path(paths_$output,"output",paste("atl-",atlas,"_pca_variance_accounted.csv",sep="")),row.names=F)
    print("Finished calculating PCA of FC")
    
    # Plot PCA results
    print("Sarting to plot PCA of FC.")
    list_plot_pca<-plot_ca(df_src=data_pca$df_fac_indiv,list_name_covar=names(list_covar_),n_dim=data_pca$n_dim)
    for (i_dim in names(list_plot_pca)){
      for (name_covar in names(list_plot_pca[[i_dim]])){
        plot<-list_plot_pca[[i_dim]][[name_covar]]
        plot<-(plot
               + ggtitle("PCA of FC"))
        ggsave(paste("atl-",atlas,"_dim-",sprintf("%02d",as.numeric(i_dim)),"-",sprintf("%02d",as.numeric(i_dim)+1),"_cov-",name_covar,"_pca_fc.eps",sep=""),plot=plot,device=cairo_ps,
               path=file.path(paths_$output,"output"),dpi=300,height=10,width=10,limitsize=F)
      }
    }
    print("Finished plotting PCA of FC")
  }
  print("Finished pca_fc().")
}


#**************************************************
# Fingerprinting ==================================
#**************************************************
fingerprint<-function(paths_=paths,
                      list_atlas_=list_atlas,
                      subset_subj_=subset_subj){
  print("Starting fingerprint().")
  nullobj<-func_createdirs(paths_)
  for (atlas in list_atlas_){
    print(paste("Calculating atlas: ",atlas,sep=""))
    df_conn<-read.csv(file.path(paths_$input,"output",paste("atl-",atlas,"_fc.csv",sep="")))
    df_edge<-df_conn[which(df_conn$ID_pnTTC==df_conn[1,"ID_pnTTC"]),]
    df_edge<-df_edge[which(df_edge$ses==df_edge[1,"ses"]),c("from","to"),]
    n_edge<-dim(df_edge)[1]
    list_node<-sort(unique(c(as.character(unique(df_edge$from)),as.character(unique(df_edge$to)))))
    n_node<-length(list_node)
    list_ses_exist <- sort(unique(df_conn$ses))
    list_id_subj_exist<-list()
    for (ses in list_ses_exist){
      df_conn_ses<-df_conn[df_conn$ses==ses,]
      list_id_subj_exist[[as.character(ses)]]<-sort(unique(df_conn_ses$ID_pnTTC))
    }

    df_conn_cbind<-data.frame(matrix(nrow=n_edge,ncol=0))
    df_ses_subj<-data.frame(matrix(nrow=0,ncol=2))
    colnames(df_ses_subj)<-c("ses","ID_pnTTC")
    #list_file_tmp<-NULL
    for (ses in list_ses_exist){
      for (id_subj in list_id_subj_exist[[ses]]){
        #print(paste("Calculating Wave: ",as.character(ses), ", Subject: ",as.character(id_subj),sep=""))
        df_conn_subj<-df_conn[which(df_conn$ID_pnTTC==id_subj),]
        df_conn_subj<-df_conn_subj[which(df_conn_subj$ses==ses),]
        df_conn_cbind<-cbind(df_conn_cbind,df_conn_subj[["z_r"]])
        df_ses_subj<-rbind(df_ses_subj,data.frame(ses=ses,ID_pnTTC=id_subj))
      }
    }
    colnames(df_conn_cbind)<-as.character(seq(ncol(df_conn_cbind)))
    rownames(df_conn_cbind)<-NULL
    print("Starting to calculate correlation of fingerprints.")
    data_fingerprint<-func_cor(input=df_conn_cbind)
    df_fingerprint<-data_fingerprint$cor_flat
    print("Finished calculating correlation of fingerprints.")
    df_fingerprint$from_ses<-df_fingerprint$from_ID_pnTTC<-df_fingerprint$to_ses<-df_fingerprint$to_ID_pnTTC<-NA
    for (i in seq(dim(df_fingerprint)[1])){
      from_id<-df_fingerprint[[i,"from"]]
      to_id<-df_fingerprint[[i,"to"]]
      df_fingerprint[[i,"from_ses"]]<-df_ses_subj[[from_id,"ses"]]
      df_fingerprint[[i,"from_ID_pnTTC"]]<-df_ses_subj[[from_id,"ID_pnTTC"]]
      df_fingerprint[[i,"to_ses"]]<-df_ses_subj[[to_id,"ses"]]
      df_fingerprint[[i,"to_ID_pnTTC"]]<-df_ses_subj[[to_id,"ID_pnTTC"]]
    }
    df_fingerprint<-df_fingerprint[c("from_ses","from_ID_pnTTC","to_ses","to_ID_pnTTC","r")]
    write.csv(df_fingerprint,file.path(paths_$output,"output",paste("atl-",atlas,"_fingerprint.csv",sep="")),row.names=F)
    
    # Prepare dataframe for fingerprint correlation plot
    df_fp_plot<-data_fingerprint$cor
    list_name_subj_ses<-paste(sprintf("%05d",df_ses_subj$ID_pnTTC),as.character(df_ses_subj$ses),sep="_")
    colnames(df_fp_plot)<-rownames(df_fp_plot)<-list_name_subj_ses
    
    # Heatmap plot of fp correlation matrix
    plot_fp_heatmap<-plot_cor_heatmap(input=df_fp_plot)
    plot_fp_heatmap<-(plot_fp_heatmap
                      + scale_fill_gradientn(colors = matlab.like2(100),name="r")
                      + ggtitle("Fingerprint correlation")
                      + theme(plot.title = element_text(hjust = 0.5),
                              axis.title=element_blank()))
    
    # Save heatmap plot
    ggsave(paste("atl-",atlas,"_fp.eps",sep=""),plot=plot_fp_heatmap,device=cairo_ps,
           path=file.path(paths_$output,"output"),dpi=300,height=10,width=10,limitsize=F)
    
  }
  print("Finished fingerprint().")
}


#**************************************************
# GTA functionalities =============================
#**************************************************

edges2igraph<-function(df_conn,df_edge,list_node,dict_roi){
  edges<-data.frame(matrix(ncol=3,nrow=dim(df_edge)[1]))
  edges[,1:2]<-df_edge[,c("from","to")]
  for (i in 1:dim(df_edge)[1]){
    edges[i,3]<-as.numeric(df_conn[intersect(which(df_conn$from==df_edge[i,"from"]),
                                             which(df_conn$to==df_edge[i,"to"])),"r"])
  }
  colnames(edges)<-c("from","to","weight")
  nodes<-data.frame(id=list_node)
  for (i in seq(length(list_node))){
    nodes[i,"label"]<-as.character(dict_roi[dict_roi$id==as.character(nodes[i,"id"]),"label"])
  }
  output <- graph.data.frame(d = edges, vertices = nodes, directed = F)
  return(output)
}


#**************************************************
# Binary GTA ======================================
#**************************************************

# Subset edges according to desired cost
subset_edge<-function(input_igraph, input_cost,n_node,n_edge){
  n_edges4cost<-as.integer(n_node*(n_node-1)/2*input_cost)
  edges2delete<-head(order(E(input_igraph)$weight),(n_edge-n_edges4cost))
  output<-delete.edges(input_igraph,edges2delete)
  return(output)
}

# Calculate binary graph metrics
gta_bin_metrics<-function(input_igraph){
  node<-names(V(input_igraph))
  metrics<-data.frame(matrix(nrow=0,ncol=3))
  colnames(metrics)<-c("node","metric","value")
  ## graph-level metrics
  # characteristic path length
  metrics<-rbind(metrics,cbind(node="graph",metric="characteristic path length",
                               value=average.path.length(input_igraph)))
  # global efficiency
  eff<-1/(shortest.paths(input_igraph))
  eff[!is.finite(eff)]<-0
  metrics<-rbind(metrics,cbind(node="graph",metric="global efficiency",
                               value=mean(eff,na.rm=TRUE)))
  # global clustering coefficient
  metrics<-rbind(metrics,cbind(node="graph",metric="global clustering coefficient",
                               value=transitivity(input_igraph)))
  # average clustering coefficient
  metrics<-rbind(metrics,cbind(node="graph",metric="average clustering coefficient",
                               value=transitivity(input_igraph,type="average")))
  # local efficiency
  # modularity
  # small-worldness
  suppressWarnings(metrics<-rbind(metrics,cbind(node="graph",metric="small-world index",
                                                value=smallworldIndex(input_igraph)$index)))
  
  ## node-level metrics
  # degree centrality
  metrics<-rbind(metrics,cbind(node=node,metric="degree centrality",
                               value=centr_degree(input_igraph)$res))
  # betweenness centrality
  metrics<-rbind(metrics,cbind(node=node,metric="betweenness centrality",
                               value=centr_betw(input_igraph)$res))
  # eigenvector centrality
  metrics<-rbind(metrics,cbind(node=node,metric="eigenvector centrality",
                               value=eigen_centrality(input_igraph)$vector))
  
  rownames(metrics)<-NULL
  return(metrics)
}


gta_bin<-function(paths_=paths,
                  list_atlas_=list_atlas,
                  list_cost_=list_cost){
  print("Starting to calculate binary GTA.")
  #data_clinical<-func_clinical_data(paths_,subset_subj_)
  nullobj<-func_createdirs(paths_)
  dict_roi<-func_dict_roi(paths_)
  
  for (atlas in list_atlas_){
    print(paste("Calculate atlas: ",atlas,sep=""))
    file_conn<-paste("atl-",atlas,"_fc.csv",sep="")
    df_conn<-read.csv(file.path(paths_$input,"output",file_conn))
    df_edge<-df_conn[which(df_conn$ID_pnTTC==df_conn[1,"ID_pnTTC"]),]
    df_edge<-df_edge[which(df_edge$ses==df_edge[1,"ses"]),c("from","to"),]
    n_edge<-dim(df_edge)[1]
    list_node<-unique(c(as.character(unique(df_edge$from)),as.character(unique(df_edge$to))))
    list_node<-list_node[order(list_node)]
    n_node<-length(list_node)
    list_ses_exist <- sort(unique(df_conn$ses))
    list_id_subj_exist<-list()
    for (ses in list_ses_exist){
      df_conn_ses<-df_conn[df_conn$ses==ses,]
      list_id_subj_exist[[as.character(ses)]]<-sort(unique(df_conn_ses$ID_pnTTC))
    }
    #df_dst<-data.frame()
    list_file_tmp<-NULL
    for (ses in list_ses_exist){
      for (id_subj in list_id_subj_exist[[ses]]){
      #for (id_subj in list_id_subj_exist[[ses]][c(1,2)]){
        print(paste("Calculating Wave: ",as.character(ses), ", Subject: ",as.character(id_subj),sep=""))
        df_conn_subj<-df_conn[which(df_conn$ID_pnTTC==id_subj),]
        df_conn_subj<-df_conn_subj[which(df_conn_subj$ses==ses),]
        igraph_subj<-edges2igraph(df_conn=df_conn_subj,df_edge=df_edge,list_node=list_node,dict_roi=dict_roi)
        
        df_metric_subj<-data.frame()
        for (cost in list_cost_){
          igraph_subj_subset<-subset_edge(igraph_subj,cost,n_node,n_edge)
          E(igraph_subj_subset)$weight<-1
          metrics<-gta_bin_metrics(igraph_subj_subset)
          metrics<-cbind(cost=cost,metrics)
          df_metric_subj<-rbind(df_metric_subj,metrics)
        }
        df_metric_subj$value<-as.numeric.factor(df_metric_subj$value)
        df_metric<-df_metric_subj[which(df_metric_subj$cost==list_cost_[1]),c("node","metric")]
        average<-data.frame()
        for (i in seq(dim(df_metric)[1])){
          average<-rbind(average,
                         cbind(cost="average",node=as.character(df_metric[i,"node"]),
                               metric=as.character(df_metric[i,"metric"]),
                               value=mean(df_metric_subj[intersect(which(df_metric_subj$node==df_metric[i,"node"]),
                                                                   which(df_metric_subj$metric==df_metric[i,"metric"])),"value"])))
        }
        df_metric_subj<-rbind(df_metric_subj,average)
        rownames(df_metric_subj)<-NULL
        df_metric_subj<-cbind(ses=ses,ID_pnTTC=id_subj,df_metric_subj)
        file_metric_tmp<-paste("TMP_atl-",atlas,"_ses-",sprintf("%02d",ses),"_sub-",sprintf("%05d",id_subj),"_gta_bin.csv",sep="")
        path_file_metric_tmp<-file.path(paths_$output,"output",file_metric_tmp)
        write.csv(df_metric_subj,path_file_metric_tmp,row.names=F)
        list_file_tmp<-c(list_file_tmp,path_file_metric_tmp)
        #df_dst<-rbind(df_dst,df_metric_subj)
      }
    }
    df_dst<-data.frame()
    for (path_file_metric_tmp in list_file_tmp){
      df_tmp<-read.csv(path_file_metric_tmp)
      df_dst<-rbind(df_dst,df_tmp)
      file.remove(path_file_metric_tmp)
      print(paste("Finished binding: ",path_file_metric_tmp,sep=""))
    }
    file_dst<-paste("atl-",atlas,"_gta_bin.csv",sep="")
    write.csv(df_dst,file.path(paths_$output,"output",file_dst),row.names=F)
  }
  print("Finished gta_bin().")
}


#**************************************************
# Weighted GTA ====================================
#**************************************************

# add metric to output list in weighted GTA
AddMetric<-function(input){
  output<-data.frame(matrix(nrow=0,ncol=3))
  if (!is.null(input$graph)){
    output_add<-cbind(node="graph",metric=input$name[[1]],value=input$graph)
    output<-rbind(output,output_add)
  }
  if (!is.null(input$node)){
    output_add<-cbind(node=names(input$node),
                      metric=input$name[[1]],value=input$node)
    output<-rbind(output,output_add)
  }
  colnames(output)<-c("node","metric","value")
  return(output)
}

WeightedMetric<-function(input_igraph){
  metrics<-data.frame(matrix(nrow=0,ncol=3))
  distance<-WeightedDistance(input_igraph)$distance
  
  metrics<-rbind(metrics,AddMetric(WeightedCharPath(input_distance=distance)))
  #metrics<-rbind(metrics,AddMetric(WeightedEccentricity(input_distance = distance)))
  metrics<-rbind(metrics,AddMetric(WeightedRadius(input_distance = distance)))
  metrics<-rbind(metrics,AddMetric(WeightedDiameter(input_distance = distance)))
  metrics<-rbind(metrics,AddMetric(WeightedGlobalEfficiency(input_distance = distance)))
  metrics<-rbind(metrics,AddMetric(WeightedClustCoef(input = input_igraph)))
  metrics<-rbind(metrics,AddMetric(WeightedTransitivity(input = input_igraph)))
  #metrics<-rbind(metrics,AddMetric(WeightedLocalEfficiency(input = input_igraph)))
  metrics<-rbind(metrics,AddMetric(WeightedModularity(input = input_igraph)))
  #metrics<-rbind(metrics,AddMetric(WeightedStrength(input = input_igraph)))
  #metrics<-rbind(metrics,AddMetric(WeightedClosenessCentrality(input_distance = distance)))
  #metrics<-rbind(metrics,AddMetric(WeightedBetweennessCentrality(input = input_igraph)))
  #metrics<-rbind(metrics,AddMetric(WeightedEigenvectorCentrality(input = input_igraph)))
  #metrics<-rbind(metrics,AddMetric(WeightedNeighborDegree(input = input_igraph)))
  metrics<-rbind(metrics,AddMetric(WeightedAssortativityCoef(input = input_igraph)))
  
  colnames(metrics)<-c("node","metric","value")
  rownames(metrics)<-NULL
  return(metrics)
}


gta_weight<-function(absolute=T,
                     threshold=NA,
                     paths_=paths,
                     list_atlas_=list_atlas,
                     list_wave_=list_wave,
                     subset_subj_=subset_subj){
  print("Starting gta_weight().")
  nullobj<-func_createdirs(paths_)
  dict_roi<-func_dict_roi(paths_)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  df_clin<-func_clinical_data_long(paths_,list_wave_)
  data_subset_clin<-func_subset_clin(df_clin,
                                     list_wave_,subset_subj_,
                                     list_covar=NULL,
                                     rem_na_clin=T)
  df_clin_subset<-data_subset_clin$df_clin

  for (atlas in list_atlas_){
    print(paste("Calculate atlas: ",atlas,sep=""))
    file_conn<-paste("atl-",atlas,"_fc.csv",sep="")
    df_conn<-read.csv(file.path(paths_$input,"output",file_conn))
    df_edge<-df_conn[which(df_conn$ID_pnTTC==df_conn[1,"ID_pnTTC"]),]
    df_edge<-df_edge[which(df_edge$ses==df_edge[1,"ses"]),c("from","to"),]
    n_edge<-dim(df_edge)[1]
    list_node<-unique(c(as.character(unique(df_edge$from)),as.character(unique(df_edge$to))))
    list_node<-list_node[order(list_node)]
    n_node<-length(list_node)
    list_ses_exist <- sort(unique(df_conn$ses))
    list_id_subj_exist<-list()
    for (ses in list_ses_exist){
      df_conn_ses<-df_conn[df_conn$ses==ses,]
      id_subj_exist<-unique(df_conn_ses$ID_pnTTC)
      id_subj_subset<-df_clin_subset[df_clin_subset$wave==ses,"ID_pnTTC"]
      id_subj_exist<-intersect(id_subj_exist,id_subj_subset)
      list_id_subj_exist[[as.character(ses)]]<-sort(id_subj_exist)
    }
    #df_dst<-data.frame()
    list_file_tmp<-NULL
    for (ses in list_ses_exist){
      for (id_subj in list_id_subj_exist[[ses]]){
        #for (id_subj in list_id_subj_exist[[ses]][c(1,2)]){
        print(paste("Calculating Wave: ",as.character(ses), ", Subject: ",as.character(id_subj),sep=""))
        df_conn_subj<-df_conn[which(df_conn$ID_pnTTC==id_subj),]
        df_conn_subj<-df_conn_subj[which(df_conn_subj$ses==ses),]
        igraph_subj<-edges2igraph(df_conn=df_conn_subj,df_edge=df_edge,list_node=list_node,dict_roi=dict_roi)
        if (absolute){
          E(igraph_subj)$weight<-abs(E(igraph_subj)$weight)
        }
        E(igraph_subj)$weight[is.na(E(igraph_subj)$weight)]<-0
        df_metric_subj<-WeightedMetric(igraph_subj)
        df_metric_subj$value<-as.numeric.factor(df_metric_subj$value)
        rownames(df_metric_subj)<-NULL
        df_metric_subj<-cbind(ses=ses,ID_pnTTC=id_subj,df_metric_subj)
        file_metric_tmp<-paste("TMP_atl-",atlas,"_ses-",sprintf("%02d",ses),"_sub-",sprintf("%05d",id_subj),"_gta_weight.csv",sep="")
        path_file_metric_tmp<-file.path(paths_$output,"output",file_metric_tmp)
        write.csv(df_metric_subj,path_file_metric_tmp,row.names=F)
        list_file_tmp<-c(list_file_tmp,path_file_metric_tmp)
        #df_dst<-rbind(df_dst,df_metric_subj)
      }
    }
    df_dst<-data.frame()
    for (path_file_metric_tmp in list_file_tmp){
      df_tmp<-read.csv(path_file_metric_tmp)
      df_dst<-rbind(df_dst,df_tmp)
      file.remove(path_file_metric_tmp)
      print(paste("Finished binding: ",path_file_metric_tmp,sep=""))
    }
    file_dst<-paste("atl-",atlas,"_gta_weight.csv",sep="")
    write.csv(df_dst,file.path(paths_$output,"output",file_dst),row.names=F)
  }
  print("Finished gta_weight().")
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

