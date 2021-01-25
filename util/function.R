#**************************************************
# Description =====================================
#**************************************************

# R script for common MRI analysis functions


#**************************************************
# Libraries =======================================
#**************************************************
library(tidyverse)
library(dplyr)
library(Hmisc)
library(FactoMineR)
library(missMDA)
library(ica)
library(parallel)
#library(DescTools)
library(pbapply)
library(ggpubr)


#**************************************************
# Network-based statistics ========================
#**************************************************
func_nbs<-function(df_fc,df_clin,df_roi,list_mod,list_plot,thr_p_cdt,progressbar){
  df_join<-join_fc_clin(df_fc,df_clin)
  # value as slope of z_r longitudinal difference against age (z(r(wave=2))-z(r(wave=1)))/delta(age)
  df_join$value<-df_join$value/df_join$diff_age
  data_gamm<-iterate_gamm(df_join,df_roi,list_mod,calc_parallel=T,calc_identical=F,
                          list_sex=list(c(1,2)),progressbar=progressbar)
  list_out_bfs<-list()
  for (model in names(list_mod)){
    list_out_model<-list()
    for (plot in names(list_plot)){
      var_exp<-list_plot[[plot]][["var_exp"]]
      df_gamm_subset<-data_gamm$df_out_gamm[data_gamm$df_out_gamm$model==model & data_gamm$df_out_gamm$term==var_exp,]
      if (nrow(df_gamm_subset)>0){
        df_gamm_sign<-df_gamm_subset[df_gamm_subset$p<thr_p_cdt*2,] # multiply with 2: two-sided to one-sided
        df_m<-df_gamm_sign[df_gamm_sign$t<0,]
        data_bfs_m<-func_bfs(df_m)
        df_f<-df_gamm_sign[df_gamm_sign$t>0,]
        data_bfs_f<-func_bfs(df_f)
        list_out_plot<-list(list("m"=data_bfs_m,"f"=data_bfs_f))
        names(list_out_plot)<-plot
        list_out_model<-c(list_out_model,list_out_plot)
      }
    }
    list_out_model<-list(list_out_model)
    names(list_out_model)<-model
    list_out_bfs<-c(list_out_bfs,list_out_model)
  }
  return(list_out_bfs)
}

func_nbs_threshold<-function(paths,data_nbs,list_max,list_mod,list_plot,thr_p_perm,
                             atlas,df_roi,df_grp){
  df_net<-df_size_net<-df_perm<-df_thr_size<-NULL
  for (model in names(data_nbs)){
    for (plot in names(data_nbs[[model]])){
      for (sex in c("m","f")){
        data_nbs_subset<-data_nbs[[model]][[plot]][[sex]]
        list_max_subset<-list_max[[model]][[plot]][[sex]]
        list_max_subset_sort<-sort(list_max_subset)
        thr_size_perm<-list_max_subset_sort[ceiling(length(list_max_subset_sort)*(1-thr_p_perm))]
        
        if (sex=="m"){
          title_sex<-"M>F"
          color_plt<-"steelblue2"
        }else{
          title_sex<-"F>M"
          color_plt<-"lightcoral"
        }
        title_plot<-list_plot[[plot]][["title"]]
        plot_permutation(paths,list_max=list_max_subset_sort,thr_size_perm,
                         atlas,model,plot,sex,title_plot,title_sex,color_plt)

        df_head<-data.frame(atlas=atlas,mod=model,plot=plot,sex=sex)
        list_network_sign<-list()
        if(length(data_nbs_subset$list_network)>0){
          for (idx_net in seq(length(data_nbs_subset$list_network))){
            network<-data_nbs_subset$list_network[[idx_net]]
            df_net<-rbind(df_net,
                          cbind(df_head,
                                data.frame(id_net=idx_net,network$df_edge)))
            df_size_net<-rbind(df_size_net,
                               cbind(df_head,
                                     data.frame(id_net=idx_net,size=network$size_net)))
            if (network$size_net>=thr_size_perm){
              list_network_sign<-c(list(network))
              plot_sex_diff_fc(paths,network$df_edge,atlas,df_roi,df_grp,
                               model,plot,sex,title_plot,title_sex,idx_net)
            }
          }
        }
        data_nbs[[model]][[plot]][[sex]][["list_network_sign_perm"]]<-list_network_sign
        data_nbs[[model]][[plot]][[sex]][["list_max_size_perm"]]<-list_max_subset
        data_nbs[[model]][[plot]][[sex]][["thr_size_perm"]]<-thr_size_perm
        df_thr_size<-rbind(df_thr_size,
                           cbind(df_head,data.frame(thr_size=thr_size_perm)))
        df_perm<-rbind(df_perm,
                       cbind(df_head,data.frame(id_perm=seq(length(list_max_subset)),
                                                max_size_net=list_max_subset)))
      }
    }
  }
  write.csv(df_net,file.path(paths$output,"output","temp",
                             paste("atl-",atlas,"_net.csv",sep="")),row.names=F)
  write.csv(df_size_net,file.path(paths$output,"output","temp",
                             paste("atl-",atlas,"_size_net.csv",sep="")),row.names=F)
  write.csv(df_thr_size,file.path(paths$output,"output","temp",
                             paste("atl-",atlas,"_thr_perm.csv",sep="")),row.names=F)
  write.csv(df_perm,file.path(paths$output,"output","temp",
                              paste("atl-",atlas,"_perm.csv",sep="")),row.names=F)
  return(data_nbs)
}


#**************************************************
# Join FC and clinical data =======================
#**************************************************
join_fc_clin<-function(df_fc,df_clin){
  df_fc$z_r[which(is.nan(df_fc$z_r))]<-0
  colnames(df_fc)[colnames(df_fc)=="z_r"]<-"value"
  colnames(df_fc)[colnames(df_fc)=="ses"]<-"wave"
  df_fc$r<-df_fc$p<-NULL
  df_clin$wave<-as.character(df_clin$wave)
  
  # Join clinical and FC data frames
  #print('Joining clinical and FC data.')
  df_join<-inner_join(df_fc,df_clin,by=c('ID_pnTTC','wave'))
  for (key in c('ID_pnTTC','wave','sex')){
    if (key %in% colnames(df_join)){
      df_join[,key]<-as.factor(df_join[,key])
    }
  }
  df_join$value<-as.numeric.factor(df_join$value)
  return(df_join)
}

join_fc_clin_cs<-function(df_fc,df_clin,wave_clin,wave_mri){
  df_fc$z_r[which(is.nan(df_fc$z_r))]<-0
  colnames(df_fc)[colnames(df_fc)=="z_r"]<-"value"
  df_fc<-df_fc[df_fc$ses==wave_mri,]
  df_fc$ses<-NULL
  #colnames(df_fc)[colnames(df_fc)=="ses"]<-"wave"
  df_fc<-df_fc[,c("ID_pnTTC","from","to","value")]
  
  df_clin<-df_clin[df_clin$wave==wave_clin,]
  df_clin$wave<-NULL
  
  # Join clinical and FC data frames
  #print('Joining clinical and FC data.')
  #df_join<-inner_join(df_fc,df_clin,by=c('ID_pnTTC','wave'))
  df_join<-inner_join(df_fc,df_clin,by='ID_pnTTC')
  for (key in c('ID_pnTTC','wave','sex')){
    if (key %in% colnames(df_join)){
      df_join[,key]<-as.factor(df_join[,key])
    }
  }
  df_join$value<-as.numeric.factor(df_join$value)
  return(df_join)
}

#**************************************************
# Prepare longitudinal FC data ====================
#**************************************************
prep_data_fc<-function(paths,atlas,key_group,include_diff=F){
  dict_roi <- func_dict_roi(paths)
  
  df_fc<-as.data.frame(fread(file.path(paths$input,"output",
                                       paste("atl-",atlas,"_fc.csv",sep=""))))
  if (!include_diff){
    df_fc<-df_fc[df_fc$ses!="2-1",]
  }
  
  # Prepare dataframe of ROIs
  list_roi<-sort(unique(c(as.character(df_fc$from),as.character(df_fc$to))))
  df_roi<-dict_roi[is.element(dict_roi$id,list_roi),c("id","label",key_group)]
  colnames(df_roi)[colnames(df_roi)==key_group]<-"group"
  
  # Prepare dataframe of group-wise FC averages
  list_group<-unique(as.character(df_roi$group))
  df_fc_temp<-df_fc
  df_fc_temp$z_r[which(is.nan(df_fc_temp$z_r))]<-0
  df_fc_temp<-inner_join(df_fc_temp,df_roi[,c("id","group")],by=c("from"="id"))
  colnames(df_fc_temp)[colnames(df_fc_temp)=="group"]<-"from_group"
  df_fc_temp<-inner_join(df_fc_temp,df_roi[,c("id","group")],by=c("to"="id"))
  colnames(df_fc_temp)[colnames(df_fc_temp)=="group"]<-"to_group"
  df_subj<-NULL
  list_subj<-sort(unique(df_fc$ID_pnTTC))
  for (id_subj in list_subj){
    list_ses<-sort(unique(df_fc[df_fc$ID_pnTTC==id_subj,"ses"]))
    #list_ses<-list_ses[list_ses!="2-1"]
    df_subj<-rbind(df_subj,data.frame(ID_pnTTC=id_subj,ses=list_ses))
  }
  df_subj$ses<-as.character(df_subj$ses)
  
  df_fc_grp<-data.frame()
  for (idx_subj_ses in seq(dim(df_subj)[1])){
    #print(paste(df_subj[idx_subj_ses,"ID_pnTTC"],df_subj[idx_subj_ses,"ses"]))
    df_fc_subset1<-df_fc_temp[df_fc_temp$ID_pnTTC==df_subj[idx_subj_ses,"ID_pnTTC"]
                              & df_fc_temp$ses==df_subj[idx_subj_ses,"ses"],]
    for (idx_grp1 in seq(length(list_group))){
      for (idx_grp2 in seq(idx_grp1,length(list_group))){
        # data in df_fc_subset2 is doubled for connections within same group,
        # but does not affect z_r average calculation
        df_fc_subset2<-rbind(df_fc_subset1[df_fc_subset1$from_group==list_group[idx_grp1]
                                           & df_fc_subset1$to_group==list_group[idx_grp2],],
                             df_fc_subset1[df_fc_subset1$from_group==list_group[idx_grp2]
                                           & df_fc_subset1$to_group==list_group[idx_grp1],])
        df_fc_grp<-rbind(df_fc_grp,
                         cbind(ID_pnTTC=df_subj[idx_subj_ses,"ID_pnTTC"],ses=df_subj[idx_subj_ses,"ses"],
                               from=list_group[idx_grp1],to=list_group[idx_grp2],
                               z_r=mean(df_fc_subset2$z_r)))
      }
    }
  }
  df_fc_grp$ID_pnTTC<-as.numeric(as.numeric.factor(df_fc_grp$ID_pnTTC))
  df_fc_grp$ses<-as.character(as.numeric.factor(df_fc_grp$ses))
  
  df_grp<-data.frame(id=list_group,label=str_to_title(gsub("_"," ",as.character(list_group))))
  return(list("df_fc"=df_fc,"df_fc_grp"=df_fc_grp,"df_roi"=df_roi,"df_grp"=df_grp))
}



#**************************************************
# Breadth-first search of connected network =======
#**************************************************
func_bfs<-function(df_edge){
  df_edge_remain<-df_edge
  list_network<-list()
  list_size<-NULL
  while (nrow(df_edge_remain)>0){ # Examine new subnetwork as long as any edge is remaining
    # Node as the origin of subnetwork
    list_node_todo<-list_node_net<-as.character(df_edge_remain[[1,"from"]])
    df_edge_net<-data.frame()
    while (length(list_node_todo)>0){ # Search as long as any connected and unexamined node exist
      node_check<-as.character(list_node_todo[1]) # Node of interest
      list_node_todo<-list_node_todo[-1]
      df_edge_new_from<-df_edge_remain[df_edge_remain$from==node_check,]
      df_edge_remain<-df_edge_remain[rownames(df_edge_remain) %nin% rownames(df_edge_new_from),]
      list_node_new_from<-as.character(df_edge_new_from[,"to"])
      df_edge_new_to<-df_edge_remain[df_edge_remain$to==node_check,]
      df_edge_remain<-df_edge_remain[rownames(df_edge_remain) %nin% rownames(df_edge_new_to),]
      list_node_new_to<-as.character(df_edge_new_to[,"from"])
      # Edge and node connected to node_check
      df_edge_new<-rbind(df_edge_new_to,df_edge_new_from)
      list_node_new<-c(list_node_new_from,list_node_new_to)
      # Add to subnetwork currently examined
      df_edge_net<-rbind(df_edge_net,df_edge_new)
      list_node_net<-c(list_node_net,list_node_new)
      # Add to nodes to be examined
      list_node_todo<-unique(c(list_node_todo,list_node_new))
    }
    size_net<-nrow(df_edge_net)
    list_size<-c(list_size,size_net)
    list_node_net<-unique(list_node_net)
    list_network<-c(list_network,list(list("df_edge"=df_edge_net,"list_node"=list_node_net,"size_net"=size_net)))
  }
  if(is.null(list_size)){
    max_size<-0
    n_network<-0
  }else{
    max_size<-max(list_size)
    n_network<-length(list_size)
  }
  output<-list("list_network"=list_network,"list_size"=list_size,"max_size"=max_size,"n_network"=n_network)
  return(output)
}


#**************************************************
# Fisher transformation of Correlation to Z =======
#**************************************************
func_fisherz<-function(rho){
  return((log((1+rho)/(1-rho)))/2)
}


#**************************************************
# Create path list ================================
#**************************************************
func_path<-function(list_path_root = c("C:/Users/atiro","D:/atiro","/home/atiroms","D:/NICT_WS"),
                    path_exp_,dir_in_,dir_out_,path_exp_full_=NULL
                    ){
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
  path_common <- file.path(path_root,"Dropbox/MRI_img/pnTTC/puberty/common")
  if (is.null(path_exp_full_)){
    path_io     <- file.path(path_root,path_exp_)
  }else{
    path_io     <- path_exp_full_
  }
  path_in     <- file.path(path_io,dir_in_)
  path_out    <- file.path(path_io,dir_out_)
  output <- list("script"=path_script,"io"=path_io,"input"=path_in,"output"=path_out,
                 "common"=path_common,"dir_in"=dir_in_,"dir_out"=dir_out_)
  return(output)
}

#paths<-func_path()


#**************************************************
# Iterate GAM/GLM over ROI paiers in FC ===========
#**************************************************

gamm_core<-function(data_src){
  df_src<-data_src$df_src
  list_mod_<-data_src$list_mod
  list_sex<-data_src$list_sex
  
  #list_sex<-sort(unique(as.numeric.factor(df_src$sex)))
  df_out_aic_add<-df_out_gamm_add<-df_out_anova_add<-data.frame()
  for (idx_mod in names(list_mod_)){
    for (idx_sex in list_sex){
      df_src_sex<-NULL
      label_sex<-NULL
      for (subidx_sex in idx_sex){
        df_src_sex<-rbind(df_src_sex,df_src[df_src$sex==subidx_sex,])
        if (is.null(label_sex)){
          label_sex<-as.character(subidx_sex)
        }else{
          label_sex<-paste(label_sex,subidx_sex,sep="_")
        }
      }
      
      df_src_sex$value<-as.numeric(df_src_sex$value)
      if (data_src$calc_parallel){
        mod<-try(gam(as.formula(list_mod_[[idx_mod]]),data=df_src_sex,method="REML",control=list(nthreads=1)), silent=F)
      }else{
        mod<-try(gam(as.formula(list_mod_[[idx_mod]]),data=df_src_sex,method="REML"), silent=F)
      }
      if (class(mod)[1]!="try-error"){
        p_table<-summary.gam(mod)$p.table
        df_out_gamm_add_add<-data.frame(term=rownames(p_table),estimate=p_table[,'Estimate'],
                                        se=p_table[,'Std. Error'],F=NA,t=p_table[,'t value'],
                                        p=p_table[,'Pr(>|t|)'])
        s_table<-summary.gam(mod)$s.table
        if(!is.null(s_table)){
          df_out_gamm_add_add<-rbind(df_out_gamm_add_add,
                                     data.frame(term=rownames(s_table),estimate=NA,se=NA,F=s_table[,'F'],
                                                t=NA,p=s_table[,'p-value']))
        }
        p_table_anova<-anova.gam(mod)$pTerms.table
        colnames(p_table_anova)<-c("df","F","p")
        df_out_gamm_add<-rbind(df_out_gamm_add,
                               cbind(sex=label_sex,model=idx_mod,df_out_gamm_add_add))
        df_out_aic_add<-rbind(df_out_aic_add,
                              data.frame(sex=label_sex,model=idx_mod,aic=mod$aic,aic_best_among_models=0))
        df_out_anova_add<-rbind(df_out_anova_add,
                                cbind(sex=label_sex,model=idx_mod,term=rownames(p_table_anova),
                                      p_table_anova))
      }
    } # Finished looping over sex
  }# Finished looping over model
  
  # Compare AICs of GAMM models
  df_out_aic_add_sex_rbind<-data.frame()
  for (idx_sex in list_sex){
    label_sex<-NULL
    for (subidx_sex in idx_sex){
      if (is.null(label_sex)){
        label_sex<-as.character(subidx_sex)
      }else{
        label_sex<-paste(label_sex,subidx_sex,sep="_")
      }
    }
    df_out_aic_add_sex<-df_out_aic_add[df_out_aic_add$sex==label_sex,]
    df_out_aic_add_sex[which(df_out_aic_add_sex$aic==min(df_out_aic_add_sex$aic)),
                       'aic_best_among_models']<-1
    df_out_aic_add_sex_rbind<-rbind(df_out_aic_add_sex_rbind,df_out_aic_add_sex)
  }
  
  # Prepare output dataframe
  df_id<-data.frame("from"=data_src$id_from,"to"=data_src$id_to,"label_from"=data_src$label_from,"label_to"=data_src$label_to)
  df_out_gamm_add<-cbind(df_id,df_out_gamm_add)
  df_out_aic_add_sex_rbind<-cbind(df_id,df_out_aic_add_sex_rbind)
  df_out_anova_add<-cbind(df_id,df_out_anova_add)
  
  return(list("df_out_gamm_add"=df_out_gamm_add,"df_out_aic_add"=df_out_aic_add_sex_rbind,
              "df_out_anova_add"=df_out_anova_add))
}

combine_gamm<-function(list_dst_sub){
  df_gamm<-df_aic<-df_anova<-data.frame()
  for (dst_sub in list_dst_sub){
    df_gamm<-rbind(df_gamm,dst_sub$df_out_gamm_add)
    df_aic<-rbind(df_aic,dst_sub$df_out_aic_add)
    df_anova<-rbind(df_anova,dst_sub$df_out_anova_add)
  }
  return(list("df_out_gamm_add"=df_gamm,"df_out_aic_add"=df_aic,
              "df_out_anova_add"=df_anova))
}

iterate_gamm<-function(df_join,df_roi,list_mod_,calc_parallel=T,calc_identical=F,list_sex=NULL,progressbar=T){
  list_roi<-df_roi$id
  if (is.null(list_sex)){
    list_sex<-sort(unique(as.numeric.factor(df_join$sex)))
  }
  
  # Prepare dataset for multi-core processing
  #print("Preparing dataset for parallel processing.")
  list_src_gamm<-list()
  if (calc_identical){
    list_id_from<-list_roi
  }else{
    list_id_from<-list_roi[-length(list_roi)]
  }
  for (id_from in list_id_from){
    df_join_from<-df_join[df_join$from==id_from,]
    label_from<-as.character(df_roi[df_roi$id==id_from,"label"])
    if (calc_identical){
      list_id_to<-list_roi[seq(which(list_roi==id_from),length(list_roi))]
    }else{
      list_id_to<-list_roi[seq(which(list_roi==id_from)+1,length(list_roi))]
    }
    for(id_to in list_id_to){
      label_to<-as.character(df_roi[df_roi$id==id_to,"label"])
      df_src<-df_join_from[df_join_from$to==id_to,]
      df_src$from<-df_src$to<-NULL
      list_src_gamm<-c(list_src_gamm,list(list("df_src"=df_src,"id_from"=id_from,"id_to"=id_to,
                                               "label_from"=label_from,"label_to"=label_to,
                                               "list_mod"=list_mod_,"list_sex"=list_sex,
                                               "calc_parallel"=calc_parallel)))
    }
  }
  
  # Parallel processing
  if (calc_parallel){
    n_cluster<-min(floor(detectCores()*3/4),length(list_src_gamm))
  }else{
    n_cluster<-2
  }
  clust<-makeCluster(n_cluster)
  #print(paste("Calculating GAM in parallel,",as.character(n_cluster),"cores.",sep=" "))
  clusterExport(clust,
                varlist=c("list_mod_","sort","gam","as.formula","summary.gam",
                          "anova.gam","as.numeric.factor"),
                envir=environment())
  if (progressbar){
    list_dst_gamm<-pblapply(list_src_gamm,gamm_core,cl=clust)
  }else{
    list_dst_gamm<-parLapply(clust,list_src_gamm,gamm_core)
  }
  stopCluster(clust)
  
  # Collect data into dataframes
  len_list<-length(list_dst_gamm)
  len_sublist<-floor(sqrt(len_list)*2)
  n_sublist<-ceil(len_list/len_sublist)
  #print(paste("Dividing results into", as.character(n_sublist), "sublists.",sep=" "))
  list_dst_gamm_sub<-list()
  for (idx_sublist in 1:n_sublist){
    #print(paste("Subgroup",as.character(idx_sublist),sep=" "))
    if (idx_sublist!=n_sublist){
      list_dst_gamm_sub<-c(list_dst_gamm_sub,
                           list(list_dst_gamm[((idx_sublist-1)*len_sublist+1):(idx_sublist*len_sublist)]))
    }else{
      list_dst_gamm_sub<-c(list_dst_gamm_sub,
                           list(list_dst_gamm[((idx_sublist-1)*len_sublist+1):len_list]))
    }
  }
  list_dst_gamm<-NULL
  gc()
  
  n_cluster<-floor(detectCores()*3/4)
  clust<-makeCluster(n_cluster)
  #print(paste("Combining within sublists,",as.character(n_cluster),"cores.",sep=" "))
  clusterExport(clust,
                varlist=NULL,
                envir=environment())
  list_dst_gamm<-parLapply(clust,list_dst_gamm_sub,combine_gamm)
  stopCluster(clust)
  
  #print("Combining sublists.")
  df_out_gamm<-df_out_aic<-df_out_anova<-NULL
  for (dst_gamm in list_dst_gamm){
    df_out_gamm<-rbind(df_out_gamm,dst_gamm$df_out_gamm_add)
    df_out_aic<-rbind(df_out_aic,dst_gamm$df_out_aic_add)
    df_out_anova<-rbind(df_out_anova,dst_gamm$df_out_anova_add)
  }
  list_dst_gamm<-NULL
  gc()
  
  rownames(df_out_gamm)<-rownames(df_out_aic)<-rownames(df_out_anova)<-NULL
  return(list("df_out_gamm"=df_out_gamm,"df_out_aic"=df_out_aic,"df_out_anova"=df_out_anova))
}


iterate_gamm_old<-function(df_join,df_roi,list_mod_){
  list_roi<-df_roi$id
  df_out_gamm<-df_out_aic<-NULL
  for (id_from in list_roi[-length(list_roi)]){
    for(id_to in list_roi[seq(which(list_roi==id_from)+1,length(list_roi))]){
      label_from<-as.character(df_roi[df_roi$id==id_from,"label"])
      label_to<-as.character(df_roi[df_roi$id==id_to,"label"])
      df_src=df_join[df_join$from==id_from & df_join$to==id_to,]
      
      print(paste("GLM/GAM: ",id_from," - ",id_to," (",label_from," - ",label_to,")",sep=""))
      df_out_aic_add<-df_out_gamm_add<-data.frame()
      for (idx_mod in names(list_mod_)){
        list_plot<-list()
        list_sex<-sort(unique(as.numeric.factor(df_src$sex)))
        for (idx_sex in list_sex){
          df_src_sex<-df_src[df_src$sex==idx_sex,]
          #mod<-gam(as.formula(list_mod_[[idx_mod]]),data=df_src_sex)
          mod<-try(gam(as.formula(list_mod_[[idx_mod]]),data=df_src_sex,method="REML"), silent=F)
          if (class(mod)[1]=="try-error"){
            print(paste("Error fiting ",idx_mod, ", sex= ",idx_sex,".",sep=''))
          }else{
            p_table<-summary.gam(mod)$p.table
            if (is.null(summary.gam(mod)$s.table)){
              df_out_gamm_add_add<-data.frame(from=id_from,to=id_to,label_from=label_from,label_to=label_to,
                                              #roi=roi,label_roi=label_roi,group=group,measure=measure,
                                              sex=idx_sex,model=idx_mod,term=rownames(p_table),
                                              estimate=p_table[,'Estimate'],se=p_table[,'Std. Error'],F=NA,
                                              t=p_table[,'t value'],p=p_table[,'Pr(>|t|)'])
              
            }else{
              s_table<-summary.gam(mod)$s.table
              df_out_gamm_add_add<-rbind(data.frame(from=id_from,to=id_to,label_from=label_from,label_to=label_to,
                                                    #roi=roi,label_roi=label_roi,group=group,measure=measure,
                                                    sex=idx_sex,model=idx_mod,term=rownames(p_table),
                                                    estimate=p_table[,'Estimate'],se=p_table[,'Std. Error'],F=NA,
                                                    t=p_table[,'t value'],p=p_table[,'Pr(>|t|)']),
                                         data.frame(from=id_from,to=id_to,label_from=label_from,label_to=label_to,
                                                    #roi=roi,label_roi=label_roi,group=group,measure=measure,
                                                    sex=idx_sex,model=idx_mod,term=rownames(s_table),
                                                    estimate=NA,se=NA,F=s_table[,'F'],
                                                    t=NA,p=s_table[,'p-value']))
            }
            df_out_gamm_add<-rbind(df_out_gamm_add,df_out_gamm_add_add)
            df_out_aic_add<-rbind(df_out_aic_add,
                                  data.frame(from=id_from,to=id_to,label_from=label_from,label_to=label_to,
                                             #roi=roi,label_roi=label_roi,group=group,measure=measure,
                                             sex=idx_sex,
                                             model=idx_mod,aic=mod$aic,aic_best_among_models=0))
          }
        }
      }
      
      # Compare AICs of GAMM models
      df_out_aic_add_sex_rbind<-data.frame()
      for (idx_sex in list_sex){
        df_out_aic_add_sex<-df_out_aic_add[df_out_aic_add$sex==idx_sex,]
        df_out_aic_add_sex[which(df_out_aic_add_sex$aic==min(df_out_aic_add_sex$aic)),
                           'aic_best_among_models']<-1
        df_out_aic_add_sex_rbind<-rbind(df_out_aic_add_sex_rbind,df_out_aic_add_sex)
      }
      
      df_out_gamm<-rbind(df_out_gamm,df_out_gamm_add)
      df_out_aic<-rbind(df_out_aic,df_out_aic_add_sex_rbind)
    }
  }
  rownames(df_out_gamm)<-rownames(df_out_aic)<-NULL
  output<-list("df_out_gamm"=df_out_gamm,"df_out_aic"==df_out_aic)
  return(output)
}


#**************************************************
# Add multiple comparison-corrected p values ======
#**************************************************
mltcomp_corr<-function(input){
  output<-data.frame("p_bonf"=p.adjust(input$p,method = "bonferroni"),
                     "p_holm"=p.adjust(input$p,method = "holm"),
                     "p_hoch"=p.adjust(input$p,method = "hochberg"),
                     "p_homm"=p.adjust(input$p,method = "hommel"),
                     "p_bh"=p.adjust(input$p,method="BH"),
                     "p_by"=p.adjust(input$p,method="BY"))
  return(output)
}

add_mltcmp<-function(df_out_gamm,df_roi,list_mod,list_plot,calc_seed_level=TRUE){
  df_plot_gamm_concat<-NULL
  for (idx_mod in names(list_mod)){
    for (idx_plot in names(list_plot)){
      var_exp<-list_plot[[idx_plot]][["var_exp"]]
      for (idx_sex in c(1,2)){
        # Subset GAMM result dataframe for plotting
        if (idx_sex==1){
          label_sex<-"m"
        }else{
          label_sex<-"f"
        }
        df_plot_gamm<-df_out_gamm[df_out_gamm$model==idx_mod 
                                  & df_out_gamm$term==var_exp
                                  & df_out_gamm$sex==idx_sex,]
        if (nrow(df_plot_gamm)>0){
          # Calculate graph-level multiple comparison-corrected p values
          df_plot_gamm<-cbind(df_plot_gamm,mltcomp_corr(df_plot_gamm))
          
          # Calculate seed-level multiple comparison-corrected p values
          if (calc_seed_level){
            for (idx_roi in as.character(df_roi$id)){
              list_row_seed<-sort(union(which(df_plot_gamm$from==idx_roi),
                                        which(df_plot_gamm$to==idx_roi)))
              df_plot_gamm_seed<-df_plot_gamm[list_row_seed,]
              df_p_seed<-mltcomp_corr(df_plot_gamm_seed)
              for (idx_edge in seq(length(list_row_seed))){# iterate over edges which starts / ends at idx_roi
                for (type_p in colnames(df_p_seed)){  # iterate over types of p values
                  # Enter corrected p to df_plot_gamm if empty or new value is smaller
                  df_plot_gamm[list_row_seed[idx_edge],
                               paste("seed",type_p,sep="_")]<-min(df_plot_gamm[list_row_seed[idx_edge],
                                                                               paste("seed",type_p,sep="_")],
                                                                  df_p_seed[idx_edge,type_p],
                                                                  na.rm=T)
                }
              }
            }
          }
        }
        df_plot_gamm_concat<-rbind(df_plot_gamm_concat,df_plot_gamm)
      }
    }
  }
  return(df_plot_gamm_concat)
}


#**************************************************
# Factor to numeric function ======================
#**************************************************
as.numeric.factor <- function(x) {
  if (class(x)=="factor"){
    return(as.numeric(levels(x))[x])
  }else{
    return(x)
  }
}


#**************************************************
# Experiment folder preparation ===================
#**************************************************
func_createdirs<-function(paths,str_proc="",copy_log=T){
  list_createdirs<-c(paths$output,
                     file.path(paths$output,"output"),
                     file.path(paths$output,"output","plot"),
                     file.path(paths$output,"output","temp"),
                     file.path(paths$output,"output","result"))
  for(d in list_createdirs){
    if (file.exists(d)){
      print(paste("Destination folder already exists:",d,sep=" "))
    }else{
      dir.create(d)
    }
  }
  if (copy_log){
    if (file.exists(file.path(paths$input,"log"))){
      if (file.exists(file.path(paths$output,"log"))){
        print("Destination log folder already exists.")
      }else{
        file.copy(file.path(paths$input,"log"),paths$output,recursive=T)
        list_log<-readLines(file.path(paths$output,"log","pipeline.log"))
        list_log<-c(list_log,paste("\t",str_proc,sep=""))
        list_log<-c(list_log,paths$dir_out)
        writeLines(list_log,file.path(paths$output,"log","pipeline.log"))
      }
    }else{
      print("Source log folder does not exist.")
    }
  }
}


#**************************************************
# Returns ROI dictionary ==========================
#**************************************************
func_dict_roi<-function(paths,
                        file_roi="ROI.csv"){
  output<-read.csv(file.path(paths$common,file_roi))
  #output<-read.csv('/home/atiroms/Dropbox/temp/ROI.csv')
  return(output)
}


#**************************************************
# Clinical data loading ===========================
#**************************************************
func_clinical_data<-function(paths,
                             subset_subj,
                             file_clinical= "CSUB.csv"
                             ){
  df_clinical <- read.csv(file.path(paths$common,file_clinical))
  for (list_subset in subset_subj){
    df_clinical <- df_clinical[which(df_clinical[,list_subset[["column"]]]==list_subset[["value"]]),]
  }
  list_id_subj<-df_clinical$ID_pnTTC
  df_id_subj<-df_clinical[,1:(which(colnames(df_clinical)=="Clinical")-1)]
  df_clinical<-df_clinical[,(-2):(-which(colnames(df_clinical)=="Clinical"))]
  n_subj<-length(list_id_subj)
  n_data_clinical<-ncol(df_clinical)-1
  
  output<-list("df_clinical"=df_clinical,"list_id_subj"=list_id_subj,"df_id_subj"=df_id_subj,
               "n_subj"=n_subj,"n_data_clinical"=n_data_clinical)
  return(output)
}


#**************************************************
# Longitudinal clinical data loading ==============
#**************************************************
print_log<-function(str_in,str_add,print_terminal=T){
  if(print_terminal){
    print(str_add)
  }
  str_out<-c(str_in,str_add)
  return(str_out)
}

func_clinical_data_long<-function(paths,list_wave,subset_subj,list_covar,rem_na_clin,
                                  file_clin= "CSUB.csv",prefix="",print_terminal=T
                                  ){
  df_src <- read.csv(file.path(paths$common,file_clin))
  
  # Vertically concatenate wave-wise clinical data
  df_clin_long <- data.frame(matrix(nrow=0,ncol=ncol(df_src)+1))
  for (wave in list_wave){
    df_tmp<-df_src['ID_pnTTC']
    df_tmp$wave<-rep(wave,nrow(df_tmp))
    df_tmp<-cbind(df_tmp,df_src[,colnames(df_src)!='ID_pnTTC'])
    df_clin_long<-rbind(df_clin_long,df_tmp)
  }
  colnames(df_clin_long)<-c('ID_pnTTC','wave',colnames(df_src)[colnames(df_src)!='ID_pnTTC'])
  #return(df_clin_long)
  list_log<-NULL
  
  df_clin_subset<-data.frame(matrix(nrow=0,ncol=ncol(df_clin_long)))
  
  # Subset clinical data according to subsetting condition
  #print('Clinical: subsetting clinical data according to specified condition.')
  list_id_subset<-list()
  for (wave in list_wave){
    str_wave<-as.character(wave)
    #print(paste('Clinical: checking wave ',str_wave,sep=''))
    df_clin_wave<-df_clin_long[df_clin_long['wave']==wave,]
    #print(paste('Clinical: ',as.character(nrow(df_clin_wave)),' source clinical data identified',sep=''))
    list_log<-print_log(list_log,paste('Clinical: ',as.character(nrow(df_clin_wave)),' source clinical data identified.',sep=''),print_terminal)
    id_intersect<-df_clin_wave[,'ID_pnTTC']
    list_id_subset_wave<-list("src"=id_intersect)
    for (key_condition in subset_subj[[str_wave]]){
      key_subset<-key_condition$key
      #value_subset<-key_condition$value
      condition_subset<-key_condition$condition
      #id_meet_cond<-df_clin_wave[df_clin_wave[key_subset]==value_subset,'ID_pnTTC']
      id_meet_cond<-df_clin_wave[eval(parse(text=paste('df_clin_wave[[key_subset]]',condition_subset,sep=''))),'ID_pnTTC']
      id_meet_cond<-id_meet_cond[!is.na(id_meet_cond)]
      id_intersect<-intersect(id_intersect,id_meet_cond)
      #print(paste('Clinical: ',as.character(length(id_meet_cond)),' subjects meeting ',key_subset, ' = ',as.character(value_subset),sep=''))
      list_log<-print_log(list_log,paste('Clinical: ',as.character(length(id_meet_cond)),' subjects satisfying ',key_subset, condition_subset,sep=''),print_terminal)
      id_meet_cond<-list(id_meet_cond)
      names(id_meet_cond)<-key_subset
      list_id_subset_wave<-c(list_id_subset_wave,id_meet_cond)
    }
    list_log<-print_log(list_log,paste('Clinical: ',as.character(length(id_intersect)),' subjects satisfying all conditions.',sep=''),print_terminal)
    df_clin_wave<-df_clin_wave[df_clin_wave$ID_pnTTC %in% id_intersect,]
    df_clin_subset<-rbind(df_clin_subset,df_clin_wave)
    id_intersect<-list('intersect'=id_intersect)
    list_id_subset_wave<-c(list_id_subset_wave,id_intersect)
    list_id_subset_wave<-list(list_id_subset_wave)
    names(list_id_subset_wave)<-str_wave
    list_id_subset<-c(list_id_subset,list_id_subset_wave)
  } # Finished looping over waves
  colnames(df_clin_subset)<-colnames(df_clin_long)
  
  # Subset unused columns of clinical data
  # Subset clinical data with missing covariate value
  list_log<-print_log(list_log,'Clinical: subsetting clinical data with missing covariate value.',print_terminal)
  df_clin_exist<-data.frame(matrix(nrow=0,ncol=length(list_covar)+2))
  list_id_exist<-list()
  for (wave in list_wave){
    str_wave<-as.character(wave)
    list_log<-print_log(list_log,paste('Clinical: checking wave ',str_wave,sep=''),print_terminal)
    df_clin_exist_wave<-df_clin_subset[df_clin_subset$wave==wave,c('ID_pnTTC','wave')]
    id_exist_intersect<-df_clin_exist_wave$ID_pnTTC
    list_log<-print_log(list_log,paste('Clinical: ',as.character(length(id_exist_intersect)),' source clinical data identified.',sep=''),print_terminal)
    list_id_exist_wave<-list('src'=id_exist_intersect)
    if (length(list_covar)>0){
      for (id_covar in seq(length(list_covar))){
        list_name_covar_src<-list_covar[[id_covar]][[str_wave]]
        name_covar_dst<-names(list_covar)[id_covar]
        
        # Choose non-NA value from list of possible source covariate columns
        n_subj_wave<-dim(df_clin_subset[df_clin_subset$wave==wave,])[1]
        df_clin_exist_wave_add<-data.frame(matrix(nrow=n_subj_wave,ncol=1))
        for (name_covar_src in list_name_covar_src){
          df_clin_exist_wave_add<-pmax(df_clin_exist_wave_add,
                                       df_clin_subset[df_clin_subset$wave==wave,name_covar_src,drop=F],
                                       na.rm=TRUE)
        }
        
        colnames(df_clin_exist_wave_add)<-name_covar_dst
        df_clin_exist_wave<-cbind(df_clin_exist_wave,df_clin_exist_wave_add)
        id_exist<-df_clin_exist_wave[!is.na(df_clin_exist_wave[name_covar_dst]),'ID_pnTTC']
        id_exist_intersect<-intersect(id_exist_intersect,id_exist)
        list_log<-print_log(list_log,paste('Clinical: ',as.character(length(id_exist)),' subjects with non-NA values of covariate: ',paste(list_name_covar_src,collapse="/"),sep=''),print_terminal)
        id_exist<-list(id_exist)
        names(id_exist)<-name_covar_dst
        list_id_exist_wave<-c(list_id_exist_wave,id_exist)
      }
      list_log<-print_log(list_log,paste('Clinical: ',as.character(length(id_exist_intersect)),' subjects with non-NA values of all convariates.',sep=''),print_terminal)
    }
    n_subj_pre<-nrow(df_clin_exist_wave)
    if (rem_na_clin){
      df_clin_exist_wave<-df_clin_exist_wave[df_clin_exist_wave$ID_pnTTC %in% id_exist_intersect,]
      n_subj_deleted<-n_subj_pre-nrow(df_clin_exist_wave)
      list_log<-print_log(list_log,paste('Clinical: ',as.character(n_subj_deleted),' subjects with NA values in any covariate >> deleted.',sep=''),print_terminal)
    }else{
      n_subj_na<-n_subj_pre-nrow(df_clin_exist_wave[df_clin_exist_wave$ID_pnTTC %in% id_exist_intersect,])
      list_log<-print_log(list_log,paste('Clinical: ',as.character(n_subj_na),' subjects with NA values in any covariate >> NOT deleted.',sep=''),print_terminal)
    }
    list_id_exist_wave<-c(list_id_exist_wave,list('intersect'=id_exist_intersect))
    df_clin_exist<-rbind(df_clin_exist,df_clin_exist_wave)
    list_id_exist_wave<-list(list_id_exist_wave)
    names(list_id_exist_wave)<-str_wave
    list_id_exist<-c(list_id_exist,list_id_exist_wave)
  }
  colnames(df_clin_exist)<-c('ID_pnTTC','wave',names(list_covar))
  rownames(df_clin_exist)<-NULL
  
  writeLines(list_log, file.path(paths$output,"output","temp",paste(prefix,"clin_long.txt",sep="_")))
  
  output<-list('df_clin'=df_clin_exist,'list_id_subset'=list_id_subset,'list_id_exist'=list_id_exist)
  return(output)
}


#**************************************************
# Calculate difference and means of clinical data =
#**************************************************
func_clinical_data_diffmean<-function(df_src,list_id_subj,list_covar){
  
  # Classify covaiates to fixed values and unfixed values
  list_covar_fix<-list_covar_change<-NULL
  for (id_covar in names(list_covar)){
    if (list_covar[[id_covar]][["1"]][1]==list_covar[[id_covar]][["2"]][1]){
      list_covar_fix<-c(list_covar_fix,id_covar)
    }else{
      list_covar_change<-c(list_covar_change,id_covar)
    }
  }
  
  # Prepare fixed covariates
  df_clin_fix<-data.frame(ID_pnTTC=list_id_subj)
  for (id_covar in list_covar_fix){
    df_clin_fix[[id_covar]]<-df_src[df_src$ID_pnTTC %in% list_id_subj & df_src$ses==1,id_covar]
  }
  
  # Prepare unfixed values
  df_clin_change<-list()
  for (ses in c(1,2)){
    df_clin_change_ses<-df_src[df_src$ses==ses & df_src$ID_pnTTC %in% list_id_subj, list_covar_change,drop=F]
    df_clin_change_ses<-list(df_clin_change_ses)
    names(df_clin_change_ses)<-as.character(ses)
    df_clin_change<-c(df_clin_change,df_clin_change_ses)
  }
  
  # Calculate difference and mean
  df_clin_diff<-df_clin_change[["2"]]-df_clin_change[["1"]]
  df_clin_mean<-(df_clin_change[["1"]]+df_clin_change[["2"]])/2
  
  # Change column names
  colnames(df_clin_change[["1"]])<-c(paste("ses1_",colnames(df_clin_change[["1"]]),sep=''))
  colnames(df_clin_change[["2"]])<-c(paste("ses2_",colnames(df_clin_change[["2"]]),sep=''))
  colnames(df_clin_diff)<-c(paste("diff_",colnames(df_clin_diff),sep=''))
  colnames(df_clin_mean)<-c(paste("mean_",colnames(df_clin_mean),sep=''))
  
  # Join dataframes
  df_join<-cbind(df_clin_fix,df_clin_change[["1"]],df_clin_change[["2"]],df_clin_diff,df_clin_mean)
  
  return(df_join)
}


#**************************************************
# Subsetting Structural data  =====================
#**************************************************
func_subset_str<-function(df_str,
                          list_measure,
                          list_str_group,
                          dict_roi,
                          key_group='group_3'){
  
  df_str_dst<-data.frame(matrix(nrow=0,ncol=ncol(df_str)))
  for (measure in list_measure){
    print(paste('Subsetting ',measure,' measures.',sep=''))
    df_str_tmp<-df_str[df_str$measure==measure,]
    list_roi<-unique(df_str_tmp$roi)
    list_roi<-list_roi[order(list_roi)]
    print(paste(as.character(length(list_roi)),' ROIs exist in the source data.',sep=''))
    list_roi_subset<-NULL
    for (roi in list_roi){
      group_roi<-as.character(dict_roi[dict_roi$id==roi,key_group])
      if (group_roi %in% list_str_group){
        list_roi_subset<-c(list_roi_subset,roi)
        df_str_tmp[df_str_tmp$roi==roi,"group"]<-group_roi
      }
    }
    print(paste(as.character(length(list_roi_subset)),' ROIs remaining.',sep=''))
    df_str_tmp<-df_str_tmp[df_str_tmp$roi %in% list_roi_subset,]
    df_str_dst<-rbind(df_str_dst,df_str_tmp)
  }
  colnames(df_str_dst)<-c(colnames(df_str),"group")
  return(df_str_dst)
}


#**************************************************
# General correlation calculation =================
#**************************************************
func_cor<-function(input,type="pearson"){
  cor <-rcorr(as.matrix(input), type=type)
  n_node<-ncol(input)
  cor_flat<-data.frame(matrix(nrow=n_node*(n_node-1)/2,ncol=4))
  colnames(cor_flat)<-c("from","to","r","p")
  k<-0
  for (i in 1:(n_node-1)){
    for (j in (i+1):n_node){
      k<-k+1
      cor_flat[k,1:4]<-c(rownames(cor$r)[i],
                         colnames(cor$r)[j],
                         cor$r[i,j],
                         cor$P[i,j])
    }
  }
  #mean_cor<-mean(cor$r,na.rm=TRUE)
  #sd_cor<-sd(cor$r,na.rm=TRUE)
  #cor_flat$z_r<-(as.numeric(cor_flat$r)-mean_cor)/sd_cor
  
  #cor_flat$z_r<-FisherZ(as.numeric(cor_flat$r))
  cor_flat$z_r<-func_fisherz(as.numeric(cor_flat$r))
  output<-list("cor"=data.frame(cor$r),"cor_flat"=cor_flat)
  return(output)
}


#**************************************************
# General PCA calculation =========================
#**************************************************
func_pca<-function(df_src,df_var=NULL,df_indiv=NULL,dim_ca=NULL,calc_corr=F){
  if (sum(is.na(df_src))>0){
    # Estimate number of dimensions
    if (is.null(dim_ca)){
      ncp_estim<-estim_ncpPCA(df_src,ncp.max=ncol(df_src))$ncp
      ncp_calc<-ncp_estim
    }else{
      ncp_estim<-estim_ncpPCA(df_src,ncp.max=dim_ca)$ncp
      if (ncp_estim==dim_ca){
        print(paste("PCA data dimension may be greater than: ",as.character(ncp_estim),sep=""))
      }
      ncp_calc<-dim_ca
    }
    # Impute data
    df_src<-imputePCA(df_src,ncp=ncp_calc)$completeObs
  }else{
    ncp_calc<-dim_ca
  }
  
  print(paste("Calculating PCA, dimension: ",as.character(ncp_calc),sep=""))
  # PCA calculation
  data_pca<-PCA(df_src,scale.unit = TRUE, ncp = ncp_calc, graph = FALSE)
  
  # Component-imaging variable matrix
  # Row: MRI variable, Column: component(factor)
  df_comp_mri<-data.frame(data_pca$var$coord)
  if(!is.null(df_var)){
    df_comp_mri<-cbind(df_var,df_comp_mri)
    colnames(df_comp_mri)<-c(colnames(df_var),sprintf("comp_%03d",1:ncp_calc))
  }else{
    colnames(df_comp_mri)<-sprintf("comp_%03d",1:ncp_calc)
  }
  rownames(df_comp_mri)<-NULL
  df_comp_mri<-rownames_to_column(df_comp_mri,"id")
  df_comp_mri$id<-as.numeric(df_comp_mri$id)
  colnames_out<-colnames(df_var)
  for (id_comp in 1:ncp_calc){
    df_comp_mri$abs<-abs(df_comp_mri[[sprintf("comp_%03d",id_comp)]])
    df_comp_mri<-df_comp_mri[order(df_comp_mri$abs,decreasing=T),]
    df_comp_mri[[sprintf("rank_%03d",id_comp)]]<-1:nrow(df_comp_mri)
    colnames_out<-c(colnames_out,sprintf("comp_%03d",id_comp), sprintf("rank_%03d",id_comp))
  }
  df_comp_mri<-df_comp_mri[order(df_comp_mri$id),colnames_out]
  
  # Component-Individual matrix
  # Row: subject, Column: (clinical variable +) component(factor)
  df_comp_subj<-data.frame(data_pca$ind$coord)
  colnames(df_comp_subj)<-sprintf("comp_%03d",1:ncp_calc)

  df_cor<-NULL
  df_cor_flat<-NULL
  if(!is.null(df_indiv)){
    df_comp_subj<-cbind(df_indiv,df_comp_subj)
    
    if (calc_corr){
      # Calculate correlation between component attribution and clinical covariate
      df_covar<-df_indiv[,-which(colnames(df_indiv) %in% c("ID_pnTTC","ses"))]
      n_covar<-ncol(df_covar)
      df_comp_subj_covar<-cbind(df_covar,df_comp_subj)
      data_cor<-func_cor(df_comp_subj_covar)
      df_cor<-data_cor$cor
      df_cor<-df_cor[(n_covar+1):nrow(df_cor),1:n_covar]
      df_cor_flat<-data_cor$cor_flat
      df_cor_flat<-df_cor_flat[df_cor_flat$from %in% colnames(df_covar) & df_cor_flat$to %in% colnames(df_comp_subj),]
      df_cor_flat<-df_cor_flat[,c("from","to","r","p")]
      colnames(df_cor_flat)<-c("covar","component","r","p")
    }
  }
  rownames(df_comp_subj)<-NULL
  
  # Matrix of variance accounted
  # Row: component(factor)
  df_vaf<-data.frame(data_pca$eig)
  colnames(df_vaf)<-c("eigenvalue","vaf","cumul_vaf")
  df_vaf$vaf<-df_vaf$vaf/100
  df_vaf$cumul_vaf<-df_vaf$cumul_vaf/100
  df_vaf$comp<-seq(1,dim(df_vaf)[1])
  df_vaf<-df_vaf[1:ncp_calc,c("comp","vaf","cumul_vaf","eigenvalue")]
  rownames(df_vaf)<-NULL
  
  return(list('df_comp_mri'=df_comp_mri,'df_comp_subj'=df_comp_subj,
              'df_vaf'=df_vaf,'dim'=ncp_calc,
              'df_comp_clin'=df_cor,'df_comp_clin_flat'=df_cor_flat))
}


#**************************************************
# General ICA calculation =========================
#**************************************************
func_ica<-function(df_src,df_var=NULL,df_indiv=NULL,dim_ca=NULL,calc_corr){
  
  if (is.null(dim_ca)){
    ncp_calc<-nrow(df_src)-1
  }else{
    ncp_calc<-dim_ca
  }
  
  print(paste("Calculating ICA, dimension: ",as.character(ncp_calc),sep=""))
  
  # Imputation using means
  if (sum(is.na(df_src))>0){
    df_src<-impute(df_src,mean)
  }
  
  df_src<-data.matrix(df_src)
  
  # ICA calculation
  #data_ica <-icafast(df_src, nc=ncp_calc,center=TRUE,maxit=100,tol=1e-6,alg="par",fun="logcosh",alpha=1)
  #data_ica <-icafast(df_src, nc=ncp_calc,center=TRUE)
  data_ica <-icaimax(df_src, nc=ncp_calc,center=TRUE)
  
  # Component-imaging variable matrix
  # Row: MRI variable, Column: component(factor)
  df_comp_mri<-data.frame(data_ica$M)
  if(!is.null(df_var)){
    df_comp_mri<-cbind(df_var,df_comp_mri)
    colnames(df_comp_mri)<-c(colnames(df_var),sprintf("comp_%03d",1:ncp_calc))
  }else{
    colnames(df_comp_mri)<-sprintf("comp_%03d",1:ncp_calc)
  }
  rownames(df_comp_mri)<-NULL
  df_comp_mri<-rownames_to_column(df_comp_mri,"id")
  df_comp_mri$id<-as.numeric(df_comp_mri$id)
  colnames_out<-colnames(df_var)
  for (id_comp in 1:ncp_calc){
    df_comp_mri$abs<-abs(df_comp_mri[[sprintf("comp_%03d",id_comp)]])
    df_comp_mri<-df_comp_mri[order(df_comp_mri$abs,decreasing=T),]
    df_comp_mri[[sprintf("rank_%03d",id_comp)]]<-1:nrow(df_comp_mri)
    colnames_out<-c(colnames_out,sprintf("comp_%03d",id_comp), sprintf("rank_%03d",id_comp))
  }
  df_comp_mri<-df_comp_mri[order(df_comp_mri$id),colnames_out]
  
  # Component-Individual matrix
  # Row: subject, Column: (clinical variable +) component(factor)
  df_comp_subj<-data.frame(data_ica$S)
  colnames(df_comp_subj)<-sprintf("comp_%03d",1:ncp_calc)
  
  df_cor<-NULL
  df_cor_flat<-NULL
  if(!is.null(df_indiv)){
    df_comp_subj<-cbind(df_indiv,df_comp_subj)
    
    if (calc_corr){
      # Calculate correlation between component attribution and clinical covariate
      df_covar<-df_indiv[,-which(colnames(df_indiv) %in% c("ID_pnTTC","ses"))]
      n_covar<-ncol(df_covar)
      df_comp_subj_covar<-cbind(df_covar,df_comp_subj)
      data_cor<-func_cor(df_comp_subj_covar)
      df_cor<-data_cor$cor
      df_cor<-df_cor[(n_covar+1):nrow(df_cor),1:n_covar]
      df_cor_flat<-data_cor$cor_flat
      df_cor_flat<-df_cor_flat[df_cor_flat$from %in% colnames(df_covar) & df_cor_flat$to %in% colnames(df_comp_subj),]
      df_cor_flat<-df_cor_flat[,c("from","to","r","p")]
      colnames(df_cor_flat)<-c("covar","component","r","p")
    }
  }
  rownames(df_comp_subj)<-NULL
  
  # Matrix of variance accounted
  # Row: component(factor)
  df_vaf<-data.frame(data_ica$vafs)
  colnames(df_vaf)<-"vaf"
  df_vaf[1,"cumul_vaf"]<-df_vaf[1,"vaf"]
  for (idx_row in 2:nrow(df_vaf)){
    df_vaf[idx_row,"cumul_vaf"]<-df_vaf[idx_row-1,"cumul_vaf"]+df_vaf[idx_row,"vaf"]
  }
  df_vaf$comp<-seq(1,dim(df_vaf)[1])
  df_vaf<-df_vaf[c("comp","vaf","cumul_vaf")]
  rownames(df_vaf)<-NULL
  
  return(list('df_comp_mri'=df_comp_mri,'df_comp_subj'=df_comp_subj,
              'df_vaf'=df_vaf,'dim'=ncp_calc,
              'df_comp_clin'=df_cor,'df_comp_clin_flat'=df_cor_flat))
}