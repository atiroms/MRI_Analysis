#**************************************************
# Description =====================================
#**************************************************
# R script to analyze relationship between structural/functional connection and clinical data.
# Inputs can be functional correlation from rsfMRI data, or Jackknife esimate of structural covariance from T1 data.


#**************************************************
# Parameters ======================================
#**************************************************

#path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"
#dir_in<-"450_fc_test"
#dir_out<-"451_gammfc_test"
#list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")
#list_atlas<-"aal116"

path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_CONN"
dir_in<-"56.2_fc"
dir_out<-"56.3_gamm_fc"
#list_atlas<-c("cnn","hoa","power264")
#list_atlas<-"cnn"
list_atlas<-"hoa"

list_wave <- c(1,2)

list_covar<-list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
                 "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
                 "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
                 "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S"),
                 "age"  =list("1"="W1_Age_at_MRI",  "2"="W2_Age_at_MRI",  "label"="Age"),
                 "sex"  =list("1"="Sex",            "2"="Sex",            "label"="Sex"))

subset_subj <- list("1"=list(list("key"="W1_T1QC","condition"="==1"),
                             list("key"="W1_rsfMRIexist","condition"="==1"),
                             list("key"="W1_Censor","condition"="<126")),
                    "2"=list(list("key"="W2_T1QC","condition"="==1"),
                             list("key"="W2_rsfMRIexist","condition"="==1"),
                             list("key"="W2_Censor","condition"="<126")))

list_mod <- list("lin"= "value ~ age + testo + s(ID_pnTTC,bs='re')")
                 #"add"= "value ~ s(age,k=3) + s(testo,k=3) + s(ID_pnTTC,bs='re')",
                 #"quad"="value ~ poly(age,2) + poly(testo,2) + s(ID_pnTTC,bs='re')")

list_plot <-list(#"a"=list("title"="Age effect","var_exp"="age"),
                 #"sa"=list("title"="Age effect","var_exp"="s(age)"),
                 #"pa1"=list("title"="Age effect","var_exp"="poly(age, 2)1"),
                 #"pa2"=list("title"="Age effect","var_exp"="poly(age, 2)2"),
                 "t"=list("title"="Testosterone effect","var_exp"="testo")
                 #"st"=list("title"="Testosterone effect","var_exp"="s(testo)")
                 #"st"=list("title"="Testosterone effect","var_exp"="s(testo)"),
                 #"pt1"=list("title"="Testosterone effect","var_exp"="poly(testo, 2)1"),
                 #"pt2"=list("title"="Testosterone effect","var_exp"="poly(testo, 2)2")
                 )


list_type_p=c("p","p_bh","seed_p_bh")
#list_type_p="p"
thr_p <- 0.05

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
library(dplyr)
library(parallel)
library(mgcv)
library(car)


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
source(file.path(paths$script,"util/function.R"))
#source(file.path(paths$script,"util/glm_function.R"))
source(file.path(paths$script,"util/plot.R"))
source(file.path(paths$script,"util/gta_function.R"))


#**************************************************
# GLM/GAM of FCs ==================================
#**************************************************
iterate_gamm<-function(df_join,df_roi,list_mod_){
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
      df_out_aic<-rbind(df_out_aic,df_out_aic_add)
    }
  }
  rownames(df_out_gamm)<-rownames(df_out_aic)<-NULL
  output<-list("df_out_gamm"=df_out_gamm,"df_out_aic"==df_out_aic)
  return(output)
}

join_fc_clin<-function(df_fc,df_clin){
  df_fc$z_r[which(is.nan(df_fc$z_r))]<-0
  colnames(df_fc)[colnames(df_fc)=="z_r"]<-"value"
  colnames(df_fc)[colnames(df_fc)=="ses"]<-"wave"
  df_fc<-df_fc[,c(-which(colnames(df_fc)=="r"),
                  -which(colnames(df_fc)=="p"))]
  
  # Join clinical and FC data frames
  print('Joining clinical and FC data.')
  df_join<-inner_join(df_fc,df_clin,by=c('ID_pnTTC','wave'))
  for (key in c('ID_pnTTC','wave','sex')){
    if (key %in% colnames(df_join)){
      df_join[,key]<-as.factor(df_join[,key])
    }
  }
  return(df_join)
}

add_mltcmp<-function(df_out_gamm,df_roi,analysis,atlas,list_mod,list_plot,calc_seed_level=TRUE){
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

plot_gamm_fc<-function(df_plot_gamm,df_roi,analysis,atlas,list_mod,list_plot,
                       list_type_p,thr_p,paths_){
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
        df_plot_gamm_subset<-df_plot_gamm[df_plot_gamm$model==idx_mod 
                                  & df_plot_gamm$term==var_exp
                                  & df_plot_gamm$sex==idx_sex,]
        if (nrow(df_plot_gamm_subset)>0){
          print(paste("GAMM output, atlas: ",atlas,", model: ",idx_mod,", plot: ",var_exp,", sex: ",label_sex,sep=""))
          # Convert GAMM rseult into igraph object
          if (!is.na(df_plot_gamm_subset[1,"estimate"])){
            df_plot_gamm_subset<-rename(df_plot_gamm_subset,"weight"="estimate")
          }else{
            df_plot_gamm_subset<-rename(df_plot_gamm_subset,"weight"="F")
          }
          
          # Convert FC dataframe into iGraph object
          list_roi<-as.character(df_roi$id)
          df_node<-data.frame(id=list_roi,stringsAsFactors = F)
          for (idx_node in seq(dim(df_node)[1])){
            df_node[idx_node,"label"]<-as.character(df_roi[df_roi$id==df_node[idx_node,"id"],"label"])
          }
          df_edge<-df_plot_gamm_subset
          df_edge$from<-as.character(df_edge$from)
          df_edge$to<-as.character(df_edge$to)
          igraph_gamm <- graph_from_data_frame(d = df_edge, vertices = df_node, directed = F)
          
          # Plot and save circular graph
          for (type_p in list_type_p){
            if(type_p %in% colnames(df_plot_gamm_subset)){
              plot<-plot_circular(igraph_in=igraph_gamm,
                                  type_p=type_p,thr_p=thr_p,
                                  limit_color=NULL)
              plot<-plot +
                ggtitle(paste("GLM/GAM sex: ",label_sex, ", model: ",idx_mod,", expvar: ",var_exp,
                              "\nanalysis: ",analysis," threshold: ",type_p,sep="")) +
                theme(plot.title = element_text(hjust = 0.5))
              ggsave(paste("atl-",atlas,"_anl-",analysis,"_mod-",idx_mod,"_plt-",var_exp,
                           "_sex-",label_sex,"_pval-",type_p,"_gamm_fc.eps",sep=""),
                     plot=plot,device=cairo_ps,path=file.path(paths_$output,"output"),
                     dpi=300,height=10,width=10,limitsize=F)
            }
          }
        }
      }
    }
  }
}

gamm_fc<-function(paths_=paths,subset_subj_=subset_subj,list_covar_=list_covar,
                  list_wave_=list_wave,list_atlas_=list_atlas,
                  #list_measure_=list_measure,list_str_group_=list_str_group,
                  list_mod_=list_mod,list_plot_=list_plot,key_group_='group_3',
                  list_type_p_=list_type_p,thr_p_=thr_p
                  ){
  print("Starting gamm_fc().")
  nullobj<-func_createdirs(paths_,copy_log=T)
  dict_roi <- func_dict_roi(paths_)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,
                                     list_covar=list_covar_,rem_na_clin=T)
  df_clin<-data_clin$df_clin
  
  
  for (atlas in list_atlas_){
    
    #****************************
    # ROI-wise FC GAMM calculation
    #****************************
    # Load ROI-wise FC data
    print(paste('Loading FC data, atlas:',atlas,sep=' '))
    df_fc<-read.csv(file.path(paths_$input,'output',paste('atl-',atlas,'_fc.csv',sep='')))
    df_join<-join_fc_clin(df_fc,df_clin)
    write.csv(df_join,file.path(paths_$output,"output",paste("atl-",atlas,"_anl-roi_src.csv",sep="")),
              row.names=F)
    
    # Calculate and save ROI-wise GAMM of FC
    print(paste('Calculating GAMM, atlas: ',atlas,sep=''))
    list_roi<-sort(unique(c(as.character(df_join$from),as.character(df_join$to))))
    df_roi<-dict_roi[is.element(dict_roi$id,list_roi),c("id","label",key_group_)]
    colnames(df_roi)[colnames(df_roi)==key_group_]<-"group"
    data_gamm<-iterate_gamm(df_join,df_roi,list_mod_)
    write.csv(data_gamm$df_out_gamm,
              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-roi_gamm.csv",sep="")),row.names = F)
    write.csv(data_gamm$df_out_aic,
              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-roi_gamm_aic.csv",sep="")),row.names = F)
    
    # Calculate multiple comparison-corrected p values
    df_plot_gamm<-add_mltcmp(data_gamm$df_out_gamm,df_roi,analysis="roi",atlas,
                             list_mod,list_plot,calc_seed_level=T)
    write.csv(df_plot_gamm,
              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-roi_gamm_plt.csv",sep="")),row.names = F)
    
    # Graphical output of ROI-wise GAMM of FC
    plot_gamm_fc(df_plot_gamm,df_roi,analysis="roi",atlas,list_mod,list_plot,
                 list_type_p_,thr_p,paths_)
    
    #****************************
    # Group-wise FC GAMM calculation
    #****************************
    # Load group-wise FC data
    print(paste('Loading group FC data, atlas:',atlas,sep=' '))
    df_fc_grp<-read.csv(file.path(paths_$input,'output',paste('atl-',atlas,'_fc_grp.csv',sep='')))
    df_join_grp<-join_fc_clin(df_fc_grp,df_clin)
    write.csv(df_join_grp,file.path(paths_$output,"output",paste("atl-",atlas,"_anl-grp_src.csv",sep="")),
              row.names=F)
    
    # Calculate and save group-wise GAMM of FC
    print(paste('Calculating GAMM, atlas: ',atlas,sep=''))
    list_roi_grp<-sort(unique(c(as.character(df_join_grp$from),as.character(df_join_grp$to))))
    #df_roi<-dict_roi[is.element(dict_roi$id,list_roi),c("id","label",key_group_)]
    df_roi_grp<-data.frame(id=list_roi_grp,label=capitalize(list_roi_grp),group="group")
    data_gamm_grp<-iterate_gamm(df_join_grp,df_roi_grp,list_mod_)
    write.csv(data_gamm_grp$df_out_gamm,
              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-grp_gamm.csv",sep="")),row.names = F)
    write.csv(data_gamm_grp$df_out_aic,
              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-grp_gamm_aic.csv",sep="")),row.names = F)
    
    # Calculate multiple comparison-corrected p values
    df_plot_gamm_grp<-add_mltcmp(data_gamm_grp$df_out_gamm,df_roi_grp,analysis="grp",atlas,
                                 list_mod,list_plot,calc_seed_level=T)
    write.csv(df_plot_gamm_grp,
              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-grp_gamm_plt.csv",sep="")),row.names = F)
    
    # Graphical output of group-wise GAMM of FC
    plot_gamm_fc(df_plot_gamm_grp,df_roi_grp,analysis="grp",atlas,list_mod,list_plot,
                 list_type_p_,thr_p,paths_)
    
    #****************************
    # Multi-scale FC GAMM calculation
    #****************************
    # Subset ROI-wise GAMM result to include only within-group connections
    df_gamm_ms<-NULL
    for (group in list_roi_grp){
      list_roi_within_grp<-as.character(df_roi[df_roi$group==group,"id"])
      df_gamm_ms_add<-data_gamm$df_out_gamm[which(is.element(as.character(data_gamm$df_out_gamm[,"from"]),list_roi_within_grp)
                                            & is.element(as.character(data_gamm$df_out_gamm[,"to"]),list_roi_within_grp)),]
      df_gamm_ms_add<-cbind(group=group,df_gamm_ms_add)
      df_gamm_ms<-rbind(df_gamm_ms,df_gamm_ms_add)
    }
    
    # Combine within-group ROI-wise GAMM results and between-group GAMM results
    df_gamm_ms<-rbind(df_gamm_ms,cbind(group="group",data_gamm_grp$df_out_gamm))
    
    # Calculate multiple comparison-corrected p values
    df_plot_gamm_ms<-add_mltcmp(df_gamm_ms,df_roi_grp,analysis="grp",atlas,list_mod,list_plot,
                                calc_seed_level=F)
    write.csv(df_plot_gamm_ms,
              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-ms_gamm_plt.csv",sep="")),row.names = F)
    
    # Split data into ROI-wise and group-wise GAMM results, graphical output
    for (group in list_roi_grp){
      df_plot_gamm_ms_split<-df_plot_gamm_ms[df_plot_gamm_ms$group==group,-1]
      df_roi_split<-df_roi[df_roi$group==group,]
      label_analysis<-paste("ms_grp-",group,sep="")
      plot_gamm_fc(df_plot_gamm_ms_split,df_roi_split,analysis=label_analysis,atlas,list_mod,list_plot,
                   list_type_p_,thr_p,paths_)
    }
    df_plot_gamm_ms_split<-df_plot_gamm_ms[df_plot_gamm_ms$group=="group",-1]
    plot_gamm_fc(df_plot_gamm_ms_split,df_roi_grp,analysis="ms_grp-group",atlas,list_mod,list_plot,
                 list_type_p_,thr_p,paths_)

  }
  print('Finished gamm_fc().')
}


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

# Core function for parallelization of fp_fc()
fp_fc_core<-function(data_zr){
  measure<-"fc"
  group_1<-data_zr$group[[1]]
  group_2<-data_zr$group[[2]]
  df_zr<-data_zr$df_zr
  df_ses_subj<-data_zr$df_ses_subj
  n_edge<-dim(df_zr)[1]
  
  # Calculate correlation matrix
  data_fingerprint<-func_cor(input=df_zr)
  df_fp_subnet<-data_fingerprint$cor_flat
  
  # Rename correlation matrix to sessions and subjects
  df_fp_subnet$from_ses<-df_fp_subnet$from_ID_pnTTC<-df_fp_subnet$to_ses<-df_fp_subnet$to_ID_pnTTC<-NA
  for (i in seq(dim(df_fp_subnet)[1])){
    from_id<-df_fp_subnet[[i,"from"]]
    to_id<-df_fp_subnet[[i,"to"]]
    df_fp_subnet[[i,"from_ses"]]<-df_ses_subj[[from_id,"ses"]]
    df_fp_subnet[[i,"from_ID_pnTTC"]]<-df_ses_subj[[from_id,"ID_pnTTC"]]
    df_fp_subnet[[i,"to_ses"]]<-df_ses_subj[[to_id,"ses"]]
    df_fp_subnet[[i,"to_ID_pnTTC"]]<-df_ses_subj[[to_id,"ID_pnTTC"]]
  }
  df_fp_subnet$measure<-measure
  df_fp_subnet$group_1<-group_1
  df_fp_subnet$group_2<-group_2
  df_fp_subnet<-df_fp_subnet[c("measure","group_1","group_2","from_ses","from_ID_pnTTC","to_ses","to_ID_pnTTC","r","z_r")]
  
  # rbind to output dataframe
  #df_fp<-rbind(df_fp,df_fp_subnet)
  
  # Prepare dataframe for fingerprint correlation plot
  df_fp_plot<-data_fingerprint$cor
  list_name_subj_ses<-paste(sprintf("%05d",df_ses_subj$ID_pnTTC),as.character(df_ses_subj$ses),sep="_")
  colnames(df_fp_plot)<-rownames(df_fp_plot)<-list_name_subj_ses
  
  # Heatmap plot of fp correlation matrix
  plot_fp_heatmap<-plot_cor_heatmap(input=df_fp_plot)
  suppressMessages(plot_fp_heatmap<-(plot_fp_heatmap
                                     + scale_fill_gradientn(colors = matlab.like2(100),name="r")
                                     + ggtitle(paste("FP Cor,",atlas,measure,group_1,group_2,sep=" "))
                                     + theme(plot.title = element_text(hjust = 0.5),
                                             axis.title=element_blank())))
  
  # Save heatmap plot
  ggsave(paste("atl-",atlas,"_msr-",measure,"_grp1-",group_1,"_grp2-",group_2,"_fp.eps",sep=""),plot=plot_fp_heatmap,device=cairo_ps,
         path=file.path(paths_$output,"output"),dpi=300,height=10,width=10,limitsize=F)
  
  return(df_fp_subnet)
}

# Main function for fingerprint computing
fp_fc<-function(paths_=paths,
                list_atlas_=list_atlas,
                key_roigroup="group_3"){
  print("Starting fp_fc().")
  nullobj<-func_createdirs(paths_)
  dict_roi<-func_dict_roi(paths_)
  dict_roi<-data.frame(id=as.character(dict_roi$id),group=as.character(dict_roi[,key_roigroup]),stringsAsFactors = F)
  
  for (atlas in list_atlas_){
    # Load connection data
    df_conn<-read.csv(file.path(paths_$input,"output",paste("atl-",atlas,"_fc.csv",sep="")))
    df_edge<-df_conn[which(df_conn$ID_pnTTC==df_conn[1,"ID_pnTTC"]),]
    df_edge<-df_edge[which(df_edge$ses==df_edge[1,"ses"]),c("from","to"),]
    df_edge$from<-as.character(df_edge$from)
    df_edge$to<-as.character(df_edge$to)
    
    # Examine existing subject IDs and sessions in connection data
    df_ses_subj<-data.frame(matrix(nrow=0,ncol=2))
    colnames(df_ses_subj)<-c("ses","ID_pnTTC")
    list_ses_exist <- sort(unique(df_conn$ses))
    for (ses in list_ses_exist){
      df_ses_subj<-rbind(df_ses_subj,
                         data.frame(ses=ses,ID_pnTTC=sort(unique(df_conn[df_conn$ses==ses,ID_pnTTC]))))
    }
    
    # Add node subgroup column to df_edge
    df_edge<-left_join(df_edge,dict_roi,by=c("from"="id"))
    colnames(df_edge)[colnames(df_edge)=="group"]<-"from_group"
    df_edge<-left_join(df_edge,dict_roi,by=c("to"="id"))
    colnames(df_edge)[colnames(df_edge)=="group"]<-"to_group"
    
    # List groups of existing nodes
    list_group<-sort(unique(c(df_edge[,"from_group"],df_edge[,"to_group"])))
    if (!("whole" %in% list_group)){
      list_group<-c("whole",list_group)
    }
    n_group<-length(list_group)
    print(paste("Atlas: ",atlas, ", ", as.character(n_group)," groups:",sep=""))
    print(list_group)
    
    # Split and combine z_r data for each subgroup of networks for parallel computing
    list_data_zr<-list()
    
    for (idx_group_1 in seq(n_group)){
      for (idx_group_2 in seq(idx_group_1,idx_group_2)){
        group_1<-list_group[idx_group_1]
        group_2<-list_group[idx_group_2]
        # 1. whole <-> whole
        # 2. (whole-group_2) <-> group_2 and group_2 <-> (whole-group_2)
        # 3. group_1 <-> group_2 and group_2 <-> group_1 (including group_1=group_2)
        if (group_1=="whole"){
          if (group_2=="whole"){
            # 1. whole <-> whole
            df_edge_group<-df_edge
          }else{
            # 2. (whole-group_2) <-> group_2 and group_2 <-> (whole-group_2)
            df_edge_group<-rbind(df_edge[df_edge$from_group!=group_2 & df_edge$to_group==group_2,],
                                 df_edge[df_edge$from_group==group_2 & df_edge$to_group!=group_2,])
          }
        }else{
          # 3. group_1 <-> group_2 and group_2 <-> group_1 (including group_1=group_2)
          df_edge_group<-rbind(df_edge[df_edge$from_group==group_1 & df_edge$to_group==group_2,],
                               df_edge[df_edge$from_group==group_2 & df_edge$to_group==group_1,])
          
        }
        
        n_edge_group<-dim(df_edge_group)[1]
        
        if (n_edge_group<5){
          print(paste("Atlas: ",atlas,", Group: ",group_1,"<->",group_2, ", Edges: ",as.character(n_edge_group)," < 5, fp calculation skipped.",sep=""))
        }else{
          # Create combined dataframe of Z-transformed correlation coefficients
          # according to pre-calculated edge and subject data
          df_conn_cbind<-data.frame(matrix(nrow=n_edge_group,ncol=0))
          
          for (idx_subj in seq(nrow(df_ses_subj))){
            df_conn_subj<-df_conn[df_conn$ID_pnTTC==df_ses_subj[idx_subj,ID_pnTTC]
                                  & df_conn$ses==df_ses_subj[idx_subj,ses],]
            df_conn_edge<-data.frame(matrix(nrow=n_edge_group,ncol=1))
            for (idx_edge in seq(n_edge_group)){
              df_conn_edge[idx_edge,,1]<-df_conn_subj[df_conn_subj$from_group==df_edge_group[idx_edge,from_group]
                                                      & df_conn_subj$to_group==df_edge_group[idx_edge,to_group],"z_r"]
              df_conn_cbind<-cbind(df_conn_cbind,df_conn_edge)
            }
          }
          colnames(df_conn_cbind)<-as.character(seq(ncol(df_conn_cbind)))
          rownames(df_conn_cbind)<-NULL
          
          list_data_zr<-c(list_data_zr,list(list("group"=c(group_1,group_2),"df_zr"=df_conn_cbind,"df_ses_subj"=df_ses_subj)))
        }
      }
    }
    
    # Parallel fingerprint correlation computing over groups of subnetworks
    n_cluster<-min(floor(detectCores()*3/4),length(list_data_zr))
    clust<-makeCluster(n_cluster)
    clusterExport(clust,
                  varlist=c("paths_","atlas","func_cor",
                            "plot_cor_heatmap","rcorr","FisherZ","rownames_to_column","gather",
                            "ggplot","aes","geom_tile","scale_fill_gradientn",
                            "matlab.like2","scale_y_discrete","scale_x_discrete",
                            "theme_light","theme","element_text","element_blank",
                            "ggtitle","ggsave"),
                  envir=environment())
    list_df_fp<-parLapply(clust,list_data_zr,fp_fc_core)
    stopCluster(clust)
    
    # Output dataframe
    df_fp<-NULL
    for (df_fp_subnet in list_df_fp){
      if (!is.null(df_fp_subnet)){
        df_fp<-rbind(df_fp,df_fp_subnet)
      }
    }
    
    # Save fingerprint correlation
    write.csv(df_fp,file.path(paths_$output,"output",paste("atl-",atlas,"_fp.csv",sep="")),row.names=F)
  }
  print("Finished fp_fc().")
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

