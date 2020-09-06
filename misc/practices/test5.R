paths_=paths
subset_subj_=model_fp_subset_subj
list_atlas_=list_atlas
list_wave_=list_wave
skip_ancova=T
list_covar_tanner_=model_fp_list_covar_tanner
list_tanner_=model_fp_list_tanner
list_mod_tanner_=model_fp_list_mod_tanner
list_graph_tanner_=model_fp_list_graph_tanner
list_strat_tanner_=model_fp_list_strat_tanner
list_covar_hormone_=model_fp_list_covar_hormone
list_hormone_=model_fp_list_hormone
list_mod_hormone_=model_fp_list_mod_hormone
list_graph_hormone_=model_fp_list_graph_hormone



print("Starting model_fp_multi().")
nullobj<-func_createdirs(paths_,str_proc="model_fp_multi()")

df_out_lm<-df_out_aic<-NULL
# Loop over clinical variables



idx_tanner<-names(list_tanner_)[[1]]



print(paste("Tanner type: ",list_tanner_[[idx_tanner]][["label"]],sep=""))
list_covar<-list_covar_tanner_
list_covar[["tanner"]]<-list_tanner_[[idx_tanner]]




list_covar_=list_covar
list_mod_=list_mod_tanner_
list_graph_=list_graph_tanner_
suffix=paste("var-",idx_tanner,sep="")



# Load and subset clinical data according to specified subsetting condition and covariate availability
print('Loading clinical data.')
data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,
                                   list_covar=list_covar_,rem_na_clin=T,prefix=suffix,print_terminal=F)
df_clin<-data_clin$df_clin
colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"

list_src_glm<-list_src_ancova<-df_src_glm_bind<-NULL



atlas<-list_atlas_[1]



# Load fingerprint data]
print(paste("Loading FP, atlas:",atlas,sep=" "))
df_fp<-as.data.frame(fread(file.path(paths_$input,"output",paste("atl-",atlas,"_fp.csv",sep=""))))

list_measure<-as.character(sort(unique(df_fp$measure)))



measure<-list_measure[1]



print(paste("Preparing dataset, atlas: ",atlas,", measure: ",measure,sep=""))
df_fp_meas<-df_fp[df_fp$measure==measure,]

# Create list of subjects who meet subsetting condition and whose MRI data exist
list_ses_exist <- sort(unique(c(df_fp_meas$from_ses,df_fp_meas$to_ses)))
list_id_subj_exist<-list()
for (ses in list_ses_exist){
  id_subj_exist_ses<-sort(unique(c(df_fp_meas[df_fp_meas$from_ses==ses,'from_ID_pnTTC'],
                                   df_fp_meas[df_fp_meas$to_ses==ses,'to_ID_pnTTC'])))
  id_subj_subset_ses<-df_clin[df_clin$ses==ses,"ID_pnTTC"]
  id_subj_exist_ses<-intersect(id_subj_exist_ses,id_subj_subset_ses)
  list_id_subj_exist[[as.character(ses)]]<-sort(id_subj_exist_ses)
}

# Identify subjects with longitudinal data
list_id_subj_exist_twice<-sort(intersect(list_id_subj_exist[["1"]],
                                         list_id_subj_exist[["2"]]))
n_id_subj_exist_twice<-length(list_id_subj_exist_twice)

# Create list of existing ROI subgroups
list_group<-sort(unique(c(as.character(df_fp_meas$group_1),as.character(df_fp_meas$group_2))))
if ("whole" %in% list_group){
  list_group<-c("whole",list_group[list_group!="whole"])
}
#list_group<-"whole      # Disable subgroup-wise analysis

n_group<-length(list_group)
for (idx_group_1 in seq(n_group)){
  group_1<-list_group[idx_group_1]
  df_fp_meas_grp1<-df_fp_meas[df_fp_meas$group_1==group_1,]
  for (idx_group_2 in seq(idx_group_1,n_group)){
    group_2<-list_group[idx_group_2]
    df_fp_meas_grp2<-df_fp_meas_grp1[df_fp_meas_grp1$group_2==group_2,]
    df_cor_fp<-data.frame(ID_pnTTC=list_id_subj_exist_twice)
    for (id_subj in list_id_subj_exist_twice){
      if (any(df_fp_meas_grp2$from_ID_pnTTC==id_subj & df_fp_meas_grp2$to_ID_pnTTC==id_subj)){
        df_cor_fp[df_cor_fp$ID_pnTTC==id_subj,"value"]<-df_fp_meas_grp2[df_fp_meas_grp2$from_ID_pnTTC==id_subj
                                                                        & df_fp_meas_grp2$to_ID_pnTTC==id_subj,
                                                                        "z_r"]
      }else{
        df_cor_fp[df_cor_fp$ID_pnTTC==id_subj,"value"]<-NA
      }
    }
    
    # Subset those without longitudinal fp correlation
    list_id_subj_nonna<-df_cor_fp[!is.na(df_cor_fp$value),"ID_pnTTC"]
    df_cor_fp<-df_cor_fp[df_cor_fp$ID_pnTTC %in% list_id_subj_nonna,]
    n_id_subj_exist_twice<-length(list_id_subj_nonna)
    #print(paste("Atlas: ",atlas,", Measure: ",measure,", Group: ",group_1," and ",group_2,
    #            ", ",as.character(n_id_subj_exist_twice)," subjects with longitudinal non-NA data.",sep=""))
    
    
    if (n_id_subj_exist_twice>0){
      # Create dataframe for GLM analysis
      df_src_glm<-func_clinical_data_join(df_src=df_clin,
                                          list_id_subj=list_id_subj_nonna,
                                          list_covar=list_covar_)
      df_src_glm<-inner_join(df_src_glm,df_cor_fp,by="ID_pnTTC")
      df_src_glm$ID_pnTTC<-as.factor(df_src_glm$ID_pnTTC)
      df_src_glm$sex<-as.factor(df_src_glm$sex)
      
      df_src_glm_bind<-rbind(df_src_glm_bind,
                             cbind(atlas=atlas,measure=measure,group_1=group_1,group_2=group_2,
                                   df_src_glm))
      
      list_src_glm<-c(list_src_glm,
                      list(list("df_src"=df_src_glm,"atlas"=atlas,"measure"=measure,
                                "group_1"=group_1,"group_2"=group_2)))
    }
    # Calculate GLM
    #out_glm<-glm_core(df_src=df_join_grp,atlas,measure,group_1,group_2,
    #                  list_mod_,list_graph_,list_covar_,paths_)
    #df_out_lm<-rbind(df_out_lm,out_glm$df_out_lm_add)
    #df_out_aic<-rbind(df_out_aic,out_glm$df_out_aic_add)
    
    # Prepare ANCOVA dataset for later parallel computing
    if (!skip_ancova){
      print(paste("Atlas: ",atlas,", Measure: ",measure,", Group: ",group_1," and ",group_2,", ANCOVA preparation.",  sep=""))
      # Create list of input dataframes for parallel ANCOVA calculation
      for (group_tanner in names(list_strat_tanner_)){
        # group by longitudinal Tanner stage
        df_join_grp_tanner<-df_join_grp
        for (ses in c(1,2)){
          list_strat_tanner_ses<-names(list_strat_tanner_[[group_tanner]][[as.character(ses)]])
          for (label_tanner in list_strat_tanner_ses){
            list_strat_tanner_ses_group<-list_strat_tanner_[[group_tanner]][[as.character(ses)]][[label_tanner]]
            #print(list_strat_tanner_ses_group)
            df_join_grp_tanner[df_join_grp_tanner[[paste('ses',as.character(ses),'_tanner',sep='')]] %in% list_strat_tanner_ses_group,
                               paste('ses',as.character(ses),'_tanner_label',sep='')]<-label_tanner
          }
        }
        df_join_grp_tanner$long_tanner<-paste(as.character(df_join_grp_tanner$ses1_tanner_label),
                                              as.character(df_join_grp_tanner$ses2_tanner_label),sep="_")
        df_join_grp_tanner$long_tanner<-as.factor(df_join_grp_tanner$long_tanner)
        
        list_sex<-list("all"=c(1,2),"male"=1,"female"=2)
        
        for (idx_sex in names(list_sex)){
          df_join_grp_tanner_sex<-df_join_grp_tanner[df_join_grp_tanner$sex %in% list_sex[[idx_sex]],]
          list_src_ancova<-c(list_src_ancova,
                             list(list("atlas"=atlas,"measure"=measure,"group_network"=c(group_1,group_2),
                                       "group_tanner"=group_tanner,"group_tanner_content"=list_strat_tanner_[[group_tanner]],
                                       "group_sex"=idx_sex,"df_src_ancova"=df_join_grp_tanner_sex)))
        }
      }
    }
  }
} # Finished looping over ROI subgroups

# Parallel GLM calculation
n_cluster<-floor(detectCores()*3/4)
print(paste("Calculating GLM/GAM in parallel, threads: ",as.character(n_cluster),sep=""))
clust<-makeCluster(n_cluster)
clusterExport(clust,
              varlist=c("list_mod_","list_graph_","list_covar_","paths_","as.numeric.factor",
                        "gam","as.formula","summary.gam","plot_gamm","ggplot",
                        "predict.gam","geom_line","aes","geom_ribbon","geom_point",
                        "theme_light","element_text","ggtitle",
                        "xlab","ylab","theme","ggsave"),
              envir=environment())
list_dst_glm<-pblapply(list_src_glm,glm_core,cl=clust)
stopCluster(clust)

df_out_lm<-df_out_aic<-NULL
for (dst_glm in list_dst_glm){
  df_out_lm<-rbind(df_out_lm,dst_glm$df_out_lm_add)
  df_out_aic<-rbind(df_out_aic,dst_glm$df_out_aic_add)
}

# GLM/GAM Data saving
rownames(df_out_lm)<-rownames(df_out_aic)<-NULL
write.csv(df_out_lm, file.path(paths_$output,"output",paste(suffix,"_fp_glm.csv",sep="")),row.names = F)
write.csv(df_out_aic,file.path(paths_$output,"output",paste(suffix,"_fp_glm_aic.csv",sep="")),row.names = F)


# group-wise GLM/GAM heatmap visualization
list_atlas=list_atlas_
list_mod=list_mod_


for (atlas in list_atlas){
  df_out_lm_subset<-df_out_lm[df_out_lm$atlas==atlas,]
  list_measure<-as.character(sort(unique(df_out_lm_subset$measure)))
  for (measure in list_measure){
    df_out_lm_subset<-df_out_lm[df_out_lm$atlas==atlas & df_out_lm$measure==measure,]
    for (model in names(list_mod)){
      for (idx_sex in c(1,2)){
        df_out_lm_subset<-df_out_lm[df_out_lm$atlas==atlas
                                    & df_out_lm$measure==measure
                                    & df_out_lm$sex==idx_sex
                                    & df_out_lm$model==model,]
        list_term<-sort(unique(as.character(df_out_lm_subset$term)))
        list_term<-list_term[list_term!="(Intercept)"]
        for (term in list_term){
          df_out_lm_subset<-df_out_lm[df_out_lm$atlas==atlas
                                      & df_out_lm$measure==measure
                                      & df_out_lm$sex==idx_sex
                                      & df_out_lm$model==model
                                      & df_out_lm$term==term,]
          list_group<-sort(unique(c(as.character(df_out_lm_subset$group_1),
                                    as.character(df_out_lm_subset$group_2))))
          list_group<-list_group[list_group!="whole"]
          n_group<-length(list_group)
          mat_static<-data.frame(matrix(nrow=n_group,ncol=n_group))
          mat_pval<-data.frame(matrix(nrow=n_group,ncol=n_group))
          colnames(mat_static)<-rownames(mat_static)<-colnames(mat_pval)<-rownames(mat_pval)<-list_group
          for (idx_group_1 in seq(n_group)){
            for (idx_group_2 in seq(idx_group_1,n_group)){
              group_1<-list_group[idx_group_1]
              group_2<-list_group[idx_group_2]
              df_out_lm_subset<-df_out_lm[df_out_lm$atlas==atlas
                                          & df_out_lm$measure==measure
                                          & df_out_lm$sex==idx_sex
                                          & df_out_lm$model==model
                                          & df_out_lm$term==term
                                          & df_out_lm$group_1==group_1
                                          & df_out_lm$group_2==group_2,c("estimate","F","p")]
              if (nrow(df_out_lm_subset)>0){
                name_static<-colnames(df_out_lm_subset[,c("estimate","F")])[!is.na(df_out_lm_subset[1,c("estimate","F")])]
                mat_static[group_1,group_2]<-mat_static[group_2,group_1]<-df_out_lm_subset[[1,name_static]]
                #mat_pval[group_1,group_2]<-mat_pval[group_2,group_1]<-df_out_lm_subset[[1,"p"]]
                pval<-df_out_lm_subset[[1,"p"]]
                if (pval<0.001){
                  mat_pval[group_1,group_2]<-mat_pval[group_2,group_1]<-"**"
                }else if (pval<0.05){
                  mat_pval[group_1,group_2]<-mat_pval[group_2,group_1]<-"*"
                }else{
                  mat_pval[group_1,group_2]<-mat_pval[group_2,group_1]<-""
                }
              }else{
                mat_static[group_1,group_2]<-mat_static[group_2,group_1]<-NA
                mat_pval[group_1,group_2]<-mat_pval[group_2,group_1]<-""
              }
            }
          }
          
          #df_pval<-rownames_to_column(mat_pval,"row")
          #df_pval<-gather(df_pval,key=column,value=p,2:ncol(df_pval))
          if (idx_sex==1){
            label_sex<-"male"
          }else{
            label_sex<-"female"
          }
          #mat_pval<-round(mat_pval,3)
          plot_stat<-plot_cor_heatmap(mat_static,mat_pval)
          suppressMessages(plot_stat<-(plot_stat
                                       + scale_fill_gradientn(colors = matlab.like2(100),
                                                              lim=c(-max(max(mat_static),-min(mat_static)),max(max(mat_static),-min(mat_static))),
                                                              name=name_static)
                                       + ggtitle(paste("GLM-FP,",atlas,measure,label_sex,model,term,sep=" "))
                                       + theme(plot.title = element_text(hjust = 0.5),
                                               axis.title=element_blank())))
          ggsave(paste("atl-",atlas,"_msr-",measure,"_sex-",label_sex,"_mod-",model,
                       "_plt-",term,"_fp_glm_heatmap.eps",sep=""),plot=plot_stat,device=cairo_ps,
                 path=file.path(paths_$output,"output","plot"),dpi=300,height=5,width=5,limitsize=F)
        }
      }
    }
  }
}
