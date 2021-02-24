source('C:/Users/NICT_WS/GitHub/MRI_Analysis/analyze/connection.R')

####

paths_=paths
list_atlas_=list_atlas
param<-param_gam_fc_diff

####

print("Starting gam_fc_diff().")
nullobj<-func_createdirs(paths_,str_proc="gam_fc_diff()",copy_log=T,list_param=param)
memory.limit(1000000)

####

atlas<-list_atlas_[1]

####

print(paste("Preparing FC data: ",atlas,sep=""))
data_fc<-prep_data_fc(paths_,atlas,param$key_group,include_diff=T,abs_nfc=param$abs_nfc)
data_fc$df_edge$id_edge<-seq(nrow(data_fc$df_edge))
data_fc$df_edge_grp$id_edge<-seq(nrow(data_fc$df_edge_grp))

####

idx_tanner<-names(param$list_tanner)[1]

####

print(paste("Atlas: ",atlas,", Tanner type: ",param$list_tanner[[idx_tanner]][["label"]],sep=""))
list_covar<-param$list_covar_tanner
list_covar[["tanner"]]<-param$list_tanner[[idx_tanner]]

####

paths<-paths_
list_wave<-param$list_wave
list_sex<-list(1,2)
subset_subj<-param$subset_subj
list_mod<-param$list_mod_tanner
list_plot<-param$list_plot_tanner
idx_var<-idx_tanner
list_p<-param$list_p
calc_parallel=F
test_mod=F
#test_mod=T

####

#print("Calculated result already exists.")
#df_gamm<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm.csv",sep=""))))
#df_gamm_grp<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_grp.csv",sep=""))))
#df_anova<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova.csv",sep=""))))
#df_anova_grp<-as.data.frame(fread(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_anova_grp.csv",sep=""))))

####

df_plot<-read.csv(file.path(paths$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_plot.csv",sep="")))

####

df_net<-df_node<-df_size_net<-NULL
#list_output<-list()
for (idx_mod in param$param_nbs$list_mod){
  for (idx_term in param$param_nbs$list_term){
    for (idx_sex in list_sex){
      df_sign<-df_plot[df_plot$p_type==param$param_nbs$p_cdt_type & df_plot$p_threshold==param$param_nbs$p_cdt_threshold
                       & df_plot$model==idx_mod & df_plot$term==idx_term & df_plot$sex==idx_sex,]
      data_bfs<-func_bfs(df_sign)
      if(length(data_bfs$list_network)>0){
        for (idx_net in seq(length(data_bfs$list_network))){
          network<-data_bfs$list_network[[idx_net]]
          #list_output<-c(list_output,
          #               list(plot_sex_diff_fc(paths,network$df_edge,atlas,wave,df_roi,df_grp,
          #                                     model,plot,sex,title_plot,title_sex,idx_net)))
          df_head<-data.frame(model=idx_mod,term=idx_term,sex=idx_sex)
          df_net<-rbind(df_net,data.frame(id_net=idx_net,network$df_edge))
          df_node_add<-inner_join(data.frame(id_net=idx_net,node=network$list_node),data_fc$df_roi,by=c("node"="id"))
          df_node_add<-dplyr::rename(df_node_add,"label_node"="label","group_node"="group")
          df_node<-rbind(df_node,cbind(df_head,df_node_add))
          df_size_net<-rbind(df_size_net,cbind(df_head,data.frame(id_net=idx_net,size=network$size_net)))
        }
      }
    }
  }
}
