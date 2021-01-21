paths_=paths
subset_subj_=sex_diff_fc_subset_subj
list_wave_=list_wave
list_atlas_=list_atlas
key_group_='group_3'
list_covar_=sex_diff_fc_list_covar
list_mod_=sex_diff_fc_list_mod
list_plot_=sex_diff_fc_list_plot
list_type_p_=sex_diff_fc_list_type_p
thr_p_=sex_diff_fc_thr_p

####

print("Starting sex_diff_fc_multi().")
nullobj<-func_createdirs(paths_,str_proc="sex_diff_fc_multi()",copy_log=T)
dict_roi <- func_dict_roi(paths_)
memory.limit(1000000)

####

atlas<-list_atlas_[1]

####

print(paste("Preparing FC data: ",atlas,sep=""))
data_fc<-prep_data_fc(paths_,atlas,key_group_)

# Prepare clinical data
data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar_,rem_na_clin=T,
                                   prefix=NULL,print_terminal=F)

# Calculate ROI-wise GAMM of FC
df_join<-join_fc_clin(data_fc$df_fc,data_clin$df_clin)
data_gamm<-iterate_gamm(df_join,data_fc$df_roi,list_mod_,calc_parallel=F,calc_identical=F,list_sex=list(c(1,2)))

####

####

df_gamm_sign<-data_gamm$df_out_gamm
#df_gamm_sign<-df_gamm_sign[df_gamm_sign$term=="sex2",]
df_gamm_sign<-df_gamm_sign[df_gamm_sign$term=="sex2" & df_gamm_sign$p<thr_p_,]
df_m<-df_gamm_sign[df_gamm_sign$t<0,]
df_f<-df_gamm_sign[df_gamm_sign$t>0,]

####

#df_edge<-df_m
data_bfs<-func_bfs(df_m)
data_bfs<-func_bfs(df_f)
data_bfs<-func_bfs(data.frame())

####

func_bfs<-function(df_edge){
  df_edge_remain<-df_edge
  list_network<-list()
  list_size<-NULL
  #df_edge_new<-df_edge_remain[1,]
  while (nrow(df_edge_remain)>0){ # for each subnetwork
    list_node_todo<-list_node_net<-df_edge_remain[[1,"from"]]
    #node_orig<-df_edge_remain[[1,"from"]]
    #list_node_new<-node_orig
    df_edge_net<-data.frame()
    while (length(list_node_todo)>0){
      node_check<-list_node_todo[1]
      df_edge_new_from<-df_edge_remain[df_edge_remain$from==node_check,]
      df_edge_remain<-df_edge_remain[rownames(df_edge_remain) %nin% rownames(df_edge_new_from),]
      list_node_new_from<-df_edge_new_from[,"to"]
      df_edge_new_to<-df_edge_remain[df_edge_remain$to==node_check,]
      df_edge_remain<-df_edge_remain[rownames(df_edge_remain) %nin% rownames(df_edge_new_to),]
      list_node_new_to<-df_edge_new_to[,"from"]
      df_edge_new<-rbind(df_edge_new_to,df_edge_new_from)
      list_node_new<-c(list_node_new_from,list_node_new_to)
      
      df_edge_net<-rbind(df_edge_net,df_edge_new)
      list_node_net<-c(list_node_net,list_node_new)
      list_node_todo<-c(list_node_todo,list_node_new)
      
      list_node_todo<-list_node_todo[-1]
    }
    
    size_net<-nrow(df_edge_net)
    list_size<-c(list_size,size_net)
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

####

####

df_roi<-data_fc$df_roi
calc_parallel=F
calc_identical=F
list_sex=list(c(1,2))

####

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

####

data_src<-list_src_gamm[[1]]

####

df_src<-data_src$df_src
list_mod_<-data_src$list_mod
list_sex<-data_src$list_sex

####

gam.control(nthreads=1)

start_time <- Sys.time()

####

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
    #if (data_src$calc_parallel){
    #  mod<-try(gam(as.formula(list_mod_[[idx_mod]]),data=df_src_sex,method="REML",control=list(nthreads=1)), silent=F)
    #}else{
    #  mod<-try(gam(as.formula(list_mod_[[idx_mod]]),data=df_src_sex,method="REML"), silent=F)
    #}
    #mod<-try(gam(as.formula(list_mod_[[idx_mod]]),data=df_src_sex,method="REML"), silent=F)
    #mod<-try(bam(as.formula(list_mod_[[idx_mod]]),data=df_src_sex,method="REML"), silent=F)
    #mod<-try(gam(as.formula(list_mod_[[idx_mod]]),data=df_src_sex,method="REML",control=list(nthreads=1)), silent=F)
    #mod<-gam(as.formula(list_mod_[[idx_mod]]),data=df_src_sex,method="REML",control=list(nthreads=1))
    mod<-gam(as.formula(list_mod_[[idx_mod]]),data=df_src_sex,method="REML")
    #if (class(mod)[1]!="try-error"){
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
    #}
  } # Finished looping over sex
}# Finished looping over model

Sys.time()-start_time

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
id_from<-data_src$id_from
id_to<-data_src$id_to
label_from<-data_src$label_from
label_to<-data_src$label_to
df_out_gamm_add<-cbind(from=id_from,to=id_to,label_from=label_from,label_to=label_to,
                       df_out_gamm_add)
df_out_aic_add_sex_rbind<-cbind(from=id_from,to=id_to,label_from=label_from,label_to=label_to,
                                df_out_aic_add_sex_rbind)

####

Sys.time()-start_time

####

####

# Calculate Group-wise GAMM of FC
df_join_grp<-join_fc_clin(data_fc$df_fc_grp,data_clin$df_clin)
data_gamm_grp<-iterate_gamm(df_join_grp,data_fc$df_grp,list_mod_,calc_parallel=F,calc_identical=T,list_sex=list(c(1,2)))

# Graphical output of ROI- and group-wise GAMM of FC
plot_gam_fc(paths_,df_gam=df_plot,df_gam_grp_sign=df_plot_grp,df_gam_grp_abs=NULL,atlas,
            list_mod,list_plot,list_type_p=list_type_p_,thr_p=thr_p_,waves=NULL,idx_var)