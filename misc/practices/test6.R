gamm_fc_list_tanner<-gamm_fc_list_tanner[1]
gamm_fc_list_hormone<-gamm_fc_list_hormone[1]
list_atlas<-list_atlas[1]

####

paths_=paths
subset_subj_=gamm_fc_subset_subj
list_wave_=list_wave
list_atlas_=list_atlas
key_group_='group_3'
list_covar_tanner_=gamm_fc_list_covar_tanner
list_tanner_=gamm_fc_list_tanner
list_mod_tanner_=gamm_fc_list_mod_tanner
list_plot_tanner_=gamm_fc_list_plot_tanner
list_covar_hormone_=gamm_fc_list_covar_hormone
list_hormone_=gamm_fc_list_hormone
list_mod_hormone_=gamm_fc_list_mod_hormone
list_plot_hormone_=gamm_fc_list_plot_hormone

####

print("Starting gamm_fc_multi().")
nullobj<-func_createdirs(paths_,str_proc="gamm_fc_multi()",copy_log=T)
dict_roi <- func_dict_roi(paths_)

####

atlas<-list_atlas[1]

####
print(paste("Preparing FC data: ",atlas,sep=""))
df_fc<-as.data.frame(fread(file.path(paths_$input,"output",
                                     paste("atl-",atlas,"_fc.csv",sep=""))))
df_fc<-df_fc[df_fc$ses!="2-1",]

# Prepare dataframe of ROIs
list_roi<-sort(unique(c(as.character(df_fc$from),as.character(df_fc$to))))
df_roi<-dict_roi[is.element(dict_roi$id,list_roi),c("id","label",key_group_)]
colnames(df_roi)[colnames(df_roi)==key_group_]<-"group"

# prepare dataframe of group-wise FC averages
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
  list_ses<-list_ses[list_ses!="2-1"]
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

####

idx_tanner<-names(list_tanner_)[1]

####

print(paste("Atlas: ",atlas,", Tanner type: ",list_tanner_[[idx_tanner]][["label"]],sep=""))
list_covar<-list_covar_tanner_
list_covar[["tanner"]]<-list_tanner_[[idx_tanner]]

####

list_mod<-list_mod_tanner_
list_plot<-list_plot_tanner_
idx_var<-idx_tanner
list_type_p_=list_type_p
thr_p_=thr_p

####

# Prepare clinical data
data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar,rem_na_clin=T,
                                   prefix=paste("var-",idx_var,sep=""),print_terminal=F)
df_clin<-data_clin$df_clin

# Join FC and clinical data
df_fc$z_r[which(is.nan(df_fc$z_r))]<-0
colnames(df_fc)[colnames(df_fc)=="z_r"]<-"value"
colnames(df_fc)[colnames(df_fc)=="ses"]<-"wave"
df_fc<-df_fc[,c(-which(colnames(df_fc)=="r"),
                -which(colnames(df_fc)=="p"))]
df_clin$wave<-as.character(df_clin$wave)
df_join<-inner_join(df_fc,df_clin,by=c('ID_pnTTC','wave'))
for (key in c('ID_pnTTC','wave','sex')){
  if (key %in% colnames(df_join)){
    df_join[,key]<-as.factor(df_join[,key])
  }
}
df_join$value<-as.numeric.factor(df_join$value)

# Calculate ROI-wise GAMM of FC
data_gamm<-iterate_gamm(df_join,df_roi,list_mod,calc_parallel=F,calc_identical=F)
df_plot<-add_mltcmp(data_gamm$df_out_gamm,df_roi,list_mod,list_plot,calc_seed_level=F)

# Join group-wise FC and clinical data
colnames(df_fc_grp)[colnames(df_fc_grp)=="z_r"]<-"value"
colnames(df_fc_grp)[colnames(df_fc_grp)=="ses"]<-"wave"
df_fc_grp$ID_pnTTC<-as.character(as.numeric.factor(df_fc_grp$ID_pnTTC))
df_fc_grp$wave<-as.character(as.numeric.factor(df_fc_grp$wave))
df_clin$ID_pnTTC<-as.character(as.numeric.factor(df_clin$ID_pnTTC))
df_clin$wave<-as.character(as.numeric.factor(df_clin$wave))
df_join_grp<-inner_join(df_fc_grp,df_clin,by=c('ID_pnTTC','wave'))
for (key in c('ID_pnTTC','wave','sex')){
  if (key %in% colnames(df_join_grp)){
    df_join_grp[,key]<-as.factor(df_join_grp[,key])
  }
}
df_join_grp$value<-as.numeric.factor(df_join_grp$value)

# Calculate Group-wise GAMM of FC
list_group<-unique(as.character(df_roi$group))
df_grp<-data.frame(id=list_group,label=str_to_title(gsub("_"," ",as.character(list_group))))
data_gamm_grp<-iterate_gamm(df_join_grp,df_grp,list_mod,calc_parallel=F,calc_identical=T)
df_plot_grp<-add_mltcmp(data_gamm_grp$df_out_gamm,df_grp,list_mod,list_plot,calc_seed_level=F)

# Save results
write.csv(data_gamm$df_out_gamm,
          file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm.csv",sep="")),row.names = F)
write.csv(data_gamm$df_out_aic,
          file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic.csv",sep="")),row.names = F)
write.csv(df_plot,
          file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_plot.csv",sep="")),row.names = F)
write.csv(data_gamm_grp$df_out_gamm,
          file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_grp.csv",sep="")),row.names = F)
write.csv(data_gamm_grp$df_out_aic,
          file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_aic_grp.csv",sep="")),row.names = F)
write.csv(df_plot_grp,
          file.path(paths_$output,"output","temp",paste("atl-",atlas,"_var-",idx_var,"_gamm_plot_grp.csv",sep="")),row.names = F)

# Graphical output of ROI- and group-wise GAMM of FC
plot_gam_fc(paths_,df_gam=df_plot,df_gam_grp_sign=df_plot_grp,df_gam_grp_abs=NULL,atlas,
            list_mod,list_plot,list_type_p=list_type_p_,thr_p=thr_p_,waves=NULL,idx_var)

####

df_join=df_join_grp
df_roi=df_grp
list_mod_=list_mod
calc_parallel=F
calc_identical=T

####

list_roi<-df_roi$id
list_sex<-sort(unique(as.numeric.factor(df_join$sex)))

####

list_dst_gamm<-list()
if (calc_identical){
  list_id_from<-list_roi
  n_edge<-length(list_roi)*(length(list_roi)-1)/2+length(list_roi)
}else{
  list_id_from<-list_roi[-length(list_roi)]
  n_edge<-length(list_roi)*(length(list_roi)-1)/2
}
interval_refresh<-max(1,floor(n_edge/300))
progressbar <- txtProgressBar(min = 1, max = n_edge, style = 3)
idx_edge<-0
for (id_from in list_id_from){
  df_join_from<-df_join[df_join$from==id_from,]
  label_from<-as.character(df_roi[df_roi$id==id_from,"label"])
  if (calc_identical){
    list_id_to<-list_roi[seq(which(list_roi==id_from),length(list_roi))]
  }else{
    list_id_to<-list_roi[seq(which(list_roi==id_from)+1,length(list_roi))]
  }
  for(id_to in list_id_to){
    idx_edge<-idx_edge+1
    label_to<-as.character(df_roi[df_roi$id==id_to,"label"])
    #print(paste("Calculating",label_from,"and",label_to,sep=" "))
    df_src<-df_join_from[df_join_from$to==id_to,]
    list_dst_gamm<-c(list_dst_gamm,
                     list(gamm_core(list("df_src"=df_src,"id_from"=id_from,"id_to"=id_to,
                                         "label_from"=label_from,"label_to"=label_to,
                                         "list_mod"=list_mod_,"list_sex"=list_sex))))
    if ((idx_edge %% interval_refresh) ==0){
      Sys.sleep(0.1)
      setTxtProgressBar(progressbar, idx_edge)
    }
  }
}
close(progressbar)

####

data_src<-list("df_src"=df_src,"id_from"=id_from,"id_to"=id_to,
               "label_from"=label_from,"label_to"=label_to,
               "list_mod"=list_mod_,"list_sex"=list_sex)

####

df_src<-data_src$df_src

list_mod_<-data_src$list_mod
list_sex<-data_src$list_sex

#list_sex<-sort(unique(as.numeric.factor(df_src$sex)))
df_out_aic_add<-df_out_gamm_add<-data.frame()
for (idx_mod in names(list_mod_)){
  for (idx_sex in list_sex){
    df_src_sex<-df_src[df_src$sex==idx_sex,]
    df_src_sex$value<-as.numeric(df_src_sex$value)
    mod<-try(gam(as.formula(list_mod_[[idx_mod]]),data=df_src_sex,method="REML"), silent=F)
    if (class(mod)[1]!="try-error"){
      p_table<-summary.gam(mod)$p.table
      if (is.null(summary.gam(mod)$s.table)){
        df_out_gamm_add_add<-data.frame(term=rownames(p_table),estimate=p_table[,'Estimate'],
                                        se=p_table[,'Std. Error'],F=NA,t=p_table[,'t value'],
                                        p=p_table[,'Pr(>|t|)'])
      }else{
        s_table<-summary.gam(mod)$s.table
        df_out_gamm_add_add<-rbind(data.frame(term=rownames(p_table),estimate=p_table[,'Estimate'],
                                              se=p_table[,'Std. Error'],F=NA,t=p_table[,'t value'],
                                              p=p_table[,'Pr(>|t|)']),
                                   data.frame(term=rownames(s_table),estimate=NA,se=NA,F=s_table[,'F'],
                                              t=NA,p=s_table[,'p-value']))
      }
      
      df_out_gamm_add<-rbind(df_out_gamm_add,
                             cbind(sex=idx_sex,model=idx_mod,df_out_gamm_add_add))
      df_out_aic_add<-rbind(df_out_aic_add,
                            data.frame(sex=idx_sex,model=idx_mod,aic=mod$aic,aic_best_among_models=0))
    }
  } # Finished looping over sex
}# Finished looping over model

# Compare AICs of GAMM models
df_out_aic_add_sex_rbind<-data.frame()
for (idx_sex in list_sex){
  df_out_aic_add_sex<-df_out_aic_add[df_out_aic_add$sex==idx_sex,]
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
