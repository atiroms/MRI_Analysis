paths_=paths
list_atlas_=list_atlas
list_waves_=gam_fc_list_waves
subset_subj_=gam_fc_subset_subj
list_covar_tanner_=gam_fc_list_covar_tanner
list_tanner_=gam_fc_list_tanner
list_mod_tanner_=gam_fc_list_mod_tanner
list_plot_tanner_=gam_fc_list_plot_tanner
list_covar_hormone_=gam_fc_list_covar_hormone
list_hormone_=gam_fc_list_hormone
list_mod_hormone_=gam_fc_list_mod_hormone
list_plot_hormone_=gam_fc_list_plot_hormone
list_type_p_=list_type_p
thr_p_=thr_p

####

print("Starting gam_fc_cs_multi()")
nullobj<-func_createdirs(paths_,str_proc="gam_fc_cs_multi()",copy_log=T)

####

label_waves<-names(list_waves_)[[1]]

####

wave_clin<-list_waves_[[label_waves]]$wave_clin
wave_mri<-list_waves_[[label_waves]]$wave_mri
waves<-list_waves_[label_waves]
#print(paste("Clinical wave: ", wave_clin,", MRI wave: ",wave_mri,sep=""))

# Prepare subject subsetting condition (MRI QC criteria) according to specified waves
subset_subj_temp<-subset_subj_[[as.character(wave_mri)]]
subset_subj_temp<-list(subset_subj_temp)
names(subset_subj_temp)<-wave_clin

####

idx_tanner<-names(list_tanner_)[[1]]

####

print(paste("Clinical wave: ", wave_clin,", MRI wave: ",wave_mri,", Tanner type: ",list_tanner_[[idx_tanner]][["label"]],sep=""))
list_covar<-list_covar_tanner_
list_covar[["tanner"]]<-list_tanner_[[idx_tanner]]
#suffix<-paste("ses-",label_waves,"_var-",idx_tanner,sep="")

####

paths_=paths_
subset_subj_=subset_subj_temp
list_covar_=list_covar
list_atlas_=list_atlas_
list_mod_=list_mod_tanner_
list_plot_=list_plot_tanner_
key_group_='group_3'
list_type_p_=list_type_p_
thr_p_=thr_p_
waves_=waves
idx_var_=idx_tanner

####

#print("Starting gam_fc_cs().")
#nullobj<-func_createdirs(paths_,str_proc="gam_fc_cs()",copy_log=T)
dict_roi <- func_dict_roi(paths_)

# Load and subset clinical data according to specified subsetting condition and covariate availability
#print('Loading clinical data.')
label_waves<-names(waves_)
wave_clin<-waves_[[1]]$wave_clin
wave_mri<-waves_[[1]]$wave_mri
data_clin<-func_clinical_data_long(paths_,list_wave=wave_clin,subset_subj_,
                                   list_covar=list_covar_,rem_na_clin=T,
                                   prefix=paste("ses-",label_waves,"_var-",idx_var_,sep=""),
                                   print_terminal=F)
df_clin<-data_clin$df_clin

####

atlas<-list_atlas_[1]

####

# Load ROI-wise FC data
#print(paste('Loading FC data, atlas:',atlas,sep=' '))
#df_fc<-read.csv(file.path(paths_$input,'output',paste('atl-',atlas,'_fc.csv',sep='')))
df_fc<-as.data.frame(fread(file.path(paths_$input,'output',paste('atl-',atlas,'_fc.csv',sep=''))))
df_join<-join_fc_clin(df_fc,df_clin,wave_clin,wave_mri)
#write.csv(df_join,file.path(paths_$output,"output","temp",
#                            paste("atl-",atlas,"_ses-",label_waves,"_var-",idx_var_,"_src.csv",sep="")),
#          row.names=F)

# Calculate and save ROI-wise GAMM of FC
print(paste('Calculating GAM, atlas: ',atlas,sep=''))
list_roi<-sort(unique(c(as.character(df_join$from),as.character(df_join$to))))
df_roi<-dict_roi[is.element(dict_roi$id,list_roi),c("id","label",key_group_)]
df_roi$id<-as.character(df_roi$id)
colnames(df_roi)[colnames(df_roi)==key_group_]<-"group"
data_gamm<-iterate_gamm(df_join,df_roi,list_mod_,calc_parallel=T,calc_identical=F)

####

calc_parallel=T
calc_identical=F

####

list_roi<-df_roi$id
list_sex<-sort(unique(as.numeric.factor(df_join$sex)))

####

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
    list_src_gamm<-c(list_src_gamm,list(list("df_src"=df_src,"id_from"=id_from,"id_to"=id_to,
                                             "label_from"=label_from,"label_to"=label_to,
                                             "list_mod"=list_mod_,"list_sex"=list_sex)))
  }
}

# Parallel processing
n_cluster<-min(floor(detectCores()*3/4),length(list_src_gamm))
#n_cluster<-min(floor(detectCores()*1/8),length(list_src_gamm))
clust<-makeCluster(n_cluster)
#print(paste("Calculating GAM in parallel,",as.character(n_cluster),"cores.",sep=" "))
clusterExport(clust,
              varlist=c("list_mod_","sort","gam","as.formula","summary.gam",
                        "as.numeric.factor"),
              envir=environment())
#list_dst_gamm<-parLapply(clust,list_src_gamm,gamm_core)
list_dst_gamm<-pblapply(list_src_gamm,gamm_core,cl=clust)
stopCluster(clust)

####

#print("Combining GAM results.")
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
df_out_gamm<-df_out_aic<-NULL
for (dst_gamm in list_dst_gamm){
  df_out_gamm<-rbind(df_out_gamm,dst_gamm$df_out_gamm_add)
  df_out_aic<-rbind(df_out_aic,dst_gamm$df_out_aic_add)
}
list_dst_gamm<-NULL
gc()

rownames(df_out_gamm)<-rownames(df_out_aic)<-NULL

####

df_gam<-add_mltcmp(data_gamm$df_out_gamm,df_roi,list_mod_,list_plot_,calc_seed_level=F)
write.csv(df_gam,
          file.path(paths_$output,"output","temp",
                    paste("atl-",atlas,"_ses-",label_waves,"_var-",idx_var_,"_gam.csv",sep=""))
          ,row.names = F)
write.csv(data_gamm$df_out_aic,
          file.path(paths_$output,"output","temp",
                    paste("atl-",atlas,"_ses-",label_waves,"_var-",idx_var_,"_gam_aic.csv",sep="")),
          row.names = F)

# Calculate group-wise average of FC
df_fc<-df_fc[df_fc$ses==wave_mri,]
df_fc<-left_join(df_fc,df_roi,by=c("from"="id"))
colnames(df_fc)[colnames(df_fc)=="group"]<-"from_grp"
df_fc<-left_join(df_fc,df_roi,by=c("to"="id"))
colnames(df_fc)[colnames(df_fc)=="group"]<-"to_grp"
list_id_subj<-sort(unique(df_fc$ID_pnTTC))
list_group<-unique(df_roi$group)
df_fc_grp<-data.frame()
for (id_subj in list_id_subj){
  df_fc_subset1<-df_fc[df_fc$ID_pnTTC==id_subj,]
  for (idx_grp1 in seq(length(list_group))){
    for (idx_grp2 in seq(idx_grp1,length(list_group))){
      label_grp1<-as.character(list_group[idx_grp1])
      label_grp2<-as.character(list_group[idx_grp2])
      df_fc_subset2<-df_fc_subset1[df_fc_subset1$from_grp==label_grp1 & df_fc_subset1$to_grp==label_grp2,]
      df_fc_grp<-rbind(df_fc_grp,data.frame("ID_pnTTC"=id_subj,"from"=label_grp1,"to"=label_grp2,
                                            "mean_z_r"=mean(df_fc_subset2$z_r),
                                            "mean_abs_z_r"=mean(abs(df_fc_subset2$z_r))))
    }
  }
}

# Join group-wise FC and clinical data
df_fc_grp$mean_z_r[which(is.nan(df_fc_grp$mean_z_r))]<-0
df_fc_grp$mean_abs_z_r[which(is.nan(df_fc_grp$mean_abs_z_r))]<-0
df_clin_wave<-df_clin[df_clin$wave==wave_clin,]
df_clin_wave$wave<-NULL
df_join<-inner_join(df_fc_grp,df_clin_wave,by='ID_pnTTC')
for (key in c('ID_pnTTC','sex')){
  if (key %in% colnames(df_join)){
    df_join[,key]<-as.factor(df_join[,key])
  }
}

# Calculate and save group-wise GAMM of FC
df_grp<-data.frame("id"=list_group,"label"=str_to_title(gsub("_"," ",list_group)))
df_join_sign<-df_join_abs<-df_join
colnames(df_join_sign)[colnames(df_join_sign)=="mean_z_r"]<-"value"
df_join_sign$mean_abs_z_r<-NULL
data_gamm_grp_sign<-iterate_gamm(df_join_sign,df_grp,list_mod_,calc_parallel=T,calc_identical=T)
df_gam_grp_sign<-add_mltcmp(data_gamm_grp_sign$df_out_gamm,df_grp,list_mod_,list_plot_,calc_seed_level=F)
write.csv(df_gam_grp_sign,
          file.path(paths_$output,"output","temp",
                    paste("atl-",atlas,"_ses-",label_waves,"_var-",idx_var_,"_gam_grp_sign.csv",sep="")),
          row.names = F)
colnames(df_join_abs)[colnames(df_join_abs)=="mean_abs_z_r"]<-"value"
df_join_abs$mean_z_r_<-NULL
data_gamm_grp_abs<-iterate_gamm(df_join_abs,df_grp,list_mod_,calc_parallel=T,calc_identical=T)
df_gam_grp_abs<-add_mltcmp(data_gamm_grp_abs$df_out_gamm,df_grp,list_mod_,list_plot_,calc_seed_level=F)
write.csv(df_gam_grp_sign,
          file.path(paths_$output,"output","temp",
                    paste("atl-",atlas,"_ses-",label_waves,"_var-",idx_var_,"_gam_grp_abs.csv",sep="")),
          row.names = F)

# Graphical output of ROI-wise and group-wise GAMM of FC
plot_gam_fc(paths_,df_gam,df_gam_grp_sign,df_gam_grp_abs,atlas,
            list_mod=list_mod_,list_plot=list_plot_,
            list_type_p=list_type_p_,thr_p=thr_p_,waves=waves_,idx_var=idx_var_)


####

nullobj<-NULL
gc()

########

gam_fc_list_waves<-gam_fc_list_waves[1]
gam_fc_list_tanner<-gam_fc_list_tanner[1]
gam_fc_list_hormone<-gam_fc_list_hormone[1]
list_atlas<-list_atlas[1]