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


print("Starting gam_fc_cs_multi()")
nullobj<-func_createdirs(paths_,str_proc="gam_fc_cs_multi()",copy_log=T)


label_waves<-names(list_waves_)[1]


wave_clin<-list_waves_[[label_waves]]$wave_clin
wave_mri<-list_waves_[[label_waves]]$wave_mri
waves<-list_waves_[label_waves]
print(paste("Clinical wave: ", wave_clin,", MRI wave: ",wave_mri,sep=""))

# Prepare subject subsetting condition (MRI QC criteria) according to specified waves
subset_subj_temp<-subset_subj_[[as.character(wave_mri)]]
subset_subj_temp<-list(subset_subj_temp)
names(subset_subj_temp)<-wave_clin

#1 Tanner stage
idx_tanner<-names(list_tanner_)[1]


print(paste("Tanner type: ",list_tanner_[[idx_tanner]][["label"]],sep=""))
list_covar<-list_covar_tanner_
list_covar[["tanner"]]<-list_tanner_[[idx_tanner]]
#suffix<-paste("ses-",label_waves,"_var-",idx_tanner,sep="")


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


print("Starting gam_fc_cs().")
#nullobj<-func_createdirs(paths_,str_proc="gam_fc_cs()",copy_log=T)
dict_roi <- func_dict_roi(paths_)

# Load and subset clinical data according to specified subsetting condition and covariate availability
print('Loading clinical data.')
label_waves<-names(waves_)
wave_clin<-waves_[[1]]$wave_clin
wave_mri<-waves_[[1]]$wave_mri
paste("ses-",label_waves,"_var-",idx_tanner,sep="")
data_clin<-func_clinical_data_long(paths_,list_wave=wave_clin,subset_subj_,
                                   list_covar=list_covar_,rem_na_clin=T,
                                   prefix=paste("ses-",label_waves,"_var-",idx_var_,sep=""),
                                   print_terminal=F)
df_clin<-data_clin$df_clin


atlas<-list_atlas_[1]


# Load ROI-wise FC data
print(paste('Loading FC data, atlas:',atlas,sep=' '))
#df_fc<-read.csv(file.path(paths_$input,'output',paste('atl-',atlas,'_fc.csv',sep='')))
df_fc<-as.data.frame(fread(file.path(paths_$input,'output',paste('atl-',atlas,'_fc.csv',sep=''))))
df_join<-join_fc_clin(df_fc,df_clin,wave_clin,wave_mri)
write.csv(df_join,file.path(paths_$output,"output","temp",
                            paste("atl-",atlas,"_ses-",label_waves,"_var-",idx_var_,"_src.csv",sep="")),
          row.names=F)

# Calculate and save ROI-wise GAMM of FC
print(paste('Calculating GAM, atlas: ',atlas,sep=''))
list_roi<-sort(unique(c(as.character(df_join$from),as.character(df_join$to))))
df_roi<-dict_roi[is.element(dict_roi$id,list_roi),c("id","label",key_group_)]
df_roi$id<-as.character(df_roi$id)
colnames(df_roi)[colnames(df_roi)==key_group_]<-"group"
data_gamm<-iterate_gamm(df_join,df_roi,list_mod_,calc_parallel=T)
df_gam<-add_mltcmp(data_gamm$df_out_gamm,df_roi,list_mod_,list_plot_,calc_seed_level=F)
write.csv(df_gam,
          file.path(paths_$output,"output","temp",
                    paste("atl-",atlas,"_ses-",label_waves,"_var-",idx_var_,"_gam.csv",sep="")),
          row.names = F)
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


list_mod=list_mod_
list_plot=list_plot_
list_type_p=list_type_p_
thr_p=thr_p_
#thr_p<-1
waves=waves_
idx_var=idx_var_


label_waves<-names(waves)
wave_clin<-waves[[1]]$wave_clin
wave_mri<-waves[[1]]$wave_mri

dict_roi<-func_dict_roi(paths_)
dict_roi<-dict_roi[dict_roi$atlas==atlas,c("id","label","group_3")]

# Create list of ROIs with blanks between groups
list_group<-unique(dict_roi$group_3)
list_roi_axis<-NULL
title_axis<-"Groups: "
for (group in list_group){
  list_roi_axis<-c(list_roi_axis,as.character(dict_roi[dict_roi$group_3==group,"label"]),"")
  title_axis<-paste(title_axis,group,", ",sep="")
}
list_roi_axis<-list_roi_axis[1:length(list_roi_axis)-1]
title_axis<-substr(title_axis,1,nchar(title_axis)-2)
list_label_group<-str_to_title(gsub("_"," ",as.character(list_group)))



idx_mod<-names(list_mod)[1]
idx_plot<-names(list_plot)[1]


idx_sex<-1


var_exp<-list_plot[[idx_plot]][["var_exp"]]

# Subset GAMM result dataframe for plotting
if (idx_sex==1){
  label_sex<-"m"
}else{
  label_sex<-"f"
}
df_gam_subset<-df_gam[df_gam$model==idx_mod & df_gam$term==var_exp & df_gam$sex==idx_sex,]
df_gam_grp_sign_subset<-df_gam_grp_sign[df_gam_grp_sign$model==idx_mod
                                        & df_gam_grp_sign$term==var_exp
                                        & df_gam_grp_sign$sex==idx_sex,]
df_gam_grp_abs_subset<-df_gam_grp_abs[df_gam_grp_abs$model==idx_mod
                                      & df_gam_grp_abs$term==var_exp
                                      & df_gam_grp_abs$sex==idx_sex,]

  #print(paste("GAMM output, atlas: ",atlas,", model: ",idx_mod,", plot: ",var_exp,", sex: ",label_sex,sep=""))
  # Convert GAMM rseult into igraph object
  if (!is.na(df_gam_subset[1,"estimate"])){
    df_gam_subset<-rename(df_gam_subset,c("estimate"="weight"))
    df_gam_grp_sign_subset<-rename(df_gam_grp_sign_subset,c("estimate"="weight"))
    df_gam_grp_abs_subset<-rename(df_gam_grp_abs_subset,c("estimate"="weight"))
    label_legend<-"beta"
  }else{
    df_gam_subset<-rename(df_gam_subset,c("F"="weight"))
    df_gam_grp_sign_subset<-rename(df_gam_grp_sign_subset,c("F"="weight"))
    df_gam_grp_abs_subset<-rename(df_gam_grp_abs_subset,c("F"="weight"))
    label_legend<-"F"
  }
  
  # Plot and save circular graph



type_p<-list_type_p[1]

list_subplot<-list()

# ROI-ROI heatmap
df_edge<-df_gam_subset
limits<-max(max(df_edge$weight),-min(df_edge$weight))
limits<-c(-limits,limits)
#df_edge<-df_edge[which(df_edge[,type_p]<thr_p),]
df_edge[which(df_edge[,type_p]>thr_p),"weight"]<-NA
df_edge<-df_edge[,c("label_from","label_to","weight")]
colnames(df_edge)<-c("row","column","r")
df_edge_inv<-data.frame(row=df_edge$column, column=df_edge$row,r=df_edge$r)
df_edge_identical<-data.frame(row=dict_roi$label,column=dict_roi$label,r=NA)
df_edge<-rbind(df_edge,df_edge_inv,df_edge_identical)

plot<-(ggplot(df_edge, aes(column, row))
       + geom_tile(aes(fill = r))
       + scale_fill_gradientn(colors = matlab.like2(100),name=label_legend,limits=limits)
       + scale_y_discrete(limits = rev(list_roi_axis))
       + scale_x_discrete(limits = list_roi_axis, position="top")
       + theme_linedraw()
       + theme(
         axis.text.x = element_text(size=1.5,angle = 90,vjust=0,hjust=0),
         axis.text.y = element_text(size=1.5),
         panel.grid.major=element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         plot.title = element_text(hjust = 0.5),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         axis.ticks=element_blank()
       )
)
list_subplot<-c(list_subplot,list(plot))

# group-group heatmap
df_is_sign<-T
for (df_gam_grp_subset in list(df_gam_grp_sign_subset,df_gam_grp_abs_subset)){
  df_edge<-df_gam_grp_subset
  limits<-max(max(df_edge$weight),-min(df_edge$weight))
  limits<-c(-limits,limits)
  df_edge[which(df_gam_grp_subset[,type_p]>thr_p),"weight"]<-NA
  df_edge$from<-as.character(df_edge$from)
  df_edge$to<-as.character(df_edge$to)
  df_edge<-df_edge[,c("label_from","label_to","weight")]
  colnames(df_edge)<-c("row","column","r")
  df_edge_inv<-df_edge[df_edge$row!=df_edge$column,]
  df_edge_inv<-data.frame(row=df_edge_inv$column, column=df_edge_inv$row,r=df_edge_inv$r)
  df_edge<-rbind(df_edge,df_edge_inv)
  
  plot<-(ggplot(df_edge, aes(column, row))
         + geom_tile(aes(fill = r))
         + scale_y_discrete(limits = rev(list_label_group))
         + scale_x_discrete(limits = list_label_group, position="top")
         + theme_linedraw()
         + theme(
           axis.text.x = element_text(size=8.5,angle = 90,vjust=0,hjust=0),
           axis.text.y = element_text(size=8.5),
           panel.grid.major=element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           plot.title = element_text(hjust = 0.5),
           axis.title.x=element_blank(),
           axis.title.y=element_blank(),
           axis.ticks=element_blank()
         )
  )
  
  if (df_is_sign){
    plot<-(plot
           + scale_fill_gradientn(colors = matlab.like2(100),
                                  name=paste("mean(",label_legend,")",sep=""),limits=limits))
  }else{
    plot<-(plot
           + scale_fill_gradientn(colors = matlab.like2(100),
                                  name=paste("mean(abs(",label_legend,"))",sep=""),limits=limits))
  }
  df_is_sign<-F
  list_subplot<-c(list_subplot,list(plot))
}
arranged_plot<-ggarrange(list_subplot[[1]],
                         ggarrange(list_subplot[[2]],list_subplot[[3]],
                                   ncol=2,
                                   labels=c("Group(signed)","Group(absolute)"),
                                   label.x=0,
                                   font.label = list(size = 10,face="plain")),
                         nrow=2,heights=c(2,1),
                         labels="ROI",
                         font.label = list(size = 10,face="plain"))

arranged_plot<-annotate_figure(arranged_plot,
                               top = text_grob(paste("GLM/GAM sex: ",label_sex,", measure: ",
                                                     idx_var,", model: ",idx_mod,
                                                     ", expvar: ",var_exp,", threshold: ",type_p,sep=""),
                                               color = "black", size = 14))

ggsave(paste("atl-",atlas,"_mod-",idx_mod,"_plt-",var_exp,
             "_sex-",label_sex,"_pval-",type_p,
             "_ses-",label_waves,"_var-",idx_var_,"_gam.png",sep=""),
       plot=arranged_plot,path=file.path(paths_$output,"output","plot"),height=13,width=10,dpi=600)




##############

gam_fc_list_waves<-gam_fc_list_waves[1]
gam_fc_list_tanner<-gam_fc_list_tanner[1]
gam_fc_list_hormone<-gam_fc_list_hormone[1]
x<-gam_fc_cs_multi()

