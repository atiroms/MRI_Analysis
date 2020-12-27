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
df_fc_grp$ID_pnTTC<-as.numeric.factor(df_fc_grp$ID_pnTTC)
df_fc_grp$wave<-as.character(as.numeric.factor(df_fc_grp$wave))
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

df_gam=df_plot
df_gam_grp_sign=df_plot_grp
df_gam_grp_abs=NULL
list_type_p=list_type_p_
thr_p=thr_p_
waves=NULL

####

dict_roi<-func_dict_roi(paths_)
dict_roi<-dict_roi[dict_roi$atlas==atlas,c("id","label","group_3")]

# Create list of ROIs with blanks between groups
list_group<-unique(as.character(dict_roi$group_3))
list_roi_axis<-NULL
title_axis<-"Groups: "
for (group in list_group){
  list_roi_axis<-c(list_roi_axis,as.character(dict_roi[dict_roi$group_3==group,"label"]),"")
  title_axis<-paste(title_axis,group,", ",sep="")
}
list_roi_axis<-list_roi_axis[1:length(list_roi_axis)-1]
title_axis<-substr(title_axis,1,nchar(title_axis)-2)
list_label_group<-str_to_title(gsub("_"," ",as.character(list_group)))

####

idx_mod<-names(list_mod)[2]
idx_plot<-names(list_plot)[2]
idx_sex<-1

####

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
if (nrow(df_gam_subset)>0){
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
  
  # Plot and save heatmap
  for (type_p in list_type_p){
    if(type_p %in% colnames(df_gam_subset)){
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
      for (df_gam_grp_subset in list(df_gam_grp_sign_subset,df_gam_grp_abs_subset)){
        if (is.null(df_gam_grp_subset)){
          list_subplot<-c(list_subplot,list(NULL))
        }else{
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
                 + scale_fill_gradientn(colors = matlab.like2(100),
                                        name=label_legend,limits=limits)
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
          list_subplot<-c(list_subplot,list(plot))
        }
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
                   "_ses-",names(waves),"_var-",idx_var,"_gam.png",sep=""),
             plot=arranged_plot,path=file.path(paths_$output,"output","plot"),height=13,width=10,dpi=600)
    }
  }
}
