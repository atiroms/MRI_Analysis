paths_=paths
list_waves_=ca_fc_list_waves
subset_subj_=ca_fc_subset_subj
list_sex_=ca_fc_list_sex
list_atlas_=list_atlas
list_covar_tanner_=ca_fc_list_covar_tanner
list_tanner_=ca_fc_list_tanner
list_covar_hormone_=ca_fc_list_covar_hormone
list_hormone_=ca_fc_list_hormone
list_dim_ca_=list_dim_ca
ratio_vis_=ratio_vis


print("Starting ca_fc_cs_multi()")
nullobj<-func_createdirs(paths_,str_proc="ca_fc_cs_multi()",copy_log=T)
# Increase memory limit
memory.limit(1000000)
df_cor<-NULL

atlas<-"aal116"



waves<-names(list_waves_)[1]
wave_mri<-list_waves_[[waves]]$wave_mri


print(paste("MRI wave: ",wave_mri,", loading PCA/ICA results.",sep=""))
df_pca_mri<-read.csv(file.path(paths_$output,"output",
                               paste("atl-",atlas,"_ses-m",wave_mri,"_fc_pca_var.csv",sep="")))


label_sex<-names(list_sex_)[1]


# PCA
dim_ca<-max(list_dim_ca_)
df_pca_mri_subset<-df_pca_mri[df_pca_mri$sex==label_sex & df_pca_mri$dim==dim_ca,]
df_pca_mri_subset$sex<-df_pca_mri_subset$dim<-NULL


df_comp_mri<-df_pca_mri_subset
method<-"pca"
ses<-wave_mri


print(paste("Generationg heatmap plot of factors, Session: ",as.character(ses),
            ", Sex: ",label_sex,", Method: ",method,", Dim: ",as.character(dim_ca),sep="")) 
dict_roi<-func_dict_roi(paths_)
dict_roi<-dict_roi[dict_roi$atlas==atlas,c("id","label","group_3")]
#dict_roi<-dict_roi[order(dict_roi$group_3),]

list_group<-unique(dict_roi$group_3)
list_roi_axis<-NULL
title_axis<-"Groups: "
for (group in list_group){
   list_roi_axis<-c(list_roi_axis,dict_roi[dict_roi$group_3==group,"label"],"")
   title_axis<-paste(title_axis,group,", ",sep="")
}
list_roi_axis<-list_roi_axis[1:length(list_roi_axis)-1]
title_axis<-substr(title_axis,1,nchar(title_axis)-2)


df_comp_mri<-inner_join(df_comp_mri,dict_roi[,c("id","label")],by=c("from"="id"))
colnames(df_comp_mri)[colnames(df_comp_mri)=="label"]<-"from_label"
df_comp_mri<-inner_join(df_comp_mri,dict_roi[,c("id","label")],by=c("to"="id"))
colnames(df_comp_mri)[colnames(df_comp_mri)=="label"]<-"to_label"

idx_comp<-1


df_edge<-df_comp_mri[,c("from_label","to_label",sprintf("comp_%03d",idx_comp))]
colnames(df_edge)<-c("row","column","r")
limits<-max(max(df_edge$r),-min(df_edge$r))
limits<-c(-limits,limits)
df_edge_inv<-data.frame(row=df_edge$column, column=df_edge$row,r=df_edge$r)
df_edge_identical<-data.frame(row=dict_roi$label,column=dict_roi$label,r=NA)
df_edge<-rbind(df_edge,df_edge_inv,df_edge_identical)


plot<-(ggplot(df_edge, aes(column, row))
       + geom_tile(aes(fill = r))
       + scale_fill_gradientn(colors = matlab.like2(100),name="z",limits=limits)
       #       + scale_y_discrete(limits = rev(dict_roi$label))
       #       + scale_x_discrete(limits = dict_roi$label, position="top")
       + scale_y_discrete(limits = rev(list_roi_axis))
       + scale_x_discrete(limits = list_roi_axis, position="top")
       + ggtitle(paste("Method: ",method,", Atlas: ",atlas,", Wave: ",as.character(ses),
                       ", Component: ",sprintf("%03d",idx_comp),"/",sprintf("%03d",dim_ca),
                       ", Sex: ",label_sex,sep=""))
       + xlab(title_axis)
       + theme_light()
       + theme(axis.text.x = element_text(size=29/log(length(list_roi_axis),2),angle = 90,vjust=0,hjust=0),
               axis.text.y = element_text(size=29/log(length(list_roi_axis),2)),
               panel.grid.major=element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               #legend.title=element_blank(),
               plot.title = element_text(hjust = 0.5),
               axis.title.y=element_blank(),
               axis.ticks=element_blank()
               )
)

plot
