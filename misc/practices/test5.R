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
df_pca_mri_grp<-read.csv(file.path(paths_$output,"output",
                                   paste("atl-",atlas,"_ses-m",wave_mri,"_fc_pca_var_grp.csv",sep="")))

label_sex<-names(list_sex_)[1]


# PCA
dim_ca<-max(list_dim_ca_)
df_pca_mri_subset<-df_pca_mri[df_pca_mri$sex==label_sex & df_pca_mri$dim==dim_ca,]
df_pca_mri_grp_subset<-df_pca_mri_grp[df_pca_mri_grp$sex==label_sex & df_pca_mri_grp$dim==dim_ca,]
df_pca_mri_subset$sex<-df_pca_mri_subset$dim<-df_pca_mri_grp_subset$sex<-df_pca_mri_grp_subset$dim<-NULL


df_comp_mri<-df_pca_mri_subset
df_comp_mri_grp<-df_pca_mri_grp_subset
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

list_plot<-list()
# ROI-ROI heatmap
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
       #+ ggtitle(paste("Method: ",method,", Atlas: ",atlas,", Wave: ",as.character(ses),
       #                ", Component: ",sprintf("%03d",idx_comp),"/",sprintf("%03d",dim_ca),
       #                ", Sex: ",label_sex,sep=""))
       #+ xlab(title_axis)
       + theme_light()
       + theme(axis.text.x = element_text(size=29/log(length(list_roi_axis),2),angle = 90,vjust=0,hjust=0),
               axis.text.y = element_text(size=29/log(length(list_roi_axis),2)),
               panel.grid.major=element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               #legend.title=element_blank(),
               plot.title = element_text(hjust = 0.5),
               #axis.title.x=element_text(size=5),
               axis.title.x=element_blank(),
               axis.title.y=element_blank(),
               axis.ticks=element_blank()
       )
)

list_plot<-c(list_plot,list(plot))

# group-group heatmap

#list_plot<-list()
for (abs_mean in c(F,T)){

   df_edge<-df_comp_mri_grp[df_comp_mri_grp$abs==abs_mean,c("from","to",sprintf("comp_%03d",idx_comp))]
   colnames(df_edge)<-c("row","column","r")
   limits<-max(max(df_edge$r),-min(df_edge$r))
   limits<-c(-limits,limits)
   df_edge_inv<-df_edge[df_edge$row!=df_edge$column,]
   df_edge_inv<-data.frame(row=df_edge_inv$column, column=df_edge_inv$row,r=df_edge_inv$r)
   df_edge<-rbind(df_edge,df_edge_inv)
   
   plot<-(ggplot(df_edge, aes(column, row))
          + geom_tile(aes(fill = r))
          + scale_fill_gradientn(colors = matlab.like2(100),name="mean z",limits=limits)
          + scale_y_discrete(limits = rev(list_group))
          + scale_x_discrete(limits = list_group, position="top")
          #+ ggtitle(paste("Method: ",method,", Atlas: ",atlas,", Wave: ",as.character(ses),
          #                ", Component: ",sprintf("%03d",idx_comp),"/",sprintf("%03d",dim_ca),
          #                ", Sex: ",label_sex,sep=""))
          #+ theme_light()
          + theme_linedraw()
          + theme(axis.text.x = element_text(size=29/log(length(list_group),2),angle = 90,vjust=0,hjust=0),
                  axis.text.y = element_text(size=29/log(length(list_group),2)),
                  panel.grid.major=element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  #legend.title=element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  #axis.title.x=element_text(size=5),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  axis.ticks=element_blank()
          )
   )
   
   list_plot<-c(list_plot,list(plot))
}


arranged_plot<-ggarrange(list_plot[[1]],
                         ggarrange(list_plot[[2]],list_plot[[3]],
                                   ncol=2,
                                   labels=c("Group(signed)","Group(absolute)"),
                                   label.x=-0.05,
                                   font.label = list(size = 10,face="plain"),
                                   common.legend = TRUE, legend = "right"),
                         nrow=2,heights=c(2,1),
                         labels="ROI",
                         font.label = list(size = 10,face="plain"))

arranged_plot<-annotate_figure(arranged_plot,
                               top = text_grob(paste("Method: ",method,", Atlas: ",atlas,", Wave: ",as.character(ses),
                                                     ", Component: ",sprintf("%03d",idx_comp),"/",sprintf("%03d",dim_ca),
                                                     ", Sex: ",label_sex,sep=""), color = "black", size = 14))

arranged_plot

ggsave("test.png",
       plot=arranged_plot,
       path=file.path(paths_$output,"output","plot"),height=13,width=10,dpi=600)

#arranged_plot<-ggarrange(list_plot[[1]],list_plot[[2]],
#                         ncol=2,labels=c("signed","absolute"),
#                         common.legend = TRUE, legend = "right")

#annotate_figure(figure,
#                top = text_grob("Visualizing Tooth Growth", color = "red", face = "bold", size = 14),
#                bottom = text_grob("Data source: \n ToothGrowth data set", color = "blue",
#                                   hjust = 1, x = 1, face = "italic", size = 10),
#                left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
#                right = text_grob(bquote("Superscript: ("*kg~NH[3]~ha^-1~yr^-1*")"), rot = 90),
#                fig.lab = "Figure 1", fig.lab.face = "bold"
#)

#ggarrange(sp,                                                 # First row with scatter plot
#          ggarrange(bxp, dp, ncol = 2, labels = c("B", "C")), # Second row with box and dot plots
#          nrow = 2, 
#          labels = "A"                                        # Labels of the scatter plot
#) 
