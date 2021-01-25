#**************************************************
# Description =====================================
#**************************************************

# R script for common visualization functionalities in MRI analysis


#**************************************************
# Libraries =======================================
#**************************************************
library(ggplot2)
library(ggraph)
library(igraph)
library(colorRamps)
library(purrr)
library(viridis)


#**************************************************
# Histogram of Permutaion =========================
#**************************************************
plot_permutation<-function(paths_,list_max,thr_size_perm,
                           atlas,model,plot,sex,title_plot,title_sex,color_plt){
  plt<-(ggplot(data.frame(max=list_max), aes(x=max))
        + geom_histogram(binwidth=5,fill=color_plt)
        + geom_vline(aes(xintercept=thr_size_perm),
                     color="grey", linetype="dashed", size=1)
        + ggtitle(paste("Atlas: ",atlas,", Model: ",model,", Plot: ",title_plot,
                        ", Sex: ",title_sex,sep=""))
        + xlab("Size")
        + ylab("Count")
        + theme_light()
        + theme(plot.title = element_text(hjust = 0.5))
  )
  ggsave(paste("atl-",atlas,"_mod-",model,"_plt-",plot,
               "_sex-",sex,"_perm.png",sep=""),
         plot=plt,path=file.path(paths_$output,"output","plot"),height=5,width=7,dpi=300)
}


#**************************************************
# Heatmap Plot of sex difference NBS ==============
#**************************************************
plot_sex_diff_fc<-function(paths_,df_edge,atlas,df_roi,df_grp,mod,plot,sex,
                           title_plot,title_sex,idx_net){
  # Create list of ROIs with blanks between groups
  list_roi_axis<-NULL
  for (group in df_grp$id){
    list_roi_axis<-c(list_roi_axis,as.character(df_roi[df_roi$group==group,"label"]),"")
  }
  list_roi_axis<-list_roi_axis[1:length(list_roi_axis)-1]
  title_axis<-paste("Groups: ",paste(df_grp$label,collapse=", "),sep="")
  
  if (!is.na(df_edge[1,"t"])){
    df_edge<-rename(df_edge,c("t"="r"),warn_missing=F)
    label_legend<-"t"
  }else{
    df_edge<-rename(df_edge,c("F"="r"),warn_missing=F)
    label_legend<-"F"
  }
  df_edge<-df_edge[,c("from","to","r")]
  df_edge_inv<-df_edge
  colnames(df_edge_inv)<-c("to","from","r")
  df_edge<-rbind(df_edge,df_edge_inv)
  rownames(df_roi)<-NULL
  df_edge_full<-NULL
  for (idx1 in seq(nrow(df_roi))){
    for (idx2 in seq(nrow(df_roi))){
      if (idx1!=idx2){
        df_edge_full<-rbind(df_edge_full,data.frame(from=df_roi[idx1,"id"],to=df_roi[idx2,"id"],
                                                  label_from=df_roi[idx1,"label"],label_to=df_roi[idx2,"label"]))
      }
    }
  }
  df_edge<-left_join(df_edge_full,df_edge,by=c("from","to"))
  df_edge<-df_edge[,c("label_from","label_to","r")]
  colnames(df_edge)<-c("row","column","r")
  limits<-max(max(df_edge$r,na.rm=T),-min(df_edge$r,na.rm=T))
  limits<-c(-limits,limits)

  plt<-(ggplot(df_edge, aes(column, row))
         + geom_tile(aes(fill = r))
         + scale_fill_gradientn(colors = matlab.like2(100),name=label_legend,limits=limits)
         + scale_y_discrete(limits = rev(list_roi_axis))
         + scale_x_discrete(limits = list_roi_axis, position="top")
         + ggtitle(paste("Atlas: ",atlas,", Model: ",mod,", Plot: ",title_plot,
                         ", Contr: ",title_sex,", #",sprintf("%02d",idx_net),sep=""))
         + xlab(title_axis)
         + theme_linedraw()
         + theme(
           axis.text.x = element_text(size=1.5,angle = 90,vjust=0,hjust=0),
           axis.text.y = element_text(size=1.5),
           panel.grid.major=element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           plot.title = element_text(hjust = 0.5),
           #axis.title.x=element_blank(),
           axis.title.y=element_blank(),
           axis.ticks=element_blank()
         )
  )
  
  ggsave(paste("atl-",atlas,"_mod-",mod,"_plt-",plot,
               "_cntr-",sex,"_idx-",sprintf("%02d",idx_net),"_net.png",sep=""),
         plot=plt,path=file.path(paths_$output,"output","plot"),height=10,width=10,dpi=600)
  
}

#**************************************************
# Heatmap Plot of GAM of FC =======================
#**************************************************

#plot_gam_fc<-function(paths_,df_comp_mri,df_comp_mri_grp,atlas,dim_ca,method,label_sex,ses){
plot_gam_fc<-function(paths_,df_gam,df_gam_grp_sign,df_gam_grp_abs,atlas,
                      list_mod,list_plot,list_type_p,thr_p,waves,idx_var){

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
            df_gam_subset<-rename(df_gam_subset,c("estimate"="weight"),warn_missing=F)
            df_gam_grp_sign_subset<-rename(df_gam_grp_sign_subset,c("estimate"="weight"),warn_missing=F)
            df_gam_grp_abs_subset<-rename(df_gam_grp_abs_subset,c("estimate"="weight"),warn_missing=F)
            label_legend<-"beta"
          }else{
            df_gam_subset<-rename(df_gam_subset,c("F"="weight"),warn_missing=F)
            df_gam_grp_sign_subset<-rename(df_gam_grp_sign_subset,c("F"="weight"),warn_missing=F)
            df_gam_grp_abs_subset<-rename(df_gam_grp_abs_subset,c("F"="weight"),warn_missing=F)
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
      }
    }
  }
}


#**************************************************
# Heatamap plot of PCA/ICA factors of FC ==========
#**************************************************
plot_ca_fc_heatmap_core<-function(data_plot){
  file<-data_plot$file
  plot<-data_plot$plot
  path<-data_plot$path
  ggsave(file,plot=plot,path=path,height=13,width=10,dpi=600)
  return(T)
}

plot_ca_fc_heatmap<-function(paths_,df_comp_mri,df_comp_mri_grp,atlas,dim_ca,method,label_sex,ses){
  print(paste("Generationg heatmap plot of factors, Session: ",as.character(ses),
              ", Sex: ",label_sex,", Method: ",method,", Dim: ",as.character(dim_ca),sep="")) 
  dict_roi<-func_dict_roi(paths_)
  dict_roi<-dict_roi[dict_roi$atlas==atlas,c("id","label","group_3")]
  dict_roi$label<-as.character(dict_roi$label)
  #dict_roi<-dict_roi[order(dict_roi$group_3),]
  
  # Create list of ROIs with blanks between groups
  list_group<-unique(dict_roi$group_3)
  list_roi_axis<-NULL
  title_axis<-"Groups: "
  for (group in list_group){
    list_roi_axis<-c(list_roi_axis,dict_roi[dict_roi$group_3==group,"label"],"")
    title_axis<-paste(title_axis,group,", ",sep="")
  }
  list_roi_axis<-list_roi_axis[1:length(list_roi_axis)-1]
  title_axis<-substr(title_axis,1,nchar(title_axis)-2)
  
  # Convert ROI ID to label
  df_comp_mri<-inner_join(df_comp_mri,dict_roi[,c("id","label")],by=c("from"="id"))
  colnames(df_comp_mri)[colnames(df_comp_mri)=="label"]<-"from_label"
  df_comp_mri<-inner_join(df_comp_mri,dict_roi[,c("id","label")],by=c("to"="id"))
  colnames(df_comp_mri)[colnames(df_comp_mri)=="label"]<-"to_label"
  
  list_plot<-list()
  for (idx_comp in 1:dim_ca){
    list_subplot<-list()
    
    # ROI-ROI heatmap
    df_edge<-df_comp_mri[,c("from_label","to_label",sprintf("comp_%03d",idx_comp))]
    colnames(df_edge)<-c("row","column","r")
    limits<-max(max(df_edge$r),-min(df_edge$r))
    limits<-c(-limits,limits)
    df_edge_inv<-data.frame(row=df_edge$column, column=df_edge$row,r=df_edge$r)
    df_edge_identical<-data.frame(row=dict_roi$label,column=dict_roi$label,r=NA)
    df_edge<-rbind(df_edge,df_edge_inv,df_edge_identical)
    df_edge$row<-as.character(df_edge$row)
    df_edge$column<-as.character(df_edge$column)
    
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
           + theme_linedraw()
           + theme(
                   #axis.text.x = element_text(size=29/log(length(list_roi_axis),2),angle = 90,vjust=0,hjust=0),
                   #axis.text.y = element_text(size=29/log(length(list_roi_axis),2)),
                   axis.text.x = element_text(size=1.5,angle = 90,vjust=0,hjust=0),
                   axis.text.y = element_text(size=1.5),
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
    list_subplot<-c(list_subplot,list(plot))
    
    # group-group heatmap
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
             #+ scale_fill_gradientn(colors = matlab.like2(100),name="mean z",limits=limits)
             + scale_y_discrete(limits = rev(list_group))
             + scale_x_discrete(limits = list_group, position="top")
             #+ ggtitle(paste("Method: ",method,", Atlas: ",atlas,", Wave: ",as.character(ses),
             #                ", Component: ",sprintf("%03d",idx_comp),"/",sprintf("%03d",dim_ca),
             #                ", Sex: ",label_sex,sep=""))
             #+ theme_light()
             + theme_linedraw()
             + theme(
                     #axis.text.x = element_text(size=29/log(length(list_group),2),angle = 90,vjust=0,hjust=0),
                     #axis.text.y = element_text(size=29/log(length(list_group),2)),
                     axis.text.x = element_text(size=8.5,angle = 90,vjust=0,hjust=0),
                     axis.text.y = element_text(size=8.5),
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
      if (abs_mean){
        plot<-(plot
               + scale_fill_gradientn(colors=viridis(100),name="mean(abs(z))"))
      }else{
        plot<-(plot
               + scale_fill_gradientn(colors = matlab.like2(100),name="mean(z)",limits=limits))
      }
      list_subplot<-c(list_subplot,list(plot))
    }
    arranged_plot<-ggarrange(list_subplot[[1]],
                             ggarrange(list_subplot[[2]],list_subplot[[3]],
                                       ncol=2,
                                       labels=c("Group(signed)","Group(absolute)"),
                                       label.x=-0.05,
                                       font.label = list(size = 10,face="plain")),
                             nrow=2,heights=c(2,1),
                             labels="ROI",
                             font.label = list(size = 10,face="plain"))
    
    arranged_plot<-annotate_figure(arranged_plot,
                                   top = text_grob(paste("Method: ",method,", Atlas: ",atlas,", Wave: ",as.character(ses),
                                                         ", Component: ",sprintf("%03d",idx_comp),"/",sprintf("%03d",dim_ca),
                                                         ", Sex: ",label_sex,sep=""), color = "black", size = 14))
    
    list_plot<-c(list_plot,list(list("file"=paste("atl-",atlas,"_method-",method,"_ses-",as.character(ses),
                                                  "_sex-",label_sex,"_dim-",sprintf("%03d",dim_ca),
                                                  "_comp-",sprintf("%03d",idx_comp),"_fc_ca.png",sep=""),
                                     "path"=file.path(paths_$output,"output","plot"),
                                     "plot"=arranged_plot)))
  }
  n_cluster<-floor(detectCores()*3/4)
  clust<-makeCluster(n_cluster)
  clusterExport(clust,
                varlist=c("ggsave"),
                envir=environment())
  nullobj<-pblapply(list_plot,plot_ca_fc_heatmap_core,cl=clust)
  stopCluster(clust)
}


#**************************************************
# Circular plot of PCA/ICA factors of FC ==========
#**************************************************

plot_ca_fc_circular_core<-function(data_plot){
  file<-data_plot$file
  plot<-data_plot$plot
  path<-data_plot$path
  ggsave(file,plot=plot,path=path,height=10,width=10,dpi=600)
  return(T)
}

plot_ca_fc_circular<-function(paths_,df_comp_mri,atlas,dim_ca,ratio_vis,method,label_sex,ses){
  print(paste("Generationg circular plot of factors, Session: ",as.character(ses),
              ", Sex: ",label_sex,", Method: ",method,", Dim: ",as.character(dim_ca),sep=""))
  df_node<-func_dict_roi(paths=paths)
  df_node<-df_node[df_node$atlas==atlas,c("id","label")]
  n_vis<-floor(nrow(df_comp_mri)*ratio_vis)
  list_plot<-list()
  for (idx_comp in 1:dim_ca){
    df_edge<-df_comp_mri[,c("from","to",sprintf("comp_%03d",idx_comp),sprintf("rank_%03d",idx_comp))]
    colnames(df_edge)<-c("from","to","weight","rank")
    
    igraph_ca<- graph_from_data_frame(d = df_edge, vertices = df_node, directed = F)
    plot<-plot_circular(igraph_in=igraph_ca,type_p="rank",thr_p=n_vis,limit_color=NULL)
    plot<-plot +
      ggtitle(paste("Method: ",method,", Atlas: ",atlas,", Wave: ",as.character(ses),
                    ", Component: ",sprintf("%03d",idx_comp),"/",sprintf("%03d",dim_ca),
                    ", Sex: ",label_sex,sep="")) +
      theme(plot.title = element_text(hjust = 0.5))
    
    list_plot<-c(list_plot,list(list("file"=paste("atl-",atlas,"_method-",method,"_ses-",as.character(ses),
                                                  "_sex-",label_sex,"_dim-",sprintf("%03d",dim_ca),
                                                  "_comp-",sprintf("%03d",idx_comp),"_fc_ca.png",sep=""),
                                     "path"=file.path(paths_$output,"output","plot"),
                                     "plot"=plot)))
    
    #ggsave(paste("atl-",atlas,"_method-",method,"_sex-",label_sex,"_dim-",sprintf("%03d",dim_ca),"_comp-",sprintf("%03d",idx_comp),
    #             "_fc_ca.png",sep=""),
    #       plot=plot,path=file.path(paths_$output,"output","plot"),height=10,width=10,dpi=600)
  }
  
  n_cluster<-floor(detectCores()*3/4)
  clust<-makeCluster(n_cluster)
  clusterExport(clust,
                varlist=c("ggsave","guide_edge_colourbar"),
                envir=environment())
  nullobj<-pblapply(list_plot,plot_ca_fc_circular_core,cl=clust)
  stopCluster(clust)
}


#**************************************************
# Circular graph ==================================
#**************************************************
plot_circular<-function(igraph_in,type_p,thr_p,limit_color=NULL){
  
  # Subset edges according to p value criteria
  df_edge<-get.data.frame(igraph_in,what="edges")
  df_edge<-df_edge[which(df_edge[,type_p]<thr_p),]
  
  # Change order of nodes for circular plot aesthetics
  df_node<-get.data.frame(igraph_in,what="vertices")
  idx_node_r<-grep("^R ",df_node$label)
  idx_node_l<-rev(grep("^L ",df_node$label))
  if (length(idx_node_r)+length(idx_node_l)>0){
    df_node<-rbind(df_node[idx_node_r,],
                   df_node[c(-idx_node_r,-idx_node_l),],
                   df_node[idx_node_l,])
  }
  
  # Add circular plot specs 
  df_node$angle <- 90 - 360 * ((1:nrow(df_node))-0.5) / nrow(df_node)
  df_node$hjust<-ifelse(df_node$angle < -90, 1, 0)
  df_node$angle<-ifelse(df_node$angle < -90, df_node$angle+180, df_node$angle)
  
  # Convert edge/node dataframes into igraph object again
  igraph_plot<-graph_from_data_frame(d = df_edge, vertices = df_node, directed = F)
  
  # Calculate color limit if not specified
  if(is.null(limit_color)){
    if (nrow(df_edge)>0){
      limit_color <- max(abs(max(df_edge$weight)),abs(min(df_edge$weight)))
      limit_color <- c(-limit_color,limit_color)
    }
  }
  
  plot<-ggraph(igraph_plot, layout = "linear",circular = T) +
    geom_node_text(aes(x = x*1.03, y=y*1.03,
                       label=label, angle = angle, hjust=hjust,vjust=0.2),
                   size=min(5,10/log(nrow(df_node))), alpha=1) +
    #               size=min(5,20/log(nrow(df_node))), alpha=1) +
    geom_node_point(aes(x=x, y=y),size=1, alpha=1,colour="grey50") +
    scale_edge_color_gradientn(colors=matlab.like2(100),limits=limit_color,na.value="grey50")+
    expand_limits(x = c(-2, 2), y = c(-2, 2))+
    #expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))+
    #ggtitle(input_title) +
    theme_void() +
    #theme(plot.title = element_text(hjust = 0.5),legend.justification=c(1,1), legend.position=c(1,1))
    theme(legend.justification=c(1,1), legend.position=c(1,1))
  if (nrow(df_edge)>0){
    plot<-plot+
      geom_edge_arc(aes(color=weight),width=1,alpha=0.5)
  }
  return(plot)
}


#**************************************************
# Plot PCA/ICA result =============================
#**************************************************
#OBSOLETE
plot_ca_old<-function(df_src,list_name_covar,n_dim){
  df_plot<-df_src
  df_plot$ses<-as.factor(df_plot$ses)
  df_plot$ID_pnTTC<-as.factor(df_plot$ID_pnTTC)
  list_plot<-list()
  for (i_dim in 1:(n_dim-1)){
    list_plot_dim<-list()
    for (name_covar in list_name_covar){
      df_plot_subset<-df_plot[,c("ses","ID_pnTTC",name_covar,sprintf("comp_%03d",i_dim),sprintf("comp_%03d",i_dim+1))]
      colnames(df_plot_subset)<-c("ses","ID_pnTTC","color","x","y")
      plot<-(ggplot(df_plot_subset)
             + aes(x=x,y=y,label=ID_pnTTC,shape=ses)
             + scale_shape_manual(name=NULL,labels=c("1st wave","2nd wave"),values=c(3,4))
             + geom_point(size=2,aes(color=color))
             + scale_color_gradientn(colors = matlab.like2(100),name=name_covar)
             #+ geom_text_repel(size=2)
             + geom_path(aes(group=ID_pnTTC),size=0.5,alpha=0.5)
             #+ ggtitle("PCA of FC")
             + xlab(sprintf("Component %03d",i_dim))
             + ylab(sprintf("Component %03d",i_dim+1))
             + theme_light()
             + theme(plot.title = element_text(hjust = 0.5))
      )
      list_plot_dim_covar<-list(plot)
      names(list_plot_dim_covar)<-name_covar
      list_plot_dim<-c(list_plot_dim,list_plot_dim_covar)
      #ggsave(paste("atl-",atlas,"_dim-",sprintf("%02d",i_dim),"-",sprintf("%02d",i_dim+1),"_cov-",name_covar,"_pca_fc.eps",sep=""),plot=plot,device=cairo_ps,
      #       path=file.path(paths_$output,"output"),dpi=300,height=10,width=10,limitsize=F)
    }
    list_plot_dim<-list(list_plot_dim)
    names(list_plot_dim)<-as.character(i_dim)
    list_plot<-c(list_plot,list_plot_dim)
  }
  return(list_plot)
}


#**************************************************
# GAMM plot =======================================
#**************************************************
# modified from voxel/plotGAM

plot_gamm<-function(plot_in,mod_gamm,df_in,spec_graph){
  
  if (is.null(plot_in)){
    plot<-ggplot()
  }else{
    plot<-plot_in
  }
  
  key_df_src<-c("value",names(mod_gamm$var.summary))
  if (!("ID_pnTTC" %in% key_df_src)){
    key_df_src<-c(key_df_src, "ID_pnTTC")
  }
  if (!("sex" %in% key_df_src)){
    key_df_src<-c(key_df_src, "sex")
  }
  df_src <- df_in[key_df_src]
  
  # add prediction line + ribbon to plot
  if (!is.null(spec_graph[["smooth"]])){
    for (name_smooth in names(spec_graph[["smooth"]])){
      spec_smooth<-spec_graph[["smooth"]][[name_smooth]]
      df_smooth <- data.frame(x = seq(min(df_src[spec_graph[["x_axis"]]]),
                                      max(df_src[spec_graph[["x_axis"]]]),
                                      length.out=200))
      names(df_smooth) <- spec_graph[["x_axis"]]
      for (i in names(df_src)) {
        if (i != spec_graph[["x_axis"]]) {
          if (any(class(df_src[i][,1])[1] == c("numeric", "integer","boolean"))) {
            df_smooth[, dim(df_smooth)[2] + 1] <- mean(df_src[i][,1])
            names(df_smooth)[dim(df_smooth)[2]] <- i
          }
          else if (any(class(df_src[i][,1])[1] == c("character", "factor","ordered"))) {
            df_smooth[, dim(df_smooth)[2] + 1] <- df_src[i][1,1]
            names(df_smooth)[dim(df_smooth)[2]] <- i
          }
        }
      }
      flag_plot<-T
      if (!is.null(spec_smooth[["fix"]])){
        for (var in names(spec_smooth[["fix"]])){
          df_smooth[[var]]<-spec_smooth[["fix"]][[var]]
          if (!(spec_smooth[["fix"]][[var]] %in% df_src[[var]])){
            flag_plot<-F
          }
        }
      }
      if (flag_plot){
        df_smooth <- cbind(df_smooth,
                           as.data.frame(predict.gam(mod_gamm, df_smooth, exclude="s(ID_pnTTC)",se.fit = TRUE)))
        plot <- (plot
                 + geom_line(aes(x=!!df_smooth[,1],y=!!df_smooth[["fit"]]),
                             color=spec_smooth[["color"]],size=0.5,alpha=spec_smooth[["alpha"]]))
        if (spec_smooth[["ribbon"]]){
          plot <- (plot
                   + geom_ribbon(aes(x=!!df_smooth[,1],
                                     ymax = !!df_smooth[["fit"]]+1.96*!!df_smooth[["se.fit"]],
                                     ymin = !!df_smooth[["fit"]]-1.96*!!df_smooth[["se.fit"]],
                                     linetype=NA),
                                 fill=spec_smooth[["color"]],alpha = 0.2*spec_smooth[["alpha"]]))
        }
      }
    }
  }
  
  # add point + path to plot
  if (!is.null(spec_graph[["point"]])){
    df_point<-data.frame()
    for (name_point in names(spec_graph[["point"]])){
      spec_point<-spec_graph[["point"]][[name_point]]
      df_point_series<-df_src
      for (var_subset in names(spec_point[["subset"]])){
        df_point_series<-df_point_series[df_point_series[[var_subset]]==spec_point[["subset"]][[var_subset]],]
      }
      if (dim(df_point_series)[1]>0){
        df_point_series$name_series<-name_point
        df_point_series$color<-spec_point[["color"]]
        df_point_series$alpha<-spec_point[["alpha"]]
        df_point<-rbind(df_point,df_point_series)
      }
    }
    plot <- (plot
             + geom_point(aes(x=df_point[[spec_graph[["x_axis"]]]],
                              y=df_point[,1]),
                          color=df_point[["color"]],shape=1,size=2,alpha=0.6*df_point[["alpha"]]))
    if (any(duplicated(df_point[["ID_pnTTC"]]))){
      plot <- (plot
               + geom_path(aes(x=df_point[[spec_graph[["x_axis"]]]],
                               y=df_point[,1],
                               group=df_point[["ID_pnTTC"]]),
                           color=df_point[["color"]],size=0.3,alpha=0.3*df_point[["alpha"]],linetype="dashed"))
    }
  }
  
  # add themes
  plot <- (plot
           + theme_light()
           + theme(plot.title = element_text(hjust = 0.5),
                   legend.position = c("top","left")))
  
  return(plot)
}


#**************************************************
# Plot correlation matrix in heatmap ==============
#**************************************************
plot_cor_heatmap<-function(input,label=NULL){
  
  input_tidy<-rownames_to_column(input,"row")
  input_tidy<-gather(input_tidy,key=column,value=r,2:ncol(input_tidy))

  # Text label overlay on heatmap
  if (!is.null(label)){
    label_tidy<-rownames_to_column(label,"row")
    label_tidy<-gather(label_tidy,key=column,value=label,2:ncol(label_tidy))
    input_tidy<-inner_join(input_tidy,label_tidy,by=c("row","column"))
  }else{
    input_tidy$label<-NA
  }
  
  plot<-(ggplot(input_tidy, aes(column, row))
         + geom_tile(aes(fill = r))
         + scale_fill_gradientn(colors = matlab.like2(100),name="r",limits=c(-1,1))
         + scale_y_discrete(limits = rev(rownames(input)))
         + scale_x_discrete(limits = colnames(input), position="top")
         + theme_linedraw()
         + theme(axis.text.x = element_text(size=29/log(ncol(input),2),angle = 90,vjust=0,hjust=0),
                 axis.text.y = element_text(size=29/log(ncol(input),2)),
                 #axis.title=element_blank(),
                 panel.grid.major=element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank()))
  
  if(!is.null(label)){
    plot<-(plot
           + geom_text(aes(label=label),color="black",fontface="bold"))
  }
  return(plot)
}


#**************************************************
# Helper functions for ggpairs() ==================
#**************************************************
# corr + heatmap plot for upper half of correlogram
custom_corr_heatmap <- function(data, mapping, ...) {
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  ct <- cor.test(x,y)
  r <- unname(ct$estimate)
  r_text <- format(r, digits=3)
  sig_text <- symnum(
    ct$p.value, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 1),
    symbols = c("***", "** ", "*  ", "   ")
  )
  text<-paste(r_text,sig_text,sep=" ")
  
  
  colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
  fill <- colFn(100)[findInterval(r, seq(-1, 1, length=100))]
  
  ggally_text(
    label = text, 
    mapping = aes(),
    xP = 0.5, yP = 0.5, 
    size = 5,
    color = I("black"),
    ...
  ) + 
    theme_void() +
    theme(panel.background = element_rect(fill=fill))
}

# small line width, xlim and ylim, loess
custom_smooth <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(color = I("black"),alpha=0.01,size=0.001,stroke = 0, shape = ".") + 
    geom_smooth(method = "lm", color = I("red"),size=0.1,linetype="dotted",fill=I("pink"), ...) +
    geom_smooth(method = "loess", color = I("blue"),size=0.1,linetype="dotted",fill=I("lightblue"), ...) +
    xlim(-1,1) +
    ylim(-1,1)
}

# xlim and ylim, theme change
custom_densityDiag <- function(data, mapping, ...){
  
  mapping <- mapping_color_to_fill(mapping)
  
  ggplot(data, mapping) +
    theme_minimal() +
    scale_y_continuous() +
    xlim(-1,1) +
    geom_density(linetype="blank", fill=I("grey"),...)
}


#**************************************************
# OBSOLETE ========================================
#**************************************************

# older version in circular style
plot_gam_fc_old<-function(df_plot_gamm,df_roi,analysis,atlas,list_mod,list_plot,
                          list_type_p,thr_p,paths_,suffix_){
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
            df_plot_gamm_subset<-rename(df_plot_gamm_subset,c("estimate"="weight"))
          }else{
            df_plot_gamm_subset<-rename(df_plot_gamm_subset,c("F"="weight"))
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
              #ggsave(paste("atl-",atlas,"_anl-",analysis,"_mod-",idx_mod,"_plt-",var_exp,
              #             "_sex-",label_sex,"_pval-",type_p,"_",suffix_,"_gamm_fc.eps",sep=""),
              #       plot=plot,device=cairo_ps,path=file.path(paths_$output,"output"),
              #       dpi=300,height=10,width=10,limitsize=F)
              ggsave(paste("atl-",atlas,"_anl-",analysis,"_mod-",idx_mod,"_plt-",var_exp,
                           "_sex-",label_sex,"_pval-",type_p,"_",suffix_,"_gamm_fc.png",sep=""),
                     plot=plot,path=file.path(paths_$output,"output","plot"),height=10,width=10,dpi=600)
            }
          }
        }
      }
    }
  }
}
