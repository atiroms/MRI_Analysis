#**************************************************
# Description =====================================
#**************************************************
# R script to analyze longitudinal clinical data.


#**************************************************
# Parameters ======================================
#**************************************************
path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/clin"
#path_exp <- "Dropbox/MRI/pnTTC/Puberty/Stats/func_XCP/test_5sub"

dir_in<-""
#dir_out<-"01_clin_long"
dir_out<-"02_clin_pair"

list_wave <- c(1,2)

#subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
#                             list("key"="W1_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)),
#                    "2"=list(list("key"="W2_T1QC","value"=1),
#                             list("key"="W2_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)))
subset_subj <- list("1"=list(),
                    "2"=list())
list_mod <- list("lin" ="value ~ age + s(ID_pnTTC,bs='re')",
                 "add" ="value ~ s(age,k=3) + s(ID_pnTTC,bs='re')",
                 "quad"="value ~ poly(age,2) + s(ID_pnTTC,bs='re')")

list_tanner<-list("max" =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage (max)"),
                  "full"=list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)"),
                  "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),
                                 "2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                 "label"="Tanner stage (gonadal)"),
                  "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),
                                 "2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
                                 "label"="Tanner stage (adrenal)"))
list_covar_tanner<-list("tanner"=list("1"="W1_Tanner_Max","2"="W2_Tanner_Max","label"="Tanner stage (max)"),
                        "age"   =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                        "sex"   =list("1"="Sex",          "2"="Sex",          "label"="Sex"))
spec_graph_tanner<-list("title"="Tanner vs Age","x_axis"="age",
                        "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                      "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                        "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                     "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1)))

list_hormone<-list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
                   "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
                   "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
                   "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S"))
list_covar_hormone<-list("hormone"=list("1"="W1_Hormone",   "2"="W2_Hormone",   "label"="Hormone"),
                         "age"    =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                         "sex"    =list("1"="Sex",          "2"="Sex",          "label"="Sex"))
spec_graph_hormone<-list("title"="Hormone vs Age","x_axis"="age",
                         "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                       "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                         "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                      "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1)))

list_pair<-list(c("gonadal","testo"),c("gonadal","dheas"),c("adrenal","testo"),c("adrenal","dheas"),
                c("max","testo"),c("max","dheas"))
list_mod_pair <- list("lin" ="value ~ hormone + s(ID_pnTTC,bs='re')",
                      "add" ="value ~ s(hormone,k=3) + s(ID_pnTTC,bs='re')",
                      "quad"="value ~ poly(hormone,2) + s(ID_pnTTC,bs='re')")
list_covar_pair<-list("tanner"=list("1"="W1_Tanner_Max","2"="W2_Tanner_Max","label"="Tanner stage (max)"),
                      "hormone"=list("1"="W1_Hormone",   "2"="W2_Hormone",   "label"="Hormone"),
                      "age"   =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                      "sex"   =list("1"="Sex",          "2"="Sex",          "label"="Sex"))
spec_graph_pair<-list("title"="Tanner vs Hormone","x_axis"="hormone",
                      "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                    "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                      "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                   "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1)))


#**************************************************
# Libraries =======================================
#**************************************************
library(mgcv)
library(dplyr)
library(ggplot2)
library(itsadug)
library(ggrepel)
library(ggpubr)
library(plyr)


#**************************************************
# Create path list ================================
#**************************************************
func_path<-function(list_path_root = c("D:/atiroms","C:/Users/atiro","/home/atiroms","C:/Users/NICT_WS"),
                    path_exp_=path_exp,
                    dir_in_=dir_in,
                    dir_out_=dir_out){
  path_root<-NA
  for(p in list_path_root){
    if(file.exists(p)){
      path_root<-p
    }
  }
  if(is.na(path_root)){
    print("Error: root path could not be found.")
  }
  path_script <- file.path(path_root,"GitHub/MRI_Analysis")
  path_common <- file.path(path_root,"DropBox/MRI_img/pnTTC/puberty/common")
  path_in     <- file.path(path_root,path_exp_,dir_in_)
  path_out    <- file.path(path_root,path_exp_,dir_out_)
  output <- list("script"=path_script,"input"=path_in,"output"=path_out,
                 "common"=path_common,"dir_in"=dir_in_,"dir_out"=dir_out_)
  return(output)
}

paths<-func_path()


#**************************************************
# Original library ================================
#**************************************************
source(file.path(paths$script,"util/function.R"))
source(file.path(paths$script,"util/plot.R"))


#**************************************************
# Plot Tanner and hormone data ====================
#**************************************************
plot_pair<-function(paths_=paths,list_wave_=list_wave,subset_subj_=subset_subj,
                    list_mod_pair_=list_mod_pair,list_pair_=list_pair,
                    list_tanner_=list_tanner,list_hormone_=list_hormone,
                    spec_graph_pair_=spec_graph_pair,list_covar_pair_=list_covar_pair
                    ){
  print("Starting plot_pair().")
  nullobj<-func_createdirs(paths_)
  df_out_lm<-data.frame()
  for (pair in list_pair_){
    tanner<-pair[1]
    hormone<-pair[2]
    print(paste('Calculating: ',list_tanner_[[tanner]][["label"]],
                " and ",list_hormone_[[hormone]][["label"]],sep=""))
    list_covar<-list_covar_pair_
    list_covar[["tanner"]]<-list_tanner_[[tanner]]
    list_covar[["hormone"]]<-list_hormone_[[hormone]]
    data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar,rem_na_clin=T)
    df_plot<-data_clin$df_clin
    df_plot$sex=as.factor(df_plot$sex)
    df_plot$wave=as.factor(df_plot$wave)
    colnames(df_plot)[colnames(df_plot)==tanner]<-"tanner"
    colnames(df_plot)[colnames(df_plot)==hormone]<-"hormone"
    df_plot<-df_plot[,c("ID_pnTTC","sex","wave","tanner","hormone")]
    write.csv(df_plot,file.path(paths_$output,"output",paste("tanner_",tanner,"_hormone-",hormone,"_src.csv",sep="")),row.names = F)
    
    # GAM fit plot
    for (idx_mod in names(list_mod_pair_)){
      plot<-NULL
      for (idx_sex in c(1,2)){
        df_plot_sex<-df_plot[df_plot$sex==idx_sex,]
        colnames(df_plot_sex)[colnames(df_plot_sex)=="tanner"]<-"value"
        mod_gamm<-gam(as.formula(list_mod_pair_[[idx_mod]]),data=df_plot_sex)
        p_table<-summary.gam(mod_gamm)$p.table
        if (is.null(summary.gam(mod_gamm)$s.table)){
          df_out_lm_add<-data.frame(tanner=tanner,hormone=hormone,sex=idx_sex,model=idx_mod,term=rownames(p_table),
                                    estimate=p_table[,'Estimate'],se=p_table[,'Std. Error'],F=NA,
                                    t=p_table[,'t value'],p=p_table[,'Pr(>|t|)'])
          
        }else{
          s_table<-summary.gam(mod_gamm)$s.table
          df_out_lm_add<-rbind(data.frame(tanner=tanner,hormone=hormone,sex=idx_sex,model=idx_mod,term=rownames(p_table),
                                          estimate=p_table[,'Estimate'],se=p_table[,'Std. Error'],F=NA,
                                          t=p_table[,'t value'],p=p_table[,'Pr(>|t|)']),
                               data.frame(tanner=tanner,hormone=hormone,sex=idx_sex,model=idx_mod,term=rownames(s_table),
                                          estimate=NA,se=NA,F=s_table[,'F'],
                                          t=NA,p=s_table[,'p-value']))
        }
        df_out_lm<-rbind(df_out_lm,df_out_lm_add)
        plot<-plot_gamm(plot_in=plot,mod_gamm,df_plot_sex,spec_graph_pair_)
      }
      label_axis_x<-list_hormone_[[hormone]][["label"]]
      label_axis_y<-list_tanner_[[tanner]][["label"]]
      plot<-(plot
             + ggtitle(paste(label_axis_y,' - ',label_axis_x,'\n',idx_mod,sep=''))
             + xlab(label_axis_x)
             + ylab(label_axis_y)
             + theme(legend.position = "none"))
      filename_plot<-paste("tanner_",tanner,"_hormone-",hormone,"_mod-",idx_mod,"_pair.eps",sep="")
      ggsave(filename_plot,plot=plot,device=cairo_ps,
             path=file.path(paths_$output,"output"),dpi=600,height=7,width=7,limitsize=F)
    }
  }
  write.csv(df_out_lm,file.path(paths_$output,"output","gamm_result.csv"),row.names=F)
  print("Finished plot_pair().")
}


#**************************************************
# Plot Longitudinal clinical data =================
#**************************************************
plot_clin<-function(paths_=paths,list_wave_=list_wave,subset_subj_=subset_subj,
                    list_mod_=list_mod,
                    list_tanner_=list_tanner,list_covar_tanner_=list_covar_tanner,
                    spec_graph_tanner_=spec_graph_tanner,
                    list_hormone_=list_hormone,list_covar_hormone_=list_covar_hormone,
                    spec_graph_hormone_=spec_graph_hormone
                    ){
  print("Starting plot_clin().")
  nullobj<-func_createdirs(paths_)

  # Longitudinal clinical plot
  for (tanner in names(list_tanner_)){
    print(paste('Calculating: ',list_tanner_[[tanner]][["label"]]))
    list_covar<-list_covar_tanner_
    list_covar[["tanner"]]<-list_tanner_[[tanner]]
    data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar,rem_na_clin=T)
    df_plot<-data_clin$df_clin
    df_plot$sex=as.factor(df_plot$sex)
    df_plot$wave=as.factor(df_plot$wave)
    colnames(df_plot)[colnames(df_plot)==tanner]<-"tanner"
    df_plot<-df_plot[,c("ID_pnTTC","age","sex","wave","tanner")]
    write.csv(df_plot,file.path(paths_$output,"output",paste("var-tanner_",tanner,"_src.csv",sep="")),row.names = F)
    
    # GAM fit plot
    for (idx_mod in names(list_mod_)){
      plot<-NULL
      for (idx_sex in c(1,2)){
        df_plot_sex<-df_plot[df_plot$sex==idx_sex,]
        colnames(df_plot_sex)[colnames(df_plot_sex)=="tanner"]<-"value"
        mod_gamm<-gam(as.formula(list_mod_[[idx_mod]]),data=df_plot_sex)
        plot<-plot_gamm(plot_in=plot,mod_gamm,df_plot_sex,spec_graph_tanner_)
      }
      label_axis_y<-list_tanner_[[tanner]][["label"]]
      plot<-(plot
             + ggtitle(paste(label_axis_y,' - Age','\n',idx_mod,sep=''))
             + xlab("Age")
             + ylab(label_axis_y)
             + theme(legend.position = "none"))
      filename_plot<-paste("var-tanner_",tanner,"_mod-",idx_mod,"_clin_long.eps",sep="")
      ggsave(filename_plot,plot=plot,device=cairo_ps,
             path=file.path(paths_$output,"output"),dpi=600,height=7,width=7,limitsize=F)
    }
  }
  
  for (hormone in names(list_hormone_)){
    print(paste('Calculating: ',list_hormone_[[hormone]][["label"]]))
    list_covar<-list_covar_hormone_
    list_covar[["hormone"]]<-list_hormone_[[hormone]]
    data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar,rem_na_clin=T)
    df_plot<-data_clin$df_clin
    df_plot$sex=as.factor(df_plot$sex)
    df_plot$wave=as.factor(df_plot$wave)
    colnames(df_plot)[colnames(df_plot)==hormone]<-"hormone"
    df_plot<-df_plot[,c("ID_pnTTC","age","sex","wave","hormone")]
    write.csv(df_plot,file.path(paths_$output,"output",paste("var-hormone_",hormone,"_src.csv",sep="")),row.names = F)
    
    # GAM fit plot
    for (idx_mod in names(list_mod_)){
      plot<-NULL
      for (idx_sex in c(1,2)){
        df_plot_sex<-df_plot[df_plot$sex==idx_sex,]
        colnames(df_plot_sex)[colnames(df_plot_sex)=="hormone"]<-"value"
        mod_gamm<-gam(as.formula(list_mod_[[idx_mod]]),data=df_plot_sex)
        plot<-plot_gamm(plot_in=plot,mod_gamm,df_plot_sex,spec_graph_hormone_)
      }
      label_axis_y<-list_hormone_[[hormone]][["label"]]
      plot<-(plot
             + ggtitle(paste(label_axis_y,' - Age','\n',idx_mod,sep=''))
             + xlab("Age")
             + ylab(label_axis_y)
             + theme(legend.position = "none"))
      filename_plot<-paste("var-hormone_",hormone,"_mod-",idx_mod,"_clin_long.eps",sep="")
      ggsave(filename_plot,plot=plot,device=cairo_ps,
             path=file.path(paths_$output,"output"),dpi=600,height=7,width=7,limitsize=F)
    }
  }
  print("Finished plot_clin().")
}


# OBSOLETE
#**************************************************
# Plot Tanner stage clinical data =================
#**************************************************
#plot_clin_tanner<-function(paths_=paths,
#                           list_wave_=list_wave,list_covar_=list_covar,
#                           subset_subj_=subset_subj
#                           ){
#  print("Starting plot_clin_tanner().")
#  nullobj<-func_createdirs(paths_)
#  
#  # Load and subset clinical data according to specified subsetting condition and covariate availability
#  print('Loading clinical data.')
#  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar_,rem_na_clin=T)
#  df_clin<-data_clin$df_clin
#  df_clin$sex=as.factor(df_clin$sex)
#  df_clin$wave=as.factor(df_clin$wave)
#  
#  write.csv(df_clin,file.path(paths_$output,"output",
#                              paste("plot_src.csv",sep="")),row.names = F)
#  
#  list_plot<-list()
#  
#  # Longitudinal clinical plot
#  plot <- (ggplot(df_clin)
#           + aes(x=age,y=tanner,color=sex,shape=wave, label=ID_pnTTC)
#           + scale_colour_manual(name=NULL,labels=c("Male","Female"),values=c("steelblue2","lightcoral"))
#           + scale_shape_manual(name=NULL,labels=c("1st wave","2nd wave"),values=c(3,4))
#           + geom_point(size=4)
#           #geom_text_repel(size=2)
#           + geom_path(aes(group=ID_pnTTC,color=sex),size=0.5,alpha=0.5)
#           + ggtitle("Tanner stage vs Age")
#           + xlab("Age (day)")
#           + ylab("Tanner Stage")
#           + theme_light()
#           + theme(plot.title = element_text(hjust = 0.5),legend.justification=c(0,1),
#                   legend.position=c(0.05,0.95),legend.direction="horizontal",
#                   panel.grid.minor=element_blank()))
#  ggsave("tanner_age_path.eps",plot=plot,device=cairo_ps,
#         path=file.path(paths_$output,"output"),dpi=300,height=7,width=7,limitsize=F)
#  plot<-list(plot)
#  names(plot)<-"tanner_age_path"
#  list_plot<-c(list_plot,plot)
#  
#  # Bar plot
#  for (wave in c(1,2)){
#    df_clin_wave<-df_clin[df_clin$wave==wave,]
#    plot <- (ggplot(as.data.frame(with(df_clin_wave, table(tanner = factor(tanner), sex))),
#                 aes(factor(tanner),y = Freq, fill = sex))
#             + geom_col(width=0.7, position = position_dodge(width=0.8))
#             + scale_fill_manual(name=NULL,labels=c("Male","Female"),
#                                 values=c("steelblue2","lightcoral"),drop=FALSE)
#             + ggtitle(paste("Tanner Stage for pn-TTC Wave:",wave,sep=" "))
#             + xlab("Tanner Stage")
#             + ylab("Number of Subjects")
#             + theme_light()
#             + theme(plot.title = element_text(hjust = 0.5),legend.justification=c(1,1),
#                     legend.position=c(0.95,0.95),panel.grid.major.x=element_blank()))
#    ggsave(paste("ses-",as.character(wave),"_tanner_bar.eps",sep=""),plot=plot,device=cairo_ps,
#           path=file.path(paths_$output,"output"),dpi=300,height=5,width=5,limitsize=F)
#    plot<-list(plot)
#    names(plot)<-paste("ses-",as.character(wave),"_tanner_bar",sep="")
#    list_plot<-c(list_plot,plot)
#  }
#  
#  # Heatmap of Tanner counts
#  list_sex<-list("all"=c(1,2),"male"=1,"female"=2)
#  for (id_sex in names(list_sex)){
#    df_clin_sex<-df_clin[df_clin$sex %in% list_sex[[id_sex]],]
#    df_heatmap<-data.frame(matrix(ncol=5,nrow=5))
#    for (tanner_ses1 in seq(5)){
#      list_id_ses1<-df_clin_sex[df_clin_sex$wave==1 & df_clin_sex$tanner==tanner_ses1,"ID_pnTTC"]
#      list_id_ses1<-list_id_ses1[!is.na(list_id_ses1)]
#      for (tanner_ses2 in seq(5)){
#        list_id_ses2<-df_clin_sex[df_clin_sex$wave==2 & df_clin_sex$tanner==tanner_ses2,"ID_pnTTC"]
#        list_id_ses2<-list_id_ses2[!is.na(list_id_ses2)]
#        n_id_intersect<-length(intersect(list_id_ses1,list_id_ses2))
#        df_heatmap[tanner_ses1,tanner_ses2]<-n_id_intersect
#      }
#    }
#    colnames(df_heatmap)<-rownames(df_heatmap)<-as.character(seq(5))
#    write.csv(df_heatmap,file.path(paths_$output,"output",
#                                   paste("sex-",id_sex,"_tanner_heatmap.csv",sep="")))
#    
#    plot<-plot_cor_heatmap(df_heatmap)
#    plot <- (plot
#             + scale_fill_viridis(name="N")
#             + ggtitle("Tanner stage")
#             + xlab("2nd wave")
#             + ylab("1st wave")
#             + theme(plot.title = element_text(hjust = 0.5),
#                     axis.text.x = element_text(size=8,angle = 0,vjust=0,hjust=0.5),
#                     axis.text.y = element_text(size=8)))
#    ggsave(paste("sex-",id_sex,"_tanner_heatmap.eps",sep=""),plot=plot,device=cairo_ps,
#           path=file.path(paths_$output,"output"),dpi=300,height=5,width=5,limitsize=F)
#    plot<-list(plot)
#    list_plot<-c(list_plot,plot)
#    names(plot)<-paste("sex-",id_sex,"_tanner_heatmap",sep="")
#  }
#  
#  print("Finished plot_clin_tanner().")
#  return(list_plot)
#}