#**************************************************
# Description =====================================
#**************************************************
# R script to analyze longitudinal clinical data.


#**************************************************
# Parameters ======================================
#**************************************************
path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/clin"
#path_exp <- "Dropbox/MRI/pnTTC/Puberty/Stats/func_XCP/test_5sub"
path_exp_full<-NULL

dir_in<-""
dir_out<-"01_clin_test"
#dir_out<-"02_clin_pair"
#dir_out<-"03_clin_long_test"


#**************************************************
# Libraries =======================================
#**************************************************
library(mgcv)
library(dplyr)
library(ggplot2)
#library(itsadug)
library(ggrepel)
library(ggpubr)
library(plyr)


#**************************************************
# Original library ================================
#**************************************************
source(file.path(getwd(),"util/function.R"))
source(file.path(getwd(),"util/plot.R"))
source(file.path(getwd(),"util/parameter.R"))
paths<-func_path(path_exp_=path_exp,dir_in_=dir_in,dir_out_=dir_out,path_exp_full_=path_exp_full)


#**************************************************
# Standardize clinical data using mixed models ====
#**************************************************
std_clin<-function(paths_=paths,param=param_std_clin){
  print("Starting std_clin().")
  nullobj<-func_createdirs(paths_,"std_clin()")
  
  list_depvar<-c(param$list_tanner,param$list_hormone)
  df_gamm_bind<-df_aic_bind<-df_anova_bind<-df_fit_bind<-data.frame()
  for (depvar in names(list_depvar)){
    print(paste('Calculating: ',list_depvar[[depvar]][["label"]]))
    
    # Prepare source dataframe
    list_covar<-param$list_covar
    list_covar[["depvar"]]<-list_depvar[[depvar]]
    df_plot<-func_clinical_data_long(paths_,param$list_wave,param$subset_subj,list_covar,
                                       rem_na_clin=T,prefix=paste("depvar-",depvar,sep=""),print_terminal=F)$df_clin
    df_plot$sex=as.factor(df_plot$sex)
    df_plot$wave=as.factor(df_plot$wave)
    df_plot<-df_plot[,c("ID_pnTTC","age","sex","wave","depvar")]
    write.csv(df_plot,file.path(paths_$output,"output","temp",paste("depvar-",depvar,"_src.csv",sep="")),row.names = F)
    colnames(df_plot)[colnames(df_plot)=="depvar"]<-"value"
    
    data_gamm<-gamm_core3(df_plot,list_mod_in=param$list_mod,list_sex_in=list(1,2),
                          calc_parallel_in=F,test_mod_in=T)
    
    df_gamm_bind<-rbind(df_gamm_bind,cbind("depvar"=depvar,data_gamm$df_gamm))
    df_aic_bind<-rbind(df_aic_bind,cbind("depvar"=depvar,data_gamm$df_aic))
    df_anova_bind<-rbind(df_anova_bind,cbind("depvar"=depvar,data_gamm$df_anova))

    for (idx_mod in names(param$list_mod)){
      plot<-NULL
      for (idx_sex in c(1,2)){
        mod_gamm<-data_gamm[["mod"]][[paste("mod-",idx_mod,"_sex-",idx_sex,sep="")]]
        df_plot_sex<-df_plot[df_plot$sex==idx_sex,]
        df_fit<-cbind(df_plot_sex,as.data.frame(predict.gam(mod_gamm, df_plot_sex[,c("ID_pnTTC","age")], exclude="s(ID_pnTTC)",se.fit = TRUE)))
        df_fit$z<-(df_fit$value-df_fit$fit)/(df_fit$se.fit*sqrt(nrow(df_fit)))
        df_fit_bind<-rbind(df_fit_bind,cbind("depvar"=depvar,"mod"=idx_mod,df_fit))
        plot<-plot_gamm(plot_in=plot,mod_gamm,df_plot_sex,param$spec_graph)
      }
      label_axis_y<-list_depvar[[depvar]][["label"]]
      plot<-(plot
             + ggtitle(paste(label_axis_y,' - Age, model: ',idx_mod,sep=''))
             + xlab("Age (day)") + ylab(label_axis_y)
             + theme(legend.position = "none"))
      ggsave(paste("depvar-",depvar,"_mod-",idx_mod,"_clin_long.eps",sep=""),plot=plot,device=cairo_ps,
             path=file.path(paths_$output,"output","plot"),
             dpi=600,height=7,width=7,limitsize=F)
    }
  }
  write.csv(df_fit_bind,file.path(paths_$output,"output","result","std_clin.csv"),row.names=F)
  write.csv(df_gamm_bind,file.path(paths_$output,"output","result","gamm.csv"),row.names=F)
  write.csv(df_aic_bind,file.path(paths_$output,"output","result","aic.csv"),row.names=F)
  write.csv(df_anova_bind,file.path(paths_$output,"output","result","anova.csv"),row.names=F)
  print("Finished std_clin().")
}


#**************************************************
# Plot Tanner and hormone longitudinal data =======
#**************************************************
plot_long<-function(paths_=paths,list_wave_=list_wave,subset_subj_=subset_subj,
                    list_tanner_=list_tanner,list_hormone_=list_hormone
                    ){
  print("Starting plot_long().")
  nullobj<-func_createdirs(paths_,"plot_long()")
  list_covar<-c(list_tanner_,list_hormone_,
                list("age"   =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                     "sex"   =list("1"="Sex",          "2"="Sex",          "label"="Sex")))
  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar,rem_na_clin=T)
  df_plot<-data_clin$df_clin
  colnames(df_plot)[colnames(df_plot)=="wave"]<-"ses"
  list_subj_long<-intersect(data_clin[["list_id_exist"]][["1"]][["intersect"]],
                            data_clin[["list_id_exist"]][["2"]][["intersect"]])
  df_plot<-func_clinical_data_join(df_src=df_plot,list_id_subj=list_subj_long,
                                   list_covar=list_covar)
  list_sex<-list("mf"=c(1,2),"m"=1,"f"=2)
  for (label_sex in names(list_sex)){
    df_plot_sex<-df_plot[df_plot$sex %in% list_sex[[label_sex]],c(-1,-2)]
    data_cor<-func_cor(df_plot_sex)
    plot<-plot_cor_heatmap(data_cor$cor)
    plot<-(plot
           + ggtitle(paste("Correlation among clinical variables, ",label_sex, sep=""))
           + xlab("Clinical variables")
           + ylab("Clinical variables")
           + theme(plot.title = element_text(hjust = 0.5),legend.position = "none"))
    filename_plot<-paste("sex-",label_sex,"_long.eps",sep="")
    ggsave(filename_plot,plot=plot,device=cairo_ps,
           path=file.path(paths_$output,"output"),dpi=600,height=7,width=7,limitsize=F)
  }
  print("Finished plot_long().")
}


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
      filename_plot<-paste("tanner-",tanner,"_hormone-",hormone,"_mod-",idx_mod,"_pair.eps",sep="")
      ggsave(filename_plot,plot=plot,device=cairo_ps,
             path=file.path(paths_$output,"output"),dpi=600,height=7,width=7,limitsize=F)
    }
  }
  write.csv(df_out_lm,file.path(paths_$output,"output","gamm_result.csv"),row.names=F)
  print("Finished plot_pair().")
}


#OBSOLETE, incorporated into std_clin
#**************************************************
# Plot Longitudinal clinical data =================
#**************************************************
#plot_clin<-function(paths_=paths,list_wave_=list_wave,subset_subj_=subset_subj,
#                    list_mod_=list_mod,
#                    list_tanner_=list_tanner,list_covar_tanner_=list_covar_tanner,
#                    spec_graph_tanner_=spec_graph_tanner,
#                    list_hormone_=list_hormone,list_covar_hormone_=list_covar_hormone,
#                    spec_graph_hormone_=spec_graph_hormone
#                    ){
#  print("Starting plot_clin().")
#  nullobj<-func_createdirs(paths_)
#
#  # Longitudinal clinical plot
#  for (tanner in names(list_tanner_)){
#    print(paste('Calculating: ',list_tanner_[[tanner]][["label"]]))
#    list_covar<-list_covar_tanner_
#    list_covar[["tanner"]]<-list_tanner_[[tanner]]
#    data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar,rem_na_clin=T)
#    df_plot<-data_clin$df_clin
#    df_plot$sex=as.factor(df_plot$sex)
#    df_plot$wave=as.factor(df_plot$wave)
#    colnames(df_plot)[colnames(df_plot)==tanner]<-"tanner"
#    df_plot<-df_plot[,c("ID_pnTTC","age","sex","wave","tanner")]
#    write.csv(df_plot,file.path(paths_$output,"output",paste("var-tanner_",tanner,"_src.csv",sep="")),row.names = F)
#    
#    # GAM fit plot
#    for (idx_mod in names(list_mod_)){
#      plot<-NULL
#      for (idx_sex in c(1,2)){
#        df_plot_sex<-df_plot[df_plot$sex==idx_sex,]
#        colnames(df_plot_sex)[colnames(df_plot_sex)=="tanner"]<-"value"
#        mod_gamm<-gam(as.formula(list_mod_[[idx_mod]]),data=df_plot_sex)
#        plot<-plot_gamm(plot_in=plot,mod_gamm,df_plot_sex,spec_graph_tanner_)
#      }
#      label_axis_y<-list_tanner_[[tanner]][["label"]]
#      plot<-(plot
#             + ggtitle(paste(label_axis_y,' - Age','\n',idx_mod,sep=''))
#             + xlab("Age")
#             + ylab(label_axis_y)
#             + theme(legend.position = "none"))
#      filename_plot<-paste("var-tanner_",tanner,"_mod-",idx_mod,"_clin_long.eps",sep="")
#      ggsave(filename_plot,plot=plot,device=cairo_ps,
#             path=file.path(paths_$output,"output"),dpi=600,height=7,width=7,limitsize=F)
#    }
#  }
#  
#  for (hormone in names(list_hormone_)){
#    print(paste('Calculating: ',list_hormone_[[hormone]][["label"]]))
#    list_covar<-list_covar_hormone_
#    list_covar[["hormone"]]<-list_hormone_[[hormone]]
#    data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar,rem_na_clin=T)
#    df_plot<-data_clin$df_clin
#    df_plot$sex=as.factor(df_plot$sex)
#    df_plot$wave=as.factor(df_plot$wave)
#    colnames(df_plot)[colnames(df_plot)==hormone]<-"hormone"
#    df_plot<-df_plot[,c("ID_pnTTC","age","sex","wave","hormone")]
#    write.csv(df_plot,file.path(paths_$output,"output",paste("var-hormone_",hormone,"_src.csv",sep="")),row.names = F)
#    
#    # GAM fit plot
#    for (idx_mod in names(list_mod_)){
#      plot<-NULL
#      for (idx_sex in c(1,2)){
#        df_plot_sex<-df_plot[df_plot$sex==idx_sex,]
#        colnames(df_plot_sex)[colnames(df_plot_sex)=="hormone"]<-"value"
#        mod_gamm<-gam(as.formula(list_mod_[[idx_mod]]),data=df_plot_sex)
#        plot<-plot_gamm(plot_in=plot,mod_gamm,df_plot_sex,spec_graph_hormone_)
#      }
#      label_axis_y<-list_hormone_[[hormone]][["label"]]
#      plot<-(plot
#             + ggtitle(paste(label_axis_y,' - Age','\n',idx_mod,sep=''))
#             + xlab("Age")
#             + ylab(label_axis_y)
#             + theme(legend.position = "none"))
#      filename_plot<-paste("var-hormone_",hormone,"_mod-",idx_mod,"_clin_long.eps",sep="")
#      ggsave(filename_plot,plot=plot,device=cairo_ps,
#             path=file.path(paths_$output,"output"),dpi=600,height=7,width=7,limitsize=F)
#    }
#  }
#  print("Finished plot_clin().")
#}


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