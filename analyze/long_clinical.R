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
dir_out<-"01_clin_plot"

list_wave <- c(1,2)

#subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
#                             list("key"="W1_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)),
#                    "2"=list(list("key"="W2_T1QC","value"=1),
#                             list("key"="W2_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)))

subset_subj <- list("1"=list(),
                    "2"=list())

list_covar<-list("tanner"=list("1"="W1_Tanner_Max",
                               "2"="W2_Tanner_Max",
                               "label"="Tanner stage"),
                 "age"=list("1"="W1_Age_at_MRI",
                            "2"="W2_Age_at_MRI",
                            "label"="Age"),
                 "sex"=list("1"="Sex",
                            "2"="Sex",
                            "label"="Sex"))

#df_clinical<-read.csv("C:/Users/atiro/Dropbox/MRI_img/pnTTC/puberty/common/CSUB.csv")


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
# Plot clinical data ==============================
#**************************************************
plot_clin<-function(paths_=paths,
                    list_wave_=list_wave,list_covar_=list_covar,
                    subset_subj_=subset_subj
                    ){
  print("Starting plot_clin().")
  nullobj<-func_createdirs(paths_)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar_,rem_na_clin=T)
  df_clin<-data_clin$df_clin
  df_clin$sex=as.factor(df_clin$sex)
  df_clin$wave=as.factor(df_clin$wave)
  
  write.csv(df_clin,file.path(paths_$output,"output",
                              paste("plot_src.csv",sep="")),row.names = F)
  
  list_plot<-list()
  
  # Longitudinal clinical plot
  plot <- (ggplot(df_clin)
           + aes(x=age,y=tanner,color=sex,shape=wave, label=ID_pnTTC)
           + scale_colour_manual(name=NULL,labels=c("Male","Female"),values=c("steelblue2","lightcoral"))
           + scale_shape_manual(name=NULL,labels=c("1st wave","2nd wave"),values=c(3,4))
           + geom_point(size=4)
           #geom_text_repel(size=2)
           + geom_path(aes(group=ID_pnTTC,color=sex),size=0.5,alpha=0.5)
           + ggtitle("Tanner stage vs Age")
           + xlab("Age (day)")
           + ylab("Tanner Stage")
           + theme_light()
           + theme(plot.title = element_text(hjust = 0.5),legend.justification=c(0,1),
                   legend.position=c(0.05,0.95),legend.direction="horizontal",
                   panel.grid.minor=element_blank()))
  ggsave("tanner_age_path.eps",plot=plot,device=cairo_ps,
         path=file.path(paths_$output,"output"),dpi=300,height=7,width=7,limitsize=F)
  plot<-list(plot)
  names(plot)<-"tanner_age_path"
  list_plot<-c(list_plot,plot)
  
  # Bar plot
  for (wave in c(1,2)){
    df_clin_wave<-df_clin[df_clin$wave==wave,]
    plot <- (ggplot(as.data.frame(with(df_clin_wave, table(tanner = factor(tanner), sex))),
                 aes(factor(tanner),y = Freq, fill = sex))
             + geom_col(width=0.7, position = position_dodge(width=0.8))
             + scale_fill_manual(name=NULL,labels=c("Male","Female"),
                                 values=c("steelblue2","lightcoral"),drop=FALSE)
             + ggtitle(paste("Tanner Stage for pn-TTC Wave:",wave,sep=" "))
             + xlab("Tanner Stage")
             + ylab("Number of Subjects")
             + theme_light()
             + theme(plot.title = element_text(hjust = 0.5),legend.justification=c(1,1),
                     legend.position=c(0.95,0.95),panel.grid.major.x=element_blank()))
    ggsave(paste("ses-",as.character(wave),"_tanner_bar.eps",sep=""),plot=plot,device=cairo_ps,
           path=file.path(paths_$output,"output"),dpi=300,height=5,width=5,limitsize=F)
    plot<-list(plot)
    names(plot)<-paste("ses-",as.character(wave),"_tanner_bar",sep="")
    list_plot<-c(list_plot,plot)
  }
  
  # Heatmap of Tanner counts
  list_sex<-list("all"=c(1,2),"male"=1,"female"=2)
  for (id_sex in names(list_sex)){
    df_clin_sex<-df_clin[df_clin$sex %in% list_sex[[id_sex]],]
    df_heatmap<-data.frame(matrix(ncol=5,nrow=5))
    for (tanner_ses1 in seq(5)){
      list_id_ses1<-df_clin_sex[df_clin_sex$wave==1 & df_clin_sex$tanner==tanner_ses1,"ID_pnTTC"]
      list_id_ses1<-list_id_ses1[!is.na(list_id_ses1)]
      for (tanner_ses2 in seq(5)){
        list_id_ses2<-df_clin_sex[df_clin_sex$wave==2 & df_clin_sex$tanner==tanner_ses2,"ID_pnTTC"]
        list_id_ses2<-list_id_ses2[!is.na(list_id_ses2)]
        n_id_intersect<-length(intersect(list_id_ses1,list_id_ses2))
        df_heatmap[tanner_ses1,tanner_ses2]<-n_id_intersect
      }
    }
    colnames(df_heatmap)<-rownames(df_heatmap)<-as.character(seq(5))
    write.csv(df_heatmap,file.path(paths_$output,"output",
                                   paste("sex-",id_sex,"_tanner_heatmap.csv",sep="")))
    
    plot<-plot_cor_heatmap(df_heatmap)
    plot <- (plot
             + scale_fill_viridis(name="N")
             + ggtitle("Tanner stage")
             + xlab("2nd wave")
             + ylab("1st wave")
             + theme(plot.title = element_text(hjust = 0.5),
                     axis.text.x = element_text(size=8,angle = 0,vjust=0,hjust=0.5),
                     axis.text.y = element_text(size=8)))
    ggsave(paste("sex-",id_sex,"_tanner_heatmap.eps",sep=""),plot=plot,device=cairo_ps,
           path=file.path(paths_$output,"output"),dpi=300,height=5,width=5,limitsize=F)
    plot<-list(plot)
    list_plot<-c(list_plot,plot)
    names(plot)<-paste("sex-",id_sex,"_tanner_heatmap",sep="")
  }
  
  print("Finished plot_clin().")
  return(list_plot)
}


#**************************************************
# GAMM ============================================
#**************************************************
#mod_clinical<-gam(tanner_max ~ s(age) + s(ID_pnTTC,bs='re'),data=df_clinical_rbind)
#summary(mod_clinical)
#plot(mod_clinical,select=2)
#plot_smooth(mod_clinical,view='age')