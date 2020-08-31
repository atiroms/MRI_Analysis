#**************************************************
# Description =====================================
#**************************************************
# R script to analyze fingerprint data.


#**************************************************
# Parameters ======================================
#**************************************************
# parameters for gta_bin() and gta_weight()
path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"
#path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/str_FS"
#path_exp <- "Dropbox/MRI/pnTTC/Puberty/Stats/func_XCP/test_5sub"
path_exp_full <-NULL

dir_in<-"422_fp_aroma"
dir_out<-"428_fp_glm_aroma"
source(file.path(getwd(),"util/parameter.R"))


# Parameters for all functions
list_wave <- c(1,2)
#list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400","shen268")
#list_atlas<-c("aal116","gordon333","power264","shen268")
#list_atlas<-"aal116"
list_atlas<-"power264"

## parameters for variance_fp()
#list_covar_variance<-list("tanner"=list("1"="W1_Tanner_Max","2"="W2_Tanner_Max","label"="Tanner stage (max)"),
#                          #"age"=list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
#                          "sex"=list("1"="Sex","2"="Sex","label"="Sex"))

#thr_pvalue <- 0.05
#n_permutation<-1000
#n_permutation<-100


#**************************************************
# Libraries =======================================
#**************************************************
library(dplyr)
library(mgcv)
library(rowr)
library(multcomp)
library(parallel)
library(DescTools)
library(ggplot2)
library(data.table)
#library(GGally)
#library(igraph)
#library(qgraph)


#**************************************************
# Original library ================================
#**************************************************
source(file.path(getwd(),"util/function.R"))
source(file.path(getwd(),"util/plot.R"))
paths<-func_path(path_exp_=path_exp,dir_in_=dir_in,dir_out_=dir_out,path_exp_full_=path_exp_full)



#**************************************************
# Variance attribution ============================
#**************************************************
variance_fp<-function(paths_=paths,
                      list_atlas_=list_atlas,list_wave_=list_wave,
                      list_covar_variance_=list_covar_variance,
                      subset_subj_=subset_subj
                      ){
  print("Starting variance_fp().")
  nullobj<-func_createdirs(paths_,str_proc="variance_fp()")
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  #list_covar_variance<-list("age"=list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
  #                          "sex"=list("1"="Sex","2"="Sex","label"="Sex"))
  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,
                                     list_covar=list_covar_variance_,rem_na_clin=T)
  df_clin<-data_clin$df_clin
  colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
  
  #df_out_lm<-df_out_aic<-NULL
  #list_src_ancova<-NULL
  
  for (atlas in list_atlas_){
    # Load fingerprint data
    df_fp<-read.csv(file.path(paths_$input,"output",paste("atl-",atlas,"_fp.csv",sep="")))
    
    list_measure<-sort(unique(as.character(df_fp$measure)))
    #df_join<-NULL
    for (measure in list_measure){
      print(paste("Atlas: ",atlas,", Measure: ",measure,sep=""))
      df_fp_meas<-df_fp[df_fp$measure==measure,]
      
      # Create list of subjects who meet subsetting condition and whose MRI data exist
      list_id_subj_exist<-list()
      for (ses in list_wave_){
        id_subj_exist_ses<-sort(unique(c(df_fp_meas[df_fp_meas$from_ses==ses,'from_ID_pnTTC'],
                                         df_fp_meas[df_fp_meas$to_ses==ses,'to_ID_pnTTC'])))
        id_subj_subset_ses<-df_clin[df_clin$ses==ses,"ID_pnTTC"]
        id_subj_exist_ses<-intersect(id_subj_exist_ses,id_subj_subset_ses)
        list_id_subj_exist[[as.character(ses)]]<-sort(id_subj_exist_ses)
      }
      
      # Identify subjects with longitudinal data
      #list_id_subj_exist_twice<-sort(intersect(list_id_subj_exist[["1"]],
      #                                         list_id_subj_exist[["2"]]))
      #n_id_subj_exist_twice<-length(list_id_subj_exist_twice)
      #df_fp_meas<-df_fp_meas[df_fp_meas$from_ID_pnTTC %in% list_id_subj_exist_twice
      #                       & df_fp_meas$to_ID_pnTTC %in% list_id_subj_exist_twice,]
      
      # Create list of groups
      list_group<-sort(unique(c(as.character(df_fp_meas$group_1),
                                as.character(df_fp_meas$group_2))))
      if ("whole" %in% list_group){
        list_group<-c("whole",list_group[list_group!="whole"])
      }
      n_group<-length(list_group)
      
      df_zr<-data.frame()
      df_stat_zr<-data.frame()
      list_groups<-NULL
      for (idx_group_1 in seq(n_group)){
        for (idx_group_2 in seq(idx_group_1,n_group)){
          group_1<-list_group[idx_group_1]
          group_2<-list_group[idx_group_2]
          list_groups<-c(list_groups,paste(group_1,group_2,sep="-"))
          
          # Prepare dataframe for variance calculation
          df_fp_group<-df_fp_meas[df_fp_meas$group_1==group_1 & df_fp_meas$group_2==group_2,]
          df_fp_group<-inner_join(df_fp_group,df_clin,by=c("from_ID_pnTTC"="ID_pnTTC","from_ses"="ses"))
          colnames(df_fp_group)[colnames(df_fp_group)=="sex"]<-"from_sex"
          colnames(df_fp_group)[colnames(df_fp_group)=="tanner"]<-"from_tanner"
          df_fp_group<-inner_join(df_fp_group,df_clin,by=c("to_ID_pnTTC"="ID_pnTTC","to_ses"="ses"))
          colnames(df_fp_group)[colnames(df_fp_group)=="sex"]<-"to_sex"
          colnames(df_fp_group)[colnames(df_fp_group)=="tanner"]<-"to_tanner"
          
          df_zr_group<-data.frame()
          df_stat_zr_group<-data.frame()
          
          # M+F Group similarity
          list_zr<-df_fp_group$z_r
          df_zr_group<-rbind(df_zr_group,data.frame(strat="group",z_r=list_zr))
          df_stat_zr_group<-rbind(df_stat_zr_group,data.frame(strat="group",mean_z_r=mean(list_zr),
                                  n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
          # M+F Group + sex similarity
          list_zr<-df_fp_group[df_fp_group$from_sex==df_fp_group$to_sex,"z_r"]
          df_zr_group<-rbind(df_zr_group,data.frame(strat="group+sex",z_r=list_zr))
          df_stat_zr_group<-rbind(df_stat_zr_group,data.frame(strat="group+sex",mean_z_r=mean(list_zr),
                                  n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
          # M+F Group + wave similarity
          list_zr<-df_fp_group[df_fp_group$from_ses==df_fp_group$to_ses,"z_r"]
          df_zr_group<-rbind(df_zr_group,data.frame(strat="group+wave",z_r=list_zr))
          df_stat_zr_group<-rbind(df_stat_zr_group,data.frame(strat="group+wave",mean_z_r=mean(list_zr),
                                                              n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
          # M+F Group (+ sex) + individual similarity
          list_zr<-df_fp_group[df_fp_group$from_ID_pnTTC==df_fp_group$to_ID_pnTTC,"z_r"]
          df_zr_group<-rbind(df_zr_group,data.frame(strat="group+individual",z_r=list_zr))
          df_stat_zr_group<-rbind(df_stat_zr_group,data.frame(strat="group+individual",mean_z_r=mean(list_zr),
                                  n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
          # M+F Calculate relative mean_z_r
          max_mean<-max(df_stat_zr_group$mean_z_r)
          df_stat_zr_group$rel_mean_z_r<-df_stat_zr_group$mean_z_r/max_mean
          df_stat_zr_group<-rbind(df_stat_zr_group,data.frame(strat=c("sex","wave","individual"),mean_z_r=NA,n_z_r=NA,sd_z_r=NA,sem_z_r=NA,
                                  rel_mean_z_r=c(df_stat_zr_group[df_stat_zr_group$strat=="group+sex","rel_mean_z_r"]-df_stat_zr_group[df_stat_zr_group$strat=="group","rel_mean_z_r"],
                                                 df_stat_zr_group[df_stat_zr_group$strat=="group+wave","rel_mean_z_r"]-df_stat_zr_group[df_stat_zr_group$strat=="group","rel_mean_z_r"],
                                                 df_stat_zr_group[df_stat_zr_group$strat=="group+individual","rel_mean_z_r"]-df_stat_zr_group[df_stat_zr_group$strat=="group+sex","rel_mean_z_r"])))
          
          # M prepare dataframe
          df_fp_male<-df_fp_group[df_fp_group$from_sex==1 & df_fp_group$to_sex==1,]
          df_zr_male<-data.frame()
          df_stat_zr_male<-data.frame()
          # M Group similarity
          list_zr<-df_fp_male$z_r
          df_zr_male<-rbind(df_zr_male,data.frame(strat="male_group",z_r=list_zr))
          df_stat_zr_male<-rbind(df_stat_zr_male,data.frame(strat="male_group",mean_z_r=mean(list_zr),
                                                            n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
          # M Group + wave similarity
          list_zr<-df_fp_male[df_fp_male$from_ses==df_fp_male$to_ses,"z_r"]
          df_zr_male<-rbind(df_zr_male,data.frame(strat="male_group+wave",z_r=list_zr))
          df_stat_zr_male<-rbind(df_stat_zr_male,data.frame(strat="male_group+wave",mean_z_r=mean(list_zr),
                                                              n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
          # M Group + individual similarity
          list_zr<-df_fp_male[df_fp_male$from_ID_pnTTC==df_fp_male$to_ID_pnTTC,"z_r"]
          df_zr_male<-rbind(df_zr_male,data.frame(strat="male_group+individual",z_r=list_zr))
          df_stat_zr_male<-rbind(df_stat_zr_male,data.frame(strat="male_group+individual",mean_z_r=mean(list_zr),
                                                            n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
          # M Group + Tanner similarity
          list_zr<-df_fp_male[df_fp_male$from_tanner==df_fp_male$to_tanner,"z_r"]
          df_zr_male<-rbind(df_zr_male,data.frame(strat="male_group+tanner",z_r=list_zr))
          df_stat_zr_male<-rbind(df_stat_zr_male,data.frame(strat="male_group+tanner",mean_z_r=mean(list_zr),
                                                            n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
          # M Group + individual + Tanner similarity
          list_zr<-df_fp_male[df_fp_male$from_ID_pnTTC==df_fp_male$to_ID_pnTTC 
                              & df_fp_male$from_tanner==df_fp_male$to_tanner,"z_r"]
          df_zr_male<-rbind(df_zr_male,data.frame(strat="male_group+individual+tanner",z_r=list_zr))
          df_stat_zr_male<-rbind(df_stat_zr_male,data.frame(strat="male_group+individual+tanner",mean_z_r=mean(list_zr),
                                                            n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
          # M Calculate relative mean_z_r
          #max_mean<-max(df_stat_zr_male$mean_z_r)
          #df_stat_zr_male$rel_mean_z_r<-df_stat_zr_male$mean_z_r/max_mean
          df_stat_zr_male$rel_mean_z_r<-df_stat_zr_male$mean_z_r/df_stat_zr_male[df_stat_zr_male$strat=="male_group+individual+tanner","mean_z_r"]
          df_stat_zr_male<-rbind(df_stat_zr_male,data.frame(strat=c("male_individual","male_wave","male_tanner","male_individual+tanner"),mean_z_r=NA,n_z_r=NA,sd_z_r=NA,sem_z_r=NA,
                                                            rel_mean_z_r=c(df_stat_zr_male[df_stat_zr_male$strat=="male_group+individual","rel_mean_z_r"]-df_stat_zr_male[df_stat_zr_male$strat=="male_group","rel_mean_z_r"],
                                                                           df_stat_zr_male[df_stat_zr_male$strat=="male_group+wave","rel_mean_z_r"]-df_stat_zr_male[df_stat_zr_male$strat=="male_group","rel_mean_z_r"],
                                                                           df_stat_zr_male[df_stat_zr_male$strat=="male_group+tanner","rel_mean_z_r"]-df_stat_zr_male[df_stat_zr_male$strat=="male_group","rel_mean_z_r"],
                                                                           df_stat_zr_male[df_stat_zr_male$strat=="male_group+individual+tanner","rel_mean_z_r"]-df_stat_zr_male[df_stat_zr_male$strat=="male_group+tanner","rel_mean_z_r"])))
          
          # F prepare dataframe
          df_fp_female<-df_fp_group[df_fp_group$from_sex==2 & df_fp_group$to_sex==2,]
          df_zr_female<-data.frame()
          df_stat_zr_female<-data.frame()
          # F Group similarity
          list_zr<-df_fp_female$z_r
          df_zr_female<-rbind(df_zr_female,data.frame(strat="female_group",z_r=list_zr))
          df_stat_zr_female<-rbind(df_stat_zr_female,data.frame(strat="female_group",mean_z_r=mean(list_zr),
                                                            n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
          # F Group + wave similarity
          list_zr<-df_fp_female[df_fp_female$from_ses==df_fp_female$to_ses,"z_r"]
          df_zr_female<-rbind(df_zr_female,data.frame(strat="female_group+wave",z_r=list_zr))
          df_stat_zr_female<-rbind(df_stat_zr_female,data.frame(strat="female_group+wave",mean_z_r=mean(list_zr),
                                                            n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
          # F Group + individual similarity
          list_zr<-df_fp_female[df_fp_female$from_ID_pnTTC==df_fp_female$to_ID_pnTTC,"z_r"]
          df_zr_female<-rbind(df_zr_female,data.frame(strat="female_group+individual",z_r=list_zr))
          df_stat_zr_female<-rbind(df_stat_zr_female,data.frame(strat="female_group+individual",mean_z_r=mean(list_zr),
                                                            n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
          # F Group + Tanner similarity
          list_zr<-df_fp_female[df_fp_female$from_tanner==df_fp_female$to_tanner,"z_r"]
          df_zr_female<-rbind(df_zr_female,data.frame(strat="female_group+tanner",z_r=list_zr))
          df_stat_zr_female<-rbind(df_stat_zr_female,data.frame(strat="female_group+tanner",mean_z_r=mean(list_zr),
                                                            n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
          # F Group + individual + Tanner similarity
          list_zr<-df_fp_female[df_fp_female$from_ID_pnTTC==df_fp_female$to_ID_pnTTC 
                              & df_fp_female$from_tanner==df_fp_female$to_tanner,"z_r"]
          df_zr_female<-rbind(df_zr_female,data.frame(strat="female_group+individual+tanner",z_r=list_zr))
          df_stat_zr_female<-rbind(df_stat_zr_female,data.frame(strat="female_group+individual+tanner",mean_z_r=mean(list_zr),
                                                            n_z_r=length(list_zr),sd_z_r=sd(list_zr),sem_z_r=sd(list_zr)/sqrt(length(list_zr))))
          # F Calculate relative mean_z_r
          #max_mean<-max(df_stat_zr_female$mean_z_r)
          #df_stat_zr_female$rel_mean_z_r<-df_stat_zr_female$mean_z_r/max_mean
          df_stat_zr_female$rel_mean_z_r<-df_stat_zr_female$mean_z_r/df_stat_zr_female[df_stat_zr_female$strat=="female_group+individual+tanner","mean_z_r"]
          df_stat_zr_female<-rbind(df_stat_zr_female,data.frame(strat=c("female_individual","female_wave","female_tanner","female_individual+tanner"),mean_z_r=NA,n_z_r=NA,sd_z_r=NA,sem_z_r=NA,
                                                            rel_mean_z_r=c(df_stat_zr_female[df_stat_zr_female$strat=="female_group+individual","rel_mean_z_r"]-df_stat_zr_female[df_stat_zr_female$strat=="female_group","rel_mean_z_r"],
                                                                           df_stat_zr_female[df_stat_zr_female$strat=="female_group+wave","rel_mean_z_r"]-df_stat_zr_female[df_stat_zr_female$strat=="female_group","rel_mean_z_r"],
                                                                           df_stat_zr_female[df_stat_zr_female$strat=="female_group+tanner","rel_mean_z_r"]-df_stat_zr_female[df_stat_zr_female$strat=="female_group","rel_mean_z_r"],
                                                                           df_stat_zr_female[df_stat_zr_female$strat=="female_group+individual+tanner","rel_mean_z_r"]-df_stat_zr_female[df_stat_zr_female$strat=="female_group+tanner","rel_mean_z_r"])))
          
          # Concatenate dataframe
          df_zr_group<-data.frame(atlas=atlas,measure=measure,group_1=group_1,group_2=group_2,
                                  rbind(df_zr_group,df_zr_male,df_zr_female))
          df_zr<-rbind(df_zr,df_zr_group)
          df_stat_zr_group<-data.frame(atlas=atlas,measure=measure,group_1=group_1,group_2=group_2,
                                       rbind(df_stat_zr_group,df_stat_zr_male,df_stat_zr_female))
          df_stat_zr<-rbind(df_stat_zr,df_stat_zr_group)
        }
      } # End of loop over groups
      
      # Save results
      write.csv(df_zr,file.path(paths_$output,"output",paste("atl-",atlas,"_msr-",measure,"_zr.csv",sep="")),row.names=F)
      write.csv(df_stat_zr,file.path(paths_$output,"output",paste("atl-",atlas,"_msr-",measure,"_stat_zr.csv",sep="")),row.names=F)
      
      # Graphical output
      df_stat_zr$groups=paste(as.character(df_stat_zr$group_1),as.character(df_stat_zr$group_2),sep="-")
      label_groups<-NULL
      for (groups in list_groups){
        label_groups<-c(label_groups,paste(capitalize(substring(groups,1,gregexpr("-",groups)[[1]][1]-1)),
                              capitalize(substring(groups,gregexpr("-",groups)[[1]][1]+1)),sep=" - "))
      }
      # M+F absolute
      df_stat_zr_plot<-df_stat_zr[df_stat_zr$strat %in% c("group","group+sex","group+wave","group+individual"),]
      plot<-(ggplot(data=df_stat_zr_plot,aes(x=groups,y=mean_z_r,fill=strat))
             + geom_bar(stat="identity",color="white",width=0.7,position=position_dodge())
             + geom_errorbar(aes(ymin=mean_z_r-sd_z_r,ymax=mean_z_r+sd_z_r),width=.1,position=position_dodge(0.7))
             + scale_x_discrete(limits = list_groups,labels = label_groups)
             + scale_y_continuous(expand=c(0,0),limits=c(0,max(df_stat_zr_plot$mean_z_r)+max(df_stat_zr_plot$sd_z_r)+0.1))
             + scale_fill_brewer(palette="Greys",direction=-1,name="Stratification")
             + ggtitle(paste("Fingerprint similarity attribution, ",atlas,sep=""))
             + xlab("Subnetworks") + ylab("Mean z(r)") + theme_classic()
             + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,vjust=1,hjust=1),
                     legend.position="top",legend.justification="left",legend.direction="vertical")
             )
      ggsave(paste("atl-",atlas,"_msr-",measure,"_sex-mf_fp_similarity.eps",sep=""),plot=plot,device=cairo_ps,
             path=file.path(paths_$output,"output"),height=7,width=14,limitsize=F)
      # M+F relative
      df_stat_zr_plot<-df_stat_zr[df_stat_zr$strat %in% c("group","sex","individual"),]
      plot<-(ggplot(data=df_stat_zr_plot,aes(x=groups,y=rel_mean_z_r,fill=factor(strat,levels=c("group","sex","individual"))))
             + geom_bar(stat="identity",color="white",width=0.5,position=position_stack(vjust=0,reverse=T))
             + scale_x_discrete(limits = list_groups,labels = label_groups)
             + scale_y_continuous(expand=c(0,0),limits=c(0,1.05))
             + scale_fill_brewer(palette="Greys",direction=-1,name="Stratification")
             + ggtitle(paste("Fingerprint relative similarity attribution, ",atlas,sep=""))
             + xlab("Subnetworks") + ylab("Relative mean z(r)") + theme_classic()
             + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,vjust=1,hjust=1),
                     legend.position="top",legend.justification="left",legend.direction="vertical")
      )
      ggsave(paste("atl-",atlas,"_msr-",measure,"_sex-mf_fp_similarity_rel.eps",sep=""),plot=plot,device=cairo_ps,
             path=file.path(paths_$output,"output"),height=7,width=14,limitsize=F)
      # M absolute
      df_stat_zr_plot<-df_stat_zr[df_stat_zr$strat %in% c("male_group","male_group+wave","male_group+individual",
                                                          "male_group+tanner","male_group+individual+tanner"),]
      plot<-(ggplot(data=df_stat_zr_plot,aes(x=groups,y=mean_z_r,fill=strat))
             + geom_bar(stat="identity",color="white",width=0.7,position=position_dodge())
             + geom_errorbar(aes(ymin=mean_z_r-sd_z_r,ymax=mean_z_r+sd_z_r),width=.1,position=position_dodge(0.7))
             + scale_x_discrete(limits = list_groups,labels = label_groups)
             + scale_y_continuous(expand=c(0,0),limits=c(0,max(df_stat_zr_plot$mean_z_r)+max(df_stat_zr_plot$sd_z_r)+0.1))
             + scale_fill_brewer(palette="Blues",direction=-1,name="Stratification")
             + ggtitle(paste("Fingerprint similarity attribution, ",atlas,sep=""))
             + xlab("Subnetworks") + ylab("Mean z(r)") + theme_classic()
             + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,vjust=1,hjust=1),
                     legend.position="top",legend.justification="left",legend.direction="vertical")
      )
      ggsave(paste("atl-",atlas,"_msr-",measure,"_sex-m_fp_similarity.eps",sep=""),plot=plot,device=cairo_ps,
             path=file.path(paths_$output,"output"),height=7,width=14,limitsize=F)
      # M relative
      df_stat_zr_plot<-df_stat_zr[df_stat_zr$strat %in% c("male_group","male_tanner","male_individual+tanner"),]
      plot<-(ggplot(data=df_stat_zr_plot,aes(x=groups,y=rel_mean_z_r,fill=factor(strat,levels=c("male_group","male_tanner","male_individual+tanner"))))
             + geom_bar(stat="identity",color="white",width=0.5,position=position_stack(vjust=0,reverse=T))
             + scale_x_discrete(limits = list_groups,labels = label_groups)
             + scale_y_continuous(expand=c(0,0),limits=c(0,1.05))
             + scale_fill_brewer(palette="Blues",direction=-1,name="Stratification")
             + ggtitle(paste("Fingerprint relative similarity attribution, ",atlas,sep=""))
             + xlab("Subnetworks") + ylab("Relative mean z(r)") + theme_classic()
             + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,vjust=1,hjust=1),
                     legend.position="top",legend.justification="left",legend.direction="vertical")
      )
      ggsave(paste("atl-",atlas,"_msr-",measure,"_sex-m_fp_similarity_rel.eps",sep=""),plot=plot,device=cairo_ps,
             path=file.path(paths_$output,"output"),height=7,width=14,limitsize=F)
      # F absolute
      df_stat_zr_plot<-df_stat_zr[df_stat_zr$strat %in% c("female_group","female_group+wave","female_group+individual",
                                                          "female_group+tanner","female_group+individual+tanner"),]
      plot<-(ggplot(data=df_stat_zr_plot,aes(x=groups,y=mean_z_r,fill=strat))
             + geom_bar(stat="identity",color="white",width=0.7,position=position_dodge())
             + geom_errorbar(aes(ymin=mean_z_r-sd_z_r,ymax=mean_z_r+sd_z_r),width=.1,position=position_dodge(0.7))
             + scale_x_discrete(limits = list_groups,labels = label_groups)
             + scale_y_continuous(expand=c(0,0),limits=c(0,max(df_stat_zr_plot$mean_z_r)+max(df_stat_zr_plot$sd_z_r)+0.1))
             + scale_fill_brewer(palette="Reds",direction=-1,name="Stratification")
             + ggtitle(paste("Fingerprint similarity attribution, ",atlas,sep=""))
             + xlab("Subnetworks") + ylab("Mean z(r)") + theme_classic()
             + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,vjust=1,hjust=1),
                     legend.position="top",legend.justification="left",legend.direction="vertical")
      )
      ggsave(paste("atl-",atlas,"_msr-",measure,"_sex-f_fp_similarity.eps",sep=""),plot=plot,device=cairo_ps,
             path=file.path(paths_$output,"output"),height=7,width=14,limitsize=F)
      # F relative
      df_stat_zr_plot<-df_stat_zr[df_stat_zr$strat %in% c("female_group","female_tanner","female_individual+tanner"),]
      plot<-(ggplot(data=df_stat_zr_plot,aes(x=groups,y=rel_mean_z_r,fill=factor(strat,levels=c("female_group","female_tanner","female_individual+tanner"))))
             + geom_bar(stat="identity",color="white",width=0.5,position=position_stack(vjust=0,reverse=T))
             + scale_x_discrete(limits = list_groups,labels = label_groups)
             + scale_y_continuous(expand=c(0,0),limits=c(0,1.05))
             + scale_fill_brewer(palette="Reds",direction=-1,name="Stratification")
             + ggtitle(paste("Fingerprint relative similarity attribution, ",atlas,sep=""))
             + xlab("Subnetworks") + ylab("Relative mean z(r)") + theme_classic()
             + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,vjust=1,hjust=1),
                     legend.position="top",legend.justification="left",legend.direction="vertical")
      )
      ggsave(paste("atl-",atlas,"_msr-",measure,"_sex-f_fp_similarity_rel.eps",sep=""),plot=plot,device=cairo_ps,
             path=file.path(paths_$output,"output"),height=7,width=14,limitsize=F)
      
    }
  }
  print("Finished variance_fp()")
}


#**************************************************
# GLM and ANCOVA of Fingerprint change ============
#**************************************************

glm_core<-function(data_src){
  df_src<-data_src$df_src
  atlas<-data_src$atlas
  measure<-data_src$measure
  group_1<-data_src$group_1
  group_2<-data_src$group_2
  
  #print(paste("Atlas: ",atlas,", Measure: ",measure,", Group: ",group_1," and ",group_2,", GLM/GAM.",  sep=""))
  df_out_aic_add<-df_out_lm_add<-data.frame()
  for (idx_mod in names(list_mod_)){
    list_plot<-list()
    list_sex<-sort(unique(as.numeric.factor(df_src$sex)))
    for (idx_sex in list_sex){
      df_src_sex<-df_src[df_src$sex==idx_sex,]
      mod<-gam(as.formula(list_mod_[[idx_mod]]),data=df_src_sex,method="REML")
      p_table<-summary.gam(mod)$p.table
      if (is.null(summary.gam(mod)$s.table)){
        df_out_lm_add_add<-data.frame(sex=idx_sex,model=idx_mod,term=rownames(p_table),
                                      estimate=p_table[,'Estimate'],se=p_table[,'Std. Error'],F=NA,
                                      t=p_table[,'t value'],p=p_table[,'Pr(>|t|)'])
        
      }else{
        s_table<-summary.gam(mod)$s.table
        df_out_lm_add_add<-rbind(data.frame(sex=idx_sex,model=idx_mod,term=rownames(p_table),
                                            estimate=p_table[,'Estimate'],se=p_table[,'Std. Error'],F=NA,
                                            t=p_table[,'t value'],p=p_table[,'Pr(>|t|)']),
                                 data.frame(sex=idx_sex,model=idx_mod,term=rownames(s_table),
                                            estimate=NA,se=NA,F=s_table[,'F'],
                                            t=NA,p=s_table[,'p-value']))
      }
      df_out_lm_add<-rbind(df_out_lm_add,df_out_lm_add_add)
      df_out_aic_add<-rbind(df_out_aic_add,
                            data.frame(sex=idx_sex,model=idx_mod,aic=mod$aic,aic_best_among_models=0))
      
      # Graphical output of GLM results
      for (idx_graph in names(list_graph_)){
        if (list_graph_[[idx_graph]][["x_axis"]] %in% colnames(mod$model)){
          # Add sex-wise lines/plots to existent plot, initialize if absent
          plot<-plot_gamm(plot_in=list_plot[[idx_graph]],mod_gamm=mod,
                          df_in=df_src_sex,
                          spec_graph=list_graph_[[idx_graph]])
          list_plot[[idx_graph]]<-plot
          
          # Output
          if (idx_sex==list_sex[length(list_sex)]){
            # Prepare x axis label
            axis_x<-list_graph_[[idx_graph]][["x_axis"]]
            for (idx_prefix in list(c("",""),c("ses1_"," 1st wave"),c("ses2_"," 2nd wave"),
                                    c("diff_"," difference"),c("mean_"," mean"))){
              for (idx_covar in names(list_covar_)){
                if (axis_x==paste(idx_prefix[1],idx_covar,sep="")){
                  label_x<-paste(list_covar_[[idx_covar]][["label"]],idx_prefix[2],sep='')
                }
              }
            }
            
            plot<-(plot
                   + ggtitle(paste(list_graph_[[idx_graph]][["title"]],atlas,measure,group_1,group_2,idx_mod,sep=" "))
                   + xlab(label_x)
                   + ylab("Fingerprint correlation")
                   + theme(legend.position = "none"))
            filename_plot<-paste("atl-",atlas,"_msr-",measure,"_grp1-",group_1,"_grp2-",group_2,"_mod-",idx_mod,
                                 "_plt-",idx_graph,"_fp_glm.eps",sep="")
            ggsave(filename_plot,plot=plot,device=cairo_ps,
                   path=file.path(paths_$output,"output","plot"),dpi=300,height=5,width=5,limitsize=F)
          }
        }
      }
    }
  }
  
  # Compare AICs of GLM models
  df_out_aic_add_sex_rbind<-data.frame()
  for (idx_sex in list_sex){
    df_out_aic_add_sex<-df_out_aic_add[df_out_aic_add$sex==idx_sex,]
    df_out_aic_add_sex[which(df_out_aic_add_sex$aic==min(df_out_aic_add_sex$aic)),
                       'aic_best_among_models']<-1
    df_out_aic_add_sex_rbind<-rbind(df_out_aic_add_sex_rbind,df_out_aic_add_sex)
  }
  
  df_out_lm_add<-cbind(atlas=atlas,measure=measure,group_1=group_1,group_2=group_2,df_out_lm_add)
  df_out_aic_add_sex_rbind<-cbind(atlas=atlas,measure=measure,group_1=group_1,group_2=group_2,
                                  df_out_aic_add_sex_rbind)
  
  return(list("df_out_lm_add"=df_out_lm_add,"df_out_aic_add"=df_out_aic_add_sex_rbind))
}

ancova_core<-function(data_input){
  atlas<-data_input$atlas
  measure<-data_input$measure
  group_1=data_input$group_network[1]
  group_2=data_input$group_network[2]
  group_tanner<-data_input$group_tanner
  group_tanner_content<-data_input$group_tanner_content
  idx_sex<-data_input$group_sex
  df_src_ancova<-data_input$df_src_ancova
  
  # Calculate ANCOVA
  if (idx_sex=="all"){
    mod_ancova<-aov(value~long_tanner+diff_age+sex,data=df_src_ancova)
    color_plot<-"seagreen"
  }else if (idx_sex=="male"){
    mod_ancova<-aov(value~long_tanner+diff_age,data=df_src_ancova)
    color_plot<-"steelblue2"
  }else if (idx_sex=="female"){
    mod_ancova<-aov(value~long_tanner+diff_age,data=df_src_ancova)
    color_plot<-"lightcoral"
  }
  df_ancova<-summary(mod_ancova)[[1]]
  df_out_ancova<-data.frame(atlas=atlas,measure=measure,group_1=group_1,group_2=group_2,tanner=group_tanner,sex=idx_sex,test="ANCOVA",
                            term=rownames(df_ancova),label=NA,
                            p=df_ancova[,'Pr(>F)'],t=NA,F=df_ancova[,'F value'],
                            fit=NA,diff=NA,sigma=NA)
  
  # Extract fitting of ANCOVA
  fit_ancova<-mod_ancova$coefficients
  list_term<-names(mod_ancova$model)
  list_term<-c("(Intercept)",list_term[list_term!="value"])
  #list_term<-list_term[list_term!="value"]
  df_out_fit<-NULL
  for (term in list_term){
    if (term %in% names(mod_ancova$xlevels)){
      list_label<-mod_ancova$xlevels[[term]]
    }else{
      list_label<-NA
    }
    
    for (label in list_label){
      if (is.na(label)){
        term_label<-term
      }else{
        term_label<-paste(term,label,sep='')
      }
      if (term_label %in% names(mod_ancova$coefficients)){
        fit<-mod_ancova$coefficients[[term_label]]
      }else{
        fit<-0
      }
      df_out_fit<-rbind(df_out_fit,
                        data.frame(atlas=atlas,measure=measure,group_1=group_1,group_2=group_2,tanner=group_tanner,sex=idx_sex,test="fit",
                                   term=term,label=label,
                                   p=NA,t=NA,F=NA,
                                   fit=fit,diff=NA,sigma=NA))
    }
  }
  
  # Graphical output of ANCOVA tanner result
  mean_fit<-df_out_fit[df_out_fit$term=="(Intercept)","fit"]
  for (term in list_term[list_term!="(Intercept)" & list_term!="long_tanner"]){
    list_label<-df_out_fit[df_out_fit$term==term,"label"]
    if (is.na(list_label[1])){
      mean_fit<-mean_fit+mean(df_src_ancova[,term]*df_out_fit[df_out_fit$term==term,"fit"])
    }else{
      list_label<-sort(unique(list_label))
      df_calc_mean<-data.frame(df_out_fit[df_out_fit$term==term & df_out_fit$label %in% list_label,c("label","fit")])
      df_src_ancova_mean<-df_src_ancova
      colnames(df_src_ancova_mean)[colnames(df_src_ancova_mean)==term]<-"label"
      df_src_ancova_mean$label<-as.character(df_src_ancova_mean$label)
      df_calc_mean<-left_join(df_src_ancova_mean,df_calc_mean,by="label")
      mean_fit<-mean_fit+mean(df_calc_mean$fit)
    }
  }
  
  df_ancova_plot<-data.frame(matrix(nrow=length(group_tanner_content[["1"]]),
                                    ncol=length(group_tanner_content[["2"]])))
  rownames(df_ancova_plot)<-names(group_tanner_content[["1"]])
  colnames(df_ancova_plot)<-names(group_tanner_content[["2"]])
  for(group_tanner_1 in names(group_tanner_content[["1"]])){
    for(group_tanner_2 in names(group_tanner_content[["2"]])){
      mean_fit_tanner<-df_out_fit[df_out_fit$term=="long_tanner" & df_out_fit$label==paste(group_tanner_1,group_tanner_2,sep="_"),"fit"]
      if (length(mean_fit_tanner)>0){
        df_ancova_plot[group_tanner_1,group_tanner_2]<-mean_fit+mean_fit_tanner
      }else{
        df_ancova_plot[group_tanner_1,group_tanner_2]<-NA
      }
    }
  }
  plot_ancova<-plot_cor_heatmap(input=df_ancova_plot)
  suppressMessages(plot_ancova<-(plot_ancova
                                 + scale_fill_gradient(low="white",high=color_plot,name="r")
                                 + ggtitle(paste("FP Cor ANCOVA,",atlas,measure,group_1,group_2,group_tanner,idx_sex,sep=" "))
                                 + xlab("2nd wave Tanner stage")
                                 + ylab("1st wave Tanner stage")
                                 + theme(plot.title = element_text(hjust = 0.5),
                                         axis.text.x = element_text(size=12,angle = 0,vjust=0,hjust=0.5),
                                         axis.text.y = element_text(size=12))))
  
  ggsave(paste("atl-",atlas,"_msr-",measure,"_grp1-",group_1,"_grp2-",group_2,"_str-",group_tanner,"_sex-",idx_sex,"_fp_ancova.eps",sep=""),plot=plot_ancova,device=cairo_ps,
         path=file.path(paths_$output,"output","plot"),dpi=300,height=5,width=5,limitsize=F)
  
  # Calculate Tukey-Kramer
  tk<-summary(glht(mod_ancova, linfct = mcp('long_tanner' = 'Tukey')))$test
  df_out_posthoc<-data.frame(atlas=atlas,measure=measure,group_1=group_1,group_2=group_2,tanner=group_tanner,
                             sex=idx_sex,test="Tukey-Kramer",
                             term="long_tanner",label=names(tk$coefficients),
                             p=tk$pvalues[1:length(tk$coefficients)],t=tk$tstat,F=NA,
                             fit=NA,diff=tk$coefficients,sigma=tk$sigma)
  df_out<-rbind(df_out_ancova,df_out_fit,df_out_posthoc)
  return(df_out)
}

heatmap_gammfp<-function(paths_,df_out_lm,list_atlas,list_mod){
  
  for (atlas in list_atlas){
    df_out_lm_subset<-df_out_lm[df_out_lm$atlas==atlas,]
    list_measure<-as.character(sort(unique(df_out_lm_subset$measure)))
    for (measure in list_measure){
      df_out_lm_subset<-df_out_lm[df_out_lm$atlas==atlas & df_out_lm$measure==measure,]
      for (model in names(list_mod)){
        for (idx_sex in c(1,2)){
          df_out_lm_subset<-df_out_lm[df_out_lm$atlas==atlas
                                      & df_out_lm$measure==measure
                                      & df_out_lm$sex==idx_sex
                                      & df_out_lm$model==model,]
          list_term<-sort(unique(as.character(df_out_lm_subset$term)))
          list_term<-list_term[list_term!="(Intercept)"]
          for (term in list_term){
            df_out_lm_subset<-df_out_lm[df_out_lm$atlas==atlas
                                        & df_out_lm$measure==measure
                                        & df_out_lm$sex==idx_sex
                                        & df_out_lm$model==model
                                        & df_out_lm$term==term,]
            list_group<-sort(unique(c(as.character(df_out_lm_subset$group_1),
                                      as.character(df_out_lm_subset$group_2))))
            list_group<-list_group[list_group!="whole"]
            n_group<-length(list_group)
            mat_static<-data.frame(matrix(nrow=n_group,ncol=n_group))
            mat_pval<-data.frame(matrix(nrow=n_group,ncol=n_group))
            colnames(mat_static)<-rownames(mat_static)<-colnames(mat_pval)<-rownames(mat_pval)<-list_group
            for (idx_group_1 in seq(n_group)){
              for (idx_group_2 in seq(idx_group_1,n_group)){
                group_1<-list_group[idx_group_1]
                group_2<-list_group[idx_group_2]
                df_out_lm_subset<-df_out_lm[df_out_lm$atlas==atlas
                                            & df_out_lm$measure==measure
                                            & df_out_lm$sex==idx_sex
                                            & df_out_lm$model==model
                                            & df_out_lm$term==term
                                            & df_out_lm$group_1==group_1
                                            & df_out_lm$group_2==group_2,c("estimate","F","p")]
                name_static<-colnames(df_out_lm_subset[,c("estimate","F")])[!is.na(df_out_lm_subset[1,c("estimate","F")])]
                mat_static[group_1,group_2]<-mat_static[group_2,group_1]<-df_out_lm_subset[[1,name_static]]
                mat_pval[group_1,group_2]<-mat_pval[group_2,group_1]<-df_out_lm_subset[[1,"p"]]
              }
            }
            
            #df_pval<-rownames_to_column(mat_pval,"row")
            #df_pval<-gather(df_pval,key=column,value=p,2:ncol(df_pval))
            if (idx_sex==1){
              label_sex<-"male"
            }else{
              label_sex<-"female"
            }
            mat_pval<-round(mat_pval,3)
            plot_stat<-plot_cor_heatmap(mat_static,mat_pval)
            suppressMessages(plot_stat<-(plot_stat
                                         + scale_fill_gradientn(colors = matlab.like2(100),
                                                                lim=c(-max(max(mat_static),-min(mat_static)),max(max(mat_static),-min(mat_static))),
                                                                name=name_static)
                                         + ggtitle(paste("GLM-FP,",atlas,measure,label_sex,model,term,sep=" "))
                                         + theme(plot.title = element_text(hjust = 0.5),
                                                 axis.title=element_blank())))
            ggsave(paste("atl-",atlas,"_msr-",measure,"_sex-",label_sex,"_mod-",model,
                         "_plt-",term,"_fp_glm_heatmap.eps",sep=""),plot=plot_stat,device=cairo_ps,
                   path=file.path(paths_$output,"output","plot"),dpi=300,height=5,width=5,limitsize=F)
          }
        }
      }
    }
  }
  
}

model_fp_multi<-function(paths_=paths,subset_subj_=model_fp_subset_subj,list_atlas_=list_atlas,
                         list_wave_=list_wave,skip_ancova=T,
                         list_covar_tanner_=model_fp_list_covar_tanner,list_tanner_=model_fp_list_tanner,
                         list_mod_tanner_=model_fp_list_mod_tanner,list_graph_tanner_=model_fp_list_graph_tanner,
                         list_strat_tanner_=model_fp_list_strat_tanner,
                         list_covar_hormone_=model_fp_list_covar_hormone,list_hormone_=model_fp_list_hormone,
                         list_mod_hormone_=model_fp_list_mod_hormone,list_graph_hormone_=model_fp_list_graph_hormone
                         ){
  print("Starting model_fp_multi().")
  nullobj<-func_createdirs(paths_,str_proc="model_fp_multi()")
  
  df_out_lm<-df_out_aic<-NULL
  # Loop over clinical variables
  #1 Tanner stage
  for (idx_tanner in names(list_tanner_)){
    print(paste("Tanner type: ",list_tanner_[[idx_tanner]][["label"]],sep=""))
    list_covar<-list_covar_tanner_
    list_covar[["tanner"]]<-list_tanner_[[idx_tanner]]
    data_model<-model_fp(paths_,list_atlas_,list_wave_,list_covar_=list_covar,
                         list_mod_=list_mod_tanner_,list_graph_=list_graph_tanner_,
                         subset_subj_,list_strat_tanner_,skip_ancova,suffix=paste("var-",idx_tanner,sep=""))
    df_out_lm<-rbind(df_out_lm,cbind(variable=idx_tanner,data_model$df_out_lm))
    df_out_aic<-rbind(df_out_aic,cbind(variable=idx_tanner,data_model$df_out_aic))
  } # Finished looping over Tanner stages
  
  #2 Hormones
  for (idx_hormone in names(list_hormone_)){
    print(paste("Hormone type: ",list_hormone_[[idx_hormone]][["label"]],sep=""))
    list_covar<-list_covar_hormone_
    list_covar[["hormone"]]<-list_hormone_[[idx_hormone]]
    data_model<-model_fp(paths_,list_atlas_,list_wave_,list_covar_=list_covar,
                         list_mod_=list_mod_hormone_,list_graph_=list_graph_hormone_,
                         subset_subj_,NULL,skip_ancova,suffix=paste("var-",idx_hormone,sep=""))
    df_out_lm<-rbind(df_out_lm,cbind(variable=idx_hormone,data_model$df_out_lm))
    df_out_aic<-rbind(df_out_aic,cbind(variable=idx_hormone,data_model$df_out_aic))
  } # Finished looping over Hormones
  
  write.csv(df_out_lm,file.path(paths_$output,"output","fp_glm.csv"),row.names = F)
  write.csv(df_out_aic,file.path(paths_$output,"output","fp_glm_aic.csv"),row.names = F)
  print("Finished model_fp_multi().")
}

model_fp<-function(paths_=paths,
                   list_atlas_=list_atlas,list_wave_=list_wave,list_covar_=list_covar,
                   list_mod_=list_mod,list_graph_=list_graph,subset_subj_=subset_subj,
                   list_strat_tanner_=list_strat_tanner,
                   skip_ancova=T,suffix=""
                   ){

  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,
                                     list_covar=list_covar_,rem_na_clin=T,prefix=suffix,print_terminal=F)
  df_clin<-data_clin$df_clin
  colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
  
  list_src_glm<-list_src_ancova<-df_src_glm_bind<-NULL
  for (atlas in list_atlas_){
    # Load fingerprint data]
    print(paste("Loading FP, atlas:",atlas,sep=" "))
    df_fp<-as.data.frame(fread(file.path(paths_$input,"output",paste("atl-",atlas,"_fp.csv",sep=""))))
    
    list_measure<-as.character(sort(unique(df_fp$measure)))

    for (measure in list_measure){
      print(paste("Preparing dataset, atlas: ",atlas,", measure: ",measure,sep=""))
      df_fp_meas<-df_fp[df_fp$measure==measure,]
    
      # Create list of subjects who meet subsetting condition and whose MRI data exist
      list_ses_exist <- sort(unique(c(df_fp_meas$from_ses,df_fp_meas$to_ses)))
      list_id_subj_exist<-list()
      for (ses in list_ses_exist){
        id_subj_exist_ses<-sort(unique(c(df_fp_meas[df_fp_meas$from_ses==ses,'from_ID_pnTTC'],
                                         df_fp_meas[df_fp_meas$to_ses==ses,'to_ID_pnTTC'])))
        id_subj_subset_ses<-df_clin[df_clin$ses==ses,"ID_pnTTC"]
        id_subj_exist_ses<-intersect(id_subj_exist_ses,id_subj_subset_ses)
        list_id_subj_exist[[as.character(ses)]]<-sort(id_subj_exist_ses)
      }
      
      # Identify subjects with longitudinal data
      list_id_subj_exist_twice<-sort(intersect(list_id_subj_exist[["1"]],
                                               list_id_subj_exist[["2"]]))
      n_id_subj_exist_twice<-length(list_id_subj_exist_twice)
      
      # Create list of existing ROI subgroups
      list_group<-sort(unique(c(as.character(df_fp_meas$group_1),as.character(df_fp_meas$group_2))))
      if ("whole" %in% list_group){
        list_group<-c("whole",list_group[list_group!="whole"])
      }
      #list_group<-"whole      # Disable subgroup-wise analysis
      
      n_group<-length(list_group)
      for (idx_group_1 in seq(n_group)){
        group_1<-list_group[idx_group_1]
        df_fp_meas_grp1<-df_fp_meas[df_fp_meas$group_1==group_1,]
        for (idx_group_2 in seq(idx_group_1,n_group)){
          group_2<-list_group[idx_group_2]
          df_fp_meas_grp2<-df_fp_meas_grp1[df_fp_meas_grp1$group_2==group_2,]
          df_cor_fp<-data.frame(ID_pnTTC=list_id_subj_exist_twice)
          for (id_subj in list_id_subj_exist_twice){
            df_cor_fp[df_cor_fp$ID_pnTTC==id_subj,"value"]<-df_fp_meas_grp2[df_fp_meas_grp2$from_ID_pnTTC==id_subj
                                                                            & df_fp_meas_grp2$to_ID_pnTTC==id_subj,
                                                                            "z_r"]
          }
          
          # Subset those without longitudinal fp correlation
          list_id_subj_nonna<-df_cor_fp[!is.na(df_cor_fp$value),"ID_pnTTC"]
          df_cor_fp<-df_cor_fp[df_cor_fp$ID_pnTTC %in% list_id_subj_nonna,]
          n_id_subj_exist_twice<-length(list_id_subj_nonna)
          #print(paste("Atlas: ",atlas,", Measure: ",measure,", Group: ",group_1," and ",group_2,
          #            ", ",as.character(n_id_subj_exist_twice)," subjects with longitudinal non-NA data.",sep=""))
          
          # Create dataframe for GLM analysis
          df_src_glm<-func_clinical_data_join(df_src=df_clin,
                                               list_id_subj=list_id_subj_nonna,
                                               list_covar=list_covar_)
          df_src_glm<-inner_join(df_src_glm,df_cor_fp,by="ID_pnTTC")
          df_src_glm$ID_pnTTC<-as.factor(df_src_glm$ID_pnTTC)
          df_src_glm$sex<-as.factor(df_src_glm$sex)
          
          df_src_glm_bind<-rbind(df_src_glm_bind,
                                 cbind(atlas=atlas,measure=measure,group_1=group_1,group_2=group_2,
                                       df_src_glm))
          
          list_src_glm<-c(list_src_glm,
                          list(list("df_src"=df_src_glm,"atlas"=atlas,"measure"=measure,
                                    "group_1"=group_1,"group_2"=group_2)))
          
          # Calculate GLM
          #out_glm<-glm_core(df_src=df_join_grp,atlas,measure,group_1,group_2,
          #                  list_mod_,list_graph_,list_covar_,paths_)
          #df_out_lm<-rbind(df_out_lm,out_glm$df_out_lm_add)
          #df_out_aic<-rbind(df_out_aic,out_glm$df_out_aic_add)
          
          # Prepare ANCOVA dataset for later parallel computing
          if (!skip_ancova){
            print(paste("Atlas: ",atlas,", Measure: ",measure,", Group: ",group_1," and ",group_2,", ANCOVA preparation.",  sep=""))
            # Create list of input dataframes for parallel ANCOVA calculation
            for (group_tanner in names(list_strat_tanner_)){
              # group by longitudinal Tanner stage
              df_join_grp_tanner<-df_join_grp
              for (ses in c(1,2)){
                list_strat_tanner_ses<-names(list_strat_tanner_[[group_tanner]][[as.character(ses)]])
                for (label_tanner in list_strat_tanner_ses){
                  list_strat_tanner_ses_group<-list_strat_tanner_[[group_tanner]][[as.character(ses)]][[label_tanner]]
                  #print(list_strat_tanner_ses_group)
                  df_join_grp_tanner[df_join_grp_tanner[[paste('ses',as.character(ses),'_tanner',sep='')]] %in% list_strat_tanner_ses_group,
                                     paste('ses',as.character(ses),'_tanner_label',sep='')]<-label_tanner
                }
              }
              df_join_grp_tanner$long_tanner<-paste(as.character(df_join_grp_tanner$ses1_tanner_label),
                                                    as.character(df_join_grp_tanner$ses2_tanner_label),sep="_")
              df_join_grp_tanner$long_tanner<-as.factor(df_join_grp_tanner$long_tanner)
              
              list_sex<-list("all"=c(1,2),"male"=1,"female"=2)
              
              for (idx_sex in names(list_sex)){
                df_join_grp_tanner_sex<-df_join_grp_tanner[df_join_grp_tanner$sex %in% list_sex[[idx_sex]],]
                list_src_ancova<-c(list_src_ancova,
                                   list(list("atlas"=atlas,"measure"=measure,"group_network"=c(group_1,group_2),
                                             "group_tanner"=group_tanner,"group_tanner_content"=list_strat_tanner_[[group_tanner]],
                                             "group_sex"=idx_sex,"df_src_ancova"=df_join_grp_tanner_sex)))
              }
            }
          }
        }
      } # Finished looping over ROI subgroups
    } # Finished looping over measure
  } # Finished looping over atlas
  
  write.csv(df_src_glm_bind,file.path(paths_$output,"output",
                                      paste(suffix,"_fp_glm_src.csv",sep="")),row.names = F)

  # Parallel GLM calculation
  n_cluster<-floor(detectCores()*3/4)
  print(paste("Calculating GLM/GAM in parallel, threads: ",as.character(n_cluster),sep=""))
  clust<-makeCluster(n_cluster)
  clusterExport(clust,
                varlist=c("list_mod_","list_graph_","list_covar_","paths_","as.numeric.factor",
                          "gam","as.formula","summary.gam","plot_gamm","ggplot",
                          "predict.gam","geom_line","aes","geom_ribbon","geom_point",
                          "theme_light","element_text","ggtitle",
                          "xlab","ylab","theme","ggsave"),
                envir=environment())
  list_dst_glm<-pblapply(list_src_glm,glm_core,cl=clust)
  stopCluster(clust)
  
  df_out_lm<-df_out_aic<-NULL
  for (dst_glm in list_dst_glm){
    df_out_lm<-rbind(df_out_lm,dst_glm$df_out_lm_add)
    df_out_aic<-rbind(df_out_aic,dst_glm$df_out_aic_add)
  }
  
  # GLM/GAM Data saving
  rownames(df_out_lm)<-rownames(df_out_aic)<-NULL
  write.csv(df_out_lm, file.path(paths_$output,"output",paste(suffix,"_fp_glm.csv",sep="")),row.names = F)
  write.csv(df_out_aic,file.path(paths_$output,"output",paste(suffix,"_fp_glm_aic.csv",sep="")),row.names = F)
  
  
  # group-wise GLM/GAM heatmap visualization
  heatmap_gammfp(paths_,df_out_lm,list_atlas_,list_mod_)
  
  # Parallel ANCOVA calculation
  if (!skip_ancova){
    print("Calculating ANCOVA in parallel.")
    n_cluster<-min(floor(detectCores()*3/4),length(list_src_ancova))
    clust<-makeCluster(n_cluster)
    clusterExport(clust,
                  varlist=c("aov","summary","glht","mcp","left_join","paths_",
                            "plot_cor_heatmap","rcorr","rownames_to_column","gather",
                            "ggplot","aes","geom_tile","scale_fill_gradientn",
                            "matlab.like2","scale_y_discrete","scale_x_discrete",
                            "theme_light","theme","element_text","element_blank",
                            "ggtitle","ggsave","scale_fill_gradient","xlab","ylab"),
                  envir=environment())
    list_df_ancova<-pblapply(list_src_ancova,ancova_core,cl=clust)
    stopCluster(clust)
    df_out_ancova<-NULL
    for (df_ancova in list_df_ancova){
      df_out_ancova<-rbind(df_out_ancova,df_ancova)
    }
    # Data saving
    rownames(df_out_ancova)<-NULL
    write.csv(df_out_ancova,file.path(paths_$output,"output","fp_ancova.csv"),row.names=F)
  }

  return(list("df_out_lm"=df_out_lm,"df_out_aic"=df_out_aic,"df_src"=df_src_glm_bind))
}


#**************************************************
# Fingerprint identification ======================
#**************************************************
identify_fp<-function(paths_=paths,list_atlas_=list_atlas,list_wave_=list_wave,
                      #list_covar_=list_covar,
                      subset_subj_=subset_subj,n_permutation_=n_permutation
                      ){
  print("Starting identify_fp().")
  nullobj<-func_createdirs(paths_)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar=NULL,rem_na_clin=F)
  df_clin<-data_clin$df_clin
  colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
  
  df_out_combined<-NULL
  
  for (atlas in list_atlas_){
    # Load fingerprint data
    df_fp<-read.csv(file.path(paths_$input,"output",paste("atl-",atlas,"_fp.csv",sep="")))
    
    list_measure<-sort(unique(df_fp$measure))
    
    for (measure in list_measure){
      print(paste("Atlas: ",atlas," Measure: ",measure,sep=""))
      df_fp_meas<-df_fp[df_fp$measure==measure,]
      
      # Create list of subjects who meet subsetting condition and whose MRI data exist
      list_ses_exist <- sort(unique(c(df_fp_meas$from_ses,df_fp_meas$to_ses)))
      list_id_subj_exist<-list()
      for (ses in list_ses_exist){
        id_subj_exist_ses<-sort(unique(c(df_fp_meas[df_fp_meas$from_ses==ses,'from_ID_pnTTC'],
                                         df_fp_meas[df_fp_meas$to_ses==ses,'to_ID_pnTTC'])))
        id_subj_subset_ses<-df_clin[df_clin$ses==ses,"ID_pnTTC"]
        id_subj_exist_ses<-intersect(id_subj_exist_ses,id_subj_subset_ses)
        list_id_subj_exist[[as.character(ses)]]<-sort(id_subj_exist_ses)
      }
      
      # Extract those with longitudinal data
      list_id_subj_exist_twice<-sort(intersect(list_id_subj_exist[["1"]],
                                               list_id_subj_exist[["2"]]))
      n_id_subj_exist_twice<-length(list_id_subj_exist_twice)
      print(paste(as.character(n_id_subj_exist_twice)," subjects with two sessions.",sep=""))
      df_fp_exist_twice<-df_fp_meas[(df_fp_meas$from_ID_pnTTC %in% list_id_subj_exist_twice) & (df_fp_meas$to_ID_pnTTC %in% list_id_subj_exist_twice),]
      df_fp_exist_twice<-df_fp_exist_twice[(df_fp_exist_twice$from_ses==1 & df_fp_exist_twice$to_ses==2),]
      
      # Output subset with longitudinal data
      write.csv(df_fp_exist_twice,file.path(paths_$output,"output",paste("atl-",atlas,"_msr-",measure,"_fp_input_subset.csv",sep="")),row.names=F)
      
      list_group<-sort(unique(as.character(df_fp_exist_twice$group)))
      if ("whole" %in% list_group){
        list_group<-c("whole",list_group[list_group!="whole"])
      }
      
      # Correlation matrix graphical output
      for (group in list_group){
        df_fp_exist_twice_plot<-df_fp_exist_twice[df_fp_exist_twice$group==group,c('from_ID_pnTTC','to_ID_pnTTC','r')]
        df_fp_exist_twice_plot<-spread(df_fp_exist_twice_plot,key=to_ID_pnTTC,value=r)
        colnames(df_fp_exist_twice_plot)[-1]<-rownames(df_fp_exist_twice_plot)<-sprintf("%05d",df_fp_exist_twice_plot$from_ID_pnTTC)
        df_fp_exist_twice_plot<-df_fp_exist_twice_plot[-1]
        plot_fp_exist_twice<-plot_cor_heatmap(input=df_fp_exist_twice_plot)
        suppressMessages(plot_fp_exist_twice<-(plot_fp_exist_twice
                                               + scale_fill_gradientn(colors = matlab.like2(100),name="r")
                                               + ggtitle(paste("Long. FP Cor,",atlas,measure,group,sep=" "))
                                               + xlab("2nd wave")
                                               + ylab("1st wave")
                                               + theme(plot.title = element_text(hjust = 0.5))))
        ggsave(paste("atl-",atlas,"_msr-",measure,"_grp-",group,"_fp_id.eps",sep=""),plot=plot_fp_exist_twice,device=cairo_ps,
                     path=file.path(paths_$output,"output"),dpi=300,height=10,width=10,limitsize=F)
      }
      
      # Calculate fingerprint identification
      df_ident<-df_perm<-NULL
      
      for (group in list_group){
        df_ident_grp<-data.frame("group"=group,"target"=list_id_subj_exist_twice)
        df_perm_grp<-data.frame("group"=group,"id_perm"=seq(1,n_permutation_))
        df_fp_exist_twice_group<-df_fp_exist_twice[df_fp_exist_twice$group==group,]
        for(ses in c(1,2)){
          if(ses==1){
            df_fp_pool<-df_fp_exist_twice_group[c("from_ses","from_ID_pnTTC","to_ses","to_ID_pnTTC","r")]
          }else{
            df_fp_pool<-df_fp_exist_twice_group[c("to_ses","to_ID_pnTTC","from_ses","from_ID_pnTTC","r")]
          }
          colnames(df_fp_pool)<-c("target_ses","target_ID_pnTTC","pool_ses","pool_ID_pnTTC","r")
          for (id_subj in list_id_subj_exist_twice){
            df_fp_subset<-df_fp_pool[df_fp_pool$target_ID_pnTTC==id_subj,]
            list_id_subj_ordered<-df_fp_subset[order(df_fp_subset$r,decreasing = TRUE,na.last=NA),'pool_ID_pnTTC']
            if (id_subj %in% list_id_subj_ordered){
              rank_similarity<-which(list_id_subj_ordered==id_subj)
            }else{
              rank_similarity<-NA
            }
            df_ident_grp[df_ident_grp$target==id_subj,paste(as.character(ses),"_targeted_rank",sep='')]<-rank_similarity
            if (is.na(rank_similarity)){
              df_ident_grp[df_ident_grp$target==id_subj,paste(as.character(ses),"_targeted_identification",sep='')]<-0
            }else{
              if (rank_similarity==1){
                df_ident_grp[df_ident_grp$target==id_subj,paste(as.character(ses),"_targeted_identification",sep='')]<-1
              }else{
                df_ident_grp[df_ident_grp$target==id_subj,paste(as.character(ses),"_targeted_identification",sep='')]<-0
              }
            }
          }
          
          # Permutation test calculation
          for (i in seq(1,n_permutation_)){
            # Create dataframe for random shuffling
            df_rand<-data.frame(pool_ID_pnTTC=list_id_subj_exist_twice,rand_ID_pnTTC=sample(list_id_subj_exist_twice,length(list_id_subj_exist_twice)))
            
            # Randomize subject ID of pool session accordint to df_rand
            df_fp_rand<-left_join(df_fp_pool,df_rand,by="pool_ID_pnTTC")
            n_identified<-0
            for (id_subj in list_id_subj_exist_twice){
              df_fp_rand_subj<-df_fp_rand[df_fp_rand$target_ID_pnTTC==id_subj,]
              list_id_subj_ordered<-df_fp_rand_subj[order(df_fp_rand_subj$r,decreasing = TRUE,na.last=NA),'rand_ID_pnTTC']
              if (id_subj %in% list_id_subj_ordered){
                rank_similarity<-which(list_id_subj_ordered==id_subj)
              }else{
                rank_similarity<-NA
              }
              if (!is.na(rank_similarity)){
                if (rank_similarity==1){
                  n_identified<-n_identified+1
                }
              }
            }
            df_perm_grp[i,paste(as.character(ses),'_n_ident',sep='')]<-n_identified
            #print(paste("Iteration: ",as.character(i),", subjects identified: ",as.character(n_identified),sep=""))
          }
        }
        df_ident<-rbind(df_ident,df_ident_grp)
        df_perm<-rbind(df_perm,df_perm_grp)
        
        n_subj<-length(list_id_subj_exist_twice)
        n_id_1_tar<-sum(df_ident_grp["1_targeted_identification"])
        n_id_2_tar<-sum(df_ident_grp["2_targeted_identification"])
        n_id<-n_id_1_tar+n_id_2_tar
        prop_id_1_tar<-n_id_1_tar/n_subj
        prop_id_2_tar<-n_id_2_tar/n_subj
        prop_id<-n_id/(n_subj*2)
        p_perm_1_tar<-sum(df_perm_grp["1_n_ident"]>sum(df_ident_grp["1_targeted_identification"]))/n_permutation_
        p_perm_2_tar<-sum(df_perm_grp["2_n_ident"]>sum(df_ident_grp["2_targeted_identification"]))/n_permutation_
        p_perm<-(p_perm_1_tar+p_perm_2_tar)/2
        
        df_out_combined<-rbind(df_out_combined,
                               data.frame(atlas=atlas,measure=measure,group=group,n_subj=n_subj,n_identified=n_id,proportion_identified=prop_id,
                                          n_identified_1_targeted=n_id_1_tar,proportion_identifeid_1_targeted=prop_id_1_tar,
                                          n_identified_2_targeted=n_id_2_tar,proportion_identified_2_targeted=prop_id_2_tar,
                                          p_permutation=p_perm,
                                          p_permutation_1_targeted=p_perm_1_tar,p_permutation_2_targeted=p_perm_2_tar))
      }
      write.csv(df_ident,file.path(paths_$output,"output",paste("atl-",atlas,"_msr-",measure,"_fp_id.csv",sep="")),row.names=F)
      write.csv(df_perm,file.path(paths_$output,"output",paste("atl-",atlas,"_msr-",measure,"_fp_perm.csv",sep="")),row.names=F)
    }
  }
  write.csv(df_out_combined,file.path(paths_$output,"output","fp_id_summary.csv"),row.names=F)
  print("Finished identify_fp().")
}


#**************************************************
# heatmap subnetwork gamm fp results ==============
#**************************************************
heatmap_gammfp_old<-function(paths_=paths){
  df_gammfp<-read.csv(file.path(paths_$input,"output","fp_glm.csv"))
  df_gammfp<-df_gammfp[df_gammfp$atlas=="shen268" & df_gammfp$sex==2
                       & df_gammfp$model=="ldm" & df_gammfp$term=="diff_tanner",]
  list_group<-sort(unique(c(as.character(df_gammfp$group_1),as.character(df_gammfp$group_2))))
  list_group<-list_group[list_group!="whole"]
  df_gammfp<-df_gammfp[df_gammfp$group_1 %in% list_group
                       & df_gammfp$group_2 %in% list_group,]
  df_plot<-data.frame(matrix(ncol=length(list_group),nrow=length(list_group)))
  colnames(df_plot)<-rownames(df_plot)<-list_group
  df_sign<-data.frame(matrix(ncol=length(list_group),nrow=length(list_group)))
  colnames(df_sign)<-rownames(df_sign)<-list_group
  for (i_row in seq(dim(df_gammfp)[1])){
    g1<-as.character(df_gammfp[i_row,"group_1"])
    g2<-as.character(df_gammfp[i_row,"group_2"])
    estimate<-as.numeric(df_gammfp[i_row,"estimate"])
    if (df_gammfp[i_row,"p"]<0.001){
      sign<-"**"
    }else if (df_gammfp[i_row,"p"]<0.05){
      sign<-"*"
    }else{
      sign<-""
    }
    
    df_plot[g1,g2]<-df_plot[g2,g1]<-estimate
    df_sign[g1,g2]<-df_sign[g2,g1]<-sign
  }
  plot<-plot_cor_heatmap(df_plot,df_sign)
  suppressMessages(plot<-(plot
                          + scale_fill_gradientn(colors = matlab.like2(100),
                                                 lim=c(-max(max(df_plot),-min(df_plot)),max(max(df_plot),-min(df_plot))),
                                                 name="beta")
                          + ggtitle("Subnetwork-wise effect of Tanner difference")
                          + theme(plot.title = element_text(hjust = 0.5),
                                  axis.title=element_blank())))
  ggsave(paste("atl-shen268_msr-fc_gammfp_subnetwork.eps",sep=""),plot=plot,device=cairo_ps,
         path=file.path(paths_$output,"output"),height=7,width=7,limitsize=F)
  
}
