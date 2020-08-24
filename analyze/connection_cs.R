#**************************************************
# Description =====================================
#**************************************************
# R script to analyze relationship between structural/functional connection and clinical data.
# Inputs can be functional correlation from rsfMRI data, or Jackknife esimate of structural covariance from T1 data.


#**************************************************
# Parameters ======================================
#**************************************************

path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"

dir_in<-"421_fc_aroma"
dir_out<-"425_fc_gam_aroma"
list_waves<-list("c1m1" =list("wave_clin"="1","wave_mri"="1"),
                 "c1m2" =list("wave_clin"="1","wave_mri"="2"),
                 #"c1m21"=list("wave_clin"="1","wave_mri"="2-1"),
                 "c2m1" =list("wave_clin"="2","wave_mri"="1"),
                 "c2m2" =list("wave_clin"="2","wave_mri"="2"))
#wave_clin<-1
#wave_mri<-2

#dir_in<-"426_fc_diff_aroma"
#dir_out<-"427.2_gam_fc_diff_aroma"
#wave_clin<-1
#wave_mri<-"2-1"

#list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400","shen268")
#list_atlas<-c("aal116","gordon333","power264","schaefer400x7","shen268")
list_atlas<-"aal116"

list_covar_tanner<-list("tanner"=list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage"),
                        "age"   =list("1"="W1_Age_at_MRI",  "2"="W2_Age_at_MRI",  "label"="Age"),
                        "sex"   =list("1"="Sex",            "2"="Sex",            "label"="Sex"))

list_tanner<-list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage (max)"),
                  "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)"),
                  "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),
                                 "2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                 "label"="Tanner stage (gonadal)"),
                  "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),
                                 "2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
                                 "label"="Tanner stage (adrenal)"))

list_covar_hormone<-list("hormone"=list("1"="W1_Hormone"   ,"2"="W2_Hormone",   "label"="Hormone"),
                         "age"    =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                         "sex"    =list("1"="Sex",          "2"="Sex",          "label"="Sex"))

list_hormone<-list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
                   "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
                   "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
                   "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S"))

subset_subj <- list("1"  =list(list("key"="W1_T1QC","condition"="==1"),
                               list("key"="W1_rsfMRIexist","condition"="==1"),
                               list("key"="W1_Censor","condition"="<126")),
                    "2"  =list(list("key"="W2_T1QC","condition"="==1"),
                               list("key"="W2_rsfMRIexist","condition"="==1"),
                               list("key"="W2_Censor","condition"="<126")),
                    "2-1"=list(list("key"="W1_T1QC","condition"="==1"),
                               list("key"="W1_rsfMRIexist","condition"="==1"),
                               list("key"="W1_Censor","condition"="<126"),
                               list("key"="W2_T1QC","condition"="==1"),
                               list("key"="W2_rsfMRIexist","condition"="==1"),
                               list("key"="W2_Censor","condition"="<126")))

list_mod_tanner <- list("l"= "value ~ age + tanner")
list_mod_hormone <- list("l"= "value ~ age + hormone")
#"a"= "value ~ s(age,k=3) + s(testo,k=3) + s(ID_pnTTC,bs='re')",
#"q"="value ~ poly(age,2) + poly(testo,2) + s(ID_pnTTC,bs='re')")

list_plot_tanner <-list("t"=list("title"="Tanner effect","var_exp"="tanner"))
list_plot_hormone <-list("h"=list("title"="Hormone effect","var_exp"="hormone"))

list_type_p=c("p","p_bh","seed_p_bh")
thr_p <- 0.05


#**************************************************
# Libraries =======================================
#**************************************************
library(mgcv)
library(data.table)


#**************************************************
# Original library ================================
#**************************************************
source(file.path(getwd(),"util/function.R"))
source(file.path(getwd(),"util/plot.R"))
paths<-func_path(path_exp_=path_exp,dir_in_=dir_in,dir_out_=dir_out)



#**************************************************
# Iterate gam_fc_cs() over clinical variables =====
# and waves =======================================
#**************************************************

paths_=paths
list_waves_=list_waves
subset_subj_=subset_subj
list_covar_tanner_=list_covar_tanner
list_tanner_=list_tanner
list_mod_tanner_=list_mod_tanner
list_plot_tanner_=list_plot_tanner
list_covar_hormone_=list_covar_hormone
list_hormone_=list_hormone
list_mod_hormone_=list_mod_hormone
list_plot_hormone_=list_plot_hormone
list_type_p_=list_type_p
thr_p_=thr_p

gam_fc_cs_multi<-function(paths_=paths,list_waves_=list_waves,subset_subj_=subset_subj,
                          list_covar_tanner_=list_covar_tanner,list_tanner_=list_tanner,
                          list_mod_tanner_=list_mod_tanner,list_plot_tanner_=list_plot_tanner,
                          list_covar_hormone_=list_covar_hormone,list_hormone_=list_hormone,
                          list_mod_hormone_=list_mod_hormone,list_plot_hormone_=list_plot_hormone,
                          list_type_p_=list_type_p,thr_p_=thr_p){
  print("Starting gam_fc_cs_multi()")
  nullobj<-func_createdirs(paths_,str_proc="gam_fc_cs_multi()",copy_log=T)
  
  for (waves in names(list_waves_)){
    wave_clin<-list_waves_[[waves]]$wave_clin
    wave_mri<-list_waves_[[waves]]$wave_mri
    print(paste("Clinical wave: ", wave_clin,", MRI wave: ",wave_mri,sep=""))
    
    # Prepare subject subsetting condition (MRI QC criteria) according to specified waves
    subset_subj_temp<-subset_subj_[[as.character(wave_mri)]]
    subset_subj_temp<-list(subset_subj_temp)
    names(subset_subj_temp)<-wave_clin

    #1 Tanner stage
    for (idx_tanner in names(list_tanner_)){
      print(paste("Tanner type: ",list_tanner_[[idx_tanner]][["label"]],sep=""))
      list_covar<-list_covar_tanner_
      list_covar[["tanner"]]<-list_tanner_[[idx_tanner]]
      suffix<-paste("_wave-",waves,"_var-",idx_tanner,sep="")
      
      nullobj<-gam_fc_cs(paths_=paths_,subset_subj_=subset_subj_temp,list_covar_=list_covar,
                         wave_clin_=wave_clin,wave_mri_=wave_mri,list_atlas_=list_atlas,
                         list_mod_=list_mod_tanner_,list_plot_=list_plot_tanner_,
                         key_group_='group_3',list_type_p_=list_type_p_,thr_p_=thr_p_,
                         suffix_=suffix)
    } # finished looping over Tanner stages
    
    
    #2 Hormones
    for (idx_hormone in names(list_hormone_)){
      print(paste("Hormone type: ",list_hormone_[[idx_hormone]][["label"]],sep=""))
      list_covar<-list_covar_hormone_
      list_covar[["hormone"]]<-list_hormone_[[idx_hormone]]
      suffix<-paste("_wave-",waves,"_var-",idx_hormone,sep="")
      
      nullobj<-gam_fc_cs(paths_=paths_,subset_subj_=subset_subj_temp,list_covar_=list_covar,
                         wave_clin_=wave_clin,wave_mri_=wave_mri,list_atlas_=list_atlas,
                         list_mod_=list_mod_hormone_,list_plot_=list_plot_hormone_,
                         key_group_='group_3',list_type_p_=list_type_p_,thr_p_=thr_p_,
                         suffix_=suffix)
    } # finished looping over Hormones
    
  } # finished looping over waves
  print("Finished gam_fc_cs_multi()")
}



#**************************************************
# Additive/Linear model of FC in cross-section ====
#**************************************************

join_fc_clin<-function(df_fc,df_clin,wave_clin,wave_mri){
  df_fc$z_r[which(is.nan(df_fc$z_r))]<-0
  colnames(df_fc)[colnames(df_fc)=="z_r"]<-"value"
  df_fc<-df_fc[df_fc$ses==wave_mri,]
  df_fc$ses<-NULL
  #colnames(df_fc)[colnames(df_fc)=="ses"]<-"wave"
  df_fc<-df_fc[,c("ID_pnTTC","from","to","value")]
  
  df_clin<-df_clin[df_clin$wave==wave_clin,]
  df_clin$wave<-NULL
  
  # Join clinical and FC data frames
  print('Joining clinical and FC data.')
  #df_join<-inner_join(df_fc,df_clin,by=c('ID_pnTTC','wave'))
  df_join<-inner_join(df_fc,df_clin,by='ID_pnTTC')
  for (key in c('ID_pnTTC','wave','sex')){
    if (key %in% colnames(df_join)){
      df_join[,key]<-as.factor(df_join[,key])
    }
  }
  return(df_join)
}

paths_=paths_
subset_subj_=subset_subj_temp
list_covar_=list_covar
wave_clin_=wave_clin
wave_mri_=wave_mri
list_atlas_=list_atlas
list_mod_=list_mod_tanner_
list_plot_=list_plot_tanner_
key_group_='group_3'
list_type_p_=list_type_p_
thr_p_=thr_p_
suffix_=suffix

gam_fc_cs<-function(paths_=paths,subset_subj_=subset_subj,list_covar_=list_covar,
                    wave_clin_=wave_clin,wave_mri_=wave_mri,list_atlas_=list_atlas,
                    list_mod_=list_mod,list_plot_=list_plot,key_group_='group_3',
                    list_type_p_=list_type_p,thr_p_=thr_p,suffix_=suffix
                    ){
  print("Starting gam_fc_cs().")
  nullobj<-func_createdirs(paths_,str_proc="gam_fc_cs()",copy_log=T)
  dict_roi <- func_dict_roi(paths_)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  data_clin<-func_clinical_data_long(paths_,list_wave=wave_clin_,subset_subj_,
                                     list_covar=list_covar_,rem_na_clin=T)
  df_clin<-data_clin$df_clin
  
  
  for (atlas in list_atlas_){
    
    # Load ROI-wise FC data
    print(paste('Loading FC data, atlas:',atlas,sep=' '))
    #df_fc<-read.csv(file.path(paths_$input,'output',paste('atl-',atlas,'_fc.csv',sep='')))
    df_fc<-fread(file.path(paths_$input,'output',paste('atl-',atlas,'_fc.csv',sep='')))
    df_join<-join_fc_clin(df_fc,df_clin,wave_clin_,wave_mri_)
    write.csv(df_join,file.path(paths_$output,"output",paste("atl-",atlas,suffix_,"_src.csv",sep="")),
              row.names=F)
    
    # Calculate and save ROI-wise GAMM of FC
    print(paste('Calculating GAM, atlas: ',atlas,sep=''))
    list_roi<-sort(unique(c(as.character(df_join$from),as.character(df_join$to))))
    df_roi<-dict_roi[is.element(dict_roi$id,list_roi),c("id","label",key_group_)]
    colnames(df_roi)[colnames(df_roi)==key_group_]<-"group"
    data_gamm<-iterate_gamm(df_join,df_roi,list_mod_)
    write.csv(data_gamm$df_out_gamm,
              file.path(paths_$output,"output",paste("atl-",atlas,suffix_,"_gam.csv",sep="")),row.names = F)
    write.csv(data_gamm$df_out_aic,
              file.path(paths_$output,"output",paste("atl-",atlas,suffix_,"_gam_aic.csv",sep="")),row.names = F)
    
    # Calculate multiple comparison-corrected p values
    df_plot_gamm<-add_mltcmp(data_gamm$df_out_gamm,df_roi,analysis="roi",atlas,
                             list_mod_,list_plot_,calc_seed_level=T)
    write.csv(df_plot_gamm,
              file.path(paths_$output,"output",paste("atl-",atlas,suffix_,"_gam_plt.csv",sep="")),row.names = F)
    
    # Graphical output of ROI-wise GAMM of FC
    plot_gam_fc(df_plot_gamm,df_roi,analysis="roi",atlas,list_mod,list_plot,
                list_type_p_,thr_p,paths_,suffix_)
    
  }
  print('Finished gam_fc_cs().')
}
