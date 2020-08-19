#**************************************************
# Description =====================================
#**************************************************
# R script to analyze relationship between structural/functional connection and clinical data.
# Inputs can be functional correlation from rsfMRI data, or Jackknife esimate of structural covariance from T1 data.


#**************************************************
# Parameters ======================================
#**************************************************

path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"

#dir_in<-"421_fc_aroma"
#dir_out<-"427.1_gam_fc_cs_aroma"
#wave_clin<-1
#wave_mri<-2

dir_in<-"426_fc_diff_aroma"
dir_out<-"427.2_gam_fc_diff_aroma"
wave_clin<-1
wave_mri<-"2-1"

#list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400","shen268")
list_atlas<-"aal116"



list_covar<-list("tanner"=list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage"),
                 "age"   =list("1"="W1_Age_at_MRI",  "2"="W2_Age_at_MRI",  "label"="Age"),
                 "sex"   =list("1"="Sex",            "2"="Sex",            "label"="Sex"))

#list_covar<-list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
#                 "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
#                 "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
#                 "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S"),
#                 "age"  =list("1"="W1_Age_at_MRI",  "2"="W2_Age_at_MRI",  "label"="Age"),
#                 "sex"  =list("1"="Sex",            "2"="Sex",            "label"="Sex"))

subset_subj <- list("1"=list(list("key"="W2_T1QC","condition"="==1"),
                             list("key"="W2_rsfMRIexist","condition"="==1"),
                             list("key"="W2_Censor","condition"="<126")))

list_mod <- list("l"= "value ~ age + tanner")
#"a"= "value ~ s(age,k=3) + s(testo,k=3) + s(ID_pnTTC,bs='re')",
#"q"="value ~ poly(age,2) + poly(testo,2) + s(ID_pnTTC,bs='re')")

list_plot <-list("t"=list("title"="Tanner effect","var_exp"="tanner"))

list_type_p=c("p","p_bh","seed_p_bh")
thr_p <- 0.05


#**************************************************
# Libraries =======================================
#**************************************************

library(mgcv)


#**************************************************
# Original library ================================
#**************************************************
source(file.path(getwd(),"util/function.R"))
source(file.path(getwd(),"util/plot.R"))
paths<-func_path(path_exp_=path_exp,dir_in_=dir_in,dir_out_=dir_out)


#**************************************************
# Additive/Linear model of FC in cross-section ====
#**************************************************

join_fc_clin<-function(df_fc,df_clin,wave_clin,wave_mri){
  df_fc$z_r[which(is.nan(df_fc$z_r))]<-0
  colnames(df_fc)[colnames(df_fc)=="z_r"]<-"value"
  df_fc<-df_fc[df_fc$ses==wave_mri,]
  df_fc$ses<-NULL
  #colnames(df_fc)[colnames(df_fc)=="ses"]<-"wave"
  df_fc<-df_fc[,c(-which(colnames(df_fc)=="r"),
                  -which(colnames(df_fc)=="p"))]
  
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


gam_fc_cs<-function(paths_=paths,subset_subj_=subset_subj,list_covar_=list_covar,
                    wave_clin_=wave_clin,wave_mri_=wave_mri,list_atlas_=list_atlas,
                    list_mod_=list_mod,list_plot_=list_plot,key_group_='group_3',
                    list_type_p_=list_type_p,thr_p_=thr_p
                    ){
  print("Starting gam_fc_cs().")
  nullobj<-func_createdirs(paths_,str_proc="gam_fc_cs()",copy_log=T)
  dict_roi <- func_dict_roi(paths_)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  data_clin<-func_clinical_data_long(paths_,list_wave=wave_clin,subset_subj_,
                                     list_covar=list_covar_,rem_na_clin=T)
  df_clin<-data_clin$df_clin
  
  
  for (atlas in list_atlas_){
    
    #****************************
    # ROI-wise FC AMM calculation
    #****************************
    # Load ROI-wise FC data
    print(paste('Loading FC data, atlas:',atlas,sep=' '))
    df_fc<-read.csv(file.path(paths_$input,'output',paste('atl-',atlas,'_fc.csv',sep='')))
    df_join<-join_fc_clin(df_fc,df_clin,wave_clin_,wave_mri_)
    write.csv(df_join,file.path(paths_$output,"output",paste("atl-",atlas,"_anl-roi_src.csv",sep="")),
              row.names=F)
    
    # Calculate and save ROI-wise GAMM of FC
    print(paste('Calculating GAM, atlas: ',atlas,sep=''))
    list_roi<-sort(unique(c(as.character(df_join$from),as.character(df_join$to))))
    df_roi<-dict_roi[is.element(dict_roi$id,list_roi),c("id","label",key_group_)]
    colnames(df_roi)[colnames(df_roi)==key_group_]<-"group"
    data_gamm<-iterate_gamm(df_join,df_roi,list_mod_)
    write.csv(data_gamm$df_out_gamm,
              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-roi_gam.csv",sep="")),row.names = F)
    write.csv(data_gamm$df_out_aic,
              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-roi_gam_aic.csv",sep="")),row.names = F)
    
    # Calculate multiple comparison-corrected p values
    df_plot_gamm<-add_mltcmp(data_gamm$df_out_gamm,df_roi,analysis="roi",atlas,
                             list_mod_,list_plot_,calc_seed_level=T)
    write.csv(df_plot_gamm,
              file.path(paths_$output,"output",paste("atl-",atlas,"_anl-roi_gam_plt.csv",sep="")),row.names = F)
    
    # Graphical output of ROI-wise GAMM of FC
    plot_gam_fc(df_plot_gamm,df_roi,analysis="roi",atlas,list_mod,list_plot,
                list_type_p_,thr_p,paths_)
    
    ##****************************
    ## Group-wise FC GAMM calculation
    ##****************************
    ## Load group-wise FC data
    #print(paste('Loading group FC data, atlas:',atlas,sep=' '))
    #df_fc_grp<-read.csv(file.path(paths_$input,'output',paste('atl-',atlas,'_fc_grp.csv',sep='')))
    #df_join_grp<-join_fc_clin(df_fc_grp,df_clin,wave_clin_,wave_mri_)
    #write.csv(df_join_grp,file.path(paths_$output,"output",paste("atl-",atlas,"_anl-grp_src.csv",sep="")),
    #          row.names=F)
    #
    ## Calculate and save group-wise GAMM of FC
    #print(paste('Calculating GAM, atlas: ',atlas,sep=''))
    #list_roi_grp<-sort(unique(c(as.character(df_join_grp$from),as.character(df_join_grp$to))))
    ##df_roi<-dict_roi[is.element(dict_roi$id,list_roi),c("id","label",key_group_)]
    #df_roi_grp<-data.frame(id=list_roi_grp,label=capitalize(list_roi_grp),group="group")
    #data_gamm_grp<-iterate_gamm(df_join_grp,df_roi_grp,list_mod_)
    #write.csv(data_gamm_grp$df_out_gamm,
    #          file.path(paths_$output,"output",paste("atl-",atlas,"_anl-grp_gam.csv",sep="")),row.names = F)
    #write.csv(data_gamm_grp$df_out_aic,
    #          file.path(paths_$output,"output",paste("atl-",atlas,"_anl-grp_gam_aic.csv",sep="")),row.names = F)
    #
    ## Calculate multiple comparison-corrected p values
    #df_plot_gamm_grp<-add_mltcmp(data_gamm_grp$df_out_gamm,df_roi_grp,analysis="grp",atlas,
    #                             list_mod_,list_plot_,calc_seed_level=T)
    #write.csv(df_plot_gamm_grp,
    #          file.path(paths_$output,"output",paste("atl-",atlas,"_anl-grp_gam_plt.csv",sep="")),row.names = F)
    #
    ## Graphical output of group-wise GAMM of FC
    #plot_gam_fc(df_plot_gamm_grp,df_roi_grp,analysis="grp",atlas,list_mod,list_plot,
    #            list_type_p_,thr_p,paths_)
    #
    ##****************************
    ## Multi-scale FC GAMM calculation
    ##****************************
    ## Subset ROI-wise GAMM result to include only within-group connections
    #df_gamm_ms<-NULL
    #for (group in list_roi_grp){
    #  list_roi_within_grp<-as.character(df_roi[df_roi$group==group,"id"])
    #  df_gamm_ms_add<-data_gamm$df_out_gamm[which(is.element(as.character(data_gamm$df_out_gamm[,"from"]),list_roi_within_grp)
    #                                              & is.element(as.character(data_gamm$df_out_gamm[,"to"]),list_roi_within_grp)),]
    #  df_gamm_ms_add<-cbind(group=group,df_gamm_ms_add)
    #  df_gamm_ms<-rbind(df_gamm_ms,df_gamm_ms_add)
    #}
    #
    ## Combine within-group ROI-wise GAMM results and between-group GAMM results
    #df_gamm_ms<-rbind(df_gamm_ms,cbind(group="group",data_gamm_grp$df_out_gamm))
    #
    ## Calculate multiple comparison-corrected p values
    #df_plot_gamm_ms<-add_mltcmp(df_gamm_ms,df_roi_grp,analysis="grp",atlas,list_mod,list_plot,
    #                            calc_seed_level=F)
    #write.csv(df_plot_gamm_ms,
    #          file.path(paths_$output,"output",paste("atl-",atlas,"_anl-ms_gam_plt.csv",sep="")),row.names = F)
    #
    ## Split data into ROI-wise and group-wise GAMM results, graphical output
    #for (group in list_roi_grp){
    #  df_plot_gamm_ms_split<-df_plot_gamm_ms[df_plot_gamm_ms$group==group,-1]
    #  df_roi_split<-df_roi[df_roi$group==group,]
    #  label_analysis<-paste("ms_grp-",group,sep="")
    #  plot_gam_fc(df_plot_gamm_ms_split,df_roi_split,analysis=label_analysis,atlas,list_mod,list_plot,
    #               list_type_p_,thr_p,paths_)
    #}
    #df_plot_gamm_ms_split<-df_plot_gamm_ms[df_plot_gamm_ms$group=="group",-1]
    #plot_gam_fc(df_plot_gamm_ms_split,df_roi_grp,analysis="ms_grp-group",atlas,list_mod,list_plot,
    #             list_type_p_,thr_p,paths_)
    #
  }
  print('Finished gam_fc_cs().')
}
