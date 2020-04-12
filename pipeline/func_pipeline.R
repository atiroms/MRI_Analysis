#**************************************************
# Description =====================================
#**************************************************
# R script to sequentially analize functional data.


#**************************************************
# Create path list ================================
#**************************************************
path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"
dir_in<-dir_out<-NULL

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
source(file.path(paths$script,"analyze/timeseries.R"))
source(file.path(paths$script,"analyze/connection.R"))
source(file.path(paths$script,"analyze/fingerprint.R"))


#**************************************************
# Parameters ======================================
#**************************************************
path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"

id_dir_start<-450
suffix_dir<-"acompcor"

#dir_summary<-"300_fp_model_summary"
dir_summary<-"500_fp_model_summary"

#list_term_summary<-c("diff_tanner","mean_tanner","s(diff_tanner)","s(mean_tanner)")
list_term_summary<-c("diff_tanner","mean_tanner","s(diff_tanner)","s(mean_tanner)")
thresh_sign<-0.05
#thresh_sign<-0.001

#list_id_dir<-list("acompcor"=202,"aroma"=212,"acompcor_gsr"=232,"aroma_gsr"=242)

#subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
#                             list("key"="W1_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)),
#                    "2"=list(list("key"="W2_T1QC","value"=1),
#                             list("key"="W2_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)))

#list_id_dir<-list("acompcor"=400,"acompcor_gsr"=410,"aroma"=420,"aroma_gsr"=430,"36p"=440)
list_id_dir<-list("aroma_gsr"=430)

#list_atlas<-c("aal116","glasser360","gordon333","power264",
#              "schaefer100","schaefer200","schaefer400","shen268")
#list_atlas<-c("aal116","gordon333","power264","shen268")
list_atlas<-c("aal116","gordon333","power264","schaefer400x7","shen268")

#list_atlas<-c("aal116","desikanKilliany","glasser360","gordon333","HarvardOxford","power264",
#              "schaefer100x7","schaefer100x17","schaefer200x7","schaefer200x17","schaefer400x7","schaefer400x17",
#              "shen268")

subset_subj <- list("1"=list(list("key"="W1_T1QC","condition"="==1"),
                             list("key"="W1_rsfMRIexist","condition"="==1"),
                             list("key"="W1_Censor","condition"="<126")),
#                             list("key"="W1_Censor","condition"="<151")),
                    "2"=list(list("key"="W2_T1QC","condition"="==1"),
                             list("key"="W2_rsfMRIexist","condition"="==1"),
                             list("key"="W2_Censor","condition"="<126")))
#                             list("key"="W2_Censor","condition"="<151")))
list_covar<-list("tanner"=list("1"="W1_Tanner_Max","2"="W2_Tanner_Max","label"="Tanner stage (max)"),
                 "age"   =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                 "sex"   =list("1"="Sex",          "2"="Sex",          "label"="Sex"))
list_mod <- list("ld" ="value ~ diff_age + diff_tanner",
                 "ldm"="value ~ diff_age + diff_tanner + mean_tanner",
                 "ad" ="value ~ s(diff_age,k=3) + s(diff_tanner,k=3)",
                 "adm"="value ~ s(diff_age,k=3) + s(mean_tanner,k=3) + s(diff_tanner,k=3)")
list_graph <-list("d(a)"=list("title"="Age diff effect","x_axis"="diff_age",
                               "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                             "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                               "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                            "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))),
                  "d(t)"=list("title"="Tanner diff effect","x_axis"="diff_tanner",
                               "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                             "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                               "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                            "Female"=list("subset"=list("sex"=2), "color"="lightcoral","alpha"=1))),
                  "m(t)"=list("title"="Tanner mean effect","x_axis"="mean_tanner",
                               "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                             "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                               "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                            "Female"=list("subset"=list("sex"=2), "color"="lightcoral","alpha"=1))))
list_strat_tanner <-list("5by5"=list("1"=list("1"=1,"2"=2,"3"=3,"4"=4,"5"=5),
                                     "2"=list("1"=1,"2"=2,"3"=3,"4"=4,"5"=5)),
                         "3by3"=list("1"=list("12"=c(1,2),"3"=3,"45"=c(4,5)),
                                     "2"=list("12"=c(1,2),"3"=3,"45"=c(4,5))),
                         "2by2"=list("1"=list("12"=c(1,2),"345"=c(3,4,5)),
                                     "2"=list("123"=c(1,2,3),"45"=c(4,5))))

list_type_tanner<-list("max"    =list("1"="W1_Tanner_Max", "2"="W2_Tanner_Max", "label"="Tanner stage (max)"),
                       "full"   =list("1"="W1_Tanner_Full","2"="W2_Tanner_Full","label"="Tanner stage (full)"),
                       "gonadal"=list("1"=c("W1_Tanner_Male_Genitals","W1_Tanner_Female_Breast"),
                                      "2"=c("W2_Tanner_Male_Genitals","W2_Tanner_Female_Breast"),
                                      "label"="Tanner stage (gonadal)"),
                       "adrenal"=list("1"=c("W1_Tanner_Male_Pubic_Hair","W1_Tanner_Female_Pubic_Hair"),
                                      "2"=c("W2_Tanner_Male_Pubic_Hair","W2_Tanner_Female_Pubic_Hair"),
                                      "label"="Tanner stage (adrenal)"))

# Parameters for gamm_multi_hormone()
#list_id_dir_fp<-list("acompcor"=202,"aroma"=212,"acompcor_gsr"=232,"aroma_gsr"=242)
#list_id_dir_fp<-list("acompcor"=302,"aroma"=312,"acompcor_gsr"=332,"aroma_gsr"=342)
list_covar_hormone<-list("hormone"=list("1"="W1_Hormone"   ,"2"="W2_Hormone",   "label"="Hormone"),
                         "age"    =list("1"="W1_Age_at_MRI","2"="W2_Age_at_MRI","label"="Age"),
                         "sex"    =list("1"="Sex",          "2"="Sex",          "label"="Sex"))
list_mod_hormone <- list("ld" ="value ~ diff_age + diff_hormone",
                         "ldm"="value ~ diff_age + diff_hormone+ mean_hormone",
                         "ad" ="value ~ s(diff_age,k=3) + s(diff_hormone,k=3)",
                         "adm"="value ~ s(diff_age,k=3) + s(mean_hormone,k=3) + s(diff_hormone,k=3)")
list_graph_hormone <-list("d(a)"=list("title"="Age diff effect","x_axis"="diff_age",
                                       "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                                     "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                                       "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                                    "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))),
                          "d(h)"=list("title"="Hormone diff effect","x_axis"="diff_hormone",
                                      "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                                    "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                                      "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                                   "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))),
                          "m(h)"=list("title"="Hormone mean effect","x_axis"="mean_hormone",
                                      "smooth"=list("Male"=list("fix"=list("sex"=1),"color"="steelblue2","alpha"=1,"ribbon"=T),
                                                    "Female"=list("fix"=list("sex"=2),"color"="lightcoral","alpha"=1,"ribbon"=T)),
                                      "point"=list("Male"=list("subset"=list("sex"=1),"color"="steelblue2","alpha"=1),
                                                   "Female"=list("subset"=list("sex"=2),"color"="lightcoral","alpha"=1))))
list_hormone<-list("testo"=list("1"="W1_Testosterone","2"="W2_Testosterone","label"="Testosterone"),
                   "corti"=list("1"="W1_Cortisol",    "2"="W2_Cortisol",    "label"="Cortisol"),
                   "dhea" =list("1"="W1_DHEA",        "2"="W2_DHEA",        "label"="DHEA"),
                   "dheas"=list("1"="W1_DHEAS",       "2"="W2_DHEAS",       "label"="DHEA-S"))

n_permutation<-1000
#n_permutation<-100


#**************************************************
# Libraries =======================================
#**************************************************


#**************************************************
# Summarize model_fp() results ====================
#**************************************************

sum_model<-function(dir_summary_=dir_summary,list_id_dir_=list_id_dir,
                          list_type_tanner_=list_type_tanner,
                          list_term_summary_=list_term_summary,thresh_sign_=thresh_sign
                          ){
  print("Starting sum_model().")
  
  df_out<-data.frame()
  flag_in_first<-T
  for (i_id_dir in seq(length(list_id_dir_))){
    for (i_type_tanner in seq(length(list_type_tanner_))){
      label_preproc<-names(list_id_dir_)[[i_id_dir]]
      label_type_tanner<-names(list_type_tanner_)[[i_type_tanner]]
      
      # Prepare directories
      id_dir_in<-as.character(list_id_dir[[i_id_dir]]+3+0.1*i_type_tanner)
      dir_in<-paste(id_dir_in,"fp_model",label_preproc,label_type_tanner,sep="_")
      paths<-func_path(dir_in=dir_in,dir_out=dir_summary_)
      if (flag_in_first){
        nullobj<-func_createdirs(paths,str_proc="model_fp()")
        flag_in_first<-F
      }
      
      path_dir_in<-file.path(paths$input,"output")
      path_dir_out<-file.path(paths$output,"output")
      
      # Load GLM/GAM results
      df_fp_glm<-read.csv(file.path(path_dir_in,"fp_glm.csv"))
      
      # All groups
      #df_sign<-df_fp_glm[df_fp_glm$term %in% list_term_summary_ & df_fp_glm$p<thresh_sign_,]
      # Whole-whole only
      df_sign<-df_fp_glm[df_fp_glm$term %in% list_term_summary_ & df_fp_glm$p<thresh_sign_
                         & df_fp_glm$group_1=="whole" & df_fp_glm$group_2=="whole",]
      
      df_out<-rbind(df_out,data.frame(preproc=label_preproc,tanner=label_type_tanner,df_fp_glm))
      
      if (dim(df_sign)[1]>0){
        #df_out<-rbind(df_out,data.frame(preproc=label_preproc,tanner=label_type_tanner,df_sign))
        for (i_row in seq(dim(df_sign)[1])){
          filename_src<-paste("atl-",df_sign[i_row,"atlas"],"_msr-",df_sign[i_row,"measure"],
                              "_grp1-",df_sign[i_row,"group_1"],"_grp2-",df_sign[i_row,"group_2"],
                              "_mod-",df_sign[i_row,"model"],"_plt-d(t)_fp_glm.eps",sep="")
          filename_dst<-paste("prc-",label_preproc,"_tnr-",label_type_tanner,"_",filename_src,sep="")
          if (!file.exists(file.path(path_dir_out,filename_dst))){
            file.copy(file.path(path_dir_in,filename_src),path_dir_out)
            file.rename(file.path(path_dir_out,filename_src),file.path(path_dir_out,filename_dst))
          }
        }
      }
    }
  }
  write.csv(df_out,file.path(path_dir_out,"sum_model_fp.csv"),row.names=F)
  print("Finished sum_model().")
}


#**************************************************
# Timeseries to GLM/GAM of Fingerprints ===========
#**************************************************
pipe_func<-function(id_dir_start_=id_dir_start,suffix_dir_=suffix_dir,list_atlas_=list_atlas,
                    list_wave_=list_wave,
                    list_covar_=list_covar,list_mod_=list_mod,list_graph_=list_graph,
                    list_strat_tanner_=list_strat_tanner,list_type_tanner_=list_type_tanner,
                    list_covar_hormone_=list_covar_hormone,list_mod_hormone_=list_mod_hormone,
                    list_graph_hormone_=list_graph_hormone,
                    list_hormone_=list_hormone,
                    subset_subj_=subset_subj,n_permutation_=n_permutation,
                    skip_ts2fc=FALSE,skip_fc2fp=FALSE,skip_fp2id=FALSE){
  
  print('Starting pipe_func().')
  
  id_dir_cnt<-id_dir_start_
  # Timeseries to functional connectivity
  if(!skip_ts2fc){
    dir_in<-paste(as.character(id_dir_cnt),"ts",suffix_dir_,sep='_')
    id_dir_cnt<-id_dir_cnt+1
    id_dir_fc<-id_dir_cnt
    dir_out<-paste(as.character(id_dir_fc),"fc",suffix_dir_,sep='_')
    paths<-func_path(dir_in_=dir_in,dir_out_=dir_out)
    nullobj<-fc(paths_=paths,list_atlas_=list_atlas_)
  }else{
    id_dir_fc<-id_dir_cnt
  }
  
  # Functional connetivity to fingerprint
  if(!skip_fc2fp){
    dir_in<-paste(as.character(id_dir_fc),"fc",suffix_dir_,sep='_')
    id_dir_cnt<-id_dir_cnt+1
    id_dir_fp<-id_dir_cnt
    dir_out<-paste(as.character(id_dir_fp),"fp",suffix_dir_,sep='_')
    paths<-func_path(dir_in_=dir_in,dir_out_=dir_out)
    nullobj<-fp_fc(paths_=paths,list_atlas_=list_atlas_)
  }else{
    id_dir_fp<-id_dir_cnt
  }
  
  # Fingerprint to identification of fingerprints
  if(!skip_fp2id){
    dir_in<-paste(as.character(id_dir_fp),"fp",suffix_dir_,sep='_')
    id_dir_cnt<-id_dir_cnt+1
    id_dir_idfp<-id_dir_cnt
    dir_out<-paste(as.character(id_dir_idfp),"fp_id",suffix_dir_,sep='_')
    paths<-func_path(dir_in_=dir_in,dir_out_=dir_out)
    nullobj<-identify_fp(paths_=paths,list_atlas_=list_atlas_,
                         list_wave_=list_wave_,
                         subset_subj_=subset_subj_,
                         n_permutation_=n_permutation_)
  }
  
  # Fingerprint to GLM / ANCOVA of fingerprint difference
  # #1 Tanner stage
  dir_in<-paste(as.character(id_dir_fp),"fp",suffix_dir_,sep='_')
  id_dir_cnt<-id_dir_cnt+1
  #id_dir_cnt<-id_dir_cnt+3
  id_dir_model_fp<-id_dir_cnt
  for (idx_type_tanner in names(list_type_tanner_)){
    id_dir_model_fp<-id_dir_model_fp+0.1
    print(paste("Tanner type: ",list_type_tanner_[[idx_type_tanner]][["label"]],sep=""))
    list_covar_[["tanner"]]<-list_type_tanner[[idx_type_tanner]]
    dir_out<-paste(as.character(id_dir_model_fp),"fp_model",suffix_dir_,idx_type_tanner,sep='_')
    paths<-func_path(dir_in_=dir_in,dir_out_=dir_out)
    nullobj<-model_fp(paths_=paths,list_atlas_=list_atlas_,
                      list_wave_=list_wave_,list_covar_=list_covar_,
                      list_mod_=list_mod_,list_graph_=list_graph_,list_strat_tanner=list_strat_tanner_,
                      #subset_subj_=subset_subj_,skip_ancova=F)
                      subset_subj_=subset_subj_,skip_ancova=T)
  }
  
  # #2 Hormone
  dir_in<-paste(as.character(id_dir_fp),"fp",suffix_dir_,sep='_')
  for (idx_hormone in names(list_hormone_)){
    id_dir_model_fp<-id_dir_model_fp+0.1
    print(paste("Hormone: ",list_hormone[[idx_hormone]][["label"]],sep=""))
    list_covar_hormone_[["hormone"]]<-list_hormone[[idx_hormone]]
    dir_out<-paste(as.character(id_dir_model_fp),"fp_model",suffix_dir_,idx_hormone,sep='_')
    paths<-func_path(dir_in_=dir_in,dir_out_=dir_out)
    nullobj<-model_fp(paths_=paths,list_atlas_=list_atlas_,
                      list_wave_=list_wave_,list_covar_=list_covar_hormone_,
                      list_mod_=list_mod_hormone_,list_graph_=list_graph_hormone_,
                      list_strat_tanner=NULL,
                      subset_subj_=subset_subj_,skip_ancova=T)
  }
  
  print('Finished pipe_func().')
}


#**************************************************
# pipe_func() with multiple starting FCs ==========
#**************************************************
pipe_func_multi<-function(list_id_dir_=list_id_dir,
                          list_atlas_=list_atlas,
                          list_wave_=list_wave,list_covar_=list_covar,
                          list_mod_=list_mod,list_graph_=list_graph,
                          list_strat_tanner_=list_strat_tanner,list_type_tanner_=list_type_tanner,
                          subset_subj_=subset_subj,n_permutation_=n_permutation,
                          skip_ts2fc=FALSE,skip_fc2fp=FALSE,skip_fp2id=TRUE){
  
  print("Starting pipe_func_multi()")
  for (suffix_dir in names(list_id_dir_)){
    id_dir_start<-list_id_dir_[[suffix_dir]]
    nullobj<-pipe_func(id_dir_start_=id_dir_start,suffix_dir_=suffix_dir,list_atlas_=list_atlas_,
                       list_wave_=list_wave_,list_covar_=list_covar_,
                       list_mod_=list_mod_,list_graph_=list_graph_,
                       list_strat_tanner_=list_strat_tanner_,list_type_tanner_=list_type_tanner_,
                       subset_subj_=subset_subj_,n_permutation_=n_permutation_,
                       skip_ts2fc=skip_ts2fc,skip_fc2fp=skip_fc2fp,skip_fp2id=skip_fp2id)
  }
  print("Finished pipe_func_multi()")
}


#**************************************************
# Iterate model_fp() over Tanner types ============
#**************************************************
# OBSOLETE
itr_model_fp<-function(id_dir_fp_=id_dir_fp,id_dir_cnt_=id_dir_fp+1,
                       suffix_dir_=suffix_dir,list_atlas_=list_atlas,
                       list_wave_=list_wave,list_covar_=list_covar,
                       list_mod_=list_mod,list_graph_=list_graph,
                       list_strat_tanner_=list_strat_tanner,list_type_tanner_=list_type_tanner,
                       subset_subj_=subset_subj){
  
  print("Starting itr_model_fp()")
  # Fingerprint to GLM / ANCOVA of fingerprint difference
  dir_in<-paste(as.character(id_dir_fp_),"fp",suffix_dir_,sep='_')
  id_dir_cnt<-id_dir_cnt_+1
  id_dir_model_fp<-id_dir_cnt
  for (idx_type_tanner in names(list_type_tanner_)){
    id_dir_model_fp<-id_dir_model_fp+0.1
    print(paste("Tanner type: ",idx_type_tanner,sep=""))
    list_covar_[["tanner"]]<-list_type_tanner[[idx_type_tanner]]
    dir_out<-paste(as.character(id_dir_model_fp),"fp_model",suffix_dir_,idx_type_tanner,sep='_')
    paths<-func_path(dir_in_=dir_in,dir_out_=dir_out)
    nullobj<-model_fp(paths_=paths,list_atlas_=list_atlas_,
                      list_wave_=list_wave_,list_covar_=list_covar_,
                      list_mod_=list_mod_,list_graph_=list_graph_,list_strat_tanner=list_strat_tanner_,
                      subset_subj_=subset_subj_)
  }
  print("Finished itr_model_fp()")
}


#**************************************************
# gamm_fp() for hormonal data =====================
#**************************************************
# OBSOLETE
gamm_multi_hormone<-function(list_id_dir_=list_id_dir_fp,
                             list_atlas_=list_atlas,
                             list_wave_=list_wave,list_covar_=list_covar_hormone,
                             list_mod_=list_mod_hormone,list_graph_=list_graph_hormone,
                             list_hormone_=list_hormone,
                             subset_subj_=subset_subj){
  
  print("Starting gamm_multi_hormone()")
  for (suffix_dir in names(list_id_dir_)){
    id_dir_fp<-list_id_dir_[[suffix_dir]]
    dir_in<-paste(as.character(id_dir_fp),"fp",suffix_dir,sep='_')
    id_dir_model_fp<-id_dir_fp+2.4
    for (idx_hormone in names(list_hormone_)){
      id_dir_model_fp<-id_dir_model_fp+0.1
      print(paste("Preproc: ",suffix_dir,", Hormone: ",list_hormone[[idx_hormone]][["label"]],sep=""))
      list_covar_[["hormone"]]<-list_hormone[[idx_hormone]]
      dir_out<-paste(as.character(id_dir_model_fp),"fp_model",suffix_dir,idx_hormone,sep='_')
      paths<-func_path(dir_in_=dir_in,dir_out_=dir_out)
      nullobj<-model_fp(paths_=paths,list_atlas_=list_atlas_,
                        list_wave_=list_wave_,list_covar_=list_covar_,
                        list_mod_=list_mod_,list_graph_=list_graph_,list_strat_tanner=NULL,
                        subset_subj_=subset_subj_,skip_ancova=T)
    }
  }
  print("Finished gamm_multi_hormone()")
}