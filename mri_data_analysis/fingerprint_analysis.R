#**************************************************
# Description =====================================
#**************************************************
# R script to analyze fingerprint data.


#**************************************************
# Parameters ======================================
#**************************************************
# parameters for gta_bin() and gta_weight()
path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP"
#path_exp <- "Dropbox/MRI/pnTTC/Puberty/Stats/func_XCP/test_5sub"

dir_in<-"55_fingerprint"
#dir_out<-"55_gta_bin"
dir_out<-"58_fp_identification"

list_wave <- c(1,2)

subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
                             list("key"="W1_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)),
                    "2"=list(list("key"="W2_T1QC","value"=1),
                             list("key"="W2_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)))

list_mod <- list("a+s+st"=
                   "value ~ s(age,k=3) + sex + s(tanner,k=3,by=sex) + s(ID_pnTTC,bs='re')")

#list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")
list_atlas<-"aal116"
#list_atlas<-"schaefer400"
#list_atlas<-c("glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")
#thr_pvalue <- 0.05
#n_permutation<-1000
n_permutation<-100


#**************************************************
# Libraries =======================================
#**************************************************
library(dplyr)
#library(ggplot2)
#library(GGally)
#library(igraph)
#library(qgraph)


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
source(file.path(paths$script,"functionality/function.R"))
#source(file.path(paths$script,"functionality/glm_function.R"))
source(file.path(paths$script,"functionality/graph.R"))
#source(file.path(paths$script,"functionality/gta_function.R"))


#**************************************************
# GAMM of Fingerprint ==============================
#**************************************************
gamm_fp<-function(paths_=paths,
                  list_atlas_=list_atlas,
                  list_wave_=list_wave,
                  list_covar_=list_covar,
                  subset_subj_=subset_subj
                  ){
  print("Starting glm_fp().")
  nullobj<-func_createdirs(paths_)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,
                                     list_covar=list_covar_,rem_na_clin=F)
  df_clin<-data_clin$df_clin
  colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
  
  for (atlas in list_atlas_){
    # Load fingerprint data
    print(paste("Loading atlas: ",atlas,sep=""))
    df_fp<-read.csv(file.path(paths_$input,"output",paste("atl-",atlas,"_fingerprint.csv",sep="")))
    
    # Create list of subjects who meet subsetting condition and whose MRI data exist
    list_ses_exist <- sort(unique(c(df_fp$from_ses,df_fp$to_ses)))
    list_id_subj_exist<-list()
    for (ses in list_ses_exist){
      id_subj_exist_ses<-sort(unique(c(df_fp[df_fp$from_ses==ses,'from_ID_pnTTC'],
                                       df_fp[df_fp$to_ses==ses,'to_ID_pnTTC'])))
      id_subj_subset_ses<-df_clin[df_clin$ses==ses,"ID_pnTTC"]
      id_subj_exist_ses<-intersect(id_subj_exist_ses,id_subj_subset_ses)
      list_id_subj_exist[[as.character(ses)]]<-sort(id_subj_exist_ses)
    }
    
    # Identify those with longitudinal data
    list_id_subj_exist_twice<-sort(intersect(list_id_subj_exist[["1"]],
                                             list_id_subj_exist[["2"]]))
    n_id_subj_exist_twice<-length(list_id_subj_exist_twice)
    print(paste(as.character(n_id_subj_exist_twice)," subjects with two sessions.",sep=""))

    # Create dataframe for GLM analysis
    
    
    
    
  }
  
  
}


#**************************************************
# Fingerprint identification ======================
#**************************************************
identify_fp<-function(paths_=paths,
                      list_atlas_=list_atlas,
                      list_wave_=list_wave,
                      #list_covar_=list_covar,
                      subset_subj_=subset_subj,
                      n_permutation_=n_permutation
                      ){
  print("Starting identify_fp().")
  nullobj<-func_createdirs(paths_)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,list_covar=NULL,rem_na_clin=F)
  df_clin<-data_clin$df_clin
  colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
  
  for (atlas in list_atlas_){
    # Load fingerprint data
    print(paste("Loading atlas: ",atlas,sep=""))
    df_fp<-read.csv(file.path(paths_$input,"output",paste("atl-",atlas,"_fingerprint.csv",sep="")))
    
    # Create list of subjects who meet subsetting condition and whose MRI data exist
    list_ses_exist <- sort(unique(c(df_fp$from_ses,df_fp$to_ses)))
    list_id_subj_exist<-list()
    for (ses in list_ses_exist){
      id_subj_exist_ses<-sort(unique(c(df_fp[df_fp$from_ses==ses,'from_ID_pnTTC'],
                                       df_fp[df_fp$to_ses==ses,'to_ID_pnTTC'])))
      id_subj_subset_ses<-df_clin[df_clin$ses==ses,"ID_pnTTC"]
      id_subj_exist_ses<-intersect(id_subj_exist_ses,id_subj_subset_ses)
      list_id_subj_exist[[as.character(ses)]]<-sort(id_subj_exist_ses)
    }
    
    # Extract those with longitudinal data
    list_id_subj_exist_twice<-sort(intersect(list_id_subj_exist[["1"]],
                                             list_id_subj_exist[["2"]]))
    n_id_subj_exist_twice<-length(list_id_subj_exist_twice)
    print(paste(as.character(n_id_subj_exist_twice)," subjects with two sessions.",sep=""))
    df_fp_exist_twice<-df_fp[(df_fp$from_ID_pnTTC %in% list_id_subj_exist_twice) & (df_fp$to_ID_pnTTC %in% list_id_subj_exist_twice),]
    df_fp_exist_twice<-df_fp_exist_twice[(df_fp_exist_twice$from_ses==1 & df_fp_exist_twice$to_ses==2),]

    # Calculate fingerprint identification
    df_ident<-data.frame("target"=list_id_subj_exist_twice)
    df_perm<-data.frame("id_perm"=seq(1,n_permutation_))
    for(ses in c(1,2)){
      if(ses==1){
        df_fp_pool<-df_fp_exist_twice[c("from_ses","from_ID_pnTTC","to_ses","to_ID_pnTTC","r")]
      }else{
        df_fp_pool<-df_fp_exist_twice[c("to_ses","to_ID_pnTTC","from_ses","from_ID_pnTTC","r")]
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
        df_ident[df_ident$target==id_subj,paste(as.character(ses),"_targeted_rank",sep='')]<-rank_similarity
        if (!is.na(rank_similarity)){
          if (rank_similarity==1){
            df_ident[df_ident$target==id_subj,paste(as.character(ses),"_targeted_identification",sep='')]<-1
          }
        }
      }
      
      print("Starting permutation test.")
      for (i in seq(1,n_permutation_)){
        df_rand<-data.frame(pool_ID_pnTTC=list_id_subj_exist_twice,rand_ID_pnTTC=sample(list_id_subj_exist_twice,length(list_id_subj_exist_twice)))
        df_fp_rand<-left_join(df_fp_subset,df_rand,by="pool_ID_pnTTC")
        n_identified<-0
        for (id_subj in list_id_subj_exist_twice){
          df_fp_rand<-df_fp_rand[df_fp_rand$target_ID_pnTTC==id_subj,]
          list_id_subj_ordered<-df_fp_subset[order(df_fp_subset$r,decreasing = TRUE,na.last=NA),'pool_ID_pnTTC']
          if (id_subj %in% list_id_subj_ordered){
            rank_similarity<-which(list_id_subj_ordered==id_subj)
          }else{
            rank_similarity<-NA
          }
          if (!is.na(rank_similarity)){
            if (rank_similarity==1){
              n_dentified<-n_identified+1
            }
          }
        }
        df_perm[i,paste(as.character(ses),'_n_ident',sep='')]<-n_identified
        #print(paste("Iteration: ",as.character(i),", subjects identified: ",as.character(n_identified),sep=""))
      }
    }
    df_ident[is.na(df_ident["1_targeted_identification"]),"1_targeted_identification"]<-0
    df_ident[is.na(df_ident["2_targeted_identification"]),"2_targeted_identification"]<-0
    write.csv(df_ident,file.path(paths_$output,"output",paste("atl-",atlas,"_fp_identification.csv",sep="")),row.names=F)
    print(paste("Number of session 2 subjects identified from session 1 target: ",as.character(sum(df_ident["1_targeted_identification"])),sep=""))
    print(paste("Number of session 1 subjects identified from session 2 target: ",as.character(sum(df_ident["2_targeted_identification"])),sep=""))
    write.csv(df_perm,file.path(paths_$output,"output",paste("atl-",atlas,"_fp_permutation.csv",sep="")),row.names=F)
    p_perm_1<-sum(df_perm["1_n_ident"]>sum(df_ident["1_targeted_identification"]))/n_permutation_
    print(paste("p value from permutation 1: ",as.character(p_perm_1)))
    p_perm_2<-sum(df_perm["2_n_ident"]>sum(df_ident["2_targeted_identification"]))/n_permutation_
    print(paste("p value from permutation 2: ",as.character(p_perm_2)))
  }
  print("Finished identify_fp().")
}
