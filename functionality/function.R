#**************************************************
# Description =====================================
#**************************************************

# R script for common MRI analysis functions


#**************************************************
# Libraries =======================================
#**************************************************
library(tidyverse)
library(dplyr)

#**************************************************
# Factor to numeric function ======================
#**************************************************
as.numeric.factor <- function(x) {
  if (class(x)=="factor"){
    return(as.numeric(levels(x))[x])
  }else{
    return(x)
  }
}


#**************************************************
# Experiment folder preparation ===================
#**************************************************
func_createdirs<-function(paths,copy_log=T){
  list_createdirs<-c(paths$output,file.path(paths$output,"output"))
  for(d in list_createdirs){
    if (!file.exists(d)){
      dir.create(d)
    }
  }
  if (copy_log){
    file.copy(file.path(paths$input,"log"),paths$output,recursive=T)
  }
}


#**************************************************
# Returns ROI dictionary ==========================
#**************************************************
func_dict_roi<-function(paths,
                        file_roi="ROI.csv"){
  output<-read.csv(file.path(paths$common,file_roi))
  return(output)
}


#**************************************************
# Clinical data loading ===========================
#**************************************************
func_clinical_data<-function(paths,
                             subset_subj,
                             file_clinical= "CSUB.csv"
                             ){
  df_clinical <- read.csv(file.path(paths$common,file_clinical))
  for (list_subset in subset_subj){
    df_clinical <- df_clinical[which(df_clinical[,list_subset[["column"]]]==list_subset[["value"]]),]
  }
  list_id_subj<-df_clinical$ID_pnTTC
  df_id_subj<-df_clinical[,1:(which(colnames(df_clinical)=="Clinical")-1)]
  df_clinical<-df_clinical[,(-2):(-which(colnames(df_clinical)=="Clinical"))]
  n_subj<-length(list_id_subj)
  n_data_clinical<-ncol(df_clinical)-1
  
  output<-list("df_clinical"=df_clinical,"list_id_subj"=list_id_subj,"df_id_subj"=df_id_subj,
               "n_subj"=n_subj,"n_data_clinical"=n_data_clinical)
  return(output)
}

#**************************************************
# Longitudinal clinical data loading ==============
#**************************************************
func_clinical_data_long<-function(paths,
                                  list_wave,
                                  file_clin= "CSUB.csv"
                                  ){
  df_src <- read.csv(file.path(paths$common,file_clin))
  df_dst <- data.frame(matrix(nrow=0,ncol=ncol(df_src)+1))
  for (wave in list_wave){
    df_tmp<-df_src['ID_pnTTC']
    df_tmp$wave<-rep(wave,nrow(df_tmp))
    df_tmp<-cbind(df_tmp,df_src[,colnames(df_src)!='ID_pnTTC'])
    df_dst<-rbind(df_dst,df_tmp)
  }
  colnames(df_dst)<-c('ID_pnTTC','wave',colnames(df_src)[colnames(df_src)!='ID_pnTTC'])
  return(df_dst)
}


#**************************************************
# Subset clinical data  ===========================
#**************************************************
# used for longitudinal analysis (GAMM)

func_subset_clin<-function(df_clin,
                           list_wave,subset_subj,
                           list_covar,
                           rem_na_clin
                           ){
  df_clin_subset<-data.frame(matrix(nrow=0,ncol=ncol(df_clin)))
  
  # Subset clinical data according to subsetting condition
  print('Subsetting clinical data frame according to specified condition.')
  list_id_subset<-list()
  for (wave in list_wave){
    str_wave<-as.character(wave)
    print(paste('Checking wave ',str_wave,sep=''))
    df_clin_wave<-df_clin[df_clin['wave']==wave,]
    print(paste('Source clinical data ',as.character(nrow(df_clin_wave)),sep=''))
    id_intersect<-df_clin_wave[,'ID_pnTTC']
    list_id_subset_wave<-list("src"=id_intersect)
    for (cond_subset in subset_subj[[str_wave]]){
      key_subset<-cond_subset$key
      value_subset<-cond_subset$value
      id_meet_cond<-df_clin_wave[df_clin_wave[key_subset]==value_subset,'ID_pnTTC']
      id_intersect<-intersect(id_intersect,id_meet_cond)
      print(paste(as.character(length(id_meet_cond)),' subjects meeting ',key_subset, ' = ',as.character(value_subset),sep=''))
      id_meet_cond<-list(id_meet_cond)
      names(id_meet_cond)<-key_subset
      list_id_subset_wave<-c(list_id_subset_wave,id_meet_cond)
    }
    print(paste(as.character(length(id_intersect)),' subjects meeting all conditions.',sep=''))
    df_clin_wave<-df_clin_wave[df_clin_wave$ID_pnTTC %in% id_intersect,]
    df_clin_subset<-rbind(df_clin_subset,df_clin_wave)
    id_intersect<-list('intersect'=id_intersect)
    list_id_subset_wave<-c(list_id_subset_wave,id_intersect)
    list_id_subset_wave<-list(list_id_subset_wave)
    names(list_id_subset_wave)<-str_wave
    list_id_subset<-c(list_id_subset,list_id_subset_wave)
  }
  colnames(df_clin_subset)<-colnames(df_clin)
  
  # Subset unused columns of clinical data
  # Subset clinical data with missing covariate value
  print('Subsetting clinical data with missing covariate value.')
  df_clin_exist<-data.frame(matrix(nrow=0,ncol=length(list_covar)+2))
  list_id_exist<-list()
  for (wave in list_wave){
    str_wave<-as.character(wave)
    print(paste('Checking wave ',str_wave,sep=''))
    df_clin_exist_wave<-df_clin_subset[df_clin_subset$wave==wave,c('ID_pnTTC','wave')]
    id_exist_intersect<-df_clin_exist_wave$ID_pnTTC
    print(paste('Source clinical data ',as.character(length(id_exist_intersect)),sep=''))
    list_id_exist_wave<-list('src'=id_exist_intersect)
    if (length(list_covar)>0){
      for (id_covar in seq(length(list_covar))){
        name_covar_src<-list_covar[[id_covar]][[str_wave]]
        name_covar_dst<-names(list_covar)[id_covar]
        df_clin_exist_wave_add<-df_clin_subset[df_clin_subset$wave==wave,name_covar_src,drop=F]
        colnames(df_clin_exist_wave_add)<-name_covar_dst
        df_clin_exist_wave<-cbind(df_clin_exist_wave,df_clin_exist_wave_add)
        id_exist<-df_clin_exist_wave[!is.na(df_clin_exist_wave[name_covar_dst]),'ID_pnTTC']
        id_exist_intersect<-intersect(id_exist_intersect,id_exist)
        print(paste(as.character(length(id_exist)),' subjects with ',name_covar_src,sep=''))
        id_exist<-list(id_exist)
        names(id_exist)<-name_covar_dst
        list_id_exist_wave<-c(list_id_exist_wave,id_exist)
      }
      print(paste(as.character(length(id_exist_intersect)),' subjects with all convariates existing.',sep=''))
    }
    if (rem_na_clin){
      n_subj_pre<-nrow(df_clin_exist_wave)
      df_clin_exist_wave<-df_clin_exist_wave[df_clin_exist_wave$ID_pnTTC %in% id_exist_intersect,]
      n_subj_deleted<-n_subj_pre-nrow(df_clin_exist_wave)
      print(paste('Deleted ',as.character(n_subj_deleted),' subjects with missing covariate.',sep=''))
    }else{
      print('Did not delete subjects with missing covariate.')
    }
    list_id_exist_wave<-c(list_id_exist_wave,list('intersect'=id_exist_intersect))
    df_clin_exist<-rbind(df_clin_exist,df_clin_exist_wave)
    list_id_exist_wave<-list(list_id_exist_wave)
    names(list_id_exist_wave)<-str_wave
    list_id_exist<-c(list_id_exist,list_id_exist_wave)
  }
  colnames(df_clin_exist)<-c('ID_pnTTC','wave',names(list_covar))
  
  output<-list('df_clin'=df_clin_exist,'list_id_subset'=list_id_subset,'list_id_exist'=list_id_exist)
  return(output)
}


#**************************************************
# Subsetting Structural data  =====================
#**************************************************
func_subset_str<-function(df_str,
                          list_measure,
                          list_str_group,
                          dict_roi){
  
  df_str_dst<-data.frame(matrix(nrow=0,ncol=ncol(df_str)))
  for (measure in list_measure){
    print(paste('Subsetting ',measure,' measures.',sep=''))
    df_str_tmp<-df_str[df_str$measure==measure,]
    list_roi<-unique(df_str_tmp$roi)
    list_roi<-list_roi[order(list_roi)]
    print(paste(as.character(length(list_roi)),' ROIs exist in the source data.',sep=''))
    list_roi_subset<-NULL
    for (roi in list_roi){
      group_roi<-dict_roi[dict_roi$id==roi,'group']
      if (group_roi %in% list_str_group){
        list_roi_subset<-c(list_roi_subset,roi)
      }
    }
    print(paste(as.character(length(list_roi_subset)),' ROIs remaining.',sep=''))
    df_str_tmp<-df_str_tmp[df_str_tmp$roi %in% list_roi_subset,]
    df_str_dst<-rbind(df_str_dst,df_str_tmp)
  }
  colnames(df_str_dst)<-colnames(df_str)
  return(df_str_dst)
}


#**************************************************
# General correlation calculation =================
#**************************************************
func_cor<-function(input){
  cor <-rcorr(as.matrix(input), type="pearson")
  n_node<-ncol(input)
  cor_flat<-data.frame(matrix(nrow=n_node*(n_node-1)/2,ncol=4))
  colnames(cor_flat)<-c("from","to","r","p")
  k<-0
  for (i in 1:(n_node-1)){
    for (j in (i+1):n_node){
      k<-k+1
      cor_flat[k,1:4]<-c(rownames(cor$r)[i],
                          colnames(cor$r)[j],
                          cor$r[i,j],
                          cor$P[i,j])
    }
  }
  output<-list("cor"=cor, "cor_flat"=cor_flat)
  return(output)
}


#**************************************************
# Multiple comparison correction of p values ======
#**************************************************

mltcomp_corr<-function(input){
  output<-data.frame("p_bonferroni"=p.adjust(input$p,method = "bonferroni"),
                     "p_holm_bonferroni"=p.adjust(input$p,method = "holm"),
                     "p_hockberg"=p.adjust(input$p,method = "holm"),
                     "p_hommel"=p.adjust(input$p,method = "hommel"),
                     "p_benjamini_hochberg"=p.adjust(input$p,method="BH"),
                     "p_benjamini_yukutieli"=p.adjust(input$p,method="BY"))
  return(output)
}
