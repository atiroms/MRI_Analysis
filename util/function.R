#**************************************************
# Description =====================================
#**************************************************

# R script for common MRI analysis functions


#**************************************************
# Libraries =======================================
#**************************************************
library(tidyverse)
library(dplyr)
library(Hmisc)

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
    if (file.exists(file.path(paths$input,"log"))){
      file.copy(file.path(paths$input,"log"),paths$output,recursive=T)
    }else{
      print("Log folder does not exist.")
    }
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
                                  subset_subj,
                                  list_covar,
                                  rem_na_clin,
                                  file_clin= "CSUB.csv"
                                  ){
  df_src <- read.csv(file.path(paths$common,file_clin))
  df_clin_long <- data.frame(matrix(nrow=0,ncol=ncol(df_src)+1))
  for (wave in list_wave){
    df_tmp<-df_src['ID_pnTTC']
    df_tmp$wave<-rep(wave,nrow(df_tmp))
    df_tmp<-cbind(df_tmp,df_src[,colnames(df_src)!='ID_pnTTC'])
    df_clin_long<-rbind(df_clin_long,df_tmp)
  }
  colnames(df_clin_long)<-c('ID_pnTTC','wave',colnames(df_src)[colnames(df_src)!='ID_pnTTC'])
  #return(df_clin_long)
  
  
  df_clin_subset<-data.frame(matrix(nrow=0,ncol=ncol(df_clin_long)))
  
  # Subset clinical data according to subsetting condition
  print('Subsetting clinical data frame according to specified condition.')
  list_id_subset<-list()
  for (wave in list_wave){
    str_wave<-as.character(wave)
    print(paste('Checking wave ',str_wave,sep=''))
    df_clin_wave<-df_clin_long[df_clin_long['wave']==wave,]
    print(paste(as.character(nrow(df_clin_wave)),' source clinical data identified',sep=''))
    id_intersect<-df_clin_wave[,'ID_pnTTC']
    list_id_subset_wave<-list("src"=id_intersect)
    for (cond_subset in subset_subj[[str_wave]]){
      key_subset<-cond_subset$key
      value_subset<-cond_subset$value
      id_meet_cond<-df_clin_wave[df_clin_wave[key_subset]==value_subset,'ID_pnTTC']
      id_meet_cond<-id_meet_cond[!is.na(id_meet_cond)]
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
  colnames(df_clin_subset)<-colnames(df_clin_long)
  
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
    print(paste(as.character(length(id_exist_intersect)),' source clinical data identified.',sep=''))
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
        print(paste(as.character(length(id_exist)),' subjects with non-NA values of covariate: ',name_covar_src,sep=''))
        id_exist<-list(id_exist)
        names(id_exist)<-name_covar_dst
        list_id_exist_wave<-c(list_id_exist_wave,id_exist)
      }
      print(paste(as.character(length(id_exist_intersect)),' subjects with non-NA values of all convariates.',sep=''))
    }
    n_subj_pre<-nrow(df_clin_exist_wave)
    if (rem_na_clin){
      df_clin_exist_wave<-df_clin_exist_wave[df_clin_exist_wave$ID_pnTTC %in% id_exist_intersect,]
      n_subj_deleted<-n_subj_pre-nrow(df_clin_exist_wave)
      print(paste(as.character(n_subj_deleted),' subjects with NA values in any covariates deleted.',sep=''))
    }else{
      n_subj_na<-n_subj_pre-nrow(df_clin_exist_wave[df_clin_exist_wave$ID_pnTTC %in% id_exist_intersect,])
      print(paste(as.character(n_subj_na),' subjects with NA values in any covariates, but NOT deleted.',sep=''))
    }
    list_id_exist_wave<-c(list_id_exist_wave,list('intersect'=id_exist_intersect))
    df_clin_exist<-rbind(df_clin_exist,df_clin_exist_wave)
    list_id_exist_wave<-list(list_id_exist_wave)
    names(list_id_exist_wave)<-str_wave
    list_id_exist<-c(list_id_exist,list_id_exist_wave)
  }
  colnames(df_clin_exist)<-c('ID_pnTTC','wave',names(list_covar))
  rownames(df_clin_exist)<-NULL
  
  output<-list('df_clin'=df_clin_exist,'list_id_subset'=list_id_subset,'list_id_exist'=list_id_exist)
  return(output)
}


#**************************************************
# Joining longitudinal clinical data ==============
#**************************************************
func_clinical_data_join<-function(df_src,list_id_subj,list_covar){
  
  # Classify covaiates to fixed values and unfixed values
  list_covar_fix<-list_covar_change<-NULL
  for (id_covar in names(list_covar)){
    if (list_covar[[id_covar]][["1"]]==list_covar[[id_covar]][["2"]]){
      list_covar_fix<-c(list_covar_fix,id_covar)
    }else{
      list_covar_change<-c(list_covar_change,id_covar)
    }
  }
  
  # Prepare fixed covariates
  df_clin_fix<-data.frame(ID_pnTTC=list_id_subj)
  for (id_covar in list_covar_fix){
    df_clin_fix[[id_covar]]<-df_src[df_src$ID_pnTTC %in% list_id_subj & df_src$ses==1,id_covar]
  }
  
  # Prepare unfixed values
  df_clin_change<-list()
  for (ses in c(1,2)){
    df_clin_change_ses<-df_src[df_src$ses==ses & df_src$ID_pnTTC %in% list_id_subj, list_covar_change]
    df_clin_change_ses<-list(df_clin_change_ses)
    names(df_clin_change_ses)<-as.character(ses)
    df_clin_change<-c(df_clin_change,df_clin_change_ses)
  }
  
  # Calculate difference and mean
  df_clin_diff<-df_clin_change[["2"]]-df_clin_change[["1"]]
  df_clin_mean<-(df_clin_change[["1"]]+df_clin_change[["2"]])/2
  
  # Change column names
  colnames(df_clin_change[["1"]])<-c(paste("ses1_",colnames(df_clin_change[["1"]]),sep=''))
  colnames(df_clin_change[["2"]])<-c(paste("ses2_",colnames(df_clin_change[["2"]]),sep=''))
  colnames(df_clin_diff)<-c(paste("diff_",colnames(df_clin_diff),sep=''))
  colnames(df_clin_mean)<-c(paste("mean_",colnames(df_clin_mean),sep=''))
  
  # Join dataframes
  df_join<-cbind(df_clin_fix,df_clin_change[["1"]],df_clin_change[["2"]],df_clin_diff,df_clin_mean)
  
  return(df_join)
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
  mean_cor<-mean(cor$r,na.rm=TRUE)
  sd_cor<-sd(cor$r,na.rm=TRUE)
  cor_flat$z_r<-(as.numeric(cor_flat$r)-mean_cor)/sd_cor
  output<-list("cor"=data.frame(cor$r),"cor_flat"=cor_flat)
  return(output)
}


#**************************************************
# General PCA calculation =========================
#**************************************************
func_pca<-function(df_src,df_var=NULL,df_indiv=NULL){
  # Estimate number of dimensions
  print("Estimating PCA dimension.")
  ncp_estimate<-estim_ncpPCA(df_src,ncp.max=ncol(df_src))$ncp
  print(paste("PCA dimension: ",as.character(ncp_estimate),sep=""))
  # Ismpute data
  df_conn<-imputePCA(df_src,ncp=ncp_estimate)$completeObs
  
  # PCA calculation
  data_pca<-PCA(df_conn,scale.unit = TRUE, ncp = ncp_estimate, graph = FALSE)
  
  df_fac_var<-data.frame(data_pca$var$coord)
  if(!is.null(df_var)){
    df_fac_var<-cbind(df_var,df_fac_var)
    colnames(df_fac_var)<-c(colnames(df_var),sprintf("dim_%02d",1:ncp_estimate))
  }else{
    colnames(df_fac_var)<-sprintf("dim_%02d",1:ncp_estimate)
  }
  rownames(df_fac_var)<-NULL
  
  df_fac_indiv<-data.frame(data_pca$ind$coord)
  if(!is.null(df_indiv)){
    df_fac_indiv<-cbind(df_indiv,df_fac_indiv)
    colnames(df_fac_indiv)<-c(colnames(df_indiv),sprintf("dim_%02d",1:ncp_estimate))
  }else{
    colnames(df_fac_indiv)<-sprintf("dim_%02d",1:ncp_estimate)
  }
  rownames(df_fac_indiv)<-NULL
  
  df_var_accounted<-data.frame(data_pca$eig)
  colnames(df_var_accounted)<-c("eigenvalue","var_accounted","cumul_var_accounted")
  df_var_accounted$var_accounted<-df_var_accounted$var_accounted/100
  df_var_accounted$cumul_var_accounted<-df_var_accounted$cumul_var_accounted/100
  df_var_accounted$dim<-seq(1,dim(df_var_accounted)[1])
  df_var_accounted<-df_var_accounted[c("dim","var_accounted","cumul_var_accounted","eigenvalue")]
  rownames(df_var_accounted)<-NULL
  
  return(list('df_fac_var'=df_fac_var,'df_fac_indiv'=df_fac_indiv,'df_var_accounted'=df_var_accounted,'n_dim'=ncp_estimate))
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
