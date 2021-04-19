#**************************************************
# Description =====================================
#**************************************************
# R script to analyze structural properties from FreeSurfer data.
# Execute DoSCA(), DoPCA(), DoICA(), DoTSNE(), DoJK() to perform each analysis.


#**************************************************
# Parameters ======================================
#**************************************************

path_exp <- "Dropbox/MRI_img/pnTTC/puberty/stats/str_FS"
path_exp_full<-NULL
dir_in <-"01_extract"
dir_out <-"03_gam"
#dir_in <-"31_meas"
#dir_out <-"32_fp"

#wave <- 1
#measure <-c("volume","thickness","area")
#list_covar<-c("W1_Tanner_Max","W1_Age_at_MRI")

#key_global_covar<-"BrainSegVolNotVent"
#key_global_covar<-"eTIV"


#subset_subj <- list(list("column"="W1_T1QC","value"=1),list("column"="Sex","value"=1))
#subset_subj <- list(list("column"="W1_T1QC","value"=1),list("column"="Sex","value"=2))
#subset_subj <- list(list("column"="W1_5sub","value"=1))
#subset_subj <- list(list("column"="W1_5sub","value"=1),list("column"="Sex","value"=1))


#**************************************************
# Libraries =======================================
#**************************************************
library(easypackages)
libraries('Hmisc','FactoMineR','ica','tidyverse','Rtsne','tidyr','parallel','data.table')


#**************************************************
# Original library ================================
#**************************************************
source(file.path(getwd(),"util/function.R"))
source(file.path(getwd(),"util/plot.R"))
source(file.path(getwd(),"util/parameter.R"))
paths<-func_path(path_exp_=path_exp,dir_in_=dir_in,dir_out_=dir_out,path_exp_full_=path_exp_full)


#**************************************************
# GAM of structural measures ======================
#**************************************************

gam_str_core<-function(paths){
  # Prepare clinical data and demean
  # QC subsetting condition must accord with MRI wave, but under the name of clinical wave
  subset_subj<-param$subset_subj[wave_mri]
  names(subset_subj)<-wave_clin
  data_clin<-func_clinical_data_long(paths,wave_clin,subset_subj,list_covar,rem_na_clin=T,prefix=paste("var-",idx_var,"_wav-",label_wave,"_src",sep=""),print_terminal=F)
  df_clin<-func_demean_clin(data_clin$df_clin,separate_sex=T)$df_clin
  fwrite(df_clin,file.path(paths$output,"output","temp",paste("var-",idx_var,"_wav-",label_wave,"_src_clin.csv",sep="")),row.names=F)
  df_clin$wave<-wave_mri # Need to meet MRI wave for later joining
  
  # Join structure and clinical data
  df_join<-inner_join(df_str,df_clin,by=c('ID_pnTTC','wave'))
  for (key in c('ID_pnTTC','wave','sex')){
    if (key %in% colnames(df_join)){
      df_join[,key]<-as.factor(df_join[,key])
    }
  }
  df_join$value<-as.numeric.factor(df_join$value)
  
  for (measure in param$list_type_measure){
    print(paste("Starting to calculate GAM of ",meas,sep=""))
    df_join_measure<-df_join[which(df_join$measure==measure & df_join$wave==wave_mri),]
    list_roi<-sort(unique(df_join_measure$roi))
    df_roi<-func_dict_roi(paths)
    df_roi<-df_roi[df_roi$id %in% list_roi,c('id','label',param$key_group)]
    colnames(df_roi)<-c('id','label','group')
    
    # Prepare parallelization cluster
    if (calc_parallel){clust<-makeCluster(floor(detectCores()*3/4))}else{clust<-makeCluster(1)}
    clusterExport(clust,varlist=c("list_mod","list_sex","calc_parallel","test_mod","as.formula","as.numeric.factor",
                                  "lm","lmer","gam","summary","anova","summary.gam","anova.gam","AIC"),
                  envir=environment())
    
    # Calculate model
    data_gamm<-iterate_gamm4_str(clust,df_join_measure,progressbar=F,test_mod=test_mod)
    stopCluster(clust)
    
    
    ####
    
    if (meas=="volume"){
      data_glm<-func_glm(df_mri=df_str_meas,data_clinical,list_covar=list_covar_,df_global_covar=df_global_covar,key_global_covar=key_global_covar_)
    }else{
      data_glm<-func_glm(df_mri=df_str_meas,data_clinical,list_covar=list_covar_)
    }
    for (i in seq(length(data_glm))){
      if (is.element("roi",colnames(data_glm[[i]]))){
        id_roicol<-which(colnames(data_glm[[i]])=="roi")
        label_roi<-NULL
        for(j in data_glm[[i]]$roi){
          label_roi<-c(label_roi,as.character(dict_roi[which(dict_roi$id==j),"label"]))
        }
        data_glm[[i]]<-cbind(data_glm[[i]][,1:id_roicol],"label_roi"=label_roi,data_glm[[i]][,(1+id_roicol):ncol(data_glm[[i]])])
      }
    }
    write.csv(data_glm$glm,file.path(paths_$output,"output",paste("w",wave,"_",meas,"_glm.csv",sep="")),row.names = F)
    write.csv(data_glm$ic,file.path(paths_$output,"output",paste("w",wave,"_",meas,"_ic.csv",sep="")),row.names = F)
    write.csv(data_glm$min_ic,file.path(paths_$output,"output",paste("w",wave,"_",meas,"_min_ic.csv",sep="")),row.names = F)
    write.csv(data_glm$vif,file.path(paths_$output,"output",paste("w",wave,"_",meas,"_vif.csv",sep="")),row.names = F)
    print(paste("    Finished calculating GLM of ",meas,sep=""))
  }
}

gam_str<-function(paths_=paths,param=param_gam_str){
  print("Starting gam_str()")
  nullobj<-func_createdirs(paths_,str_proc="gam_str()",copy_log=T,list_param=param)

  df_str<-as.data.frame(fread(file.path(paths_$input,"output","fs_measure.csv")))
  #df_str$value[which(is.nan(df_str$value))]<-0
  #df_global_covar<-df_str[which(df_str$measure=="global" & df_str$wave==wave & df_str$roi==param$key_global_covar),]
  
  for (label_wave in names(param$list_wave)){
    wave_clin<-param$list_wave[[label_wave]]$clin
    wave_mri<-param$list_wave[[label_wave]]$mri
    # Loop over Tanner stages
    for (idx_tanner in names(param$list_tanner)){
      list_covar<-param$list_covar_tanner
      list_covar[["tanner"]]<-param$list_tanner[[idx_tanner]]
      gam_fc_cs_core(paths_,df_str,param,list_sex=list(1,2),list_covar,
                     list_mod=param$list_mod_tanner,list_term=param$list_term_tanner,idx_var=idx_tanner,
                     calc_parallel=T,test_mod=F)
    }
    # Loop over hormones
    for (idx_hormone in names(param$list_hormone)){
      list_covar<-param$list_covar_hormone
      list_covar[["hormone"]]<-param$list_hormone[[idx_hormone]]
      gam_fc_cs_core(paths_,atlas,param,list_sex=list(1,2),list_covar,
                     list_mod=param$list_mod_hormone,list_term=param$list_term_hormone,idx_var=idx_hormone,
                     calc_parallel=T,test_mod=F)
    } 
    

    

  }
  
  


  print("Finished gam_str().")
}

#**************************************************
# Fingerprinting ==================================
#**************************************************

# Core function for parallelization of fp_str()
fp_str_core<-function(data_str){
  measure<-data_str$measure
  group<-data_str$group
  df_str<-data_str$df_str
  df_ses_subj<-data_str$df_ses_subj
  n_roi<-dim(df_str)[1]
  
  # Calculate correlation matrix
  data_fingerprint<-func_cor(input=df_str)
  df_fp_subnet<-data_fingerprint$cor_flat
  
  # Rename correlation matrix to sessions and subjects
  df_fp_subnet$from_ses<-df_fp_subnet$from_ID_pnTTC<-df_fp_subnet$to_ses<-df_fp_subnet$to_ID_pnTTC<-NA
  for (i in seq(dim(df_fp_subnet)[1])){
    from_id<-df_fp_subnet[[i,"from"]]
    to_id<-df_fp_subnet[[i,"to"]]
    df_fp_subnet[[i,"from_ses"]]<-df_ses_subj[[from_id,"ses"]]
    df_fp_subnet[[i,"from_ID_pnTTC"]]<-df_ses_subj[[from_id,"ID_pnTTC"]]
    df_fp_subnet[[i,"to_ses"]]<-df_ses_subj[[to_id,"ses"]]
    df_fp_subnet[[i,"to_ID_pnTTC"]]<-df_ses_subj[[to_id,"ID_pnTTC"]]
  }
  df_fp_subnet$measure<-measure
  df_fp_subnet$group<-group
  df_fp_subnet<-df_fp_subnet[c("measure","group","from_ses","from_ID_pnTTC","to_ses","to_ID_pnTTC","r")]
  
  # Prepare dataframe for fingerprint correlation plot
  df_fp_plot<-data_fingerprint$cor
  list_name_subj_ses<-paste(sprintf("%05d",df_ses_subj$ID_pnTTC),as.character(df_ses_subj$ses),sep="_")
  colnames(df_fp_plot)<-rownames(df_fp_plot)<-list_name_subj_ses
  
  # Heatmap plot of fp correlation matrix
  plot_fp_heatmap<-plot_cor_heatmap(input=df_fp_plot)
  suppressMessages(plot_fp_heatmap<-(plot_fp_heatmap
                                     + scale_fill_gradientn(colors = matlab.like2(100),name="r")
                                     + ggtitle(paste("FP Cor,",atlas,measure,group,sep=" "))
                                     + theme(plot.title = element_text(hjust = 0.5),
                                             axis.title=element_blank())))
  
  # Save heatmap plot
  ggsave(paste("atl-",atlas,"_msr-",measure,"_grp-",group,"_fp.eps",sep=""),plot=plot_fp_heatmap,device=cairo_ps,
         path=file.path(paths_$output,"output"),dpi=300,height=10,width=10,limitsize=F)
  
  return(df_fp_subnet)
}

# Main function for fingerprint computing
fp_str<-function(paths_=paths,
                 atlas="dk",
                 key_roigroup="group_3"){
  print("Starting fp_str().")
  nullobj<-func_createdirs(paths_)
  dict_roi<-func_dict_roi(paths_)
  dict_roi<-data.frame(id=as.character(dict_roi$id),group=as.character(dict_roi[,key_roigroup]),stringsAsFactors = F)
  
  # Load structure data
  df_str<-read.csv(file.path(paths_$input,"output","fs_measure.csv"))
  list_measure<-unique(as.character(df_str$measure))
  
  list_data_str<-list()
  for (measure in list_measure){
    df_str_measure<-df_str[df_str$measure==measure,]
    
    # Examine existing subject IDs and sessions in connection data
    list_ses_exist <- sort(unique(df_str_measure$ses))
    list_id_subj_exist<-list()
    for (ses in list_ses_exist){
      df_str_ses<-df_str_measure[df_str_measure$ses==ses,]
      list_id_subj_exist[[as.character(ses)]]<-sort(unique(df_str_ses$ID_pnTTC))
    }
    
    # Examine existing ROI and their groups
    df_roi<-df_str_measure[df_str_measure$ID_pnTTC==df_str_measure[1,"ID_pnTTC"],]
    df_roi<-df_roi[df_roi$ses==df_roi[1,"ses"],"roi",drop=F]
    df_roi$roi<-as.character(df_roi$roi)
    df_roi<-df_roi[order(df_roi$roi),,drop=F]
    df_roi<-left_join(df_roi,dict_roi,by=c("roi"="id"))
    
    # List groups of existing rois
    list_group<-sort(unique(c(df_roi$group)))
    if (!("whole" %in% list_group)){
      if (length(list_group)>1){
        list_group<-c("whole",list_group)
      }
    }
    print(paste("Measure: ",measure, ", ", as.character(length(list_group))," groups:",sep=""))
    print(list_group)
    
    # Split and combine structure data for each subgroup of rois for later parallel computing

    for (group in list_group){
      # Create dataframe of structural measures within each group
      if (group=="whole"){
        df_roi_group<-df_roi
      }else{
        df_roi_group<-df_roi[df_roi$group==group,]
      }
      n_roi_group<-dim(df_roi_group)[1]
      
      if (n_roi_group<5){
        print(paste("Measure: ",measure,", group: ",group, ", nodes: ",as.character(n_roi_group)," < 5, fp calculation skipped.",sep=""))
      }else{
        # Create combined dataframe of structural measures
        # Also create dataframe of sessions and subjects
        df_str_cbind<-data.frame(matrix(nrow=n_roi_group,ncol=0))
        df_ses_subj<-data.frame(matrix(nrow=0,ncol=2))
        colnames(df_ses_subj)<-c("ses","ID_pnTTC")
        for (ses in list_ses_exist){
          for (id_subj in list_id_subj_exist[[ses]]){
            df_str_subj<-df_str_measure[df_str_measure$ID_pnTTC==id_subj & df_str_measure$ses==ses,]
            df_str_subj<-df_str_subj[df_str_subj$roi %in% df_roi_group$roi,]
            df_str_cbind<-cbind(df_str_cbind,df_str_subj[["value"]])
            df_ses_subj<-rbind(df_ses_subj,data.frame(ses=ses,ID_pnTTC=id_subj))
          }
        }
        colnames(df_str_cbind)<-as.character(seq(ncol(df_str_cbind)))
        rownames(df_str_cbind)<-NULL
        
        list_data_str<-c(list_data_str,list(list("measure"=measure,"group"=group,
                                                 "df_str"=df_str_cbind,"df_ses_subj"=df_ses_subj)))
      }
    }
  }
  
  # Parallel fingerprint correlation computing over groups of subnetworks
  print("Calculating structure fingerprint correlation in parallel.")
  n_cluster<-min(floor(detectCores()*3/4),length(list_data_str))
  clust<-makeCluster(n_cluster)
  clusterExport(clust,
                varlist=c("paths_","func_cor","atlas",
                          "plot_cor_heatmap","rcorr","rownames_to_column","gather",
                          "ggplot","aes","geom_tile","scale_fill_gradientn",
                          "matlab.like2","scale_y_discrete","scale_x_discrete",
                          "theme_light","theme","element_text","element_blank",
                          "ggtitle","ggsave"),
                envir=environment())
  list_df_fp<-parLapply(clust,list_data_str,fp_str_core)
  stopCluster(clust)
  
  # Output dataframe
  df_fp<-NULL
  for (df_fp_subnet in list_df_fp){
    if (!is.null(df_fp_subnet)){
      df_fp<-rbind(df_fp,df_fp_subnet)
    }
  }
  
  # Save fingerprint correlation
  write.csv(df_fp,file.path(paths_$output,"output",paste("atl-",atlas,"_fp.csv",sep="")),row.names=F)
  print("Finished fp_str().")
}


#**************************************************
# GLM of structural measures ======================
#**************************************************
glm_str<-function(paths_=paths,subset_subj_=subset_subj,list_covar_=list_covar,file_input_=file_input,
                  wave_=wave,key_global_covar_=key_global_covar
                  ){
  print("Starting glm_str()")
  data_clinical<-func_clinical_data(paths_,subset_subj_)
  nullobj<-func_createdirs(paths_,copy_log=T)
  dict_roi<-func_dict_roi(paths_)
  df_str<-read.csv(file.path(paths_$input,"output",file_input_))
  df_str$value[which(is.nan(df_str$value))]<-0
  df_global_covar<-df_str[which(df_str$measure=="global" & df_str$wave==wave & df_str$roi==key_global_covar_),]
  for (meas in measure){
    print(paste("    Starting to calculate GLM of ",meas,sep=""))
    df_str_meas<-df_str[which(df_str$measure==meas & df_str$wave==wave_),]
    if (meas=="volume"){
      data_glm<-func_glm(df_mri=df_str_meas,data_clinical,list_covar=list_covar_,df_global_covar=df_global_covar,key_global_covar=key_global_covar_)
    }else{
      data_glm<-func_glm(df_mri=df_str_meas,data_clinical,list_covar=list_covar_)
    }
    for (i in seq(length(data_glm))){
      if (is.element("roi",colnames(data_glm[[i]]))){
        id_roicol<-which(colnames(data_glm[[i]])=="roi")
        label_roi<-NULL
        for(j in data_glm[[i]]$roi){
          label_roi<-c(label_roi,as.character(dict_roi[which(dict_roi$id==j),"label"]))
        }
        data_glm[[i]]<-cbind(data_glm[[i]][,1:id_roicol],"label_roi"=label_roi,data_glm[[i]][,(1+id_roicol):ncol(data_glm[[i]])])
      }
    }
    write.csv(data_glm$glm,file.path(paths_$output,"output",paste("w",wave,"_",meas,"_glm.csv",sep="")),row.names = F)
    write.csv(data_glm$ic,file.path(paths_$output,"output",paste("w",wave,"_",meas,"_ic.csv",sep="")),row.names = F)
    write.csv(data_glm$min_ic,file.path(paths_$output,"output",paste("w",wave,"_",meas,"_min_ic.csv",sep="")),row.names = F)
    write.csv(data_glm$vif,file.path(paths_$output,"output",paste("w",wave,"_",meas,"_vif.csv",sep="")),row.names = F)
    print(paste("    Finished calculating GLM of ",meas,sep=""))
  }
  print("Finished glm_str().")
  return(data_glm)
  
}


#**************************************************
# OBSOLETE ========================================
#**************************************************

##### Parameters ####
#p_uncorrected<-0.001
#p_corrected<-0.05
#
#n_components<-10
##n_components<-30
##n_components<-5
#tsne_dims<-2
#tsne_perplexity<-30
#tsne_max_itr<-1000
#
##### Data Loading ####
#
#source(file.path(script_dir,"Functionalities/LoadClinicalData.R"))
#HeatmapPlot(clinical_data,"Clinical Data","Clinical Measure",
#            colnames(clinical_data)[-1],scale_data = T)
#structural_data<-read.csv(file.path(input_dir,structural_file))
#structural_data$flag<-F
#for (i in subject_id){
#  structural_data[which(structural_data$ID_pnTTC==i),"flag"]<-T
#}
#structural_data<-structural_data[which(structural_data$flag),-(ncol(structural_data))]
#colnames(structural_data)[-1]<-ConvertID(colnames(structural_data)[-1],roi_data,input_roi_type,"ID_long")
#n_rois<-ncol(structural_data)-1
#HeatmapPlot(structural_data,
#            "Scaled ROI Measures",
#            "ROI",
#            ConvertID(colnames(structural_data)[-1],roi_data,"ID_long","label_proper"),
#            scale_data=T)
#
#
##### General Linear Model Analysis ####
#
#DoGLM<-function(){
#  dirname<-ExpDir("GLM")
#  structural_data_tidy<-gather(structural_data,key=ROI,value=value,-ID_pnTTC)
#  structural_data_tidy$ROI_label<-ConvertID(structural_data_tidy$ROI,roi_data,"ID_long","label_proper")
#  glm<-CommonGLM(structural_data_tidy,covariate_label,global_covariate=F,dirname,"GLM_Structure.csv")
#  models_expvars<-glm[which(glm[,"ROI"]==glm[1,"ROI"]),
#                      c("model","exp_var")]
#  glm_ordered<-NULL
#  for (i in 1:nrow(models_expvars)){
#    id_obs<-which(glm[,"model"]==models_expvars[i,"model"])
#    id_obs<-intersect(id_obs,which(glm[,"exp_var"]==models_expvars[i,"exp_var"]))
#    glm_subset<-glm[id_obs,]
#    glm_subset<-MultCompCorr(glm_subset)
#    glm_ordered<-rbind(glm_ordered,glm_subset)
#  }
#  write.csv(glm_ordered, file.path(dirname,"GLM_ordered.csv"),row.names=F)
#  output<-list(glm,glm_ordered)
#  names(output)<-c("GLM","GLM_ordered")
#  return(output)
#}
#
#
##### Structural Correlation Analysis ####
#
#DoSCA<-function(){
#  dirname<-ExpDir("SCA")
#  corr<-CalcCorr(structural_data[-1], dirname,"SCA")
#  graph<-Corr2Graph(corr)
#  fig1<-corr[[3]]
#  fig2<-CircularPlot(graph,
#                     pvalue_type="p_Benjamini_Hochberg",
#                     input_title="Structural Covariance for All Subjects")
#  return(list(corr,fig1,fig2))
#}
#
#
##### Principal Component Analysis ####
#
#DoPCA<-function(){
#  dirname<-ExpDir("PCA")
#  data<-structural_data[-1]
#  indexcolumn<-structural_data[1]
#  pca <-PCA(data,scale.unit = TRUE, ncp = n_components, graph = FALSE)
#  varfactor<-data.frame(pca$var$coord)
#  varfactor<-cbind(colnames(data),ConvertID(colnames(data),roi_data,"ID_long","label_proper"),varfactor)
#  colnames(varfactor)<-c("ROI_ID","ROI_proper",sprintf("Dim_%02d",1:n_components))
#  rownames(varfactor)<-NULL
#  indfactor<-data.frame(pca$ind$coord)
#  colnames(indfactor)<-sprintf("Dim_%02d",1:n_components)
#  indfactor<-cbind(indexcolumn,indfactor)
#  varianceaccounted<-pca$eig
#  
#  write.csv(varfactor, file.path(dirname,"VariableFactor.csv"),row.names=F)
#  write.csv(indfactor, file.path(dirname,"IndividualFactor.csv"),row.names=F)
#  write.csv(varianceaccounted, file.path(dirname,"VarianceAccounted.csv"))
#  clincorr<-MeasClinicalCorr(indfactor,dirname)
#  
#  return(list(varfactor,indfactor,varianceaccounted,clincorr))
#}
#
#
##### Independent Component Analysis ####
#
#DoICA<-function(){
#  dirname<-ExpDir("ICA")
#  data<-structural_data[-1]
#  indexcolumn<-structural_data[1]
#  ica <-icafast(data, nc=n_components,center=TRUE,maxit=100,tol=1e-6,alg="par",fun="logcosh",alpha=1)
#  varfactor<-data.frame(ica$M)
#  varfactor<-cbind(colnames(data),ConvertID(colnames(data),roi_data,"ID_long","label_proper"),varfactor)
#  colnames(varfactor)<-c("ROI_ID","ROI_proper",sprintf("Dim_%02d",1:n_components))
#  rownames(varfactor)<-NULL
#  indfactor<-data.frame(ica$S)
#  colnames(indfactor)<-sprintf("Dim_%02d",1:n_components)
#  indfactor<-cbind(indexcolumn,indfactor)
#  varianceaccounted<-ica$vafs
#  
#  write.csv(varfactor, file.path(dirname,"VariableFactor.csv"),row.names=F)
#  write.csv(indfactor, file.path(dirname,"IndividualFactor.csv"),row.names=F)
#  write.csv(varianceaccounted, file.path(dirname,"VarianceAccounted.csv"))
#  clincorr<-MeasClinicalCorr(indfactor,dirname)
#  
#  return(list(varfactor,indfactor,varianceaccounted,clincorr))
#}
#
#
##### t-SNE Analysis ####
#
#DoTSNE<-function(){
#  dirname<-ExpDir("tSNE")
#  data<-structural_data[-1]
#  indexcolumn<-structural_data[1]
#  #  tsne <- Rtsne(data[-1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
#  tsne<-Rtsne(data, dims = tsne_dims, perplexity=tsne_perplexity, verbose=TRUE, max_iter = tsne_max_itr)
#  indfactor<-data.frame(tsne$Y)
#  colnames(indfactor)<-sprintf("Dim_%02d",1:tsne_dims)
#  #  pairs(indfactor)
#  indfactor<-cbind(indexcolumn,indfactor)
#  graph<-ComponentPlot(indfactor[-1],"t-SNE")
#  write.csv(indfactor, file.path(dirname,"Coordinates.csv"),row.names=F)
#  return(list(tsne,graph))
#}
#
#
##### Jackknife Estimate of Structural Covariance ####
#
#DoJK<-function(){
#  dirname<-ExpDir("JK")
#  overall_corr<-CalcCorr(structural_data[-1],plot=F,save=F)[[2]]
#  jkstats<-data.frame(matrix(nrow=0,ncol=6))
#  jkstats<-rbind(jkstats,cbind(metric=rep("overall_corr",nrow(overall_corr)),overall_corr[,-6]))
#  jkpe<-jkpv<-jkz<-data.frame(matrix(ncol=6,nrow=0))
#  for (i in 1:n_subject){
#    jk<-CalcCorr(structural_data[-i,-1],plot=F,save=F)[[2]]
#    jkpe<-rbind(jkpe,cbind(ID_pnTTC=rep(subject_id[i],nrow(jk)),jk[,-6]))
#    jkpv<-rbind(jkpv,cbind(ID_pnTTC=rep(subject_id[i],nrow(jk)),jk[,1:4],
#                           n_subject*overall_corr$r-(n_subject-1)*jk$r))
#  }
#  for (j in 1:nrow(overall_corr)){
#    mean_jkpe<-mean(as.numeric(jkpe[intersect(which(jkpe$from==overall_corr[j,"from"]),
#                                              which(jkpe$to==overall_corr[j,"to"])),"r"]))
#    sd_jkpe<-sd(as.numeric(jkpe[intersect(which(jkpe$from==overall_corr[j,"from"]),
#                                          which(jkpe$to==overall_corr[j,"to"])),"r"]))
#    jkstats<-rbind(jkstats,cbind(metric="mean_jkpe",overall_corr[j,1:4],r=mean_jkpe))
#    jkstats<-rbind(jkstats,cbind(metric="sd_jkpe",overall_corr[j,1:4],r=sd_jkpe))
#    for (i in 1:n_subject){
#      jkz<-rbind(jkz,
#                 cbind(ID_pnTTC=subject_id[i], overall_corr[j,1:4],
#                       r=(-1)*(jkpe[intersect(intersect(which(jkpe$ID_pnTTC==subject_id[i]),
#                                                        which(jkpe$from==overall_corr[j,"from"])),
#                                              which(jkpe$to==overall_corr[j,"to"])),"r"]-mean_jkpe)/sd_jkpe))
#    }
#  }
#  jkz<-jkz[order(jkz$ID_pnTTC),]
#  jkstats<-jkstats[c(which(jkstats$metric=="overall_corr"),
#                     which(jkstats$metric=="mean_jkpe"),
#                     which(jkstats$metric=="sd_jkpe")),]
#  colnames(jkstats)<-c("metric","from","from_label","to","to_label","value")
#  colnames(jkpe)<-colnames(jkpv)<-colnames(jkz)<-c("ID_pnTTC","from","from_label","to","to_label","r")
#  write.csv(jkstats, file.path(dirname,"JKStats.csv"),row.names=F)
#  write.csv(jkpe, file.path(dirname,"JKPE.csv"),row.names=F)
#  write.csv(jkpv, file.path(dirname,"JKPV.csv"),row.names=F)
#  write.csv(jkz, file.path(dirname,"JKZ.csv"),row.names=F)
#  return(list(jkstats,jkpe,jkpv,jkz))
#}
#
#
##### Laterality Index Calculation ####
#
#DoLI<-function(){
#  dirname<-ExpDir("LI")
#  structural_data_tidy<-gather(structural_data,key=ROI,value=value,-ID_pnTTC)
#  li<-CommonLI(structural_data_tidy,"ROI",dirname,"LI_Structure.csv")
##  clincorr<-MeasClinicalCorr(output[c(-1,-2),],dirname)
#  li_tidy<-li[,c("ID_pnTTC","ROI","L_ROI_ID","R_ROI_ID","Laterality_Index")]
#  colnames(li_tidy)[5]<-"value"
#  glm<-CommonGLM(li_tidy,covariate_label,F,dirname,"GLM_LI_Structure.csv")
#  models_expvars<-glm[which(glm[,"ROI"]==glm[1,"ROI"]),
#                      c("model","exp_var")]
#  glm_ordered<-NULL
#  for (i in 1:nrow(models_expvars)){
#    id_obs<-which(glm[,"model"]==models_expvars[i,"model"])
#    id_obs<-intersect(id_obs,which(glm[,"exp_var"]==models_expvars[i,"exp_var"]))
#    glm_subset<-glm[id_obs,]
#    glm_subset<-cbind(glm_subset,MultCompCorr(glm_subset))
#    glm_ordered<-rbind(glm_ordered,glm_subset)
#  }
#  write.csv(glm_ordered, file.path(dirname,"GLM_LI_Structure_ordered.csv"),row.names=F)
#  output<-list(li,glm,glm_ordered)
#  names(output)<-c("Laterality_Index","GLM_of_LI","GLM_of_LI_ordered")
#  return(output)
#}
#