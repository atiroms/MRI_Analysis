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

dir_in<-"103_fp_acompcor"
dir_out<-"104_fp_id_acompcor"
#dir_out<-"105_fp_model_acompcor"

#dir_in<-"113_fp_aroma"
#dir_out<-"115_fp_model_aroma"

#dir_out<-"55_gta_bin"
#dir_out<-"56_fp_identification"
#dir_out<-"58_glm_ancova_acompcor"

list_wave <- c(1,2)

subset_subj <- list("1"=list(list("key"="W1_T1QC","value"=1),
                             list("key"="W1_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)),
                    "2"=list(list("key"="W2_T1QC","value"=1),
                             list("key"="W2_T1QC_new_mild_rsfMRIexist_motionQC3","value"=1)))

list_covar<-list("tanner"=list("1"="W1_Tanner_Max",
                               "2"="W2_Tanner_Max",
                               "label"="Tanner stage"),
                 "age"=list("1"="W1_Age_at_MRI",
                            "2"="W2_Age_at_MRI",
                            "label"="Age"),
                 "sex"=list("1"="Sex",
                            "2"="Sex",
                            "label"="Sex"))

list_mod <- list("lin_diff_a"=
                   "value ~ diff_age + sex",
                 "lin_diff_t"=
                   "value ~ sex + sex:diff_tanner",
                 "lin_diff_at"=
                   "value ~ diff_age + sex + sex:diff_tanner",
                 "add_diff_a"=
                   "value ~ s(diff_age,k=3) + sex",
                 "add_diff_t"=
                   "value ~ sex + s(diff_tanner,k=3,by=sex)",
                 "add_diff_at"=
                   "value ~ s(diff_age,k=3) + sex + s(diff_tanner,k=3,by=sex)"
                 )

#list_mod <- list("additive"=
#                   "value ~ s(age,k=3) + sex + s(tanner,k=3,by=sex)",
#                 "additive_mixed"=
#                   "value ~ s(age,k=3) + sex + s(tanner,k=3,by=sex) + s(ID_pnTTC,bs='re')",
#                 "linear"=
#                   "value ~ age + sex + sex:tanner",
#                 "linear_mixed"=
#                   "value ~ age + sex + sex:tanner + s(ID_pnTTC,bs='re')")

list_graph <-list("a"=list("title"="Effect of age difference",
                           "x_axis"="diff_age",
                           "smooth"=list("Male"=list("fix"=list("sex"=1),
                                                     "color"="steelblue2","alpha"=1,"ribbon"=T),
                                         "Female"=list("fix"=list("sex"=2),
                                                       "color"="lightcoral","alpha"=1,"ribbon"=T)),
                           "point"=list("Male"=list("subset"=list("sex"=1),
                                                    "color"="steelblue2","alpha"=1),
                                        "Female"=list("subset"=list("sex"=2),
                                                      "color"="lightcoral","alpha"=1))),
                  "st"=list("title"="Effect of Tanner stage difference",
                            "x_axis"="diff_tanner",
                            "smooth"=list("Male"=list("fix"=list("sex"=1),
                                                      "color"="steelblue2","alpha"=1,"ribbon"=T),
                                          "Female"=list("fix"=list("sex"=2),
                                                        "color"="lightcoral","alpha"=1,"ribbon"=T)),
                            "point"=list("Male"=list("subset"=list("sex"=1),
                                                     "color"="steelblue2","alpha"=1),
                                         "Female"=list("subset"=list("sex"=2),
                                                       "color"="lightcoral","alpha"=1))),
                  "sat"=list("title"="Age difference-Tanner stage difference interaction",
                             "x_axis"="diff_age",
                             "smooth"=list("Male delta TS = -1"=list("fix"=list("sex"=1,"diff_tanner"=-1),
                                                              "color"="Steelblue2","alpha"=0.4,"ribbon"=F),
                                           "Male delta TS = 1"=list("fix"=list("sex"=1,"diff_tanner"=1),
                                                              "color"="steelblue2","alpha"=0.7,"ribbon"=F),
                                           "Male delta TS = 3"=list("fix"=list("sex"=1,"diff_tanner"=3),
                                                              "color"="steelblue2","alpha"=1,"ribbon"=F),
                                           "Female delta TS = -1"=list("fix"=list("sex"=2,"diff_tanner"=-1),
                                                                "color"="lightcoral","alpha"=0.4,"ribbon"=F),
                                           "Female delta TS = 1"=list("fix"=list("sex"=2,"diff_tanner"=1),
                                                                "color"="lightcoral","alpha"=0.7,"ribbon"=F),
                                           "Female delta TS = 3"=list("fix"=list("sex"=2,"diff_tanner"=3),
                                                                "color"="lightcoral","alpha"=1,"ribbon"=F)),
                             "point"=list("Male"=list("subset"=list("sex"=1),
                                                      "color"="steelblue2","alpha"=1),
                                          "Female"=list("subset"=list("sex"=2),
                                                        "color"="lightcoral","alpha"=1))))


list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")
#list_atlas<-"aal116"
#list_atlas<-"schaefer400"
#list_atlas<-c("glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")
#thr_pvalue <- 0.05
n_permutation<-1000
#n_permutation<-100


#**************************************************
# Libraries =======================================
#**************************************************
library(dplyr)
library(mgcv)
library(rowr)
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
source(file.path(paths$script,"util/function.R"))
#source(file.path(paths$script,util/glm_function.R"))
source(file.path(paths$script,"util/plot.R"))
#source(file.path(paths$script,"util/gta_function.R"))


#**************************************************
# GLM and ANCOVA of Fingerprint change ============
#**************************************************

paths_=paths
list_atlas_=list_atlas
list_wave_=list_wave
list_covar_=list_covar
list_mod_=list_mod
list_graph_=list_graph
subset_subj_=subset_subj

model_fp<-function(paths_=paths,
                   list_atlas_=list_atlas,
                   list_wave_=list_wave,
                   list_covar_=list_covar,
                   list_mod_=list_mod,
                   list_graph_=list_graph,
                   subset_subj_=subset_subj
                   ){
  print("Starting glm_ancova_fp().")
  nullobj<-func_createdirs(paths_)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,
                                     list_covar=list_covar_,rem_na_clin=T)
  df_clin<-data_clin$df_clin
  colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
  
  df_out_lm<-data.frame(matrix(nrow=0,ncol=8))
  colnames(df_out_lm)<-c("atlas","model","term","estimate","se","F","t","p")
  df_out_aic<-data.frame(matrix(nrow=0,ncol=4))
  colnames(df_out_aic)<-c("atlas","model","aic","aic_best_among_models")
  
  df_out_ancova<-data.frame(matrix(ncol=10,nrow=0))
  colnames(df_out_ancova)<-c("atlas","sex","test","term","comparison","p","F","diff","ci_l","ci_u")
  
  for (atlas in list_atlas_){
    # Load fingerprint data
    print(paste("Loading atlas: ",atlas,sep=""))
    df_fp<-read.csv(file.path(paths_$input,"output",paste("atl-",atlas,"_fp.csv",sep="")))
    
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

    # Collect longitudinal fp correlation data
    df_cor_fp<-data.frame(ID_pnTTC=list_id_subj_exist_twice)
    for (id_subj in list_id_subj_exist_twice){
      df_cor_fp[df_cor_fp$ID_pnTTC==id_subj,"value"]<-df_fp[df_fp$from_ID_pnTTC==id_subj &df_fp$to_ID_pnTTC==id_subj,"r"]
    }
    
    # Subset those without longitudinal fp correlation
    list_id_subj_exist_twice<-df_cor_fp[!is.na(df_cor_fp$value),"ID_pnTTC"]
    df_cor_fp<-df_cor_fp[df_cor_fp$ID_pnTTC %in% list_id_subj_exist_twice,]
    n_id_subj_exist_twice<-length(list_id_subj_exist_twice)
    print(paste(as.character(n_id_subj_exist_twice)," subjects with non-NA data for two sessions.",sep=""))
    
    # Create dataframe for GLM analysis
    df_join<-func_clinical_data_join(df_src=df_clin,
                                     list_id_subj=list_id_subj_exist_twice,
                                     list_covar=list_covar_)
    
    df_join<-inner_join(df_join,df_cor_fp,by="ID_pnTTC")
    df_join$ID_pnTTC<-as.factor(df_join$ID_pnTTC)
    df_join$sex<-as.factor(df_join$sex)
    write.csv(df_join,file.path(paths_$output,"output",
                                     paste("atl-",atlas,"_fp_glm_ancova_src.csv",sep="")),row.names = F)
    
    # Calculate GLM
    print('Calculating GLM.')
    list_mod_gamm<-list()
    df_out_aic_atlas<-data.frame()
    for (mod in names(list_mod_)){
      print(paste("GLM of model: ", mod, sep=""))
      list_mod_gamm[[mod]]<-gam(as.formula(list_mod_[[mod]]),data=df_join)
      p_table<-summary.gam(list_mod_gamm[[mod]])$p.table
      if (is.null(summary.gam(list_mod_gamm[[mod]])$s.table)){
        df_out_lm_add<-data.frame(atlas=atlas,model=mod,term=rownames(p_table),
                                  estimate=p_table[,'Estimate'],se=p_table[,'Std. Error'],F=NA,
                                  t=p_table[,'t value'],p=p_table[,'Pr(>|t|)'])
        
      }else{
        s_table<-summary.gam(list_mod_gamm[[mod]])$s.table
        df_out_lm_add<-rbind(data.frame(atlas=atlas,model=mod,term=rownames(p_table),
                                        estimate=p_table[,'Estimate'],se=p_table[,'Std. Error'],F=NA,
                                        t=p_table[,'t value'],p=p_table[,'Pr(>|t|)']),
                             data.frame(atlas=atlas,model=mod,term=rownames(s_table),
                                        estimate=NA,se=NA,F=s_table[,'F'],
                                        t=NA,p=s_table[,'p-value']))
      }
      df_out_lm<-rbind(df_out_lm,df_out_lm_add)
      df_out_aic_atlas<-rbind(df_out_aic_atlas,
                             data.frame(atlas=atlas,model=mod,aic=list_mod_gamm[[mod]]$aic,
                                        aic_best_among_models=0))
      
      for (idx_graph in names(list_graph_)){
        if (list_graph_[[idx_graph]][["x_axis"]] %in% colnames(list_mod_gamm[[mod]]$model)){
          plot<-plot_gamm(mod_gamm=list_mod_gamm[[mod]],
                          df_join,
                          spec_graph=list_graph_[[idx_graph]])
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
                 + ggtitle(list_graph_[[idx_graph]][["title"]])
                 + xlab(label_x)
                 + ylab("Fingerprint correlation")
                 + theme(legend.position = "none"))
          filename_plot<-paste("atl-",atlas,"_mod-",mod,"_plt-",idx_graph,"_fp_glm.eps",sep="")
          ggsave(filename_plot,plot=plot,device=cairo_ps,
                 path=file.path(paths_$output,"output"),dpi=300,height=5,width=5,limitsize=F)
        }
      }
    }
    
    # Compare AICs of GLM models
    df_out_aic_atlas[which(df_out_aic_atlas$aic==min(df_out_aic_atlas$aic)),'aic_best_among_models']<-1
    df_out_aic<-rbind(df_out_aic,df_out_aic_atlas)
    
    # Calculate ANCOVA
    print('Calculating ANCOVA.')
    df_join$long_tanner<-paste(as.character(df_join$ses1_tanner),
                               as.character(df_join$ses2_tanner),sep="_")
    df_join$long_tanner<-as.factor(df_join$long_tanner)
    list_sex<-list("all"=c(1,2),"male"=1,"female"=2)

    for (id_sex in names(list_sex)){
      print(paste("ANCOVA of sex: ",id_sex,sep=''))
      df_join_sex<-df_join[df_join$sex %in% list_sex[[id_sex]],]
      df_heatmap<-data.frame(matrix(ncol=5,nrow=5))
      df_id<-data.frame(matrix(ncol=0,nrow=0))
      for (tanner_ses1 in seq(5)){
        for (tanner_ses2 in seq(5)){
          list_id_intersect<-df_join_sex[which(df_join_sex$ses1_tanner==tanner_ses1 & df_join_sex$ses2_tanner==tanner_ses2),"ID_pnTTC"]
          list_id_intersect<-sort(list_id_intersect[!is.na(list_id_intersect)])
          n_id_intersect<-length(list_id_intersect)
          df_heatmap[tanner_ses1,tanner_ses2]<-n_id_intersect
          df_id<-cbind.fill(df_id,list_id_intersect,fill=NA)
          colnames(df_id)[dim(df_id)[2]]<-paste("tanner",as.character(tanner_ses1),
                                                as.character(tanner_ses2),sep="_")
        }
      }
      colnames(df_heatmap)<-rownames(df_heatmap)<-as.character(seq(5))
      df_id$init<-NULL
      write.csv(df_heatmap,file.path(paths_$output,"output",
                                     paste("atl-",atlas,"_sex-",
                                           id_sex,"_tanner_heatmap.csv",sep="")))
      write.csv(df_id,file.path(paths_$output,"output",
                                     paste("atl-",atlas,"_sex-",
                                           id_sex,"_tanner_id.csv",sep="")),row.names=F)
      if (id_sex=="all"){
        mod_ancova<-aov(value~long_tanner+diff_age+sex,data=df_join_sex)
      }else{
        mod_ancova<-aov(value~long_tanner+diff_age,data=df_join_sex)
      }
      df_ancova<-summary(mod_ancova)[[1]]
      df_out_ancova_add<-data.frame(atlas=atlas,sex=id_sex,test="ANCOVA",term=rownames(df_ancova),comparison=NA,
                                    p=df_ancova[,'Pr(>F)'],F=df_ancova[,'F value'],diff=NA,ci_l=NA,ci_u=NA)
      df_out_ancova<-rbind(df_out_ancova,df_out_ancova_add)
      
      # Calculate Tukey-Kramer
      df_tc<-TukeyHSD(mod_ancova,which='long_tanner')[[1]]
      df_out_ancova_add<-data.frame(atlas=atlas,sex=id_sex,test="Tukey-Kramer",term="long_tanner",comparison=rownames(df_tc),
                                    p=df_tc[,'p adj'],F=NA,diff=df_tc[,'diff'],ci_l=df_tc[,'lwr'],ci_u=df_tc[,'upr'])
      df_out_ancova<-rbind(df_out_ancova,df_out_ancova_add)
    }
  }
  rownames(df_out_lm)<-rownames(df_out_aic)<-rownames(df_out_ancova)<-NULL
  write.csv(df_out_lm, file.path(paths_$output,"output","fp_glm.csv"),row.names = F)
  write.csv(df_out_aic,file.path(paths_$output,"output","fp_glm_aic.csv"),row.names = F)
  write.csv(df_out_ancova,file.path(paths_$output,"output","fp_ancova.csv"),row.names=F)
  print("Finished glm_ancova_fp()")
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
  
  df_out_combined<-data.frame(matrix(nrow=0,ncol=11))
  colnames(df_out_combined)<-c("atlas","n_subj","n_identified","proportion_identified",
                               "n_identified_1_targeted","proportion_identified_1_targeted",
                               "n_identified_2_targeted","proportion_identified_2_targeted",
                               "p_permutation","p_permutation_1_targeted","p_permutation_2_targeted")
  
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
    
    # Output subset with longitudinal data
    write.csv(df_fp_exist_twice,file.path(paths_$output,"output",paste("atl-",atlas,"_fp_input_subset.csv",sep="")),row.names=F)
    df_fp_exist_twice_plot<-df_fp_exist_twice[c('from_ID_pnTTC','to_ID_pnTTC','r')]
    df_fp_exist_twice_plot<-spread(df_fp_exist_twice_plot,key=to_ID_pnTTC,value=r)
    colnames(df_fp_exist_twice_plot)[-1]<-rownames(df_fp_exist_twice_plot)<-sprintf("%05d",df_fp_exist_twice_plot$from_ID_pnTTC)
    df_fp_exist_twice_plot<-df_fp_exist_twice_plot[-1]
    plot_fp_exist_twice<-plot_cor_heatmap(input=df_fp_exist_twice_plot)
    plot_fp_exist_twice<-(plot_fp_exist_twice
                          + scale_fill_gradientn(colors = matlab.like2(100),name="r")
                          + ggtitle("Fingerprint correlation of two waves")
                          + xlab("2nd wave")
                          + ylab("1st wave")
                          + theme(plot.title = element_text(hjust = 0.5)))
    ggsave(paste("atl-",atlas,"_fp_identification.eps",sep=""),plot=plot_fp_exist_twice,device=cairo_ps,
           path=file.path(paths_$output,"output"),dpi=300,height=10,width=10,limitsize=F)
    
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
    
    n_subj<-length(list_id_subj_exist_twice)
    n_id_1_tar<-sum(df_ident["1_targeted_identification"])
    n_id_2_tar<-sum(df_ident["2_targeted_identification"])
    n_id<-n_id_1_tar+n_id_2_tar
    prop_id_1_tar<-n_id_1_tar/n_subj
    prop_id_2_tar<-n_id_2_tar/n_subj
    prop_id<-n_id/(n_subj*2)
    p_perm_1_tar<-sum(df_perm["1_n_ident"]>sum(df_ident["1_targeted_identification"]))/n_permutation_
    p_perm_2_tar<-sum(df_perm["2_n_ident"]>sum(df_ident["2_targeted_identification"]))/n_permutation_
    p_perm<-(p_perm_1_tar+p_perm_2_tar)/2
    
    df_out_combined<-rbind(df_out_combined,
                           data.frame(atlas=atlas,n_subj=n_subj,n_identified=n_id,proportion_identified=prop_id,
                                      n_identified_1_targeted=n_id_1_tar,proportion_identifeid_1_targeted=prop_id_1_tar,
                                      n_identified_2_targeted=n_id_2_tar,proportion_identified_2_targeted=prop_id_2_tar,
                                      p_permutation=p_perm,
                                      p_permutation_1_targeted=p_perm_1_tar,p_permutation_2_targeted=p_perm_2_tar))

    write.csv(df_ident,file.path(paths_$output,"output",paste("atl-",atlas,"_fp_identification.csv",sep="")),row.names=F)
    print(paste("Number of session 2 subjects identified from session 1 target: ",as.character(sum(df_ident["1_targeted_identification"])),sep=""))
    print(paste("Number of session 1 subjects identified from session 2 target: ",as.character(sum(df_ident["2_targeted_identification"])),sep=""))
    write.csv(df_perm,file.path(paths_$output,"output",paste("atl-",atlas,"_fp_permutation.csv",sep="")),row.names=F)
  }
  write.csv(df_out_combined,file.path(paths_$output,"output","fp_identification_summary.csv"),row.names=F)
  print("Finished identify_fp().")
}
