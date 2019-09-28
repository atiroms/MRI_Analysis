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

dir_in<-"203_fp_acompcor"
dir_out<-"206_fp_model_acompcor_test"

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

#list_mod <- list("lin_diff_t"=
#                   "value ~ sex + sex:diff_tanner",
#                 "lin_diff_at"=
#                   "value ~ diff_age + sex + sex:diff_tanner",
#                 "lin_diff_at_mean_t"=
#                   "value ~ diff_age + sex + sex:mean_tanner + sex:diff_tanner",
#                 "lin_diff_a_ses_t"=
#                   "value ~ diff_age + sex + sex:ses1_tanner + sex:ses2_tanner",
#                 "add_diff_t"=
#                   "value ~ sex + s(diff_tanner,k=3,by=sex)",
#                 "add_diff_at"=
#                   "value ~ s(diff_age,k=3) + sex + s(diff_tanner,k=3,by=sex)",
#                 "add_diff_at_mean_t"=
#                   "value ~ s(diff_age,k=3) + sex + s(mean_tanner,k=3,by=sex) + s(diff_tanner,k=3,by=sex)",
#                 "add_diff_a_ses_t"=
#                   "value ~ s(diff_age,k=3) + sex + s(ses1_tanner,k=3,by=sex) + s(ses2_tanner,k=3,by=sex)")

list_mod <- list("lin_diff_t"=
                   "value ~ diff_tanner",
                 "lin_diff_at"=
                   "value ~ diff_age + diff_tanner",
                 "lin_diff_at_mean_t"=
                   "value ~ diff_age + mean_tanner + diff_tanner",
                 "lin_diff_a_ses_t"=
                   "value ~ diff_age + ses1_tanner + ses2_tanner",
                 "add_diff_t"=
                   "value ~ s(diff_tanner,k=3)",
                 "add_diff_at"=
                   "value ~ s(diff_age,k=3) + s(diff_tanner,k=3)",
                 "add_diff_at_mean_t"=
                   "value ~ s(diff_age,k=3) + s(mean_tanner,k=3) + s(diff_tanner,k=3)",
                 "add_diff_a_ses_t"=
                   "value ~ s(diff_age,k=3) + s(ses1_tanner,k=3) + s(ses2_tanner,k=3)")

list_graph <-list("a"=list("title"="Age diff effect",
                           "x_axis"="diff_age",
                           "smooth"=list("Male"=list("fix"=list("sex"=1),
                                                     "color"="steelblue2","alpha"=1,"ribbon"=T),
                                         "Female"=list("fix"=list("sex"=2),
                                                       "color"="lightcoral","alpha"=1,"ribbon"=T)),
                           "point"=list("Male"=list("subset"=list("sex"=1),
                                                    "color"="steelblue2","alpha"=1),
                                        "Female"=list("subset"=list("sex"=2),
                                                      "color"="lightcoral","alpha"=1))),
                  "tdiff"=list("title"="Tanner diff effect",
                            "x_axis"="diff_tanner",
                            "smooth"=list("Male"=list("fix"=list("sex"=1),
                                                      "color"="steelblue2","alpha"=1,"ribbon"=T),
                                          "Female"=list("fix"=list("sex"=2),
                                                        "color"="lightcoral","alpha"=1,"ribbon"=T)),
                            "point"=list("Male"=list("subset"=list("sex"=1),
                                                     "color"="steelblue2","alpha"=1),
                                         "Female"=list("subset"=list("sex"=2),
                                                       "color"="lightcoral","alpha"=1))),
                  "tmean"=list("title"="Tanner mean effect",
                            "x_axis"="mean_tanner",
                            "smooth"=list("Male"=list("fix"=list("sex"=1),
                                                      "color"="steelblue2","alpha"=1,"ribbon"=T),
                                          "Female"=list("fix"=list("sex"=2),
                                                        "color"="lightcoral","alpha"=1,"ribbon"=T)),
                            "point"=list("Male"=list("subset"=list("sex"=1),
                                                     "color"="steelblue2","alpha"=1),
                                         "Female"=list("subset"=list("sex"=2),
                                                       "color"="lightcoral","alpha"=1))),
                  "t1"=list("title"="1st Tanner effect",
                            "x_axis"="ses1_tanner",
                            "smooth"=list("Male"=list("fix"=list("sex"=1),
                                                      "color"="steelblue2","alpha"=1,"ribbon"=T),
                                          "Female"=list("fix"=list("sex"=2),
                                                        "color"="lightcoral","alpha"=1,"ribbon"=T)),
                            "point"=list("Male"=list("subset"=list("sex"=1),
                                                     "color"="steelblue2","alpha"=1),
                                         "Female"=list("subset"=list("sex"=2),
                                                       "color"="lightcoral","alpha"=1))),
                  "t2"=list("title"="2nd Tanner effect",
                            "x_axis"="ses2_tanner",
                            "smooth"=list("Male"=list("fix"=list("sex"=1),
                                                      "color"="steelblue2","alpha"=1,"ribbon"=T),
                                          "Female"=list("fix"=list("sex"=2),
                                                        "color"="lightcoral","alpha"=1,"ribbon"=T)),
                            "point"=list("Male"=list("subset"=list("sex"=1),
                                                     "color"="steelblue2","alpha"=1),
                                         "Female"=list("subset"=list("sex"=2),
                                                       "color"="lightcoral","alpha"=1))))

list_tanner <-list("5by5"=
                     list("1"=list("1"=1,"2"=2,"3"=3,"4"=4,"5"=5),
                          "2"=list("1"=1,"2"=2,"3"=3,"4"=4,"5"=5)),
                   "3by3"=
                     list("1"=list("12"=c(1,2),"3"=3,"45"=c(4,5)),
                          "2"=list("12"=c(1,2),"3"=3,"45"=c(4,5))),
                   "2by2"=
                     list("1"=list("12"=c(1,2),"345"=c(3,4,5)),
                          "2"=list("123"=c(1,2,3),"45"=c(4,5))))



#list_atlas<-c("aal116","glasser360","gordon333","power264","schaefer100","schaefer200","schaefer400")
list_atlas<-"aal116"
#list_atlas<-"schaefer400"
#list_atlas<-"dk"
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
library(multcomp)
library(parallel)
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

glm_core<-function(df_src,atlas,measure,group,list_mod_,list_graph_,list_covar_,paths_){
  print(paste("Atlas: ",atlas,", Measure: ",measure,", Group: ",group,", GLM/GAM.",  sep=""))
  df_out_aic_add<-df_out_lm_add<-data.frame()
  for (idx_mod in names(list_mod_)){
    list_plot<-list()
    list_sex<-sort(unique(as.numeric.factor(df_src$sex)))
    for (idx_sex in list_sex){
      df_src_sex<-df_src[df_src$sex==idx_sex,]
      mod<-gam(as.formula(list_mod_[[idx_mod]]),data=df_src_sex)
      p_table<-summary.gam(mod)$p.table
      if (is.null(summary.gam(mod)$s.table)){
        df_out_lm_add_add<-data.frame(atlas=atlas,measure=measure,group=group,sex=idx_sex,model=idx_mod,term=rownames(p_table),
                                  estimate=p_table[,'Estimate'],se=p_table[,'Std. Error'],F=NA,
                                  t=p_table[,'t value'],p=p_table[,'Pr(>|t|)'])
        
      }else{
        s_table<-summary.gam(mod)$s.table
        df_out_lm_add_add<-rbind(data.frame(atlas=atlas,measure=measure,group=group,sex=idx_sex,model=idx_mod,term=rownames(p_table),
                                        estimate=p_table[,'Estimate'],se=p_table[,'Std. Error'],F=NA,
                                        t=p_table[,'t value'],p=p_table[,'Pr(>|t|)']),
                             data.frame(atlas=atlas,measure=measure,group=group,sex=idx_sex,model=idx_mod,term=rownames(s_table),
                                        estimate=NA,se=NA,F=s_table[,'F'],
                                        t=NA,p=s_table[,'p-value']))
      }
      df_out_lm_add<-rbind(df_out_lm_add,df_out_lm_add_add)
      df_out_aic_add<-rbind(df_out_aic_add,
                            data.frame(atlas=atlas,measure=measure,group=group,sex=idx_sex,
                                       model=idx_mod,aic=mod$aic,aic_best_among_models=0))
      
      # Graphical output of GLM results
      for (idx_graph in names(list_graph_)){
        if (list_graph_[[idx_graph]][["x_axis"]] %in% colnames(mod$model)){
          # Add sex-wise lines/plots to existent plot, initialize if absent
          plot<-plot_gamm(plot_in=list_plot[[idx_graph]],mod_gamm=mod,
                          df_join_measure_roi=df_src_sex,
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
                   + ggtitle(paste(list_graph_[[idx_graph]][["title"]],atlas,measure,group,idx_mod,sep=" "))
                   + xlab(label_x)
                   + ylab("Fingerprint correlation")
                   + theme(legend.position = "none"))
            filename_plot<-paste("atl-",atlas,"_msr-",measure,"_grp-",group,"_mod-",idx_mod,
                                 "_plt-",idx_graph,"_fp_glm.eps",sep="")
            ggsave(filename_plot,plot=plot,device=cairo_ps,
                   path=file.path(paths_$output,"output"),dpi=300,height=5,width=5,limitsize=F)
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
  
  return(list("df_out_lm_add"=df_out_lm_add,"df_out_aic_add"=df_out_aic_add_sex_rbind))
}

ancova_core<-function(data_input){
  atlas<-data_input$atlas
  measure<-data_input$measure
  group=data_input$group_network
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
  df_out_ancova<-data.frame(atlas=atlas,measure=measure,group=group,tanner=group_tanner,sex=idx_sex,test="ANCOVA",
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
                        data.frame(atlas=atlas,measure=measure,group=group,tanner=group_tanner,sex=idx_sex,test="fit",
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
      df_calc_mean<-data.frame(df_out_fit[df_out_fit$term==term & df_out_fit$label==list_label,c("label","fit")])
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
                                 + ggtitle(paste("FP Cor Model,",atlas,measure,group,group_tanner,idx_sex,sep=" "))
                                 + xlab("2nd wave")
                                 + ylab("1st wave")
                                 + theme(plot.title = element_text(hjust = 0.5),
                                         axis.text.x = element_text(size=8,angle = 0,vjust=0,hjust=0.5),
                                         axis.text.y = element_text(size=8))))
  
  ggsave(paste("atl-",atlas,"_msr-",measure,"_grp-",group,"_tan-",group_tanner,"_sex-",idx_sex,"_fp_ancova.eps",sep=""),plot=plot_ancova,device=cairo_ps,
         path=file.path(paths_$output,"output"),dpi=300,height=5,width=5,limitsize=F)
  
  # Calculate Tukey-Kramer
  tk<-summary(glht(mod_ancova, linfct = mcp('long_tanner' = 'Tukey')))$test
  df_out_posthoc<-data.frame(atlas=atlas,measure=measure,group=group,tanner=group_tanner,
                             sex=idx_sex,test="Tukey-Kramer",
                             term="long_tanner",label=names(tk$coefficients),
                             p=tk$pvalues[1:length(tk$coefficients)],t=tk$tstat,F=NA,
                             fit=NA,diff=tk$coefficients,sigma=tk$sigma)
  df_out<-rbind(df_out_ancova,df_out_fit,df_out_posthoc)
  return(df_out)
}

model_fp<-function(paths_=paths,
                   list_atlas_=list_atlas,
                   list_wave_=list_wave,
                   list_covar_=list_covar,
                   list_mod_=list_mod,
                   list_graph_=list_graph,
                   list_tanner_=list_tanner,
                   subset_subj_=subset_subj
                   ){
  print("Starting model_fp().")
  nullobj<-func_createdirs(paths_)
  
  # Load and subset clinical data according to specified subsetting condition and covariate availability
  print('Loading clinical data.')
  data_clin<-func_clinical_data_long(paths_,list_wave_,subset_subj_,
                                     list_covar=list_covar_,rem_na_clin=T)
  df_clin<-data_clin$df_clin
  colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
  
  df_out_lm<-df_out_aic<-NULL
  list_src_ancova<-NULL
  
  for (atlas in list_atlas_){
    print(paste("Atlas: ",atlas,sep=""))
    
    # Load fingerprint data
    df_fp<-read.csv(file.path(paths_$input,"output",paste("atl-",atlas,"_fp.csv",sep="")))
    
    list_measure<-sort(unique(df_fp$measure))
    df_join<-NULL
    for (measure in list_measure){
      print(paste("Atlas: ",atlas,", Measure: ",measure,sep=""))
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
      
      list_group<-sort(unique(as.character(df_fp_meas$group)))
      if ("whole" %in% list_group){
        list_group<-c("whole",list_group[list_group!="whole"])
      }
                    
      for (group in list_group){
        # Collect longitudinal fp correlation data
        df_cor_fp<-data.frame(ID_pnTTC=list_id_subj_exist_twice)
        for (id_subj in list_id_subj_exist_twice){
          df_cor_fp[df_cor_fp$ID_pnTTC==id_subj,"value"]<-df_fp_meas[df_fp_meas$from_ID_pnTTC==id_subj & df_fp_meas$to_ID_pnTTC==id_subj &df_fp_meas$group==group,"r"]
        }
        
        # Subset those without longitudinal fp correlation
        list_id_subj_nonna<-df_cor_fp[!is.na(df_cor_fp$value),"ID_pnTTC"]
        df_cor_fp<-df_cor_fp[df_cor_fp$ID_pnTTC %in% list_id_subj_nonna,]
        n_id_subj_exist_twice<-length(list_id_subj_nonna)
        print(paste("Atlas: ",atlas,", Measure: ",measure,", Group: ",group,", ",as.character(n_id_subj_exist_twice)," subjects with longitudinal non-NA data.",sep=""))
        
        # Create dataframe for GLM analysis
        df_join_grp<-func_clinical_data_join(df_src=df_clin,
                                         list_id_subj=list_id_subj_nonna,
                                         list_covar=list_covar_)
        df_join_grp<-inner_join(df_join_grp,df_cor_fp,by="ID_pnTTC")
        df_join_grp$ID_pnTTC<-as.factor(df_join_grp$ID_pnTTC)
        df_join_grp$sex<-as.factor(df_join_grp$sex)
        df_join_grp<-cbind(measure=measure,group=group,df_join_grp)
        
        df_join<-rbind(df_join,df_join_grp)
        
        # Calculate GLM
        out_glm<-glm_core(df_src=df_join_grp,atlas,measure,group,
                          list_mod_,list_graph_,list_covar_,paths_)
        df_out_lm<-rbind(df_out_lm,out_glm$df_out_lm_add)
        df_out_aic<-rbind(df_out_aic,out_glm$df_out_aic_add)
        
        # Prepare ANCOVA calculation for later parallel computing
        print(paste("Atlas: ",atlas,", Measure: ",measure,", Group: ",group,", ANCOVA preparation.",  sep=""))
        # Create list of input dataframes for parallel ANCOVA calculation
        for (group_tanner in names(list_tanner_)){
          # group by longitudinal Tanner stage
          df_join_grp_tanner<-df_join_grp
          for (ses in c(1,2)){
            list_tanner_ses<-names(list_tanner_[[group_tanner]][[as.character(ses)]])
            for (label_tanner in list_tanner_ses){
              list_tanner_ses_group<-list_tanner_[[group_tanner]][[as.character(ses)]][[label_tanner]]
              #print(list_tanner_ses_group)
              df_join_grp_tanner[df_join_grp_tanner[[paste('ses',as.character(ses),'_tanner',sep='')]] %in% list_tanner_ses_group,
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
                               list(list("atlas"=atlas,
                                         "measure"=measure,
                                         "group_network"=group,
                                         "group_tanner"=group_tanner,
                                         "group_tanner_content"=list_tanner_[[group_tanner]],
                                         "group_sex"=idx_sex,
                                         "df_src_ancova"=df_join_grp_tanner_sex)))
          }
        }
      }
    }
    write.csv(df_join,file.path(paths_$output,"output",
                                paste("atl-",atlas,"_fp_model_src.csv",sep="")),row.names = F)
  }
  
  # Data saving
  rownames(df_out_lm)<-rownames(df_out_aic)<-NULL
  write.csv(df_out_lm, file.path(paths_$output,"output","fp_glm.csv"),row.names = F)
  write.csv(df_out_aic,file.path(paths_$output,"output","fp_glm_aic.csv"),row.names = F)
  
  # Parallel ANCOVA calculation
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
  list_df_ancova<-parLapply(clust,list_src_ancova,ancova_core)
  stopCluster(clust)
  df_out_ancova<-NULL
  for (df_ancova in list_df_ancova){
    df_out_ancova<-rbind(df_out_ancova,df_ancova)
  }
  
  # Data saving
  rownames(df_out_ancova)<-NULL
  write.csv(df_out_ancova,file.path(paths_$output,"output","fp_ancova.csv"),row.names=F)
  print("Finished model_fp()")
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
