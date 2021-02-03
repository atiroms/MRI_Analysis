#source('D:/atiro/GitHub/MRI_Analysis/analyze/connection_cs.R')
source('C:/Users/NICT_WS/GitHub/MRI_Analysis/analyze/connection_cs.R')

####

paths_=paths
list_atlas_=list_atlas
key_group_='group_3'
subset_subj_=sex_diff_fc_cs_subset_subj
list_covar_=sex_diff_fc_cs_list_covar
list_mod_cs_=sex_diff_fc_cs_list_mod_cs
list_mod_diff_=sex_diff_fc_cs_list_mod_diff
list_mod_long_=sex_diff_fc_cs_list_mod_long
list_plot_=sex_diff_fc_cs_list_plot
thr_p_cdt_=sex_diff_fc_cs_thr_p_cdt
thr_p_perm_=sex_diff_fc_cs_thr_p_perm
n_perm_=sex_diff_fc_cs_n_perm

####

print("Starting sex_diff_fc_cs()")
nullobj<-func_createdirs(paths_,str_proc="sex_diff_fc_cs()",copy_log=T)
# Increase memory limit
memory.limit(1000000)

####

#atlas<-list_atlas[1]
atlas<-"aal116"

####

print(paste("Preparing data: ",atlas,sep=""))
data_fc<-prep_data_fc(paths_,atlas,key_group_,include_diff=T,include_grp=F)
data_clin<-func_clinical_data_long(paths_,list_wave=c("1","2"),subset_subj_,list_covar_,rem_na_clin=T,
                                   prefix=paste("atl-",atlas,sep=""),print_terminal=F)

# Longitudinal change analysis
df_fc_diff<-data_fc$df_fc[data_fc$df_fc$ses=="2-1",]
# List of subjects meeting QC criteria and has non-NA covariates in both waves
list_id_subj<-sort(intersect(data_clin$list_id_exist[[1]]$intersect,data_clin$list_id_exist[[2]]$intersect))
df_clin<-data_clin$df_clin
colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
df_clin_diff<-func_clinical_data_diffmean(df_clin,list_id_subj,list_covar_)
df_clin_diff$wave<-"2-1"
#df_clin_diff<-func_demean_clin(df_clin_diff,thr_cont=10)$df_clin
mean_mean_age<-mean(df_clin_diff$mean_age)
df_clin_diff$mean_age<-df_clin_diff$mean_age-mean_mean_age

# Network-based statistics
#func_nbs(paths=paths_,atlas=atlas,wave="diff",
#         df_fc=df_fc_diff,df_clin=df_clin_diff,
#         list_mod=list_mod_diff_,list_plot=list_plot_,list_sex=list(c(1,2)),
#         df_roi=data_fc$df_roi,df_edge=data_fc$df_edge,df_grp=data_fc$df_grp,
#         thr_p_cdt=thr_p_cdt_,n_perm=n_perm_,thr_p_perm=thr_p_perm_,
#         calc_parallel=T,test_mod=F)

####

#test_mod<-func_nbs(paths=paths_,atlas=atlas,wave="diff",
#         df_fc=df_fc_diff,df_clin=df_clin_diff,
#         list_mod=list_mod_diff_,list_plot=list_plot_,list_sex=list(c(1,2)),
#         df_roi=data_fc$df_roi,df_edge=data_fc$df_edge,df_grp=data_fc$df_grp,
#         thr_p_cdt=thr_p_cdt_,n_perm=n_perm_,thr_p_perm=thr_p_perm_,
#         calc_parallel=T,test_mod=T)
#summary(test_mod[[1]])
#summary(test_mod[[2]])


####

paths=paths_
wave="diff"
df_fc=df_fc_diff
df_clin=df_clin_diff
list_mod=list_mod_diff_
list_plot=list_plot_
list_sex=list(c(1,2))
df_roi=data_fc$df_roi
df_edge=data_fc$df_edge
df_grp=data_fc$df_grp
thr_p_cdt=thr_p_cdt_
n_perm=n_perm_
thr_p_perm=thr_p_perm_
calc_parallel=T
test_mod=F

####
print(paste("Calculating model, atlas: ",atlas,", wave: ",wave,sep=""))
clust<-makeCluster(floor(detectCores()*3/4))
clusterExport(clust,
              varlist=c("list_mod","list_sex","calc_parallel","test_mod","sort","gam","as.formula","summary.gam",
                        "anova.gam","as.numeric.factor"),
              envir=environment())
data_nbs<-func_nbs_core(clust=clust,df_fc=df_fc,df_clin=df_clin,
                        df_roi=df_roi,df_edge=df_edge,list_mod=list_mod,
                        thr_p_cdt=thr_p_cdt,list_plot=list_plot,
                        progressbar=F,output_gamm=F,calc_slope=T,test_mod=test_mod)
if(test_mod){
  return(data_nbs)
}else{
  data_nbs<-data_nbs$data_nbs
  
  # Permutation test
  print(paste("Calculating permutation, atlas: ",atlas,", wave: ",wave,sep=""))
  set.seed(0)
  list_max<-list()
  pb<-txtProgressBar(min=0,max=n_perm,style=3,width=50)
  for (idx_perm in seq(n_perm)){
    df_clin_perm<-df_clin
    df_clin_perm$sex<-sample(df_clin_perm$sex)
    data_nbs_perm<-func_nbs_core(clust=clust,df_fc=df_fc,df_clin=df_clin_perm,
                                 df_roi=df_roi,df_edge=df_edge,list_mod=list_mod,
                                 thr_p_cdt=thr_p_cdt,list_plot=list_plot,
                                 progressbar=F,output_gamm=F,calc_slope=T)$data_nbs
    for (model in names(data_nbs_perm)){
      for (plot in names(data_nbs_perm[[model]])){
        if(idx_perm==1){
          list_max_sex<-list("m"=data_nbs_perm[[model]][[plot]][["m"]][["max_size"]],
                             "f"=data_nbs_perm[[model]][[plot]][["f"]][["max_size"]])
          list_max[[model]][[plot]]<-list_max_sex
        }else{
          for (sex in c("m","f")){
            list_max[[model]][[plot]][[sex]]<-c(list_max[[model]][[plot]][[sex]],
                                                data_nbs_perm[[model]][[plot]][[sex]][["max_size"]])
          }
        }
      }
    }
    setTxtProgressBar(pb,idx_perm)
  }
  close(pb)
  
  # Summarize permutation and threshold subgraphs
  print(paste("Preparing output, atlas: ",atlas,", wave: ",wave,sep=""))
  df_net<-df_size_net<-df_perm<-df_thr_size<-NULL
  list_output<-list()
  for (model in names(data_nbs)){
    for (plot in names(data_nbs[[model]])){
      for (sex in c("m","f")){
        data_nbs_subset<-data_nbs[[model]][[plot]][[sex]]
        list_max_subset<-list_max[[model]][[plot]][[sex]]
        list_max_subset_sort<-sort(list_max_subset)
        thr_size_perm<-list_max_subset_sort[ceiling(length(list_max_subset_sort)*(1-thr_p_perm))]
        
        if (sex=="m"){
          title_sex<-"M>F"
          color_plt<-"steelblue2"
        }else{
          title_sex<-"F>M"
          color_plt<-"lightcoral"
        }
        title_plot<-list_plot[[plot]][["title"]]
        #list_output<-c(list_output,
        #               list(plot_permutation(paths,list_max=list_max_subset_sort,thr_size_perm,
        #                                     atlas,wave,model,plot,sex,title_plot,title_sex,color_plt)))
        plot_permutation(paths,list_max=list_max_subset_sort,thr_size_perm,
                         atlas,wave,model,plot,sex,title_plot,title_sex,color_plt)
        df_head<-data.frame(atlas=atlas,wave=wave,mod=model,plot=plot,sex=sex)
        list_network_sign<-list()
        if(length(data_nbs_subset$list_network)>0){
          for (idx_net in seq(length(data_nbs_subset$list_network))){
            network<-data_nbs_subset$list_network[[idx_net]]
            list_output<-c(list_output,
                           list(plot_sex_diff_fc(paths,network$df_edge,atlas,wave,df_roi,df_grp,
                                                 model,plot,sex,title_plot,title_sex,idx_net)))
            df_net<-rbind(df_net,cbind(df_head, data.frame(id_net=idx_net,network$df_edge)))
            df_size_net<-rbind(df_size_net,cbind(df_head,data.frame(id_net=idx_net,size=network$size_net)))
            if (network$size_net>=thr_size_perm){
              list_network_sign<-c(list(network))
            }
          }
        }
        data_nbs[[model]][[plot]][[sex]][["list_network_sign_perm"]]<-list_network_sign
        data_nbs[[model]][[plot]][[sex]][["list_max_size_perm"]]<-list_max_subset
        data_nbs[[model]][[plot]][[sex]][["thr_size_perm"]]<-thr_size_perm
        df_thr_size<-rbind(df_thr_size,
                           cbind(df_head,data.frame(thr_size=thr_size_perm)))
        df_perm<-rbind(df_perm,
                       cbind(df_head,data.frame(id_perm=seq(length(list_max_subset)),
                                                max_size_net=list_max_subset)))
      }
    }
  }
  plot_parallel(clust,list_output)
  stopCluster(clust)
  write.csv(df_net,file.path(paths$output,"output","temp",
                             paste("atl-",atlas,"_wave-",wave,"_net.csv",sep="")),row.names=F)
  write.csv(df_size_net,file.path(paths$output,"output","temp",
                                  paste("atl-",atlas,"_wave-",wave,"_size_net.csv",sep="")),row.names=F)
  write.csv(df_thr_size,file.path(paths$output,"output","temp",
                                  paste("atl-",atlas,"_wave-",wave,"_thr_perm.csv",sep="")),row.names=F)
  write.csv(df_perm,file.path(paths$output,"output","temp",
                              paste("atl-",atlas,"_wave-",wave,"_perm.csv",sep="")),row.names=F)
}

####
####

plot_long<-(ggplot(df_src,aes(x=age,y=value,color=sex))
            + geom_point()
            + scale_color_manual(values = c("steelblue2", "lightcoral"),labels=c("M","F"))
            + xlab("age(d)")
            + ylab("z(r)")
            + geom_path(aes(x=age,y=value,group=ID_pnTTC),
                        color="grey",size=0.3,linetype="dashed")
            + geom_smooth(method="lm",se=F,fullrange=T)
            + theme_light()
            
)
plot_long

df_src$abs_value<-abs(df_src$value)
plot_long_abs<-(ggplot(df_src,aes(x=age,y=abs_value,color=sex))
            + geom_point()
            + scale_color_manual(values = c("steelblue2", "lightcoral"),labels=c("M","F"))
            + xlab("age(d)")
            + ylab("z(r)")
            + geom_path(aes(x=age,y=abs_value,group=ID_pnTTC),
                        color="grey",size=0.3,linetype="dashed")
            + geom_smooth(method="lm",se=F,fullrange=T)
            + theme_light()
            
)
plot_long_abs

id_long<-sort(as.numeric(intersect(unique(as.character(df_src[df_src$wave==1,"ID_pnTTC"])),
                        unique(as.character(df_src[df_src$wave==2,"ID_pnTTC"])))))


####

# Longitudinal change analysis
df_fc_diff<-data_fc$df_fc[data_fc$df_fc$ses=="2-1",]
# List of subjects meeting QC criteria and has non-NA covariates in both waves
list_id_subj<-sort(intersect(data_clin$list_id_exist[[1]]$intersect,data_clin$list_id_exist[[2]]$intersect))
df_clin<-data_clin$df_clin
colnames(df_clin)[colnames(df_clin)=="wave"]<-"ses"
df_clin_diffmean<-func_clinical_data_diffmean(df_clin,list_id_subj,list_covar_)
df_clin_diffmean$wave<-"2-1"

# Network-based statistics
print(paste("Calculating model: ",atlas,sep=""))
data_nbs<-func_nbs(df_fc=df_fc_diff,df_clin=df_clin_diffmean,
                   df_roi=data_fc$df_roi,list_mod=list_mod_diff_,
                   thr_p_cdt=thr_p_cdt_,list_plot=list_plot_,
                   progressbar=F,output_gamm=T,calc_slope=T)
data_nbs<-data_nbs$data_nbs

####

df_temp_s<-data_nbs[["int"]][["s"]][["f"]][["list_network"]][[1]][["df_edge"]]
df_temp_s<-df_temp_s[order(df_temp_s$p),]
df_temp_sxm<-data_nbs$int[["sxm(a)"]]$m$list_network[[1]][["df_edge"]]
df_temp_sxm<-df_temp_sxm[order(df_temp_sxm$p),]

df_join<-join_fc_clin(df_fc_diff,df_clin_diffmean)
df_join$value<-df_join$value/df_join$diff_age

list_plot_mean<-list()
for (idx_row in seq(nrow(df_temp_sxm))){
  df_join_subset<-df_join[df_join$from==df_temp_sxm[idx_row,"from"] & df_join$to==df_temp_sxm[idx_row,"to"],]
  plot_mean<-(ggplot(df_join_subset,aes(x=mean_age,y=value,color=sex))
              + geom_point()
              + scale_color_manual(values = c("steelblue2", "lightcoral"),labels=c("M","F"))
              + xlab("mean(age)(d)")
              + ylab("delta(z(r))")
              + geom_smooth(method="lm",se=F,fullrange=T)
              + theme_light()
  )
  list_plot_mean<-c(list_plot_mean,list(plot_mean))
}

#df_join_subset<-df_join[df_join$from=="aal116_00011" & df_join$to=="aal116_00044",]
#df_join_subset_m<-df_join_subset[df_join_subset$sex=="1",]
#df_join_subset_f<-df_join_subset[df_join_subset$sex=="2",]



plot_ses1<-(ggplot(df_join_subset,aes(x=ses1_age,y=value,color=sex))
            + geom_point()
            + scale_color_manual(values = c("steelblue2", "lightcoral"))
            + geom_smooth(method="lm")
)
plot_ses1
plot_ses2<-(ggplot(df_join_subset,aes(x=ses2_age,y=value,color=sex))
            + geom_point()
            + scale_color_manual(values = c("steelblue2", "lightcoral"))
            + geom_smooth(method="lm")
)
plot_ses2

plot_diff<-(ggplot(df_join_subset,aes(x=diff_age,y=value,color=sex))
            + geom_point()
            + scale_color_manual(values = c("steelblue2", "lightcoral"))
            + geom_smooth(method="lm")
)
plot_diff

####



####

data_gamm<-data_nbs$data_gamm
df_gamm<-data_gamm$df_out_gamm
df_gamm_subset<-df_gamm[df_gamm$model=="int" & df_gamm$from=="aal116_00011" & df_gamm$to=="aal116_00044",]
rownames(df_gamm_subset)<-NULL



mean(df_join_subset[df_join_subset$sex=="2","value"])-mean(df_join_subset[df_join_subset$sex=="1","value"])
df_gamm_subset[2,"estimate"]-df_gamm_subset[1,"estimate"]

####

mod1<-gam(value ~ sex + mean_age + sex*mean_age,data=df_join_subset,method="REML",control=list(nthreads=1))
summary(mod1)

mod2<-gam(value ~ sex*mean_age,data=df_join_subset,method="REML",control=list(nthreads=1))
summary(mod2)

mod_m<-gam(value ~ mean_age,data=df_join_subset_m,method="REML",control=list(nthreads=1))
summary(mod_m)

mod_f<-gam(value ~ mean_age,data=df_join_subset_f,method="REML",control=list(nthreads=1))
summary(mod_f)

