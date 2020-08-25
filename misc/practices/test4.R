paths_=paths
list_waves_=list_waves
subset_subj_=subset_subj
list_atlas_=list_atlas
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


print("Starting gam_fc_cs_multi()")
nullobj<-func_createdirs(paths_,str_proc="gam_fc_cs_multi()",copy_log=T)


#waves<-names(list_waves_)[[1]]
waves<-names(list_waves_)[[3]]


wave_clin<-list_waves_[[waves]]$wave_clin
wave_mri<-list_waves_[[waves]]$wave_mri
print(paste("Clinical wave: ", wave_clin,", MRI wave: ",wave_mri,sep=""))

# Prepare subject subsetting condition (MRI QC criteria) according to specified waves
subset_subj_temp<-subset_subj_[[as.character(wave_mri)]]
subset_subj_temp<-list(subset_subj_temp)
names(subset_subj_temp)<-wave_clin

idx_tanner<-"max"

print(paste("Tanner type: ",list_tanner_[[idx_tanner]][["label"]],sep=""))
list_covar<-list_covar_tanner_
list_covar[["tanner"]]<-list_tanner_[[idx_tanner]]
suffix<-paste("wave-",waves,"_var-",idx_tanner,sep="")


paths_=paths_
subset_subj_=subset_subj_temp
list_covar_=list_covar
wave_clin_=wave_clin
wave_mri_=wave_mri
list_atlas_=list_atlas_
list_mod_=list_mod_tanner_
list_plot_=list_plot_tanner_
key_group_='group_3'
list_type_p_=list_type_p_
thr_p_=thr_p_
suffix_=suffix


print("Starting gam_fc_cs().")
nullobj<-func_createdirs(paths_,str_proc="gam_fc_cs()",copy_log=T)
dict_roi <- func_dict_roi(paths_)

# Load and subset clinical data according to specified subsetting condition and covariate availability
print('Loading clinical data.')
data_clin<-func_clinical_data_long(paths_,list_wave=wave_clin_,subset_subj_,
                                   list_covar=list_covar_,rem_na_clin=T,suffix=suffix_)
df_clin<-data_clin$df_clin

#atlas<-"power264"
atlas<-"aal116"


# Load ROI-wise FC data
print(paste('Loading FC data, atlas:',atlas,sep=' '))
#df_fc<-read.csv(file.path(paths_$input,'output',paste('atl-',atlas,'_fc.csv',sep='')))
df_fc<-fread(file.path(paths_$input,'output',paste('atl-',atlas,'_fc.csv',sep='')))
df_join<-join_fc_clin(df_fc,df_clin,wave_clin_,wave_mri_)
#write.csv(df_join,file.path(paths_$output,"output",paste("atl-",atlas,suffix_,"_src.csv",sep="")),
#          row.names=F)

# Calculate and save ROI-wise GAMM of FC
print(paste('Calculating GAM, atlas: ',atlas,sep=''))
list_roi<-sort(unique(c(as.character(df_join$from),as.character(df_join$to))))
df_roi<-dict_roi[is.element(dict_roi$id,list_roi),c("id","label",key_group_)]
colnames(df_roi)[colnames(df_roi)==key_group_]<-"group"


list_roi<-df_roi$id

# Prepare dataset for multi-core processing
print("Preparing dataset for parallel processing.")
list_src_gamm<-list()
for (id_from in list_roi[-length(list_roi)]){
  df_join_from<-df_join[df_join$from==id_from,]
  label_from<-as.character(df_roi[df_roi$id==id_from,"label"])
  for(id_to in list_roi[seq(which(list_roi==id_from)+1,length(list_roi))]){
    label_to<-as.character(df_roi[df_roi$id==id_to,"label"])
    df_src<-df_join_from[df_join_from$to==id_to,]
    list_src_gamm<-c(list_src_gamm,list(list("df_src"=df_src,
                                             "id_from"=id_from,"id_to"=id_to,
                                             "label_from"=label_from,"label_to"=label_to)))
  }
}

# Parallel processing
n_cluster<-min(floor(detectCores()*3/4),length(list_src_gamm))
clust<-makeCluster(n_cluster)
print(paste("Parallel processing,",as.character(n_cluster),"cores.",sep=" "))
clusterExport(clust,
              varlist=c("list_mod_","sort","gam","as.formula","summary.gam",
                        "as.numeric.factor"),
              envir=environment())
list_dst_gamm<-parLapply(clust,list_src_gamm,gamm_core)
stopCluster(clust)

# Collect data into dataframes
len_list<-length(list_dst_gamm)
len_sublist<-floor(sqrt(len_list)*2)
n_sublist<-ceil(len_list/len_sublist)
print(paste("Dividing results into", as.character(n_sublist), "sublists.",sep=" "))
list_dst_gamm_sub<-list()
for (idx_sublist in 1:n_sublist){
  #print(paste("Subgroup",as.character(idx_sublist),sep=" "))
  if (idx_sublist!=n_sublist){
    list_dst_gamm_sub<-c(list_dst_gamm_sub,
                         list(list_dst_gamm[((idx_sublist-1)*len_sublist+1):(idx_sublist*len_sublist)]))
  }else{
    list_dst_gamm_sub<-c(list_dst_gamm_sub,
                         list(list_dst_gamm[((idx_sublist-1)*len_sublist+1):len_list]))
  }
}
list_dst_gamm<-NULL
gc()

n_cluster<-floor(detectCores()*3/4)
clust<-makeCluster(n_cluster)
print(paste("Combining within sublists,",as.character(n_cluster),"cores.",sep=" "))
clusterExport(clust,
              varlist=NULL,
              envir=environment())
list_dst_gamm<-parLapply(clust,list_dst_gamm_sub,combine_gamm)
stopCluster(clust)

print("Combining sublists.")
df_out_gamm<-df_out_aic<-NULL
for (dst_gamm in list_dst_gamm){
  df_out_gamm<-rbind(df_out_gamm,dst_gamm$df_out_gamm_add)
  df_out_aic<-rbind(df_out_aic,dst_gamm$df_out_aic_add)
}
list_dst_gamm<-NULL
gc()

rownames(df_out_gamm)<-rownames(df_out_aic)<-NULL



data_gamm<-list("df_out_gamm"=df_out_gamm,"df_out_aic"=df_out_aic)


write.csv(data_gamm$df_out_gamm,
          file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_gam.csv",sep="")),row.names = F)
write.csv(data_gamm$df_out_aic,
          file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_gam_aic.csv",sep="")),row.names = F)

# Calculate multiple comparison-corrected p values
df_plot_gamm<-add_mltcmp(data_gamm$df_out_gamm,df_roi,analysis="roi",atlas,
                         list_mod_,list_plot_,calc_seed_level=T)
write.csv(df_plot_gamm,
          file.path(paths_$output,"output",paste("atl-",atlas,"_",suffix_,"_gam_plt.csv",sep="")),row.names = F)




analysis="roi"
list_mod=list_mod_
list_plot=list_plot_
list_type_p=list_type_p_
thr_p=thr_p_


idx_mod<-names(list_mod)[[1]]
idx_plot<-names(list_plot)[[1]]
var_exp<-list_plot[[idx_plot]][["var_exp"]]
idx_sex<-1


# Subset GAMM result dataframe for plotting
if (idx_sex==1){
  label_sex<-"m"
}else{
  label_sex<-"f"
}
df_plot_gamm_subset<-df_plot_gamm[df_plot_gamm$model==idx_mod 
                                  & df_plot_gamm$term==var_exp
                                  & df_plot_gamm$sex==idx_sex,]


print(paste("GAMM output, atlas: ",atlas,", model: ",idx_mod,", plot: ",var_exp,", sex: ",label_sex,sep=""))
# Convert GAMM rseult into igraph object
if (!is.na(df_plot_gamm_subset[1,"estimate"])){
  df_plot_gamm_subset<-rename(df_plot_gamm_subset,"weight"="estimate")
}else{
  df_plot_gamm_subset<-rename(df_plot_gamm_subset,"weight"="F")
}

# Convert FC dataframe into iGraph object
list_roi<-as.character(df_roi$id)
df_node<-data.frame(id=list_roi,stringsAsFactors = F)
for (idx_node in seq(dim(df_node)[1])){
  df_node[idx_node,"label"]<-as.character(df_roi[df_roi$id==df_node[idx_node,"id"],"label"])
}
df_edge<-df_plot_gamm_subset
df_edge$from<-as.character(df_edge$from)
df_edge$to<-as.character(df_edge$to)
igraph_gamm <- graph_from_data_frame(d = df_edge, vertices = df_node, directed = F)


type_p<-list_type_p[1]


plot<-plot_circular(igraph_in=igraph_gamm,
                    type_p=type_p,thr_p=thr_p,
                    limit_color=NULL)
plot<-plot +
  ggtitle(paste("GLM/GAM sex: ",label_sex, ", model: ",idx_mod,", expvar: ",var_exp,
                "\nanalysis: ",analysis," threshold: ",type_p,sep="")) +
  theme(plot.title = element_text(hjust = 0.5))


write_rds(plot,file.path(paths_$output,"output","save_plot.Rds"))

ggsave(paste("atl-",atlas,"_anl-",analysis,"_mod-",idx_mod,"_plt-",var_exp,
             "_sex-",label_sex,"_pval-",type_p,"_",suffix_,"_gamm_fc.eps",sep=""),
       plot=plot,device=cairo_ps,path=file.path(paths_$output,"output"),
       dpi=300,height=10,width=10,limitsize=F)

library(ggplot2)
library(colorRamps)
library(ggraph)
library(readr)
saved_plot<-read_rds(file.path("C:/Users/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP/425_fc_gam_cs_aroma/output","save_plot.Rds"))

paths_<-list("output"="C:/Users/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP/425_fc_gam_cs_aroma")
ggsave("test_gamm_fc.eps",
       plot=saved_plot,device=cairo_ps,path=file.path(paths_$output,"output"))

ggsave("test_gamm_fc.png",
       plot=saved_plot,path=file.path(paths_$output,"output"))

