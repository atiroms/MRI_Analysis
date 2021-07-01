# speed up TFNBS
t_start<-Sys.time()
data_tfnbs<-func_iterate_tfnbs(paths,df_gamm,df_anova,df_deltah_in=NULL,data_fc,plot_result=T,return_nbs=T,
                               atlas,param,list_mod,list_term,idx_var,label_wave)
print(Sys.time()-t_start)


####
####
list_plot<-list()
df_tfnbs<-data.frame(); df_max<-data.frame(); df_deltah<-data.frame()
####
idx_mod<-param$param_nbs$list_mod[[1]]
####
set_term<-param$param_nbs$list_term[[1]]
####
idx_term_detect<-set_term$term_detect[[1]]
var_exp_detect<-list_term[[idx_term_detect]]$var_exp
####
idx_sex<-param$list_sex[[1]]
####
if (idx_sex==1){label_sex<-"m"}else if (idx_sex==2){label_sex<-"f"}else{label_sex<-"mf"}
df_stat<-df_gamm[df_gamm$model==idx_mod & df_gamm$term==var_exp_detect & df_gamm$sex==idx_sex,]
if (nrow(df_stat)==0){
  df_anova_subset<-df_anova[df_anova$model==idx_mod & df_anova$term==var_exp_detect & df_anova$sex==idx_sex,]
  if (nrow(df_anova_subset)>0){
    # In case the term does not exist in df_gamm, plot using df_anova instead
    df_stat<-df_anova_subset
  }
}
if (nrow(df_stat)>0){
  if (is.na(df_stat[1,"F"])){ # GAMM result
    list_df_stat<-list("both"=df_stat,"pos"=df_stat[df_stat$estimate>0,],"neg"=df_stat[df_stat$estimate<0,])
  }else{ # ANCOVA result
    list_df_stat<-list("both"=df_stat,"pos"=data.frame(),"neg"=data.frame())
  }
}else{
  list_df_stat<-list("both"=data.frame(),"pos"=data.frame(),"neg"=data.frame())
}
####
type_sign<-names(list_df_stat)[[1]]
####
df_stat_temp<-list_df_stat[[type_sign]]
####
if (is.null(df_deltah_in)){
  delta_h_in<-NULL
}else{
  delta_h_in<-df_deltah_in[df_deltah_in$model==idx_mod & df_deltah_in$term==var_exp_detect
                           & df_deltah_in$sex==idx_sex & df_deltah_in$sign==type_sign,"delta_h"]
}

####

data_tfnbs<-func_tfnbs2(df_stat_temp,delta_h_in,param,clust)

data_tfnbs<-func_tfnbs(df_stat_temp,delta_h_in,param)


####
####

func_tfnbs_core<-function(src_tfnbs){
  df_stat_sign<-src_tfnbs$df_stat_sign
  thresh_h=src_tfnbs$thresh_h
  data_bfs<-func_bfs(df_stat_sign[!is.na(df_stat_sign$stat),])
  for (subnet in data_bfs$list_network){
    nbs_increment<-((subnet$size_net)**(param$param_tfnbs$e))*(thresh_h**(param$param_tfnbs$h))
    df_edge<-subnet$df_edge
    df_edge[df_edge$stat>=0,"stat"]<-nbs_increment
    df_edge[df_edge$stat<0,"stat"]<-(-1)*nbs_increment
    df_out<-left_join(df_stat_sign[,c("from","to")],df_edge,by=c("from","to"))
    df_out[is.na(df_out$stat),"stat"]<-0
    list_out=df_out$stat
  }
  return(list_out)
}

if (calc_parallel){clust<-makeCluster(floor(detectCores()*3/4))}else{clust<-makeCluster(1)}
clusterExport(clust,varlist=c("param","func_bfs","%nin%","left_join"),
              envir=environment())


func_tfnbs2<-function(df_stat,delta_h_in,param,clust){
  if (is.na(df_stat[1,"F"])){ # GAMM result
    df_stat<-df_stat[,c("from","to","t")]
    df_stat<-dplyr::rename(df_stat,"stat"="t")
  }else{ # ANCOVA result
    df_stat<-df_stat[,c("from","to","F")]
    df_stat<-dplyr::rename(df_stat,"stat"="F")
  }
  df_edge<-df_stat[,c("from","to")]
  max_thresh_h<-max(abs(df_stat$stat))
  if (is.null(delta_h_in)){
    delta_h<-max_thresh_h/param$param_tfnbs$n_thresh_h
  }else{
    delta_h<-delta_h_in
    #print(paste("calculate with delta_h:",as.character(delta_h),sep=''))
  }
  list_src_tfnbs<-list()
  for (thresh_h in seq(0,max_thresh_h,delta_h)){
    df_stat_sign<-df_stat
    df_stat_sign[abs(df_stat_sign$stat)<thresh_h,"stat"]<-NA
    list_src_tfnbs<-c(list_src_tfnbs,list(list("df_stat_sign"=df_stat_sign,"thresh_h"=thresh_h)))
  }

  list_dst_tfnbs<-parLapply(clust,list_src_tfnbs,func_tfnbs_core)
  df_out<-data.frame(df_stat[,c("from","to")],nbs=Reduce(`+`, list_dst_tfnbs))

  return(list("df_tfnbs"=df_out,"max_nbs"=max(abs(df_out$nbs)),"delta_h"=delta_h))
}

stopCluster(clust)
