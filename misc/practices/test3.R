

#setMKLthreads(1)
#getMKLthreads()

####

source('D:/atiro/GitHub/MRI_Analysis/analyze/connection_cs.R')
#source('C:/Users/NICT_WS/GitHub/MRI_Analysis/analyze/connection_cs.R')

####

paths_=paths
list_atlas_=list_atlas
key_group_='group_3'
subset_subj_=sex_diff_fc_subset_subj
list_covar_=sex_diff_fc_list_covar
list_mod_diff_=sex_diff_fc_list_mod_diff
list_mod_long_=sex_diff_fc_list_mod_long
list_mod_cs_=sex_diff_fc_list_mod_cs
list_plot_=sex_diff_fc_list_plot
thr_p_cdt_=sex_diff_fc_thr_p_cdt
thr_p_perm_=sex_diff_fc_thr_p_perm
n_perm_=sex_diff_fc_n_perm

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
df_clin_diff<-data_clin$df_clin
colnames(df_clin_diff)[colnames(df_clin_diff)=="wave"]<-"ses"
df_clin_diff<-func_clinical_data_diffmean(df_clin_diff,list_id_subj,list_covar_)
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


df_fc_long<-data_fc$df_fc[data_fc$df_fc$ses %in% c(1,2),]
df_clin_long<-data_clin$df_clin
#func_nbs(paths=paths_,atlas=atlas,wave="long",
#         df_fc=df_fc_long,df_clin=df_clin_long,
#         list_mod=list_mod_long_,calc_slope=F,list_plot=list_plot_,list_sex=list(c(1,2)),
#         df_roi=data_fc$df_roi,df_edge=data_fc$df_edge,df_grp=data_fc$df_grp,
#         thr_p_cdt=thr_p_cdt_,n_perm=n_perm_,thr_p_perm=thr_p_perm_,
#         calc_parallel=T,test_mod=F)

####

paths=paths_
atlas=atlas
#wave="long"
#df_fc=df_fc_long
#df_clin=df_clin_long
#list_mod=list_mod_long_
#calc_slope=F
wave="diff"
df_fc=df_fc_diff
df_clin=df_clin_diff
list_mod=list_mod_diff_
calc_slope=T
list_plot=list_plot_
list_sex=list(c(1,2))
df_roi=data_fc$df_roi
df_edge=data_fc$df_edge
df_grp=data_fc$df_grp
thr_p_cdt=thr_p_cdt_
n_perm=n_perm_
thr_p_perm=thr_p_perm_
calc_parallel=T
#calc_parallel=F
test_mod=F
#test_mod=T

####

print(paste("Calculating model, atlas: ",atlas,", wave: ",wave,sep=""))
if (calc_parallel){
  clust<-makeCluster(floor(detectCores()*3/4))
}else{
  clust<-makeCluster(1)
}
clusterExport(clust,
              varlist=c("list_mod","list_sex","calc_parallel","test_mod","sort","gam","as.formula","summary.gam",
                        "anova.gam","as.numeric.factor"),
              envir=environment())
#data_nbs<-func_nbs_core(clust=clust,df_fc=df_fc,df_clin=df_clin,
#                        df_roi=df_roi,df_edge=df_edge,list_mod=list_mod,
#                        thr_p_cdt=thr_p_cdt,list_plot=list_plot,
#                        progressbar=F,output_gamm=F,calc_slope=calc_slpoe,test_mod=test_mod)

####

progressbar=F
output_gamm=F

####

df_join<-join_fc_clin(df_fc,df_clin)
if(calc_slope){
  # value as slope of z_r longitudinal difference against age (z(r(wave=2))-z(r(wave=1)))/delta(age)
  df_join$value<-df_join$value/df_join$diff_age
}
df_edge$label_from<-df_edge$label_to<-NULL
df_edge$id_edge<-seq(nrow(df_edge))
#data_gamm<-iterate_gamm3(clust,df_join,df_edge,progressbar=progressbar,test_mod=test_mod)


####
####

df_join<-inner_join(df_join,df_edge,by=c("from","to"))
list_src_gamm<-split(df_join,df_join$id_edge)

if(test_mod){
  list_src_gamm<-list_src_gamm[1]
}

if (progressbar){
  list_dst_gamm<-pblapply(list_src_gamm,gamm_core3,cl=clust)
}else{
  list_dst_gamm<-parLapply(clust,list_src_gamm,gamm_core3)
}
df_gamm<-rbindlist(ListExtract(list_dst_gamm,"df_gamm"))
df_aic<-rbindlist(ListExtract(list_dst_gamm,"df_aic"))
df_anova<-rbindlist(ListExtract(list_dst_gamm,"df_anova"))
rownames(df_gamm)<-rownames(df_aic)<-rownames(df_anova)<-NULL





####

df_join$id_edge<-factor(df_join$id_edge)

df_join_1<-df_join[df_join$id_edge==1,]
df_join_2<-df_join[df_join$id_edge==2,]
df_join_12<-df_join[df_join$id_edge %in% c(1,2),]
df_join_10<-df_join[df_join$id_edge %in% seq(10),]
df_join_100<-df_join[df_join$id_edge %in% seq(100),]

####

func_model<-function(df_src,formula_model){
  t_start<-Sys.time()
  #mod<-gam(as.formula(formula_model),data=df_src,method="REML")
  #mod<-bam(as.formula(formula_model),data=df_src,method="REML")
  mod<-bam(as.formula(formula_model),data=df_src,method="fREML",nthreads=1)
  t_elapsed<-Sys.time()-t_start
  return(list("mod"=mod,"t_elapsed"=t_elapsed))
}

mod<-func_model(df_join_1,"value ~ sex + age + s(ID_pnTTC,bs='re')")
mod<-func_model(df_join_2,"value ~ sex + age + s(ID_pnTTC,bs='re')")
mod<-func_model(df_join_12,"value ~ id_edge + sex:id_edge + age:id_edge + s(ID_pnTTC,id_edge,bs='re')")
mod<-func_model(df_join_12,"value ~ id_edge + sex:id_edge + age:id_edge + te(ID_pnTTC,id_edge,bs='re')")
mod<-func_model(df_join_10,"value ~ id_edge + sex:id_edge + age:id_edge + s(ID_pnTTC,id_edge,bs='re')")

summary(mod$mod)
print(mod$t_elapsed)
#plot(mod)


####

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

