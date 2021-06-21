# check ordered factor regression result

df_fwep<-fread("D:/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP/423.3_fc_gam_diff_aroma_test11/output/temp/atl-ho112_var-gonadal_wav-2-1_perm_fwep.csv")
df_fwep<-df_fwep[df_fwep$p_fwe<0.05,]


df_pred<-fread("D:/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP/423.3_fc_gam_diff_aroma_test11/output/temp/atl-ho112_var-gonadal_wav-2-1_gamm_pred.csv")

df_edge<-fread("D:/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/stats/func_XCP/423.3_fc_gam_diff_aroma_test11/output/temp/atl-ho112_var-gonadal_wav-2-1_bfs_edge.csv")
View(df_edge)


df_edge_subset<-df_edge[df_edge$sign=="pos" & df_edge$id_net==1 & df_edge$p_threshold==0.01 & df_edge$sex==1 & df_edge$term=="ses2_tanner.L",]

df_plot<-data.frame()
for (idx_row in seq(nrow(df_edge_subset))){
  node_from<-as.character(df_edge_subset[idx_row,"from"])
  node_to<-as.character(df_edge_subset[idx_row,"to"])
  df_plot<-rbind(df_plot,
                 data.frame(id_edge=idx_row,
                 df_pred[df_pred$from==node_from & df_pred$to==node_to & df_pred$ses1_tanner==1 & df_pred$sex==1,]))
}


plt <- (ggplot(data=df_plot)
         + geom_point(aes(x=ses2_tanner,y=prediction),color="grey50",shape=1,size=2)
         + geom_path(aes(x=ses2_tanner,y=prediction,group=id_edge),
                     color="grey50",size=0.3,alpha=0.5,linetype="dashed")
         #+ xlab(list_term[[idx_term]][["title"]])
         + ylab("Predicted z(r)")
         + theme_light()
         + theme(plot.title = element_text(hjust = 0.5)))

plt

