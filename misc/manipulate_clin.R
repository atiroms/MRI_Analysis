library(easypackages)
libraries('data.table',"dplyr","plyr","ggplot2","colorRamps","viridis","tidyverse")

path_src<-'D:/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/stats/clin/04_clin_combine'
path_common<-'D:/NICT_WS/Dropbox/MRI_img/pnTTC/puberty/common'


#**************************************************
# Factor to numeric function ======================
#**************************************************
as.numeric.factor <- function(x) {
  if (class(x)[1] %in% c("factor","ordered")){
    return(as.numeric(levels(x))[x])
  }else{
    return(x)
  }
}

#**************************************************
# Pick up required clinical data ==================
#**************************************************
df_clin<-data.frame(fread(file.path(path_src,'01_item','spss_data03.csv'),encoding='UTF-8'))
df_tanner<-data.frame(fread(file.path(path_src,'04_max_tanner','wavemax_tanner02.csv'),encoding='UTF-8'))
df_id<-data.frame(ID_pnTTC=seq(max(df_clin$ID_pnTTC,na.rm=T)))
df_clin<-dplyr::left_join(df_id,df_clin,by='ID_pnTTC')

# Calculate Handedness 1=right,2=left,3=both
# if one is missing use the other one
# if one is 3, the result is 3
# if one is 1 and the other is 2, the result is 3
df_temp<-df_clin[,c('ID_pnTTC','W1_Handedness','W2_Handedness')]
for (id_row in seq(nrow(df_temp))){
  if (is.na(df_temp[id_row,2])){
    if (is.na(df_temp[id_row,3])){
      df_temp[id_row,'Handedness']<-NA
    }else{
      df_temp[id_row,'Handedness']<-df_temp[id_row,3]
    }
  }else{
    if (is.na(df_temp[id_row,3])){
      df_temp[id_row,'Handedness']<-df_temp[id_row,2]
    }else{
      if (df_temp[id_row,2]==df_temp[id_row,3]){
        df_temp[id_row,'Handedness']<-df_temp[id_row,2]
      }else{
        df_temp[id_row,'Handedness']<-3
      }
    }
  }
}
df_temp$W1_Handedness<-df_temp$W2_Handedness<-NULL

df_clin$W1_Tanner_Female_Breast<-df_clin$W1_Tanner_Female_Pubic_Hair<-df_clin$W1_Tanner_Male_Genitals<-df_clin$W1_Tanner_Male_Pubic_Hair<-NULL
df_clin$W2_Tanner_Female_Breast<-df_clin$W2_Tanner_Female_Pubic_Hair<-df_clin$W2_Tanner_Male_Genitals<-df_clin$W2_Tanner_Male_Pubic_Hair<-NULL

df_clin<-dplyr::left_join(df_clin,df_tanner,by='ID_pnTTC')
df_clin<-dplyr::left_join(df_clin,df_temp,by='ID_pnTTC')

list_colname<-c('ID_pnTTC','Sex','W1_Handedness','W2_Handedness','Handedness',
                colnames(df_tanner)[2:length(colnames(df_tanner))],
                'W1_WISC_PC_Score','W1_WISC_IF_Score','W2_WISC_PC_Score','W2_WISC_IF_Score','W2_WISC_DS_Score','W2_WISC_CD_Score',
                'W1_CBCL_G1_Score','W1_CBCL_G2_Score','W1_CBCL_G3_Score','W1_CBCL_G4_Score','W1_CBCL_G5_Score',
                'W1_CBCL_G6_Score','W1_CBCL_G7_Score','W1_CBCL_G8_Score','W1_CBCL_Int_Score','W1_CBCL_Ext_Score',
                'W2_CBCL_G1_Score','W2_CBCL_G2_Score','W2_CBCL_G3_Score','W2_CBCL_G4_Score','W2_CBCL_G5_Score',
                'W2_CBCL_G6_Score','W2_CBCL_G7_Score','W2_CBCL_G8_Score','W2_CBCL_Int_Score','W2_CBCL_Ext_Score',
                'W1_SDQ_ES','W1_SDQ_CP','W1_SDQ_HI','W1_SDQ_PP','W1_SDQ_PB','W1_SDQ_TD',
                'W2_SDQ_ES','W2_SDQ_CP','W2_SDQ_HI','W2_SDQ_PP','W2_SDQ_PB','W2_SDQ_TD',
                'W1_Height','W1_Weight','W2_Height','W2_Weight')
df_clin<-df_clin[,list_colname]

fwrite(df_clin,file.path(path_src,'05_clin_combine','clin_combine01.csv'))


#**************************************************
# Calculate Tanner as past maximum ================
#**************************************************
#list_type_tanner<-c('Tanner_Female_Breast','Tanner_Female_Pubic_Hair','Tanner_Male_Genitals','Tanner_Male_Pubic_Hair')
df_tanner<-data.frame(fread(file.path(path_src,'03_excel','clin_combine01.csv'),encoding='UTF-8'))
df_tanner<-df_tanner[!is.na(df_tanner$Sex),]
df_tanner$Age<-as.numeric(df_tanner$Age)
for (idx_row in seq(nrow(df_tanner))){
  df_tanner[idx_row,'gonadal']<-as.numeric(max(df_tanner[idx_row,'Tanner_Female_Breast'],df_tanner[idx_row,'Tanner_Male_Genitals'],na.rm=T))
  df_tanner[idx_row,'adrenal']<-as.numeric(max(df_tanner[idx_row,'Tanner_Female_Pubic_Hair'],df_tanner[idx_row,'Tanner_Male_Pubic_Hair'],na.rm=T))
}
df_tanner$Tanner_Female_Breast<-df_tanner$Tanner_Female_Pubic_Hair<-df_tanner$Tanner_Male_Genitals<-df_tanner$Tanner_Male_Pubic_Hair<-NULL

list_id<-seq(max(df_tanner$ID_pnTTC))
df_dst<-data.frame()
for (id_subj in list_id){
  df_tanner_subj<-df_tanner[df_tanner$ID_pnTTC==id_subj,]
  df_mri_subj<-df_tanner_subj[df_tanner_subj$Source=='M',]
  df_tanner_subj<-df_tanner_subj[df_tanner_subj$Source!='M',]
  for (wave in df_mri_subj$Wave){
    age<-df_mri_subj[df_mri_subj$Wave==wave,'Age']
    df_tanner_subj_wave<-df_tanner_subj[df_tanner_subj$Wave==wave,]
    df_temp<-df_mri_subj[df_mri_subj$Wave==wave,]
    df_temp$Source<-df_temp$Timing<-df_temp$Age<-NULL
    df_replace<-df_tanner_subj_wave[df_tanner_subj_wave$Source=='E' & df_tanner_subj_wave$Timing=='M',c("gonadal","adrenal")]
    if (nrow(df_replace)>0){
      df_temp[1,c("gonadal","adrenal")]<-df_replace
    }
    for (type_tanner in c("gonadal","adrenal")){
      # max of all past data
      list_tanner<-as.numeric(df_tanner_subj[df_tanner_subj$Age<=age,type_tanner])
      if (any(!is.na(list_tanner))){
        max_tanner<-max(list_tanner,na.rm=T)
      }else{
        max_tanner<-NA
      }
      df_temp[1,paste('max',type_tanner,sep='_')]<-max_tanner
      
      list_tanner<-as.numeric(df_tanner_subj_wave[df_tanner_subj_wave$Age<=age,type_tanner])
      if (any(!is.na(list_tanner))){
        max_tanner<-max(list_tanner,na.rm=T)
      }else{
        max_tanner<-NA
      }
      df_temp[1,paste('wavemax',type_tanner,sep='_')]<-max_tanner
    }
    df_temp[1,'diff_max_gonadal']<-df_temp[1,'max_gonadal']-df_temp[1,'gonadal']
    df_temp[1,'diff_max_adrenal']<-df_temp[1,'max_adrenal']-df_temp[1,'adrenal']
    df_temp[1,'diff_wavemax_gonadal']<-df_temp[1,'wavemax_gonadal']-df_temp[1,'gonadal']
    df_temp[1,'diff_wavemax_adrenal']<-df_temp[1,'wavemax_adrenal']-df_temp[1,'adrenal']
    df_temp<-df_temp[,c('ID_pnTTC','Sex','Wave',
                        'gonadal','max_gonadal','wavemax_gonadal','diff_max_gonadal','diff_wavemax_gonadal',
                        'adrenal','max_adrenal','wavemax_adrenal','diff_max_adrenal','diff_wavemax_adrenal')]
    df_dst<-rbind(df_dst,df_temp)
  }
}
fwrite(df_dst,file.path(path_src,'04_max_tanner','max_tanner01.csv'))

# Organize in CSUB.csv format
df_src<-data.frame(fread(file.path(path_src,'04_max_tanner','max_tanner01.csv'),encoding='UTF-8'))
list_id<-seq(max(df_src$ID_pnTTC))
df_dst_max<-df_dst_wavemax<-data.frame(matrix(nrow=length(list_id),ncol=13))
list_colname<-'ID_pnTTC'
for (wave in c(1,2)){
  list_colname<-c(list_colname,paste('W',as.character(wave),'_',c('Tanner_Female_Breast','Tanner_Female_Pubic_Hair','Tanner_Male_Genitals','Tanner_Male_Pubic_Hair','Tanner_Max','Tanner_Full'),sep=''))
}
colnames(df_dst_max)<-colnames(df_dst_wavemax)<-list_colname
df_dst_max$ID_pnTTC<-df_dst_wavemax$ID_pnTTC<-list_id
for (id_row in seq(nrow(df_src))){
  id_subj<-df_src[id_row,'ID_pnTTC']
  wave<-df_src[id_row,'Wave']
  sex<-df_src[id_row,'Sex']
  if (sex==1){
    colname_dst<-paste("W",as.character(wave),'_',c('Tanner_Male_Genitals','Tanner_Male_Pubic_Hair','Tanner_Max','Tanner_Full'),sep='')
  }else{
    colname_dst<-paste("W",as.character(wave),'_',c('Tanner_Female_Breast','Tanner_Female_Pubic_Hair','Tanner_Max','Tanner_Full'),sep='')
  }
  df_src_row<-df_src[id_row,c("max_gonadal","max_adrenal")]
  if (any(!is.na(df_src_row))){
    df_src_row[,c('max','full')]<-c(max(df_src_row,na.rm=T),max(df_src_row,na.rm=F))
  }else{
    df_src_row[,c('max','full')]<-c(NA,NA)
  }
  df_dst_max[df_dst_max$ID_pnTTC==id_subj,colname_dst]<-df_src_row
  df_src_row<-df_src[id_row,c("wavemax_gonadal","wavemax_adrenal")]
  if (any(!is.na(df_src_row))){
    df_src_row[,c('max','full')]<-c(max(df_src_row,na.rm=T),max(df_src_row,na.rm=F))
  }else{
    df_src_row[,c('max','full')]<-c(NA,NA)
  }
  df_dst_wavemax[df_dst_wavemax$ID_pnTTC==id_subj,colname_dst]<-df_src_row
}

fwrite(df_dst_max,file.path(path_src,'04_max_tanner','max_tanner02.csv'))
fwrite(df_dst_wavemax,file.path(path_src,'04_max_tanner','wavemax_tanner02.csv'))


#**************************************************
# Tanner examination ==============================
#**************************************************
df_spss<-data.frame(fread(file.path(path_src,'03_excel','spss_data03.csv'),encoding='UTF-8'))
df_tanner_w1<-data.frame(fread(file.path(path_src,'03_excel','w1_excel_data.csv'),encoding='UTF-8'))
df_tanner_w2<-data.frame(fread(file.path(path_src,'03_excel','w2_excel_data.csv'),encoding='UTF-8'))
df_date<-data.frame(fread(file.path(path_src,'02_date','date_data01.csv'),encoding='UTF-8'))
df_csub<-data.frame(fread(file.path(path_common,'CSUB.csv')))

df_combine<-data.frame()

df_spss<-df_spss[!is.na(df_spss$ID_pnTTC),
                 c('ID_pnTTC','Sex','W1_Tanner_Female_Breast','W1_Tanner_Female_Pubic_Hair','W1_Tanner_Male_Genitals','W1_Tanner_Male_Pubic_Hair',
                   'W2_Tanner_Female_Breast','W2_Tanner_Female_Pubic_Hair','W2_Tanner_Male_Genitals','W2_Tanner_Male_Pubic_Hair')]
df_spss<-dplyr::left_join(df_spss,df_date,by='ID_pnTTC')

for (wave in c(1,2)){
  colnames_src<-paste('W',as.character(wave),'_',c('Age_Month','Tanner_Female_Breast','Tanner_Female_Pubic_Hair','Tanner_Male_Genitals','Tanner_Male_Pubic_Hair'),sep='')
  colnames_dst<-c('Age','Tanner_Female_Breast','Tanner_Female_Pubic_Hair','Tanner_Male_Genitals','Tanner_Male_Pubic_Hair')
  df_temp<-df_spss[,c('ID_pnTTC',colnames_src)]
  colnames(df_temp)<-c('ID_pnTTC',colnames_dst)
  df_temp<-data.frame('Source'='S','Wave'=wave,'Timing'='Q',df_temp)
  df_combine<-rbind(df_combine,df_temp)
  colnames_src<-paste('W',as.character(wave),'_MRI_Age_Month',sep='')
  df_temp<-df_spss[,c('ID_pnTTC',colnames_src)]
  colnames(df_temp)<-c('ID_pnTTC','Age')
  df_temp<-data.frame('Source'='M','Wave'=wave,'Timing'='M',df_temp)
  df_combine<-plyr::rbind.fill(df_combine,df_temp)
}

df_tanner_w1[df_tanner_w1=='na']<-NA
df_tanner_w1[df_tanner_w1==9999]<-NA
df_tanner_w2[df_tanner_w2=='na']<-NA
df_tanner_w2[df_tanner_w2==9999]<-NA
for (label_timing in c('Q','M','S')){
  colnames_src<-paste(label_timing,c('Age','77','78','79','80'),sep='_')
  colnames_dst<-c('ID_pnTTC','Age','Tanner_Female_Breast','Tanner_Female_Pubic_Hair','Tanner_Male_Genitals','Tanner_Male_Pubic_Hair')
  df_temp<-df_tanner_w1[,c('ID_pnTTC',colnames_src)]
  colnames(df_temp)<-colnames_dst
  df_temp<-data.frame('Source'='E','Wave'=1,'Timing'=label_timing,df_temp)
  df_combine<-rbind(df_combine,df_temp)
}
for (label_timing in c('M','S')){
  colnames_src<-paste(label_timing,c('Age','7','8','9','10'),sep='_')
  colnames_dst<-c('ID_pnTTC','Age','Tanner_Female_Breast','Tanner_Female_Pubic_Hair','Tanner_Male_Genitals','Tanner_Male_Pubic_Hair')
  df_temp<-df_tanner_w2[,c('ID_pnTTC',colnames_src)]
  colnames(df_temp)<-colnames_dst
  df_temp<-data.frame('Source'='E','Wave'=2,'Timing'=label_timing,df_temp)
  df_combine<-rbind(df_combine,df_temp)
}

for (colname in c('Tanner_Female_Breast','Tanner_Female_Pubic_Hair','Tanner_Male_Genitals','Tanner_Male_Pubic_Hair')){
  df_combine[!is.na(df_combine[colname]) & df_combine[colname]==6,colname]<-NA
}
df_combine<-df_combine[!is.na(df_combine$Age),]
df_sex<-df_csub[c('ID_pnTTC','Sex')]
df_combine<-dplyr::left_join(df_combine,df_sex,by='ID_pnTTC')

colnames_dst<-colnames(df_combine)
colnames_dst<-c('ID_pnTTC','Sex',colnames_dst[(colnames_dst!='ID_pnTTC') & colnames_dst!='Sex'])
df_combine<-df_combine[,colnames_dst]
df_combine<-df_combine[order(df_combine$Age),]
df_combine<-df_combine[order(df_combine$ID_pnTTC),]
fwrite(df_combine,file.path(path_src,'03_excel','clin_combine01.csv'))


#**************************************************
# Check Tanner ====================================
#**************************************************
df_tanner<-data.frame(fread(file.path(path_src,'03_excel','clin_combine01.csv'),encoding='UTF-8'))
df_tanner<-df_tanner[!is.na(df_tanner$Sex),]
df_tanner$Age<-as.numeric(df_tanner$Age)
for (idx_row in seq(nrow(df_tanner))){
  df_tanner[idx_row,'gonadal']<-as.numeric(max(df_tanner[idx_row,'Tanner_Female_Breast'],df_tanner[idx_row,'Tanner_Male_Genitals'],na.rm=T))
  df_tanner[idx_row,'adrenal']<-as.numeric(max(df_tanner[idx_row,'Tanner_Female_Pubic_Hair'],df_tanner[idx_row,'Tanner_Male_Pubic_Hair'],na.rm=T))
}
df_tanner$Tanner_Female_Breast<-df_tanner$Tanner_Female_Pubic_Hair<-df_tanner$Tanner_Male_Genitals<-df_tanner$Tanner_Male_Pubic_Hair<-NULL

list_id<-sort(unique(df_tanner$ID_pnTTC))
df_check<-data.frame()
df_inconsis<-data.frame()
df_mri_diff<-data.frame()
list_id_inconsis<-list_id_decrease<-NULL
for (id_subj in list_id){
  df_tanner_subj<-df_tanner[df_tanner$ID_pnTTC==id_subj,]
  df_mri_subj<-df_tanner_subj[df_tanner_subj$Source=='M',]
  df_tanner_subj<-df_tanner_subj[df_tanner_subj$Source!='M',]
  # Check if same timing data with different value exists
  list_age_dupl<-sort(unique(df_tanner_subj[duplicated(df_tanner_subj$Age),'Age']))
  inconsis<-F
  for (age_dupl in list_age_dupl){
    df_tanner_dupl<-df_tanner_subj[df_tanner_subj$Age==age_dupl,]
    inconsis_gonadal<-length(unique(df_tanner_dupl$gonadal))>1
    inconsis_adrenal<-length(unique(df_tanner_dupl$gonadal))>1
    if (inconsis_gonadal|inconsis_adrenal){
      inconsis<-T
      df_inconsis<-rbind(df_inconsis,df_tanner_dupl)
    }
  }
  if (inconsis){
    list_id_inconsis<-c(list_id_inconsis,id_subj)
  }
  # Check if decreasing Tanner data exists
  list_gonadal<-df_tanner_subj$gonadal
  list_adrenal<-df_tanner_subj$adrenal
  list_gonadal<-list_gonadal[!is.na(list_gonadal)]
  list_adrenal<-list_adrenal[!is.na(list_adrenal)]
  list_gonadal<-diff(list_gonadal)
  list_adrenal<-diff(list_adrenal)
  decrease<-F
  if (min(list_gonadal)<0|min(list_adrenal)<0){
    list_id_decrease<-c(list_id_decrease,id_subj)
    decrease<-T
  }
  
  # Check if MRI age (calculated from day) differs from Excel MRI age
  mri_diff<-F
  if (nrow(df_mri_subj)>0){
    for (wave in c(1,2)){
      age_mri<-df_mri_subj[df_mri_subj$Wave==wave,'Age']
      age_excel<-df_tanner_subj[df_tanner_subj$Wave==wave & df_tanner_subj$Timing=='M','Age']
      if (length(age_mri)>0 & length(age_excel)>0){
        if (age_mri!=age_excel){
          df_mri_diff<-rbind(df_mri_diff,df_mri_subj[df_mri_subj$Wave==wave,])
          mri_diff<-T
        }
      }
    }
  }
  
  df_check_add<-data.frame(ID_pnTTC=id_subj,Tanner_inconsistency=inconsis,Tanner_decrease=decrease,MRI_Age_diff=mri_diff)
  df_check<-rbind(df_check,df_check_add)
}
df_check_out<-df_check[df_check$Tanner_inconsistency | df_check$Tanner_decrease | df_check$MRI_Age_diff,]
fwrite(df_check_out,file.path(path_src,'03_excel','check01.csv'))

  
#**************************************************
# Visualize Tanner ================================
#**************************************************
df_tanner<-data.frame(fread(file.path(path_src,'03_excel','clin_combine01.csv'),encoding='UTF-8'))
df_tanner<-df_tanner[!is.na(df_tanner$Sex),]
df_tanner$label<-paste(df_tanner$Source,as.character(df_tanner$Wave),df_tanner$Timing,sep='')
df_tanner$Age<-as.numeric(df_tanner$Age)
df_tanner$ID_pnTTC<-factor(df_tanner$ID_pnTTC)
df_mri<-df_tanner[df_tanner$Source=='M',]
df_tanner<-df_tanner[df_tanner$Source!='M',]

for (idx_row in seq(nrow(df_tanner))){
  df_tanner[idx_row,'gonadal']<-as.numeric(max(df_tanner[idx_row,'Tanner_Female_Breast'],df_tanner[idx_row,'Tanner_Male_Genitals'],na.rm=T))
  df_tanner[idx_row,'adrenal']<-as.numeric(max(df_tanner[idx_row,'Tanner_Female_Pubic_Hair'],df_tanner[idx_row,'Tanner_Male_Pubic_Hair'],na.rm=T))
}

list_id<-sort(unique(df_tanner$ID_pnTTC))
list_id_plot<-list_id[1:50]

df_plot_tanner<-df_tanner
df_plot_tanner<-df_plot_tanner[df_plot_tanner$Sex==1,]
df_plot_tanner<-df_plot_tanner[df_plot_tanner$ID_pnTTC %in% list_id_plot,]

df_plot_mri<-df_mri
df_plot_mri<-df_plot_mri[df_plot_mri$Sex==1,]
df_plot_mri<-df_plot_mri[df_plot_mri$ID_pnTTC %in% list_id_plot,]

plt<-(ggplot()
      +geom_point(data=df_plot_tanner,aes(x=Age,y=ID_pnTTC,color=gonadal,fill=gonadal),shape=25,position=position_nudge(y=0.2),alpha=1,size=3)
      +geom_point(data=df_plot_tanner,aes(x=Age,y=ID_pnTTC,color=adrenal,fill=adrenal),shape=24,position=position_nudge(y=-0.2),alpha=1,size=3)
      +geom_point(data=df_plot_mri,aes(x=Age,y=ID_pnTTC),shape=3,color='black',size=3)
      +scale_color_gradientn(colors=matlab.like2(5),lim=c(1,5),name='Tanner')
      +scale_fill_gradientn(colors=matlab.like2(5),lim=c(1,5),name='Tanner')
      +geom_path(data=df_plot_tanner,aes(x=Age,y=ID_pnTTC,group=ID_pnTTC),color='black',alpha=0.1)
      +geom_text(data=df_plot_tanner,mapping=aes(x=Age, y=ID_pnTTC,label=label),
                 size=3,vjust=-1.6, hjust=0.5)
      +xlab('age(month)')
      +ylab('subject ID')
      +theme_classic())
plt



#**************************************************
# Item selection ==================================
#**************************************************

file_ib<-'191027TTC_itembank_labelling.csv'
file_spss<-'spss_data.csv'

df_ib<-data.frame(fread(file.path(path_src,file_ib),encoding='UTF-8'))
df_spss<-data.frame(fread(file.path(path_src,file_spss)))

df_item<-data.frame(item_id=colnames(df_spss))
df_item<-dplyr::left_join(df_item,df_ib,by=c("item_id"=colnames(df_ib)[2]))
fwrite(df_item,file.path(path_src,'items01.csv'))

# Open 'items01.csv' with Excel, manually delete duplicate rows, change column names, Save as 'items03.csv'
df_item<-data.frame(fread(file.path(path_parent,'items03.csv'),encoding='UTF-8'))
df_spss<-data.frame(fread(file.path(path_parent,'spss_data02.csv')))

# Subset df_spss and rename columns
df_item<-df_item[!is.na(df_item$id_label_new),]
df_item<-df_item[df_item$id_label_new!="",]
df_spss<-df_spss[,df_item$id_label]
colnames(df_spss)<-df_item$id_label_new
df_spss[df_spss=='.']<-NA
df_csub<-data.frame(fread(file.path(path_common,'CSUB.csv')))
df_csub<-df_csub[,c("ID_TTC","ID_pnTTC")]
list_colname<-c("ID_TTC","ID_pnTTC",colnames(df_spss)[2:ncol(df_spss)])
df_spss<-dplyr::left_join(df_spss,df_csub,by="ID_TTC")
df_spss<-df_spss[,list_colname]
fwrite(df_spss,file.path(path_parent,'spss_data03.csv'))


#**************************************************
# Date/Age calculation ============================
#**************************************************

# Integrate birth/exam date data
df_spss<-data.frame(fread(file.path(path_parent,'spss_data03.csv')))
df_csub<-data.frame(fread(file.path(path_common,'CSUB.csv')))
df_date_spss<-df_spss[,c('ID_TTC','ID_pnTTC','Birth_Year','Birth_Month','W1_Year','W1_Month','W1_Age_Month','W2_Year','W2_Month','W2_Age_Month')]
df_date_csub<-df_csub[,c('ID_pnTTC','W1_Age_at_MRI','W2_Age_at_MRI')]
df_date<-dplyr::left_join(df_date_spss,df_date_csub,by='ID_pnTTC')
df_date$Birth_date<-sprintf('%04d-%02d-01',df_date$Birth_Year,df_date$Birth_Month)
df_date$W1_MRI_date<-as.Date(df_date$Birth_date)+df_date$W1_Age_at_MRI
df_date$W2_MRI_date<-as.Date(df_date$Birth_date)+df_date$W2_Age_at_MRI
df_date$W1_MRI_Year<-as.numeric(format(df_date$W1_MRI_date,'%Y'))
df_date$W1_MRI_Month<-as.numeric(format(df_date$W1_MRI_date,'%m'))
df_date$W2_MRI_Year<-as.numeric(format(df_date$W2_MRI_date,'%Y'))
df_date$W2_MRI_Month<-as.numeric(format(df_date$W2_MRI_date,'%m'))
df_date$W1_MRI_Age_Month<-12*(df_date$W1_MRI_Year-df_date$Birth_Year)+(df_date$W1_MRI_Month-df_date$Birth_Month)
df_date$W2_MRI_Age_Month<-12*(df_date$W2_MRI_Year-df_date$Birth_Year)+(df_date$W2_MRI_Month-df_date$Birth_Month)
fwrite(df_date,file.path(path_parent,'date_data01.csv'))

# Extract age from date data
df_date<-data.frame(fread(file.path(path_parent,'date_data01.csv')))
df_date<-df_date[,c('ID_TTC','ID_pnTTC','W1_Age_Month','W1_MRI_Age_Month','W2_Age_Month','W2_MRI_Age_Month')]